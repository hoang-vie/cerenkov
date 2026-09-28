#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TPad.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

// Core-wall truth analysis for the general-v1.6 schema.
//
// IMPORTANT:
//   * This macro does NOT use BoundaryDebug.
//   * It reads the compact per-photon first-Core-wall summary already stored in
//     Photon: ReachedCoreWall, FirstWallTime_ns, FirstWallX/Y/Z_mm.
//   * No Geant4 geometry, optical process, or ROOT schema change is required.
//   * Face assignment is derived offline from the stored wall XYZ and RunConfig
//     CoreX/Y/Z. Ambiguous edge/corner hits are counted and excluded from the
//     six-face map instead of being assigned arbitrarily.

namespace corewall {

struct TrackKey {
  int run = 0;
  int event = 0;
  int track = 0;
  bool operator==(const TrackKey& o) const {
    return run == o.run && event == o.event && track == o.track;
  }
};

struct TrackKeyHash {
  std::size_t operator()(const TrackKey& k) const {
    // Stable integer mixing; only used for an in-memory lookup table.
    uint64_t x = uint32_t(k.run);
    x = (x * 0x9e3779b185ebca87ULL) ^ uint32_t(k.event);
    x = (x * 0xc2b2ae3d27d4eb4fULL) ^ uint32_t(k.track);
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    return static_cast<std::size_t>(x);
  }
};

struct EventKey {
  int run = 0;
  int event = 0;
  bool operator==(const EventKey& o) const {
    return run == o.run && event == o.event;
  }
};

struct EventKeyHash {
  std::size_t operator()(const EventKey& k) const {
    uint64_t x = (uint64_t(uint32_t(k.run)) << 32) | uint32_t(k.event);
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    return static_cast<std::size_t>(x);
  }
};

struct TrackInfo {
  int parent = 0;
  int root = 0;
  int pdg = 0;
  float vertexTimeNs = 0.f;
};

struct CoreGeometry {
  double hx = 0.0;
  double hy = 0.0;
  double hz = 0.0;
  bool valid = false;
};

struct FaceHit {
  int face = -1;          // canonical project convention: 0=+X,1=-X,2=+Z,3=-Z,4=+Y,5=-Y
  double u = 0.0;
  double v = 0.0;
  bool ambiguous = false; // true for edge/corner: >=2 coordinates lie on a Core surface
  bool offSurface = false;
};

struct RangeAccumulator {
  bool set = false;
  double lo = 0.0;
  double hi = 0.0;
  void Add(double x) {
    if (!std::isfinite(x)) return;
    if (!set) { lo = hi = x; set = true; return; }
    lo = std::min(lo, x);
    hi = std::max(hi, x);
  }
};

bool PdgMatch(int pdg, int requested, bool combineCharge) {
  if (requested == 0) return true;
  return combineCharge ? std::abs(pdg) == std::abs(requested) : pdg == requested;
}

TString ParticleLabel(int pdg, bool combineCharge) {
  if (pdg == 0) return "all particles";
  if (combineCharge) {
    if (std::abs(pdg) == 11) return "e^{#pm}";
    if (std::abs(pdg) == 13) return "#mu^{#pm}";
    if (std::abs(pdg) == 211) return "#pi^{#pm}";
    if (std::abs(pdg) == 321) return "K^{#pm}";
    if (std::abs(pdg) == 2212) return "p/#bar{p}";
    return Form("|PDG|=%d", std::abs(pdg));
  }
  if (pdg == 11) return "e^{-}";
  if (pdg == -11) return "e^{+}";
  if (pdg == 13) return "#mu^{-}";
  if (pdg == -13) return "#mu^{+}";
  if (pdg == 211) return "#pi^{+}";
  if (pdg == -211) return "#pi^{-}";
  if (pdg == 321) return "K^{+}";
  if (pdg == -321) return "K^{-}";
  if (pdg == 2212) return "p";
  if (pdg == -2212) return "#bar{p}";
  return Form("PDG=%d", pdg);
}

bool IsAllowedSourceMode(const std::string& mode) {
  return mode == "all" || mode == "primaryEmitter" ||
         mode == "emitter" || mode == "ancestor";
}

bool IsAllowedTimeRef(const std::string& ref) {
  return ref == "none" || ref == "absolute" || ref == "photonBirth" ||
         ref == "emitterBirth" || ref == "rootPrimaryBirth" ||
         ref == "topTrigger";
}

bool SelectPhotonSource(
    const std::unordered_map<TrackKey, TrackInfo, TrackKeyHash>& tracks,
    int runID, int eventID, int parentTrackID, int sourceTrackID,
    int selectedPDG, bool combineCharge, const std::string& sourceMode) {
  if (sourceMode == "all" || selectedPDG == 0) return true;

  if (sourceMode == "ancestor") {
    auto sourceIt = tracks.find({runID, eventID, sourceTrackID});
    if (sourceIt == tracks.end()) return false;
    auto rootIt = tracks.find({runID, eventID, sourceIt->second.root});
    return rootIt != tracks.end() &&
           PdgMatch(rootIt->second.pdg, selectedPDG, combineCharge);
  }

  auto parentIt = tracks.find({runID, eventID, parentTrackID});
  if (parentIt == tracks.end()) return false; // optical parent / incomplete detail

  if (sourceMode == "emitter")
    return PdgMatch(parentIt->second.pdg, selectedPDG, combineCharge);

  // primaryEmitter
  return parentIt->second.parent == 0 &&
         PdgMatch(parentIt->second.pdg, selectedPDG, combineCharge);
}

int ImmediateEmitterPDG(
    const std::unordered_map<TrackKey, TrackInfo, TrackKeyHash>& tracks,
    int runID, int eventID, int parentTrackID) {
  auto it = tracks.find({runID, eventID, parentTrackID});
  return (it == tracks.end()) ? 0 : it->second.pdg;
}

bool TimeValue(
    const std::string& ref,
    const std::unordered_map<TrackKey, TrackInfo, TrackKeyHash>& tracks,
    const std::unordered_map<EventKey, float, EventKeyHash>& topTriggerTime,
    int runID, int eventID, int parentTrackID, int sourceTrackID,
    double wallTimeNs, double photonBirthTimeNs,
    double& out) {
  if (ref == "none" || ref == "absolute") {
    out = wallTimeNs;
    return std::isfinite(out);
  }
  if (ref == "photonBirth") {
    out = wallTimeNs - photonBirthTimeNs;
    return std::isfinite(out);
  }
  if (ref == "emitterBirth") {
    auto it = tracks.find({runID, eventID, parentTrackID});
    if (it == tracks.end()) return false;
    out = wallTimeNs - it->second.vertexTimeNs;
    return std::isfinite(out);
  }
  if (ref == "rootPrimaryBirth") {
    auto sourceIt = tracks.find({runID, eventID, sourceTrackID});
    if (sourceIt == tracks.end()) return false;
    auto rootIt = tracks.find({runID, eventID, sourceIt->second.root});
    if (rootIt == tracks.end()) return false;
    out = wallTimeNs - rootIt->second.vertexTimeNs;
    return std::isfinite(out);
  }
  if (ref == "topTrigger") {
    auto it = topTriggerTime.find({runID, eventID});
    if (it == topTriggerTime.end()) return false;
    out = wallTimeNs - it->second;
    return std::isfinite(out);
  }
  return false;
}

FaceHit ClassifyFace(double x, double y, double z,
                     const CoreGeometry& g, double tol) {
  FaceHit h;
  if (!g.valid) { h.offSurface = true; return h; }

  const double rx = std::abs(std::abs(x) - g.hx);
  const double ry = std::abs(std::abs(y) - g.hy);
  const double rz = std::abs(std::abs(z) - g.hz);

  const bool onX = rx <= tol;
  const bool onY = ry <= tol;
  const bool onZ = rz <= tol;
  const int nOn = int(onX) + int(onY) + int(onZ);
  if (nOn == 0) { h.offSurface = true; return h; }
  if (nOn > 1) { h.ambiguous = true; return h; }

  if (onX) {
    if (x >= 0) { h.face = 0; h.u =  y; h.v = z; } // +X
    else        { h.face = 1; h.u = -y; h.v = z; } // -X
  } else if (onZ) {
    if (z >= 0) { h.face = 2; h.u =  x; h.v =  y; } // +Z
    else        { h.face = 3; h.u =  x; h.v = -y; } // -Z
  } else {
    if (y >= 0) { h.face = 4; h.u =  x; h.v = -z; } // +Y
    else        { h.face = 5; h.u =  x; h.v =  z; } // -Y
  }
  return h;
}

const char* FaceName(int face) {
  static const char* names[6] = {"+X", "-X", "+Z", "-Z", "+Y", "-Y"};
  return (face >= 0 && face < 6) ? names[face] : "?";
}

void FaceRanges(int face, const CoreGeometry& g,
                double& uMin, double& uMax, double& vMin, double& vMax) {
  if (face == 0 || face == 1) { // +/-X: u=y, v=z
    uMin = -g.hy; uMax = g.hy; vMin = -g.hz; vMax = g.hz;
  } else if (face == 2 || face == 3) { // +/-Z: u=x, v=y
    uMin = -g.hx; uMax = g.hx; vMin = -g.hy; vMax = g.hy;
  } else { // +/-Y: u=x, v=z
    uMin = -g.hx; uMax = g.hx; vMin = -g.hz; vMax = g.hz;
  }
}

std::pair<double,double> SafeHistRange(const RangeAccumulator& r,
                                       double fallbackLo, double fallbackHi) {
  if (!r.set) return {fallbackLo, fallbackHi};
  if (r.hi <= r.lo) {
    const double pad = std::max(1.0, std::abs(r.lo) * 0.05);
    return {r.lo - pad, r.hi + pad};
  }
  const double span = r.hi - r.lo;
  return {r.lo - 0.05 * span, r.hi + 0.05 * span};
}

void DrawOrigami(TCanvas& c, TH2D* hFace[6], const TString& title) {
  c.Clear();
  c.SetFillColor(0);
  c.cd();

  // Robust six-face "cross" layout.  All pad coordinates are strictly inside
  // [0,1], avoiding the ROOT TPad error seen in Step-1.
  //
  //            +Y
  //            +Z
  //       -X   -Y   +X
  //            -Z
  //
  // The histograms all share the same linear z scale.  We intentionally draw
  // without six separate palette bars, which otherwise consume most of the
  // useful area in small pads.
  const int slotFace[6] = {4, 2, 1, 5, 0, 3}; // +Y,+Z,-X,-Y,+X,-Z
  const char* padName[6] = {"pTop","pFront","pLeft","pCenter","pRight","pBack"};

  const double xy[6][4] = {
    {0.39,0.73,0.61,0.95}, // +Y
    {0.39,0.51,0.61,0.73}, // +Z
    {0.17,0.29,0.39,0.51}, // -X
    {0.39,0.29,0.61,0.51}, // -Y
    {0.61,0.29,0.83,0.51}, // +X
    {0.39,0.07,0.61,0.29}  // -Z
  };

  double commonMax = 0.0;
  for (int f=0; f<6; ++f) commonMax = std::max(commonMax, hFace[f]->GetMaximum());

  for (int s=0; s<6; ++s) {
    auto* p = new TPad(padName[s], padName[s], xy[s][0], xy[s][1], xy[s][2], xy[s][3]);
    p->SetFillColor(0);
    p->SetFrameFillColor(0);
    p->SetRightMargin(0.05);
    p->SetLeftMargin(0.18);
    p->SetBottomMargin(0.18);
    p->SetTopMargin(0.15);
    p->Draw();
    p->cd();

    TH2D* h2 = hFace[slotFace[s]];
    h2->SetStats(0);
    h2->SetMinimum(0.0);
    if (commonMax > 0) h2->SetMaximum(commonMax);
    h2->GetXaxis()->SetTitleSize(0.055);
    h2->GetYaxis()->SetTitleSize(0.055);
    h2->GetXaxis()->SetLabelSize(0.045);
    h2->GetYaxis()->SetLabelSize(0.045);
    h2->GetXaxis()->SetTitleOffset(1.05);
    h2->GetYaxis()->SetTitleOffset(1.20);
    h2->SetTitleSize(0.060);
    h2->Draw("COL");

    c.cd();
  }

  TLatex tx;
  tx.SetNDC();
  tx.SetTextAlign(22);
  tx.SetTextSize(0.034);
  tx.DrawLatex(0.50,0.985,title);
  tx.SetTextSize(0.022);
  tx.DrawLatex(0.50,0.025,
               Form("Common linear colour scale across all faces; max bin = %.0f photons",
                    commonMax));
  c.Update();
}

} // namespace corewall

void analyze_corewall(const char* inputName,
                      int selectedEventID = -1,
                      int selectedPDG = 0,
                      bool combineCharge = true,
                      const char* sourceModeC = "all",
                      const char* timeRefC = "none",
                      double timeMin_ns = -1e99,
                      double timeMax_ns =  1e99,
                      double wallTolerance_mm = 1e-3) {
  using namespace corewall;

  const std::string sourceMode = sourceModeC ? sourceModeC : "all";
  const std::string timeRef = timeRefC ? timeRefC : "none";
  if (!IsAllowedSourceMode(sourceMode)) {
    std::cerr << "sourceMode must be all, primaryEmitter, emitter, or ancestor\n";
    return;
  }
  if (!IsAllowedTimeRef(timeRef)) {
    std::cerr << "timeRef must be none, absolute, photonBirth, emitterBirth, rootPrimaryBirth, or topTrigger\n";
    return;
  }
  if (!(timeMin_ns < timeMax_ns)) {
    std::cerr << "timeMin_ns must be smaller than timeMax_ns\n";
    return;
  }
  if (!(wallTolerance_mm > 0)) {
    std::cerr << "wallTolerance_mm must be > 0\n";
    return;
  }

  TFile* f = TFile::Open(inputName, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << inputName << "\n";
    return;
  }

  auto* tPho = dynamic_cast<TTree*>(f->Get("Photon"));
  auto* tTrack = dynamic_cast<TTree*>(f->Get("Track"));
  auto* tCfg = dynamic_cast<TTree*>(f->Get("RunConfig"));
  auto* tTrig = dynamic_cast<TTree*>(f->Get("TriggerEvent"));
  if (!tPho || !tTrack || !tCfg) {
    std::cerr << "Need Photon, Track and RunConfig from the general schema.\n";
    f->Close();
    return;
  }

  // Geometry/provenance: current project convention is one simulation run per ROOT file.
  int cfgRun = 0, boundaryDebug = -1;
  float coreX = 0.f, coreY = 0.f, coreZ = 0.f;
  tCfg->SetBranchAddress("RunID", &cfgRun);
  tCfg->SetBranchAddress("CoreX_mm", &coreX);
  tCfg->SetBranchAddress("CoreY_mm", &coreY);
  tCfg->SetBranchAddress("CoreZ_mm", &coreZ);
  if (tCfg->GetBranch("BoundaryDebug")) tCfg->SetBranchAddress("BoundaryDebug", &boundaryDebug);
  if (tCfg->GetEntries() < 1) {
    std::cerr << "RunConfig is empty.\n";
    f->Close();
    return;
  }
  tCfg->GetEntry(0);
  CoreGeometry geom{0.5*coreX, 0.5*coreY, 0.5*coreZ,
                    coreX>0 && coreY>0 && coreZ>0};
  if (!geom.valid) {
    std::cerr << "Invalid Core dimensions in RunConfig: "
              << coreX << " x " << coreY << " x " << coreZ << " mm\n";
    f->Close();
    return;
  }
  if (tCfg->GetEntries() > 1) {
    std::cerr << "WARNING: RunConfig has " << tCfg->GetEntries()
              << " rows; this Step-1 analyzer uses the first Core geometry.\n";
  }

  // Build non-optical truth lookup for emitter/ancestor selection and timing references.
  std::unordered_map<TrackKey, TrackInfo, TrackKeyHash> tracks;
  int trRun=0, trEvent=0, trID=0, trParent=0, trRoot=0, trPDG=0;
  float trVertexTime=0.f;
  tTrack->SetBranchAddress("RunID", &trRun);
  tTrack->SetBranchAddress("EventID", &trEvent);
  tTrack->SetBranchAddress("TrackID", &trID);
  tTrack->SetBranchAddress("ParentID", &trParent);
  tTrack->SetBranchAddress("RootPrimaryTrackID", &trRoot);
  tTrack->SetBranchAddress("PDG", &trPDG);
  tTrack->SetBranchAddress("VertexTime_ns", &trVertexTime);
  tracks.reserve(static_cast<std::size_t>(tTrack->GetEntries()*1.2 + 16));
  for (Long64_t i=0; i<tTrack->GetEntries(); ++i) {
    tTrack->GetEntry(i);
    tracks[{trRun,trEvent,trID}] = {trParent,trRoot,trPDG,trVertexTime};
  }

  std::unordered_map<EventKey, float, EventKeyHash> topTriggerTime;
  if (tTrig) {
    int tgRun=0, tgEvent=0;
    float tTop=0.f;
    tTrig->SetBranchAddress("RunID", &tgRun);
    tTrig->SetBranchAddress("EventID", &tgEvent);
    tTrig->SetBranchAddress("TimeTop_ns", &tTop);
    for (Long64_t i=0; i<tTrig->GetEntries(); ++i) {
      tTrig->GetEntry(i);
      if (std::isfinite(tTop) && tTop >= 0) topTriggerTime[{tgRun,tgEvent}] = tTop;
    }
  }

  // Photon branch bindings: all are present in general-v1.6 compact Photon output.
  int pRun=0, pEvent=0, photonID=0, parentTrackID=0, sourceTrackID=0;
  int reached=0;
  float wavelength=0.f, birthTime=0.f, wallTime=0.f, wx=0.f, wy=0.f, wz=0.f;
  tPho->SetBranchAddress("RunID", &pRun);
  tPho->SetBranchAddress("EventID", &pEvent);
  tPho->SetBranchAddress("PhotonID", &photonID);
  tPho->SetBranchAddress("ParentTrackID", &parentTrackID);
  tPho->SetBranchAddress("SourceTrackID", &sourceTrackID);
  tPho->SetBranchAddress("Wavelength_nm", &wavelength);
  tPho->SetBranchAddress("BirthTime_ns", &birthTime);
  tPho->SetBranchAddress("ReachedCoreWall", &reached);
  tPho->SetBranchAddress("FirstWallTime_ns", &wallTime);
  tPho->SetBranchAddress("FirstWallX_mm", &wx);
  tPho->SetBranchAddress("FirstWallY_mm", &wy);
  tPho->SetBranchAddress("FirstWallZ_mm", &wz);

  const bool useTimeCut = (timeRef != "none") &&
                          (timeMin_ns > -0.9e99 || timeMax_ns < 0.9e99);

  long long nPhotonRows = 0;
  long long nSourceSelected = 0;
  long long nReached = 0;
  long long nTimeRefMissing = 0;
  long long nTimeRejected = 0;
  long long nAfterTimeSelection = 0;
  long long nAcceptedMap = 0;
  long long nNegativeFlightTime = 0;
  long long nAmbiguous = 0;
  long long nOffSurface = 0;
  long long faceCounts[6] = {0,0,0,0,0,0};
  std::map<int,long long> emitterPdgCounts;
  RangeAccumulator rAbs, rFlight, rEmitter, rRoot, rTrigger;

  auto basicSelected = [&]() -> bool {
    if (selectedEventID >= 0 && pEvent != selectedEventID) return false;
    return SelectPhotonSource(tracks,pRun,pEvent,parentTrackID,sourceTrackID,
                              selectedPDG,combineCharge,sourceMode);
  };

  // Pass 1: QA/counting/ranges. No large per-photon vector is kept in memory.
  for (Long64_t i=0; i<tPho->GetEntries(); ++i) {
    tPho->GetEntry(i);
    ++nPhotonRows;
    if (!basicSelected()) continue;
    ++nSourceSelected;
    if (!reached) continue;
    ++nReached;

    double tSel = wallTime;
    bool haveTimeRef = TimeValue(timeRef,tracks,topTriggerTime,pRun,pEvent,
                                 parentTrackID,sourceTrackID,wallTime,birthTime,tSel);
    if (!haveTimeRef) {
      if (useTimeCut) { ++nTimeRefMissing; continue; }
    } else if (useTimeCut && (tSel < timeMin_ns || tSel >= timeMax_ns)) {
      ++nTimeRejected;
      continue;
    }
    ++nAfterTimeSelection;

    const int emitterPDG = ImmediateEmitterPDG(tracks,pRun,pEvent,parentTrackID);
    emitterPdgCounts[emitterPDG]++;

    const double flightDt = double(wallTime)-double(birthTime);
    if (flightDt < -1e-6) ++nNegativeFlightTime;
    rAbs.Add(wallTime);
    rFlight.Add(flightDt);
    double tmp=0;
    if (TimeValue("emitterBirth",tracks,topTriggerTime,pRun,pEvent,parentTrackID,sourceTrackID,wallTime,birthTime,tmp)) rEmitter.Add(tmp);
    if (TimeValue("rootPrimaryBirth",tracks,topTriggerTime,pRun,pEvent,parentTrackID,sourceTrackID,wallTime,birthTime,tmp)) rRoot.Add(tmp);
    if (TimeValue("topTrigger",tracks,topTriggerTime,pRun,pEvent,parentTrackID,sourceTrackID,wallTime,birthTime,tmp)) rTrigger.Add(tmp);

    const FaceHit hit = ClassifyFace(wx,wy,wz,geom,wallTolerance_mm);
    if (hit.offSurface) { ++nOffSurface; continue; }
    if (hit.ambiguous) { ++nAmbiguous; continue; }
    ++nAcceptedMap;
    faceCounts[hit.face]++;
  }

  // Histograms use actual Core dimensions; there is no hard-coded 100 mm cube assumption.
  TH2D* hFace[6] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
  for (int face=0; face<6; ++face) {
    double u0,u1,v0,v1;
    FaceRanges(face,geom,u0,u1,v0,v1);
    const int nU = std::max(20, int(std::round((u1-u0)/1.0))); // default display bin ~1 mm
    const int nV = std::max(20, int(std::round((v1-v0)/1.0)));
    hFace[face] = new TH2D(Form("hCoreWallFace%d",face),
                           Form("Core wall %s;Local U (mm);Local V (mm)",FaceName(face)),
                           nU,u0,u1,nV,v0,v1);
    hFace[face]->SetStats(0);
  }

  TH1D hFaceCount("hCoreWallFaceCount","First Core-wall hits by face;Face;Photons",6,-0.5,5.5);
  for (int fce=0; fce<6; ++fce) {
    hFaceCount.GetXaxis()->SetBinLabel(fce+1,FaceName(fce));
    hFaceCount.SetBinContent(fce+1,faceCounts[fce]);
  }
  hFaceCount.SetStats(0);

  TH1D hWavelength("hCoreWallWavelength","First Core-wall photon wavelength;Wavelength (nm);Photons",140,200,900);

  auto absRange = SafeHistRange(rAbs,0,100);
  auto flightRange = SafeHistRange(rFlight,0,20);
  auto emitRange = SafeHistRange(rEmitter,0,20);
  auto rootRange = SafeHistRange(rRoot,0,100);
  auto trigRange = SafeHistRange(rTrigger,-10,50);
  TH1D hAbs("hCoreWallAbsoluteTime","First Core-wall global time;Global time (ns);Photons",160,absRange.first,absRange.second);
  TH1D hFlight("hCoreWallPhotonFlightTime","Photon birth to first Core wall;t_{wall}-t_{#gamma,birth} (ns);Photons",160,flightRange.first,flightRange.second);
  TH1D hEmitter("hCoreWallEmitterDelay","Emitter birth to first Core wall;t_{wall}-t_{emitter,birth} (ns);Photons",160,emitRange.first,emitRange.second);
  TH1D hRoot("hCoreWallRootPrimaryDelay","Root-primary birth to first Core wall;t_{wall}-t_{root,birth} (ns);Photons",160,rootRange.first,rootRange.second);
  TH1D hTrigger("hCoreWallTopTriggerDelay","Top-trigger to first Core wall;t_{wall}-t_{top} (ns);Photons",160,trigRange.first,trigRange.second);

  std::vector<std::pair<int,long long>> pdgVec(emitterPdgCounts.begin(),emitterPdgCounts.end());
  std::sort(pdgVec.begin(),pdgVec.end(),[](const auto& a,const auto& b){return a.second>b.second;});
  TH1D hEmitterPDG("hCoreWallEmitterPDG","Immediate emitter composition at first Core wall;Emitter PDG;Photons",
                   std::max<int>(1,pdgVec.size()),0,std::max<int>(1,pdgVec.size()));
  hEmitterPDG.SetStats(0);
  for (std::size_t i=0; i<pdgVec.size(); ++i) {
    hEmitterPDG.SetBinContent(i+1,pdgVec[i].second);
    hEmitterPDG.GetXaxis()->SetBinLabel(i+1,Form("%d",pdgVec[i].first));
  }

  // Pass 2: fill plots with the exact same selection used in pass 1.
  for (Long64_t i=0; i<tPho->GetEntries(); ++i) {
    tPho->GetEntry(i);
    if (!basicSelected() || !reached) continue;

    double tSel = wallTime;
    bool haveTimeRef = TimeValue(timeRef,tracks,topTriggerTime,pRun,pEvent,
                                 parentTrackID,sourceTrackID,wallTime,birthTime,tSel);
    if (useTimeCut) {
      if (!haveTimeRef || tSel < timeMin_ns || tSel >= timeMax_ns) continue;
    }

    const FaceHit hit = ClassifyFace(wx,wy,wz,geom,wallTolerance_mm);
    if (hit.offSurface || hit.ambiguous) continue;
    hFace[hit.face]->Fill(hit.u,hit.v);
    hWavelength.Fill(wavelength);
    hAbs.Fill(wallTime);
    hFlight.Fill(double(wallTime)-double(birthTime));
    double tmp=0;
    if (TimeValue("emitterBirth",tracks,topTriggerTime,pRun,pEvent,parentTrackID,sourceTrackID,wallTime,birthTime,tmp)) hEmitter.Fill(tmp);
    if (TimeValue("rootPrimaryBirth",tracks,topTriggerTime,pRun,pEvent,parentTrackID,sourceTrackID,wallTime,birthTime,tmp)) hRoot.Fill(tmp);
    if (TimeValue("topTrigger",tracks,topTriggerTime,pRun,pEvent,parentTrackID,sourceTrackID,wallTime,birthTime,tmp)) hTrigger.Fill(tmp);
  }

  TString base(inputName);
  Ssiz_t slash = base.Last('/');
  if (slash != kNPOS) base.Remove(0,slash+1);
  base.ReplaceAll(".root","");
  TString evtTag = selectedEventID >= 0 ? Form("evt%d",selectedEventID) : "evtALL";
  TString srcTag = Form("PDG%d_%s",selectedPDG,sourceMode.c_str());
  TString timeTag = timeRef.c_str();
  if (useTimeCut) {
    auto cleanNumber = [](double x) {
      TString s = Form("%.6g",x);
      s.ReplaceAll("-","m");
      s.ReplaceAll("+","p");
      s.ReplaceAll(".","p");
      return s;
    };
    timeTag += Form("_%s_to_%s_ns",cleanNumber(timeMin_ns).Data(),cleanNumber(timeMax_ns).Data());
  }
  TString pdfOut = Form("corewall_%s_%s_%s_%s.pdf",base.Data(),evtTag.Data(),srcTag.Data(),timeTag.Data());
  TString rootOut = Form("corewall_%s_%s_%s_%s.root",base.Data(),evtTag.Data(),srcTag.Data(),timeTag.Data());

  TCanvas c("cCoreWall","Core-wall truth analysis",1400,900);
  c.SetFillColor(0);
  c.Print(pdfOut+"[");

  c.Clear();
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.06,0.94,"CORE-WALL TRUTH ANALYSIS -- compact Photon first-wall data");
  tx.DrawLatex(0.06,0.88,Form("Input: %s",inputName));
  tx.DrawLatex(0.06,0.83,Form("Core: %.3f x %.3f x %.3f mm | wall tolerance %.6g mm",coreX,coreY,coreZ,wallTolerance_mm));
  tx.DrawLatex(0.06,0.78,Form("Event: %s | source: %s | selected PDG: %s",evtTag.Data(),sourceMode.c_str(),ParticleLabel(selectedPDG,combineCharge).Data()));
  if (useTimeCut)
    tx.DrawLatex(0.06,0.73,Form("Time cut: %s in [%.6g, %.6g) ns",timeRef.c_str(),timeMin_ns,timeMax_ns));
  else
    tx.DrawLatex(0.06,0.73,Form("Time cut: none (timeRef argument = %s)",timeRef.c_str()));
  tx.DrawLatex(0.06,0.67,Form("Photon rows: %lld | source-selected: %lld | reached Core wall: %lld",nPhotonRows,nSourceSelected,nReached));
  tx.DrawLatex(0.06,0.62,Form("After time selection: %lld | mapped: %lld | ambiguous edge/corner: %lld | off-surface: %lld",
                                  nAfterTimeSelection,nAcceptedMap,nAmbiguous,nOffSurface));
  tx.DrawLatex(0.06,0.57,Form("Time-rejected: %lld | missing requested time reference: %lld | negative (t_{wall}-t_{birth}): %lld",
                                  nTimeRejected,nTimeRefMissing,nNegativeFlightTime));
  tx.DrawLatex(0.06,0.51,Form("RunConfig BoundaryDebug = %d (analysis does not read BoundaryDebug)",boundaryDebug));
  tx.DrawLatex(0.06,0.45,"Face convention: 0=+X, 1=-X, 2=+Z, 3=-Z, 4=+Y, 5=-Y");
  tx.DrawLatex(0.06,0.39,"Edge/corner hits are intentionally excluded from the map instead of arbitrary face assignment.");
  tx.DrawLatex(0.06,0.33,"This macro analyzes the Photon rows that exist in the file; thinAll/triggered detail may be trigger-selected.");
  c.Print(pdfOut);

  DrawOrigami(c,hFace,Form("First Core-wall truth map -- %s -- %s",evtTag.Data(),ParticleLabel(selectedPDG,combineCharge).Data()));
  c.Print(pdfOut);

  c.Clear(); c.Divide(2,2);
  c.cd(1); hFaceCount.Draw("HIST");
  c.cd(2); hWavelength.Draw("HIST");
  c.cd(3); hEmitterPDG.Draw("HIST");
  c.cd(4); hFlight.Draw("HIST");
  c.Print(pdfOut);

  c.Clear(); c.Divide(2,2);
  c.cd(1); hAbs.Draw("HIST");
  c.cd(2); hEmitter.Draw("HIST");
  c.cd(3); hRoot.Draw("HIST");
  c.cd(4); hTrigger.Draw("HIST");
  c.Print(pdfOut);

  c.Print(pdfOut+"]");

  TFile fOut(rootOut,"RECREATE");
  for (auto* h : hFace) h->Write();
  hFaceCount.Write();
  hWavelength.Write();
  hEmitterPDG.Write();
  hAbs.Write();
  hFlight.Write();
  hEmitter.Write();
  hRoot.Write();
  hTrigger.Write();
  fOut.Close();

  std::cout << "\n=== CORE-WALL STEP-1 SUMMARY ===\n";
  std::cout << "Input: " << inputName << "\n";
  std::cout << "Core [mm]: " << coreX << " x " << coreY << " x " << coreZ << "\n";
  std::cout << "BoundaryDebug in input: " << boundaryDebug << " (NOT used by this analysis)\n";
  std::cout << "Photon rows=" << nPhotonRows
            << " sourceSelected=" << nSourceSelected
            << " reachedWall=" << nReached
            << " afterTimeSelection=" << nAfterTimeSelection
            << " mapped=" << nAcceptedMap
            << " ambiguous=" << nAmbiguous
            << " offSurface=" << nOffSurface << "\n";
  const long long mapClosure = nAcceptedMap + nAmbiguous + nOffSurface;
  std::cout << "Mapping closure: mapped+ambiguous+offSurface=" << mapClosure
            << " / afterTimeSelection=" << nAfterTimeSelection
            << (mapClosure==nAfterTimeSelection ? " [OK]" : " [MISMATCH]") << "\n";
  std::cout << "Negative photon flight times (t_wall-t_birth < -1e-6 ns): "
            << nNegativeFlightTime << "\n";
  std::cout << "Face counts:";
  for (int face=0; face<6; ++face) {
    const double frac = nAcceptedMap ? 100.0*double(faceCounts[face])/double(nAcceptedMap) : 0.0;
    std::cout << " " << FaceName(face) << "=" << faceCounts[face]
              << "(" << std::fixed << std::setprecision(3) << frac << "%)";
  }
  std::cout << "\n";
  std::cout << "Wrote " << pdfOut << "\n";
  std::cout << "Wrote " << rootOut << "\n";
  const bool qaReady = nReached>0 && nOffSurface==0 && mapClosure==nAfterTimeSelection && nNegativeFlightTime==0;
  std::cout << "CORE-WALL STEP-1.1: " << (qaReady ? "READY_FOR_REVIEW" : "CHECK_QA") << "\n";

  f->Close();
}