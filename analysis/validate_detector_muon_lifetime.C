#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLine.h>
#include <TSystem.h>
#include <TString.h>
#include <TPaveText.h>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <utility>

namespace {
int findCoreVolumeID(TFile& f) {
  auto* vm = dynamic_cast<TTree*>(f.Get("VolumeMap"));
  if (!vm) return -999999;

  int volumeID = 0;
  char volumeName[256] = {};
  vm->SetBranchAddress("VolumeID", &volumeID);
  vm->SetBranchAddress("VolumeName", volumeName);

  for (Long64_t i = 0; i < vm->GetEntries(); ++i) {
    vm->GetEntry(i);
    const std::string n(volumeName);
    if (n == "Core" || n.find("Core") != std::string::npos) return volumeID;
  }
  return -999999;
}

struct MuInfo {
  double vertexTime_ns = 0.0;
  double properTime_ns = std::numeric_limits<double>::quiet_NaN();
  double endTime_ns = std::numeric_limits<double>::quiet_NaN();
  double endY_mm = std::numeric_limits<double>::quiet_NaN();
  bool coreDecay = false;
};
}

void validate_detector_muon_lifetime(
    const char* inputName="validation_detector_muon_lifetime_run0.root") {

  constexpr double tauRef_us = 2.1969811; // reference value used only for comparison
  constexpr double michelKEMax_MeV = 52.3; // approximate kinetic endpoint guide

  TFile f(inputName, "READ");
  if (f.IsZombie()) {
    std::cerr << "Cannot open " << inputName << "\n";
    return;
  }

  auto* tr = dynamic_cast<TTree*>(f.Get("Track"));
  if (!tr) {
    std::cerr << "Missing Track tree\n";
    return;
  }

  const int coreID = findCoreVolumeID(f);
  if (coreID == -999999) {
    std::cerr << "Could not resolve Core VolumeID from VolumeMap\n";
    return;
  }

  int event=0, track=0, parent=0, pdg=0, endVolumeID=0;
  float vertexKE=0.f, vertexTime=0.f, endTime=0.f, properNs=0.f, endY=0.f;
  char endProc[256] = {};

  tr->SetBranchAddress("EventID",&event);
  tr->SetBranchAddress("TrackID",&track);
  tr->SetBranchAddress("ParentID",&parent);
  tr->SetBranchAddress("PDG",&pdg);
  tr->SetBranchAddress("VertexKE_MeV",&vertexKE);
  tr->SetBranchAddress("VertexTime_ns",&vertexTime);
  tr->SetBranchAddress("EndTime_ns",&endTime);
  tr->SetBranchAddress("EndProperTime_ns",&properNs);
  tr->SetBranchAddress("EndVolumeID",&endVolumeID);
  tr->SetBranchAddress("EndY_mm",&endY);
  tr->SetBranchAddress("EndProcessName",endProc);

  std::map<std::pair<int,int>, MuInfo> primaryMuPlus;

  // Pass 1: identify primary mu+ and the clean detector-geometry lifetime sample.
  for (Long64_t i=0; i<tr->GetEntries(); ++i) {
    tr->GetEntry(i);
    if (parent != 0 || pdg != -13) continue;

    MuInfo m;
    m.vertexTime_ns = vertexTime;
    m.properTime_ns = properNs;
    m.endTime_ns = endTime;
    m.endY_mm = endY;
    m.coreDecay = (endVolumeID == coreID &&
                   std::string(endProc).find("Decay") != std::string::npos);
    primaryMuPlus[{event,track}] = m;
  }

  TH1D hProper("hDetectorMuProper",
               "Detector-geometry stopped #mu^{+}: proper lifetime;"
               "Proper decay time (#mus);Decays",
               160,0,12);

  TH1D hDelayed("hDetectorMuDelayed",
                "Detector-like delayed timing from Michel e^{+};"
                "t_{e^{+} birth}-t_{#mu birth} (#mus);Events",
                160,0,12);

  TH1D hMichel("hDetectorMuMichel",
               "Michel e^{+} kinetic energy in detector geometry;"
               "Kinetic energy (MeV);e^{+}",
               160,0,55);

  TH1D hStopY("hDetectorMuStopY",
              "Primary #mu^{+} decay position in Water Core;"
              "End y (mm);Decays",
              120,-60,60);

  long long nPrimary = primaryMuPlus.size();
  long long nCoreDecay = 0;
  long long nMichelMatched = 0;
  double maxMichelKE = 0.0;

  for (const auto& kv : primaryMuPlus) {
    const auto& m = kv.second;
    if (!m.coreDecay) continue;
    ++nCoreDecay;
    if (std::isfinite(m.properTime_ns) && m.properTime_ns >= 0)
      hProper.Fill(m.properTime_ns/1000.0);
    if (std::isfinite(m.endY_mm))
      hStopY.Fill(m.endY_mm);
  }

  // Pass 2: select direct Michel positrons from those primary mu+.
  for (Long64_t i=0; i<tr->GetEntries(); ++i) {
    tr->GetEntry(i);

    if (pdg != -11) continue; // e+
    const auto it = primaryMuPlus.find({event,parent});
    if (it == primaryMuPlus.end()) continue;
    const MuInfo& mu = it->second;
    if (!mu.coreDecay) continue;

    ++nMichelMatched;
    hMichel.Fill(vertexKE);
    if (vertexKE > maxMichelKE) maxMichelKE = vertexKE;

    const double dt_us = (double(vertexTime) - mu.vertexTime_ns)/1000.0;
    if (std::isfinite(dt_us) && dt_us >= 0)
      hDelayed.Fill(dt_us);
  }

  if (hProper.GetEntries() < 1000) {
    std::cerr << "Too few clean Core decays for a lifetime fit: "
              << hProper.GetEntries() << "\n";
    return;
  }

  TF1 fProper("fProper","[0]*exp(-x/[1])",0.1,10.0);
  fProper.SetParameters(hProper.GetMaximum(), tauRef_us);
  hProper.Fit(&fProper,"Q0","",0.1,10.0);

  TF1 fDelayed("fDelayed","[0]*exp(-x/[1])",0.1,10.0);
  fDelayed.SetParameters(hDelayed.GetMaximum(), tauRef_us);
  if (hDelayed.GetEntries() >= 1000)
    hDelayed.Fit(&fDelayed,"Q0","",0.1,10.0);

  const double tauProper = fProper.GetParameter(1);
  const double tauProperErr = fProper.GetParError(1);
  const double pullProper = tauProperErr>0 ? (tauProper-tauRef_us)/tauProperErr : 999;

  double tauDelayed = std::numeric_limits<double>::quiet_NaN();
  double tauDelayedErr = std::numeric_limits<double>::quiet_NaN();
  double pullDelayed = std::numeric_limits<double>::quiet_NaN();
  if (hDelayed.GetEntries() >= 1000) {
    tauDelayed = fDelayed.GetParameter(1);
    tauDelayedErr = fDelayed.GetParError(1);
    if (tauDelayedErr>0) pullDelayed=(tauDelayed-tauRef_us)/tauDelayedErr;
  }

  std::cout << "\nDETECTOR-GEOMETRY mu+ LIFETIME VALIDATION\n";
  std::cout << "Input file                  : " << inputName << "\n";
  std::cout << "Primary mu+                 : " << nPrimary << "\n";
  std::cout << "Core decays                 : " << nCoreDecay << "\n";
  std::cout << "Core-decay fraction         : "
            << (nPrimary ? 100.0*double(nCoreDecay)/nPrimary : 0.0) << " %\n";
  std::cout << "Matched Michel e+           : " << nMichelMatched << "\n";
  std::cout << "Reference tau               : " << tauRef_us << " us\n";
  std::cout << std::setprecision(8)
            << "Truth proper-time tau       : " << tauProper
            << " +/- " << tauProperErr << " us\n"
            << "Truth proper-time pull      : " << pullProper << "\n";

  if (hDelayed.GetEntries() >= 1000) {
    std::cout << "Michel-delay tau            : " << tauDelayed
              << " +/- " << tauDelayedErr << " us\n"
              << "Michel-delay pull           : " << pullDelayed << "\n";
  } else {
    std::cout << "Michel-delay tau            : insufficient entries\n";
  }

  std::cout << "Max Michel e+ KE            : " << maxMichelKE << " MeV\n";

  const bool passCore = nPrimary>0 && double(nCoreDecay)/nPrimary > 0.98;
  const bool passProper = std::fabs(pullProper) < 3.0;
  const bool passMichelEndpoint = maxMichelKE > 50.0 && maxMichelKE < 53.5;
  const bool passDelayed = (hDelayed.GetEntries() < 1000) || std::fabs(pullDelayed) < 3.0;
  const bool allPass = passCore && passProper && passMichelEndpoint && passDelayed;

  std::cout << "Core stopping               : " << (passCore ? "PASS":"FAIL") << "\n";
  std::cout << "Proper lifetime             : " << (passProper ? "PASS":"FAIL") << "\n";
  std::cout << "Michel endpoint             : " << (passMichelEndpoint ? "PASS":"FAIL") << "\n";
  std::cout << "Michel delayed timing       : " << (passDelayed ? "PASS":"FAIL") << "\n";
  std::cout << "DETECTOR-GEOMETRY MUON LIFETIME VALIDATION: "
            << (allPass ? "PASS":"FAIL") << "\n";

  gSystem->mkdir("validation_plots", kTRUE);

  TCanvas c("cDetectorMuonLifetime",
            "Detector-geometry muon lifetime validation", 1300,900);
  c.Divide(2,2);

  c.cd(1);
  hProper.Draw("HIST");
  fProper.SetLineWidth(2);
  fProper.Draw("SAME");

  c.cd(2);
  hDelayed.Draw("HIST");
  if (hDelayed.GetEntries() >= 1000) {
    fDelayed.SetLineWidth(2);
    fDelayed.Draw("SAME");
  }

  c.cd(3);
  hMichel.Draw("HIST");
  TLine endpoint(michelKEMax_MeV,0,michelKEMax_MeV,hMichel.GetMaximum()*1.05);
  endpoint.SetLineStyle(2);
  endpoint.Draw("SAME");

  c.cd(4);
  hStopY.SetStats(true);
  hStopY.Draw("HIST");

  c.SaveAs("validation_plots/detector_muon_lifetime.png");
  c.SaveAs("validation_plots/detector_muon_lifetime.pdf");
}
