#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

uint64_t Key(int eventID, int trackID) {
  return (uint64_t(uint32_t(eventID)) << 32) | uint32_t(trackID);
}

struct TrackInfo {
  int parent = 0;
  int pdg = 0;
};

struct EventSummary {
  long long generated = 0;
  long long wall = 0;
  long long wallPrimaryMuon = 0;
  long long wallElectronEmitter = 0;
  long long wallOther = 0;

  double wallTimeSum = 0.0;
  double wallTimeMin = std::numeric_limits<double>::infinity();
  double wallTimeMax = -std::numeric_limits<double>::infinity();

  double flightSum = 0.0;
  double flightMin = std::numeric_limits<double>::infinity();
  double flightMax = -std::numeric_limits<double>::infinity();
  long long flightN = 0;
};

bool IsPrimaryMuon(const std::unordered_map<uint64_t,TrackInfo>& tracks,
                   int eventID, int parentTrackID) {
  auto it = tracks.find(Key(eventID,parentTrackID));
  return it != tracks.end() && it->second.parent == 0 && std::abs(it->second.pdg) == 13;
}

bool IsElectronEmitter(const std::unordered_map<uint64_t,TrackInfo>& tracks,
                       int eventID, int parentTrackID) {
  auto it = tracks.find(Key(eventID,parentTrackID));
  return it != tracks.end() && std::abs(it->second.pdg) == 11;
}

double SafeMean(double sum, long long n) {
  return n > 0 ? sum/double(n) : std::numeric_limits<double>::quiet_NaN();
}

} // namespace

void scan_corewall_events(
    const char* inputName="validation_corewall_muon_topology_run0.root",
    int topN=10) {

  TFile* f = TFile::Open(inputName,"READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << inputName << "\n";
    return;
  }

  auto* tTrack = dynamic_cast<TTree*>(f->Get("Track"));
  auto* tPho   = dynamic_cast<TTree*>(f->Get("Photon"));
  auto* tCfg   = dynamic_cast<TTree*>(f->Get("RunConfig"));
  if (!tTrack || !tPho || !tCfg) {
    std::cerr << "Missing Track, Photon, or RunConfig\n";
    f->Close();
    return;
  }

  int boundaryDebug=-1, cerenkovEnabled=-1;
  char revision[256] = {};
  char schema[128] = {};
  tCfg->SetBranchAddress("BoundaryDebug",&boundaryDebug);
  tCfg->SetBranchAddress("CerenkovEnabled",&cerenkovEnabled);
  if (tCfg->GetBranch("ProjectRevision")) tCfg->SetBranchAddress("ProjectRevision",revision);
  if (tCfg->GetBranch("SchemaVersion")) tCfg->SetBranchAddress("SchemaVersion",schema);
  tCfg->GetEntry(0);

  std::unordered_map<uint64_t,TrackInfo> tracks;
  int te=0,tid=0,parent=0,pdg=0;
  tTrack->SetBranchAddress("EventID",&te);
  tTrack->SetBranchAddress("TrackID",&tid);
  tTrack->SetBranchAddress("ParentID",&parent);
  tTrack->SetBranchAddress("PDG",&pdg);
  for (Long64_t i=0;i<tTrack->GetEntries();++i) {
    tTrack->GetEntry(i);
    tracks[Key(te,tid)] = {parent,pdg};
  }

  std::unordered_map<int,EventSummary> events;
  int pe=0,parentPho=0,reached=0;
  float wallT=0,birthT=0;
  tPho->SetBranchAddress("EventID",&pe);
  tPho->SetBranchAddress("ParentTrackID",&parentPho);
  tPho->SetBranchAddress("ReachedCoreWall",&reached);
  tPho->SetBranchAddress("FirstWallTime_ns",&wallT);
  tPho->SetBranchAddress("BirthTime_ns",&birthT);

  long long negativeFlight=0;
  for (Long64_t i=0;i<tPho->GetEntries();++i) {
    tPho->GetEntry(i);
    auto& e = events[pe];
    ++e.generated;
    if (!reached) continue;

    ++e.wall;
    e.wallTimeSum += wallT;
    e.wallTimeMin = std::min(e.wallTimeMin,double(wallT));
    e.wallTimeMax = std::max(e.wallTimeMax,double(wallT));

    const double flight = double(wallT)-double(birthT);
    if (flight < -1e-6) ++negativeFlight;
    else {
      e.flightSum += flight;
      e.flightMin = std::min(e.flightMin,flight);
      e.flightMax = std::max(e.flightMax,flight);
      ++e.flightN;
    }

    if (IsPrimaryMuon(tracks,pe,parentPho)) ++e.wallPrimaryMuon;
    else if (IsElectronEmitter(tracks,pe,parentPho)) ++e.wallElectronEmitter;
    else ++e.wallOther;
  }

  std::vector<std::pair<int,EventSummary>> rows(events.begin(),events.end());
  std::sort(rows.begin(),rows.end(),
            [](const auto& a,const auto& b){ return a.first < b.first; });

  std::cout << "\n=== CORE-WALL STEP-3A EVENT SCAN ===\n";
  std::cout << "Input: " << inputName << "\n";
  std::cout << "ProjectRevision: " << (revision[0]?revision:"?") << "\n";
  std::cout << "SchemaVersion: " << (schema[0]?schema:"?") << "\n";
  std::cout << "CerenkovEnabled=" << cerenkovEnabled
            << " BoundaryDebug=" << boundaryDebug << "\n";
  std::cout << "Events with Photon rows=" << rows.size()
            << " negative wall-birth times=" << negativeFlight << "\n\n";

  std::cout << std::setw(6)  << "Evt"
            << std::setw(11) << "Ngen"
            << std::setw(11) << "Nwall"
            << std::setw(11) << "muPrim"
            << std::setw(11) << "eEmit"
            << std::setw(11) << "other"
            << std::setw(12) << "<twall>"
            << std::setw(12) << "<flight>"
            << "\n";

  for (const auto& kv : rows) {
    const int evt = kv.first;
    const auto& e = kv.second;
    std::cout << std::setw(6)  << evt
              << std::setw(11) << e.generated
              << std::setw(11) << e.wall
              << std::setw(11) << e.wallPrimaryMuon
              << std::setw(11) << e.wallElectronEmitter
              << std::setw(11) << e.wallOther
              << std::setw(12) << std::fixed << std::setprecision(4)
              << SafeMean(e.wallTimeSum,e.wall)
              << std::setw(12) << SafeMean(e.flightSum,e.flightN)
              << "\n";
  }

  auto byWall = rows;
  std::sort(byWall.begin(),byWall.end(),
            [](const auto& a,const auto& b){ return a.second.wall > b.second.wall; });

  std::cout << "\nTop " << std::min<int>(topN,byWall.size())
            << " events by Core-wall photon count:\n";
  for (int i=0;i<std::min<int>(topN,byWall.size());++i) {
    const auto& [evt,e] = byWall[i];
    const double muFrac = e.wall ? 100.0*double(e.wallPrimaryMuon)/double(e.wall) : 0.0;
    const double eFrac  = e.wall ? 100.0*double(e.wallElectronEmitter)/double(e.wall) : 0.0;
    std::cout << "  event " << evt
              << ": Nwall=" << e.wall
              << ", primary-muon=" << std::setprecision(2) << muFrac << "%"
              << ", e-emitter=" << eFrac << "%"
              << ", wallTime=[" << std::setprecision(4)
              << e.wallTimeMin << "," << e.wallTimeMax << "] ns"
              << ", flight=[" << e.flightMin << "," << e.flightMax << "] ns"
              << "\n";
  }

  if (!rows.empty()) {
    // Recommend a representative event nearest the median Nwall.
    auto byCount = rows;
    std::sort(byCount.begin(),byCount.end(),
              [](const auto& a,const auto& b){ return a.second.wall < b.second.wall; });
    const auto& rep = byCount[byCount.size()/2];
    std::cout << "\nRecommended representative EventID for Step 3B: "
              << rep.first << " (Nwall=" << rep.second.wall << ")\n";
  }

  const bool pass = !rows.empty() && cerenkovEnabled==1 && boundaryDebug==0 &&
                    negativeFlight==0;
  std::cout << "CORE-WALL STEP-3A EVENT SCAN: "
            << (pass ? "READY_FOR_STEP3B" : "CHECK_QA") << "\n";

  f->Close();
}
