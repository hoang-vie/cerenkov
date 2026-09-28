#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <cmath>

void inspect_output(const char* inputName, long long expectedGenerated=-1, int expectedChannels=-1) {
  TFile f(inputName, "READ");
  if (f.IsZombie()) { std::cerr << "Cannot open " << inputName << "\n"; return; }

  const std::vector<const char*> names = {
    "TriggerEvent","ScintDeposit","Track","TrackVolume","Photon",
    "SensorIncident","DetectedSignal","ChannelEvent","RunConfig","ChannelMap",
    "VolumeMap","RunStatistics","BoundaryDebug","ReadoutPlaneConfig","PrimaryTruth",
    "CoreBoundaryCrossing"
  };

  std::cout << "\n=== ROOT OUTPUT STRUCTURE ===\n";
  for (auto name : names) {
    auto* t = dynamic_cast<TTree*>(f.Get(name));
    std::cout << "  " << name << ": " << (t ? std::to_string(t->GetEntries()) : "MISSING") << " rows\n";
  }

  auto* cfg = dynamic_cast<TTree*>(f.Get("RunConfig"));
  auto* stat = dynamic_cast<TTree*>(f.Get("RunStatistics"));
  auto* trig = dynamic_cast<TTree*>(f.Get("TriggerEvent"));
  auto* prim = dynamic_cast<TTree*>(f.Get("PrimaryTruth"));
  auto* cmap = dynamic_cast<TTree*>(f.Get("ChannelMap"));
  auto* inc = dynamic_cast<TTree*>(f.Get("SensorIncident"));
  auto* det = dynamic_cast<TTree*>(f.Get("DetectedSignal"));
  auto* coreBoundary = dynamic_cast<TTree*>(f.Get("CoreBoundaryCrossing"));

  long long generated=-1, triggered=-1, detailed=-1;
  double cryTime=0.0;
  if (stat && stat->GetEntries()>0) {
    int generatedRow=0, triggeredRow=0, detailedRow=0, workerID=0;
    double cryTimeRow=0.0;
    if (stat->GetBranch("WorkerThreadID")) stat->SetBranchAddress("WorkerThreadID", &workerID);
    stat->SetBranchAddress("GeneratedEvents", &generatedRow);
    stat->SetBranchAddress("TriggeredEvents", &triggeredRow);
    stat->SetBranchAddress("DetailedEvents", &detailedRow);
    if (stat->GetBranch("CRYSimulatedTime_s")) stat->SetBranchAddress("CRYSimulatedTime_s", &cryTimeRow);
    generated=triggered=detailed=0;
    for (Long64_t i=0;i<stat->GetEntries();++i) {
      stat->GetEntry(i);
      generated += generatedRow;
      triggered += triggeredRow;
      detailed += detailedRow;
      cryTime += cryTimeRow;
    }
    std::cout << "\nRunStatistics total from " << stat->GetEntries() << " worker row(s): generated=" << generated
              << " triggered=" << triggered << " detailed=" << detailed
              << " CRY_time_s=" << cryTime << "\n";
  }

  char schema[256] = {};
  char mode[256] = {};
  char sensorMode[256] = {};
  char coreBoundaryMode[256] = {};
  char coreBoundaryDetail[256] = {};
  char projectRevision[256] = {};
  char gitBranch[256] = {};
  int gitDirty=-1;
  int configuredChannels=-1, configuredThreads=-1;
  if (cfg && cfg->GetEntries()>0) {
    if (cfg->GetBranch("SchemaVersion")) cfg->SetBranchAddress("SchemaVersion", schema);
    if (cfg->GetBranch("OutputEventMode")) cfg->SetBranchAddress("OutputEventMode", mode);
    if (cfg->GetBranch("SensorMode")) cfg->SetBranchAddress("SensorMode", sensorMode);
    if (cfg->GetBranch("CoreBoundaryMode")) cfg->SetBranchAddress("CoreBoundaryMode", coreBoundaryMode);
    if (cfg->GetBranch("CoreBoundaryDetail")) cfg->SetBranchAddress("CoreBoundaryDetail", coreBoundaryDetail);
    if (cfg->GetBranch("ProjectRevision")) cfg->SetBranchAddress("ProjectRevision", projectRevision);
    if (cfg->GetBranch("GitBranch")) cfg->SetBranchAddress("GitBranch", gitBranch);
    if (cfg->GetBranch("GitDirty")) cfg->SetBranchAddress("GitDirty", &gitDirty);
    if (cfg->GetBranch("ChannelCount")) cfg->SetBranchAddress("ChannelCount", &configuredChannels);
    if (cfg->GetBranch("NThreads")) cfg->SetBranchAddress("NThreads", &configuredThreads);
    cfg->GetEntry(0);
    std::cout << "RunConfig: schema=" << (schema[0]?schema:"?")
              << " outputMode=" << (mode[0]?mode:"?")
              << " sensorMode=" << (sensorMode[0]?sensorMode:"?")
              << " coreBoundary=" << (coreBoundaryMode[0]?coreBoundaryMode:"?")
              << "/" << (coreBoundaryDetail[0]?coreBoundaryDetail:"?")
              << " channels=" << configuredChannels
              << " threads=" << configuredThreads;
    if (projectRevision[0]) {
      std::cout << " revision=" << projectRevision;
      if (gitBranch[0]) std::cout << " branch=" << gitBranch;
      if (gitDirty >= 0) std::cout << " dirty=" << gitDirty;
    }
    std::cout << "\n";
  }

  bool ok=true;
  for (auto name : names) {
    if (!f.Get(name)) {
      std::cout << "[FAIL] Missing required ntuple " << name << "\n";
      ok=false;
    }
  }
  if (cfg && cfg->GetEntries()!=1) {
    std::cout << "[FAIL] RunConfig rows=" << cfg->GetEntries() << " expected exactly 1\n"; ok=false;
  } else if (cfg) {
    std::cout << "[ OK ] RunConfig has exactly one metadata row\n";
  }
  if (expectedGenerated >= 0 && generated != expectedGenerated) {
    std::cout << "[FAIL] GeneratedEvents=" << generated << " expected=" << expectedGenerated << "\n"; ok=false;
  } else if (expectedGenerated >= 0) {
    std::cout << "[ OK ] GeneratedEvents matches expected value " << expectedGenerated << "\n";
  }
  if (expectedChannels >= 0 && configuredChannels != expectedChannels) {
    std::cout << "[FAIL] RunConfig.ChannelCount=" << configuredChannels << " expected=" << expectedChannels << "\n"; ok=false;
  } else if (expectedChannels >= 0) {
    std::cout << "[ OK ] ChannelCount matches expected value " << expectedChannels << "\n";
  }
  if (stat && configuredThreads > 0 && stat->GetEntries() != configuredThreads) {
    std::cout << "[FAIL] RunStatistics worker rows=" << stat->GetEntries()
              << " but RunConfig.NThreads=" << configuredThreads << "\n"; ok=false;
  } else if (stat && configuredThreads > 0) {
    std::cout << "[ OK ] RunStatistics has one row per configured worker thread\n";
  }
  if (cmap && configuredChannels>=0 && cmap->GetEntries()!=configuredChannels) {
    std::cout << "[FAIL] ChannelMap rows != RunConfig.ChannelCount\n"; ok=false;
  } else if (cmap && configuredChannels>=0) {
    std::cout << "[ OK ] ChannelMap row count matches ChannelCount\n";
  }

  const std::string modeStr = mode;
  const std::string sensorModeStr = sensorMode;
  if (!modeStr.empty() && generated>=0 && (modeStr=="thinAll" || modeStr=="fullAll") && trig) {
    if (trig->GetEntries()!=generated) {
      std::cout << "[FAIL] TriggerEvent rows != GeneratedEvents in " << modeStr << "\n"; ok=false;
    } else std::cout << "[ OK ] TriggerEvent contains all generated events\n";
  }
  if (!modeStr.empty() && triggered>=0 && modeStr=="triggered" && trig) {
    if (trig->GetEntries()!=triggered) {
      std::cout << "[FAIL] TriggerEvent rows != TriggeredEvents in triggered mode\n"; ok=false;
    } else std::cout << "[ OK ] TriggerEvent contains only triggered events\n";
  }

  if (!modeStr.empty() && generated>=0 && (modeStr=="thinAll" || modeStr=="fullAll") && trig && prim) {
    int np=0; Long64_t sumNP=0;
    trig->SetBranchAddress("NPrimaryNonOptical", &np);
    for (Long64_t i=0;i<trig->GetEntries();++i) { trig->GetEntry(i); sumNP += np; }
    if (prim->GetEntries()!=sumNP) {
      std::cout << "[FAIL] PrimaryTruth rows=" << prim->GetEntries()
                << " but sum TriggerEvent.NPrimaryNonOptical=" << sumNP << "\n"; ok=false;
    } else std::cout << "[ OK ] PrimaryTruth count matches non-optical event primary multiplicity\n";
  }

  if (inc && det) {
    if (det->GetEntries() > inc->GetEntries()) {
      std::cout << "[FAIL] DetectedSignal rows exceed SensorIncident rows\n"; ok=false;
    } else std::cout << "[ OK ] NDetected <= NIncident\n";
  }

  if (sensorModeStr=="ConcentratorInput" && det && det->GetEntries()!=0) {
    std::cout << "[FAIL] ConcentratorInput should not create DetectedSignal rows\n"; ok=false;
  } else if (sensorModeStr=="ConcentratorInput" && det) {
    std::cout << "[ OK ] ConcentratorInput has zero DetectedSignal rows\n";
  }

  const std::string coreBoundaryModeStr = coreBoundaryMode;
  if (coreBoundaryModeStr == "off" && coreBoundary && coreBoundary->GetEntries()!=0) {
    std::cout << "[FAIL] CoreBoundaryMode=off but CoreBoundaryCrossing has "
              << coreBoundary->GetEntries() << " row(s)\n"; ok=false;
  } else if (coreBoundaryModeStr == "off" && coreBoundary) {
    std::cout << "[ OK ] CoreBoundaryMode=off has zero CoreBoundaryCrossing rows\n";
  }

  std::cout << (ok ? "\nOUTPUT SANITY: PASS\n" : "\nOUTPUT SANITY: FAIL -- inspect messages above\n");
}
