#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TString.h>
#include <TLegend.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

namespace {
struct Point {
  double energyMeV;
  long long nPrimary = 0;
  long long nDecay = 0;
  long long nCoreDecay = 0;
  long long nCoreAnyEnd = 0;
  double meanEndY = 0.0;
};

int findCoreVolumeID(TFile& f) {
  auto* vm = dynamic_cast<TTree*>(f.Get("VolumeMap"));
  if (!vm) return -999999;

  int volumeID = 0;
  char volumeName[256] = {};
  vm->SetBranchAddress("VolumeID", &volumeID);
  vm->SetBranchAddress("VolumeName", volumeName);

  for (Long64_t i = 0; i < vm->GetEntries(); ++i) {
    vm->GetEntry(i);
    std::string name(volumeName);
    if (name == "Core" || name.find("Core") != std::string::npos)
      return volumeID;
  }
  return -999999;
}

Point inspectOne(const char* filename, double energyMeV) {
  Point p;
  p.energyMeV = energyMeV;

  TFile f(filename, "READ");
  if (f.IsZombie()) {
    std::cerr << "Cannot open " << filename << "\n";
    return p;
  }

  auto* tr = dynamic_cast<TTree*>(f.Get("Track"));
  if (!tr) {
    std::cerr << "Missing Track in " << filename << "\n";
    return p;
  }

  const int coreID = findCoreVolumeID(f);
  if (coreID == -999999) {
    std::cerr << "Could not resolve Core VolumeID in " << filename << "\n";
    return p;
  }

  int pdg = 0, parent = 0, endVolumeID = 0;
  float vertexKE = 0.f, endY = 0.f, properNs = 0.f;
  char endProc[256] = {};

  tr->SetBranchAddress("PDG", &pdg);
  tr->SetBranchAddress("ParentID", &parent);
  tr->SetBranchAddress("VertexKE_MeV", &vertexKE);
  tr->SetBranchAddress("EndVolumeID", &endVolumeID);
  tr->SetBranchAddress("EndY_mm", &endY);
  tr->SetBranchAddress("EndProperTime_ns", &properNs);
  tr->SetBranchAddress("EndProcessName", endProc);

  double sumEndY = 0.0;
  long long nEndY = 0;

  for (Long64_t i = 0; i < tr->GetEntries(); ++i) {
    tr->GetEntry(i);

    // PDG convention: mu+ = -13. Only gun primaries.
    if (parent != 0 || pdg != -13) continue;
    ++p.nPrimary;

    const bool endedInCore = (endVolumeID == coreID);
    const bool decayed = std::string(endProc).find("Decay") != std::string::npos;

    if (endedInCore) {
      ++p.nCoreAnyEnd;
      if (std::isfinite(endY)) {
        sumEndY += endY;
        ++nEndY;
      }
    }
    if (decayed) ++p.nDecay;
    if (endedInCore && decayed) ++p.nCoreDecay;
  }

  if (nEndY > 0) p.meanEndY = sumEndY / nEndY;
  return p;
}
}

void validate_detector_muon_stopping_scan(
    const char* prefix="validation_detector_muon_stopping_scan") {

  const std::vector<double> energies = {10,20,30,40,50,60,80,100,150,200};
  std::vector<Point> pts;
  pts.reserve(energies.size());

  std::cout << "\nDETECTOR-GEOMETRY STOPPED mu+ ENERGY SCAN\n";
  std::cout << "Core-stop proxy = primary mu+ ends by Decay with EndVolumeID == Core\n\n";
  std::cout << std::setw(8) << "KE[MeV]"
            << std::setw(11) << "Nprim"
            << std::setw(11) << "Ndecay"
            << std::setw(12) << "CoreEnd"
            << std::setw(13) << "CoreDecay"
            << std::setw(12) << "Frac[%]"
            << std::setw(14) << "<EndY>[mm]" << "\n";

  int bestIndex = -1;
  double bestFrac = -1.0;

  for (size_t i = 0; i < energies.size(); ++i) {
    TString fn = Form("%s_run%zu.root", prefix, i);
    Point p = inspectOne(fn.Data(), energies[i]);
    pts.push_back(p);

    const double frac = p.nPrimary ? 100.0 * double(p.nCoreDecay) / p.nPrimary : 0.0;
    if (frac > bestFrac) {
      bestFrac = frac;
      bestIndex = int(i);
    }

    std::cout << std::fixed << std::setprecision(1)
              << std::setw(8) << p.energyMeV
              << std::setw(11) << p.nPrimary
              << std::setw(11) << p.nDecay
              << std::setw(12) << p.nCoreAnyEnd
              << std::setw(13) << p.nCoreDecay
              << std::setw(12) << std::setprecision(2) << frac
              << std::setw(14) << std::setprecision(2) << p.meanEndY
              << "\n";
  }

  gSystem->mkdir("validation_plots", kTRUE);
  TH1D h("hCoreStopFraction",
         "Stopped #mu^{+} in Water Core vs input kinetic energy;"
         "Input kinetic energy (MeV);Core-decay fraction",
         energies.size(), 0.5, energies.size()+0.5);

  for (size_t i = 0; i < pts.size(); ++i) {
    const double frac = pts[i].nPrimary ?
      double(pts[i].nCoreDecay) / pts[i].nPrimary : 0.0;
    h.SetBinContent(i+1, frac);
    h.GetXaxis()->SetBinLabel(i+1, Form("%.0f", energies[i]));
  }
  h.SetMinimum(0.0);
  h.SetMaximum(1.05);
  h.SetStats(false);

  TCanvas c("cDetectorMuonStopping",
            "Detector muon stopping scan", 1000, 650);
  h.Draw("HIST TEXT");
  c.SaveAs("validation_plots/detector_muon_stopping_scan.png");
  c.SaveAs("validation_plots/detector_muon_stopping_scan.pdf");

  if (bestIndex >= 0) {
    std::cout << "\nBest tested point: "
              << energies[bestIndex] << " MeV, Core-decay fraction = "
              << std::setprecision(2) << bestFrac << "%\n";
    std::cout << "Use this only as the center of the NEXT finer scan; "
                 "do not freeze the lifetime energy yet.\n";
  }

  std::cout << "\nInterpretation:\n"
            << "  - This is NOT yet the lifetime measurement.\n"
            << "  - It finds a mu+ energy window that stops/decays inside the real Water Core.\n"
            << "  - After choosing a stable energy window, run a high-stat detector-level\n"
            << "    lifetime sample and fit EndProperTime_ns plus delayed daughter timing.\n";
}
