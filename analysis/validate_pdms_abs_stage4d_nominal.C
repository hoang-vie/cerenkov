#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

void validate_pdms_abs_stage4d_nominal() {

  constexpr Long64_t N = 1000000;

  // PDMS occupies x = 53..54 mm.
  // Source is at x = 53.0001 mm.
  constexpr double path_mm = 0.9999;

  // Generic10m table, converted by absLengthUnit = m.
  constexpr double L_mm = 10000.0;

  TFile f(
    "validation_pdms_abs_stage4d_nominal_run0.root",
    "READ"
  );

  auto* photon =
      dynamic_cast<TTree*>(f.Get("Photon"));

  auto* incident =
      dynamic_cast<TTree*>(f.Get("SensorIncident"));

  auto* detected =
      dynamic_cast<TTree*>(f.Get("DetectedSignal"));

  auto* rc =
      dynamic_cast<TTree*>(f.Get("RunConfig"));

  if (!photon || !incident || !detected || !rc) {
    std::cout
        << "Missing required tree : FAIL\n";
    return;
  }

  const Long64_t nPhoton =
      photon->GetEntries();

  const Long64_t nIncident =
      incident->GetEntries();

  const Long64_t nDetected =
      detected->GetEntries();

  const Long64_t safety =
      photon->GetEntries(
        "KilledBySafetyCut!=0"
      );

  const Long64_t absorbed =
      N - nIncident;

  const double expectedSurv =
      std::exp(-path_mm/L_mm);

  const double expectedAbsProb =
      1.0 - expectedSurv;

  const double expectedAbs =
      N*expectedAbsProb;

  const double sigma =
      std::sqrt(
        N*
        expectedAbsProb*
        (1.0-expectedAbsProb)
      );

  const double z =
      sigma>0.0
      ? (double(absorbed)-expectedAbs)/sigma
      : 0.0;

  char schema[128]={0};
  char profile[128]={0};
  char dataset[512]={0};
  char status[512]={0};

  Float_t validMin=-999;
  Float_t validMax=-999;
  Float_t measuredMin=-999;
  Float_t measuredMax=-999;

  rc->SetBranchAddress(
      "SchemaVersion",
      schema
  );

  rc->SetBranchAddress(
      "CouplingProfile",
      profile
  );

  rc->SetBranchAddress(
      "CouplingProfileAbsDataset",
      dataset
  );

  rc->SetBranchAddress(
      "CouplingProfileAbsValidMin_nm",
      &validMin
  );

  rc->SetBranchAddress(
      "CouplingProfileAbsValidMax_nm",
      &validMax
  );

  rc->SetBranchAddress(
      "CouplingProfileAbsMeasuredMin_nm",
      &measuredMin
  );

  rc->SetBranchAddress(
      "CouplingProfileAbsMeasuredMax_nm",
      &measuredMax
  );

  rc->SetBranchAddress(
      "CouplingProfileAbsStatus",
      status
  );

  rc->GetEntry(0);

  const bool metadataOK =
      rc->GetListOfBranches()->GetEntries()==154 &&
      std::string(schema)=="general-v1.8" &&
      std::string(profile)=="PDMS_Generic" &&
      std::string(dataset)==
        "OpticalGrease_PDMS_ABS_Generic10m.txt" &&
      std::abs(validMin-200.0f)<1e-4f &&
      std::abs(validMax-900.0f)<1e-4f &&
      std::abs(measuredMin+1.0f)<1e-4f &&
      std::abs(measuredMax+1.0f)<1e-4f &&
      std::string(status)==
        "Generic10m_Placeholder_NotVendorMeasurement";

  const bool pass =
      nPhoton==N &&
      nIncident<=N &&
      nDetected==0 &&
      safety==0 &&
      metadataOK &&
      std::abs(z)<5.0;

  std::cout
      << "\n=== PDMS STAGE 4D-B "
      << "GENERIC10M RUNTIME VALIDATION ===\n";

  std::cout
      << std::fixed
      << std::setprecision(6)

      << "path=" << path_mm << " mm"
      << " L=" << L_mm << " mm\n"

      << "expectedSurvival="
      << 100.0*expectedSurv << "%\n"

      << "observedSurvival="
      << 100.0*double(nIncident)/double(N)
      << "%\n"

      << "expectedAbsorbed="
      << expectedAbs << "\n"

      << "observedAbsorbed="
      << absorbed << "\n"

      << "z=" << z << "\n"

      << "Photon=" << nPhoton
      << " SensorIncident=" << nIncident
      << " DetectedSignal=" << nDetected
      << " safety=" << safety
      << "\n"

      << "metadata="
      << (metadataOK ? "OK" : "BAD")
      << "\n";

  std::cout
      << "PDMS STAGE 4D-B GENERIC10M "
      << "RUNTIME VALIDATION: "
      << (pass ? "PASS" : "FAIL")
      << "\n";
}
