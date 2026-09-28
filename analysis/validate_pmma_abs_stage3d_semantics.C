#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void validate_pmma_abs_stage3d_semantics(
    const char* prefix="validation_pmma_abs_stage3d_model_semantics") {

  const std::vector<double> wavelengths = {
      250.0, 299.0, 300.0, 500.0, 700.0
  };

  const std::vector<Long64_t> expectedOOR = {
      1, 1, 0, 0, 1
  };

  bool allPass = true;

  std::cout
      << "\n=== PMMA STAGE 3D MEASURED/MODEL SEMANTICS VALIDATION ===\n";

  for (int run=0; run<5; ++run) {

    std::ostringstream fn;
    fn << prefix << "_run" << run << ".root";

    TFile f(fn.str().c_str(), "READ");

    if (f.IsZombie()) {
      std::cout
          << "run" << run
          << " missing ROOT file : FAIL\n";
      allPass = false;
      continue;
    }

    auto* photon = dynamic_cast<TTree*>(f.Get("Photon"));
    auto* rc = dynamic_cast<TTree*>(f.Get("RunConfig"));

    if (!photon || !rc || rc->GetEntries()!=1) {
      std::cout
          << "run" << run
          << " missing Photon/RunConfig : FAIL\n";
      allPass = false;
      continue;
    }

    const Long64_t nPhoton = photon->GetEntries();
    const Long64_t oor =
        photon->GetEntries("AcrylicAbsOutOfRange!=0");

    char schema[128] = {};
    char revision[256] = {};
    char registry[256] = {};
    char coreProfile[128] = {};
    char dataset[512] = {};
    char status[512] = {};
    char policy[1024] = {};

    Float_t supportMin = -1;
    Float_t supportMax = -1;
    Float_t measuredMin = -1;
    Float_t measuredMax = -1;

    rc->SetBranchAddress("SchemaVersion", schema);
    rc->SetBranchAddress("ProjectRevision", revision);
    rc->SetBranchAddress("MaterialProfileRegistryVersion", registry);

    rc->SetBranchAddress("CoreProfile", coreProfile);
    rc->SetBranchAddress("CoreProfileAbsDataset", dataset);

    rc->SetBranchAddress(
        "CoreProfileAbsValidMin_nm", &supportMin);
    rc->SetBranchAddress(
        "CoreProfileAbsValidMax_nm", &supportMax);

    rc->SetBranchAddress(
        "CoreProfileAbsMeasuredMin_nm", &measuredMin);
    rc->SetBranchAddress(
        "CoreProfileAbsMeasuredMax_nm", &measuredMax);

    rc->SetBranchAddress("CoreProfileAbsStatus", status);
    rc->SetBranchAddress("CoreProfileOutOfRangePolicy", policy);

    rc->GetEntry(0);

    const std::string policyStr(policy);

    const bool schemaOK =
        std::string(schema) == "general-v1.8" &&
        rc->GetListOfBranches()->GetEntries() == 154;

    const bool provenanceOK =
        std::string(revision) ==
          "hotfix6_physics_audit5_pmma_stage3d_measured_support_split" &&
        std::string(registry) ==
          "curated-bulk-v2_absMeasuredRange" &&
        std::string(coreProfile) == "PMMA_SNO_A" &&
        std::string(dataset) ==
          "PMMA_SNO_UVT_SampleA_ABS_Measured300to500_Model200to900_Candidate.txt" &&
        std::string(status) ==
          "RuntimeValidationOverride_via_pmmaAbsFile";

    const bool rangeOK =
        std::abs(supportMin - 200.0f) < 1e-4f &&
        std::abs(supportMax - 900.0f) < 1e-4f &&
        std::abs(measuredMin - 300.0f) < 1e-4f &&
        std::abs(measuredMax - 500.0f) < 1e-4f;

    const bool policyOK =
        policyStr.find("ExplicitAbsTableOverride") != std::string::npos &&
        policyStr.find("MeasuredValidityFromProfile") != std::string::npos &&
        policyStr.find("NumericalSupportFromDataset") != std::string::npos &&
        policyStr.find("OutsideSupport_Geant4EndpointClamp") != std::string::npos;

    const bool oorOK =
        nPhoton == 1 &&
        oor == expectedOOR[run];

    const bool pass =
        schemaOK &&
        provenanceOK &&
        rangeOK &&
        policyOK &&
        oorOK;

    allPass &= pass;

    std::cout
        << "run" << run
        << " " << wavelengths[run] << " nm"
        << " support=[" << supportMin << "," << supportMax << "]"
        << " measured=[" << measuredMin << "," << measuredMax << "]"
        << " OOR=" << oor
        << " schema=" << (schemaOK ? "OK" : "BAD")
        << " provenance=" << (provenanceOK ? "OK" : "BAD")
        << " policy=" << (policyOK ? "OK" : "BAD")
        << " : " << (pass ? "PASS" : "FAIL")
        << "\n";
  }

  std::cout
      << "PMMA STAGE 3D MEASURED/MODEL SEMANTICS VALIDATION: "
      << (allPass ? "PASS" : "FAIL")
      << "\n";
}
