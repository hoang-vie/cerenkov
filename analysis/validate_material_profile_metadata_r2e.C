#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>

namespace {
bool Near(float a, float b, double tol=1e-5) {
  return std::abs(double(a)-double(b)) <= tol;
}
}

void validate_material_profile_metadata_r2e(
    const char* inputName="validation_material_profile_metadata_r2e_run0.root") {
  TFile* f=TFile::Open(inputName,"READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << inputName << "\n";
    return;
  }
  auto* t=dynamic_cast<TTree*>(f->Get("RunConfig"));
  if (!t || t->GetEntries()!=1) {
    std::cerr << "RunConfig missing or row count != 1\n";
    f->Close();
    return;
  }

  const char* required[] = {
    "SchemaVersion","ProjectRevision","MaterialProfileRegistryVersion",
    "LegacyRunConfigCompatibility","CoreProfile","CoreProfileG4Material",
    "CoreProfileRindexDataset","CoreProfileAbsDataset","CoreProfileRayleighDataset",
    "CoreProfileRindexValidMin_nm","CoreProfileRindexValidMax_nm",
    "CoreProfileAbsValidMin_nm","CoreProfileAbsValidMax_nm",
    "CoreProfileRayleighValidMin_nm","CoreProfileRayleighValidMax_nm",
    "CoreProfileRindexStatus","CoreProfileAbsStatus","CoreProfileRayleighStatus",
    "CoreProfileOutOfRangePolicy","UseOuterShell","OuterShellThickness_mm",
    "OuterShellProfile","OuterShellProfileG4Material","OuterShellProfileRindexDataset",
    "OuterShellProfileAbsDataset","OuterShellProfileRayleighDataset",
    "OuterShellProfileRindexValidMin_nm","OuterShellProfileRindexValidMax_nm",
    "OuterShellProfileAbsValidMin_nm","OuterShellProfileAbsValidMax_nm",
    "OuterShellProfileRayleighValidMin_nm","OuterShellProfileRayleighValidMax_nm",
    "OuterShellProfileRindexStatus","OuterShellProfileAbsStatus",
    "OuterShellProfileRayleighStatus","OuterShellProfileOutOfRangePolicy",
    "CouplingProfile","CouplingProfileG4Material","CouplingProfileRindexDataset",
    "CouplingProfileAbsDataset","CouplingProfileRayleighDataset",
    "CouplingProfileRindexValidMin_nm","CouplingProfileRindexValidMax_nm",
    "CouplingProfileAbsValidMin_nm","CouplingProfileAbsValidMax_nm",
    "CouplingProfileRayleighValidMin_nm","CouplingProfileRayleighValidMax_nm",
    "CouplingProfileRindexStatus","CouplingProfileAbsStatus",
    "CouplingProfileRayleighStatus","CouplingProfileOutOfRangePolicy",
    // legacy compatibility comparisons
    "AcrylicThickness_mm","UseAcrylicShell","AcrylicRindexDataset","AcrylicAbsDataset",
    "WaterRindexDataset","WaterAbsDataset","WaterRayleighDataset",
    "CouplingRindexDataset","CouplingAbsDataset",
    "AcrylicAbsValidMin_nm","AcrylicAbsValidMax_nm",
    "AcrylicRindexValidMin_nm","AcrylicRindexValidMax_nm",
    "CouplingRindexValidMin_nm","CouplingRindexValidMax_nm",
    "CouplingAbsValidMin_nm","CouplingAbsValidMax_nm",
    "WaterRindexValidMin_nm","WaterRindexValidMax_nm",
    "WaterAbsValidMin_nm","WaterAbsValidMax_nm",
    "WaterRayleighValidMin_nm","WaterRayleighValidMax_nm"
  };
  bool pass=true;
  for (const char* b : required) {
    if (!t->GetBranch(b)) {
      std::cerr << "Missing branch: " << b << "\n";
      pass=false;
    }
  }
  if (!pass) {
    std::cout << "R2E MATERIAL PROFILE METADATA VALIDATION: FAIL\n";
    f->Close();
    return;
  }

  // Geant4 ROOT analysis stores CreateNtupleSColumn fields as C-style
  // Char_t branches in this output.  Bind writable char buffers rather than
  // std::string*: ROOT otherwise reports a Char_t/std::string type mismatch.
  constexpr std::size_t kText = 1024;
  char schema[kText]={},revision[kText]={},registry[kText]={},compat[kText]={};
  char coreProfile[kText]={},coreMat[kText]={},coreRI[kText]={},coreAbs[kText]={},coreRay[kText]={};
  char coreRIStatus[kText]={},coreAbsStatus[kText]={},coreRayStatus[kText]={},corePolicy[kText]={};
  char shellProfile[kText]={},shellMat[kText]={},shellRI[kText]={},shellAbs[kText]={},shellRay[kText]={};
  char shellRIStatus[kText]={},shellAbsStatus[kText]={},shellRayStatus[kText]={},shellPolicy[kText]={};
  char coupProfile[kText]={},coupMat[kText]={},coupRI[kText]={},coupAbs[kText]={},coupRay[kText]={};
  char coupRIStatus[kText]={},coupAbsStatus[kText]={},coupRayStatus[kText]={},coupPolicy[kText]={};
  char legacyAcrRI[kText]={},legacyAcrAbs[kText]={},legacyWaterRI[kText]={},legacyWaterAbs[kText]={},legacyWaterRay[kText]={};
  char legacyCoupRI[kText]={},legacyCoupAbs[kText]={};
  int useShell=0,legacyUseShell=0;
  float shellT=0,legacyShellT=0;
  float coreRIMin=0,coreRIMax=0,coreAbsMin=0,coreAbsMax=0,coreRayMin=0,coreRayMax=0;
  float shellRIMin=0,shellRIMax=0,shellAbsMin=0,shellAbsMax=0,shellRayMin=0,shellRayMax=0;
  float coupRIMin=0,coupRIMax=0,coupAbsMin=0,coupAbsMax=0,coupRayMin=0,coupRayMax=0;
  float legacyAcrAbsMin=0,legacyAcrAbsMax=0,legacyAcrRIMin=0,legacyAcrRIMax=0;
  float legacyCoupRIMin=0,legacyCoupRIMax=0,legacyCoupAbsMin=0,legacyCoupAbsMax=0;
  float legacyWaterRIMin=0,legacyWaterRIMax=0,legacyWaterAbsMin=0,legacyWaterAbsMax=0,legacyWaterRayMin=0,legacyWaterRayMax=0;

#define BS(name,var) t->SetBranchAddress(name,var)
#define BN(name,var) t->SetBranchAddress(name,&var)
  BS("SchemaVersion",schema); BS("ProjectRevision",revision);
  BS("MaterialProfileRegistryVersion",registry); BS("LegacyRunConfigCompatibility",compat);
  BS("CoreProfile",coreProfile); BS("CoreProfileG4Material",coreMat);
  BS("CoreProfileRindexDataset",coreRI); BS("CoreProfileAbsDataset",coreAbs); BS("CoreProfileRayleighDataset",coreRay);
  BN("CoreProfileRindexValidMin_nm",coreRIMin); BN("CoreProfileRindexValidMax_nm",coreRIMax);
  BN("CoreProfileAbsValidMin_nm",coreAbsMin); BN("CoreProfileAbsValidMax_nm",coreAbsMax);
  BN("CoreProfileRayleighValidMin_nm",coreRayMin); BN("CoreProfileRayleighValidMax_nm",coreRayMax);
  BS("CoreProfileRindexStatus",coreRIStatus); BS("CoreProfileAbsStatus",coreAbsStatus);
  BS("CoreProfileRayleighStatus",coreRayStatus); BS("CoreProfileOutOfRangePolicy",corePolicy);
  BN("UseOuterShell",useShell); BN("OuterShellThickness_mm",shellT);
  BS("OuterShellProfile",shellProfile); BS("OuterShellProfileG4Material",shellMat);
  BS("OuterShellProfileRindexDataset",shellRI); BS("OuterShellProfileAbsDataset",shellAbs); BS("OuterShellProfileRayleighDataset",shellRay);
  BN("OuterShellProfileRindexValidMin_nm",shellRIMin); BN("OuterShellProfileRindexValidMax_nm",shellRIMax);
  BN("OuterShellProfileAbsValidMin_nm",shellAbsMin); BN("OuterShellProfileAbsValidMax_nm",shellAbsMax);
  BN("OuterShellProfileRayleighValidMin_nm",shellRayMin); BN("OuterShellProfileRayleighValidMax_nm",shellRayMax);
  BS("OuterShellProfileRindexStatus",shellRIStatus); BS("OuterShellProfileAbsStatus",shellAbsStatus);
  BS("OuterShellProfileRayleighStatus",shellRayStatus); BS("OuterShellProfileOutOfRangePolicy",shellPolicy);
  BS("CouplingProfile",coupProfile); BS("CouplingProfileG4Material",coupMat);
  BS("CouplingProfileRindexDataset",coupRI); BS("CouplingProfileAbsDataset",coupAbs); BS("CouplingProfileRayleighDataset",coupRay);
  BN("CouplingProfileRindexValidMin_nm",coupRIMin); BN("CouplingProfileRindexValidMax_nm",coupRIMax);
  BN("CouplingProfileAbsValidMin_nm",coupAbsMin); BN("CouplingProfileAbsValidMax_nm",coupAbsMax);
  BN("CouplingProfileRayleighValidMin_nm",coupRayMin); BN("CouplingProfileRayleighValidMax_nm",coupRayMax);
  BS("CouplingProfileRindexStatus",coupRIStatus); BS("CouplingProfileAbsStatus",coupAbsStatus);
  BS("CouplingProfileRayleighStatus",coupRayStatus); BS("CouplingProfileOutOfRangePolicy",coupPolicy);

  BN("AcrylicThickness_mm",legacyShellT); BN("UseAcrylicShell",legacyUseShell);
  BS("AcrylicRindexDataset",legacyAcrRI); BS("AcrylicAbsDataset",legacyAcrAbs);
  BS("WaterRindexDataset",legacyWaterRI); BS("WaterAbsDataset",legacyWaterAbs); BS("WaterRayleighDataset",legacyWaterRay);
  BS("CouplingRindexDataset",legacyCoupRI); BS("CouplingAbsDataset",legacyCoupAbs);
  BN("AcrylicAbsValidMin_nm",legacyAcrAbsMin); BN("AcrylicAbsValidMax_nm",legacyAcrAbsMax);
  BN("AcrylicRindexValidMin_nm",legacyAcrRIMin); BN("AcrylicRindexValidMax_nm",legacyAcrRIMax);
  BN("CouplingRindexValidMin_nm",legacyCoupRIMin); BN("CouplingRindexValidMax_nm",legacyCoupRIMax);
  BN("CouplingAbsValidMin_nm",legacyCoupAbsMin); BN("CouplingAbsValidMax_nm",legacyCoupAbsMax);
  BN("WaterRindexValidMin_nm",legacyWaterRIMin); BN("WaterRindexValidMax_nm",legacyWaterRIMax);
  BN("WaterAbsValidMin_nm",legacyWaterAbsMin); BN("WaterAbsValidMax_nm",legacyWaterAbsMax);
  BN("WaterRayleighValidMin_nm",legacyWaterRayMin); BN("WaterRayleighValidMax_nm",legacyWaterRayMax);
#undef BS
#undef BN

  t->GetEntry(0);
  auto eqs=[](const char* a,const char* b){ return std::strcmp(a,b)==0; };
  auto eqp=[](const char* a,const char* b){ return std::strcmp(a,b)==0; };
  auto nonempty=[](const char* a){ return a[0] != '\0'; };

  pass &= eqs(schema,"general-v1.7");
  pass &= eqs(revision,"hotfix6_physics_audit5_r2e_profile_metadata");
  pass &= eqs(registry,"curated-bulk-v1");
  pass &= std::strstr(compat,"columns 0-98 preserved") != nullptr;

  pass &= eqs(coreProfile,"Water_Default") && eqs(coreMat,"G4_WATER");
  pass &= eqp(coreRI,legacyWaterRI) && eqp(coreAbs,legacyWaterAbs) && eqp(coreRay,legacyWaterRay);
  pass &= Near(coreRIMin,legacyWaterRIMin) && Near(coreRIMax,legacyWaterRIMax);
  pass &= Near(coreAbsMin,legacyWaterAbsMin) && Near(coreAbsMax,legacyWaterAbsMax);
  pass &= Near(coreRayMin,legacyWaterRayMin) && Near(coreRayMax,legacyWaterRayMax);

  pass &= useShell==1 && legacyUseShell==1 && Near(shellT,3.0) && Near(shellT,legacyShellT);
  pass &= eqs(shellProfile,"PMMA_SNO_A") && eqs(shellMat,"G4_PLEXIGLASS");
  pass &= eqp(shellRI,legacyAcrRI) && eqp(shellAbs,legacyAcrAbs) && shellRay[0]=='\0';
  pass &= Near(shellRIMin,legacyAcrRIMin) && Near(shellRIMax,legacyAcrRIMax);
  pass &= Near(shellAbsMin,legacyAcrAbsMin) && Near(shellAbsMax,legacyAcrAbsMax);
  pass &= Near(shellRayMin,-1.0) && Near(shellRayMax,-1.0);

  pass &= eqs(coupProfile,"PDMS_Generic") && eqs(coupMat,"OpticalGrease_PDMS_Generic");
  pass &= eqp(coupRI,legacyCoupRI) && eqp(coupAbs,legacyCoupAbs) && coupRay[0]=='\0';
  pass &= Near(coupRIMin,legacyCoupRIMin) && Near(coupRIMax,legacyCoupRIMax);
  pass &= Near(coupAbsMin,legacyCoupAbsMin) && Near(coupAbsMax,legacyCoupAbsMax);
  pass &= Near(coupRayMin,-1.0) && Near(coupRayMax,-1.0);

  pass &= nonempty(coreRIStatus) && nonempty(coreAbsStatus);
  pass &= nonempty(coreRayStatus) && nonempty(corePolicy);
  pass &= nonempty(shellRIStatus) && nonempty(shellAbsStatus);
  pass &= nonempty(shellRayStatus) && nonempty(shellPolicy);
  pass &= nonempty(coupRIStatus) && nonempty(coupAbsStatus);
  pass &= nonempty(coupRayStatus) && nonempty(coupPolicy);

  std::cout << "SchemaVersion=" << schema << "\n";
  std::cout << "ProjectRevision=" << revision << "\n";
  std::cout << "Core=" << coreProfile << " / " << coreMat << "\n";
  std::cout << "OuterShell enabled=" << useShell << " t=" << shellT << " mm profile="
            << shellProfile << " / " << shellMat << "\n";
  std::cout << "Coupling=" << coupProfile << " / " << coupMat << "\n";
  std::cout << "R2E MATERIAL PROFILE METADATA VALIDATION: " << (pass?"PASS":"FAIL") << "\n";
  f->Close();
}
