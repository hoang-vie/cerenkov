#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double HC = 1239.8419843320026;

std::vector<std::pair<double,double>>
LoadTable(const char* fn) {

  std::ifstream in(fn);
  std::vector<std::pair<double,double>> out;
  std::string line;

  while (std::getline(in,line)) {
    if (line.empty() || line[0]=='#') continue;

    std::istringstream ss(line);

    double wl=0,L=0;

    if (ss>>wl>>L)
      out.push_back({wl,L});
  }

  return out;
}

double InterpEnergy(
    const std::vector<std::pair<double,double>>& t,
    double wl) {

  std::vector<std::pair<double,double>> e;

  for (const auto& p:t)
    e.push_back({HC/p.first,p.second});

  std::sort(e.begin(),e.end());

  const double E=HC/wl;

  if (E<=e.front().first) return e.front().second;
  if (E>=e.back().first)  return e.back().second;

  auto it=std::lower_bound(
      e.begin(),e.end(),
      std::make_pair(E,-1e300));

  const auto hi=*it;
  const auto lo=*(it-1);

  const double f=
      (E-lo.first)/(hi.first-lo.first);

  return lo.second+f*(hi.second-lo.second);
}

}

void validate_pmma_stage3f_production_default() {

  constexpr Long64_t N=50000;
  constexpr double path=3.0;

  const double wl[3]={280.,400.,600.};
  const Long64_t expectedOOR[3]={N,0,N};

  const auto table=LoadTable(
      "../data/"
      "PMMA_SNO_UVT_SampleA_ABS_Measured300to500_Model200to900.txt"
  );

  bool allPass=true;

  std::cout
      <<"\n=== PMMA STAGE 3F PRODUCTION DEFAULT VALIDATION ===\n";

  for (int run=0;run<3;++run) {

    const double expectedL=
        InterpEnergy(table,wl[run]);

    const double expectedSurv=
        std::exp(-path/expectedL);

    TString fn=
        Form(
          "validation_pmma_stage3f_production_default_run%d.root",
          run
        );

    TFile f(fn,"READ");

    auto* p=dynamic_cast<TTree*>(f.Get("Photon"));
    auto* rc=dynamic_cast<TTree*>(f.Get("RunConfig"));

    if (!p || !rc) {
      std::cout
          <<"run"<<run
          <<" missing Photon/RunConfig : FAIL\n";
      allPass=false;
      continue;
    }

    Float_t endL=0;
    Int_t oor=0;
    Int_t killed=0;

    p->SetBranchAddress("EndTrackLength_mm",&endL);
    p->SetBranchAddress("AcrylicAbsOutOfRange",&oor);
    p->SetBranchAddress("KilledBySafetyCut",&killed);

    Long64_t absorbed=0;
    Long64_t oorCount=0;
    Long64_t safety=0;

    for (Long64_t i=0;i<p->GetEntries();++i) {

      p->GetEntry(i);

      if (oor) ++oorCount;
      if (killed) ++safety;

      if (std::isfinite(endL) &&
          endL>=0.0 &&
          endL<path-1e-4)
        ++absorbed;
    }

    const double observedSurv=
        1.0-double(absorbed)/double(N);

    const double pAbs=
        1.0-expectedSurv;

    const double sigma=
        std::sqrt(
          N*pAbs*(1.0-pAbs));

    const double expectedAbs=
        N*pAbs;

    const double z=
        sigma>1e-12
        ? (double(absorbed)-expectedAbs)/sigma
        : 0.0;

    char schema[128]={0};
    char revision[256]={0};
    char registry[256]={0};
    char profile[128]={0};
    char dataset[512]={0};
    char status[512]={0};
    char policy[1024]={0};

    Float_t validMin=-1,validMax=-1;
    Float_t measuredMin=-1,measuredMax=-1;

    rc->SetBranchAddress("SchemaVersion",schema);
    rc->SetBranchAddress("ProjectRevision",revision);
    rc->SetBranchAddress(
        "MaterialProfileRegistryVersion",registry);

    rc->SetBranchAddress("CoreProfile",profile);
    rc->SetBranchAddress("CoreProfileAbsDataset",dataset);

    rc->SetBranchAddress(
        "CoreProfileAbsValidMin_nm",&validMin);
    rc->SetBranchAddress(
        "CoreProfileAbsValidMax_nm",&validMax);

    rc->SetBranchAddress(
        "CoreProfileAbsMeasuredMin_nm",&measuredMin);
    rc->SetBranchAddress(
        "CoreProfileAbsMeasuredMax_nm",&measuredMax);

    rc->SetBranchAddress(
        "CoreProfileAbsStatus",status);

    rc->SetBranchAddress(
        "CoreProfileOutOfRangePolicy",policy);

    rc->GetEntry(0);

    const bool metadataOK=
        rc->GetListOfBranches()->GetEntries()==154 &&
        std::string(schema)=="general-v1.8" &&
        std::string(revision)==
          "hotfix6_physics_audit5_pmma_stage3f_production_freeze" &&
        std::string(registry)==
          "curated-bulk-v2_absMeasuredRange" &&
        std::string(profile)=="PMMA_SNO_A" &&
        std::string(dataset)==
          "PMMA_SNO_UVT_SampleA_ABS_Measured300to500_Model200to900.txt" &&
        std::abs(validMin-200.0f)<1e-4f &&
        std::abs(validMax-900.0f)<1e-4f &&
        std::abs(measuredMin-300.0f)<1e-4f &&
        std::abs(measuredMax-500.0f)<1e-4f &&
        std::string(status)==
          "SNO_UVT_SampleA_Measured300to500_plus_explicit_transport_model_200to900";

    const std::string pol(policy);

    const bool policyOK=
        pol.find("Measured300to500")!=std::string::npos &&
        pol.find("UV250to300_AnchoredExp300to320")!=std::string::npos &&
        pol.find("UV200to250_Hold250")!=std::string::npos &&
        pol.find("Visible500to900_Hold500")!=std::string::npos &&
        pol.find("OutsideSupport_Geant4EndpointClamp")!=std::string::npos;

    const bool countOK=
        p->GetEntries()==N;

    const bool oorOK=
        oorCount==expectedOOR[run];

    const bool safetyOK=
        safety==0;

    const bool physicsOK=
        std::abs(z)<5.0;

    const bool pass=
        metadataOK &&
        policyOK &&
        countOK &&
        oorOK &&
        safetyOK &&
        physicsOK;

    allPass &= pass;

    std::cout
      <<std::fixed<<std::setprecision(6)
      <<"run"<<run
      <<" "<<wl[run]<<" nm"
      <<" L="<<expectedL<<" mm"
      <<" expectedSurv="<<100.0*expectedSurv<<"%"
      <<" observedSurv="<<100.0*observedSurv<<"%"
      <<" z="<<z
      <<" OOR="<<oorCount
      <<" metadata="<<(metadataOK?"OK":"BAD")
      <<" policy="<<(policyOK?"OK":"BAD")
      <<" safety="<<safety
      <<" : "<<(pass?"PASS":"FAIL")
      <<"\n";
  }

  std::cout
      <<"PMMA STAGE 3F PRODUCTION DEFAULT VALIDATION: "
      <<(allPass?"PASS":"FAIL")
      <<"\n";
}
