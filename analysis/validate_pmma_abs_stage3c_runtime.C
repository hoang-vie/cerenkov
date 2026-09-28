#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

std::vector<std::pair<double,double>> LoadTable(const char* path) {
  std::ifstream in(path);
  std::vector<std::pair<double,double>> out;
  std::string line;
  while (std::getline(in,line)) {
    if (line.empty() || line[0]=='#') continue;
    std::istringstream ss(line);
    double wl=0, L=0;
    if (ss>>wl>>L) out.push_back({wl,L});
  }
  std::sort(out.begin(),out.end());
  return out;
}

double InterpEnergy(const std::vector<std::pair<double,double>>& t, double wl_nm) {
  constexpr double hc=1239.8419843320026; // eV nm
  if (t.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::vector<std::pair<double,double>> e;
  e.reserve(t.size());
  for (const auto& p:t) e.push_back({hc/p.first,p.second});
  std::sort(e.begin(),e.end());
  const double x=hc/wl_nm;
  if (x<=e.front().first) return e.front().second;
  if (x>=e.back().first) return e.back().second;
  auto it=std::lower_bound(e.begin(),e.end(),std::make_pair(x,-1e300),
      [](const auto& a,const auto& b){ return a.first<b.first; });
  if (it==e.begin()) return it->second;
  const auto hi=*it, lo=*(it-1);
  const double f=(x-lo.first)/(hi.first-lo.first);
  return lo.second+f*(hi.second-lo.second);
}

double ConditionalExpCDF(double x, double L, double censor) {
  if (x<=0) return 0.0;
  if (x>=censor) return 1.0;
  const double den=1.0-std::exp(-censor/L);
  return (1.0-std::exp(-x/L))/den;
}

}

void validate_pmma_abs_stage3c_runtime(
    const char* prefix="validation_pmma_abs_stage3c_runtime") {

  constexpr double censor_mm=1000.0;
  constexpr Long64_t expectedEvents=20000;
  const std::vector<double> wavelengths={300.0,310.0,350.0,400.0,500.0};
  const auto candidate=LoadTable(
      "../data/PMMA_SNO_UVT_SampleA_ABS_TableI_RT_plus_TableII_300to500_Candidate.txt");
  if (candidate.empty()) {
    std::cerr<<"Cannot load PMMA candidate absorption table\n";
    return;
  }

  bool allPass=true;
  std::cout<<"\n=== PMMA SNO DENSIFIED ABS RUNTIME VALIDATION ===\n";

  for (int run=0; run<int(wavelengths.size()); ++run) {
    const double wl=wavelengths[run];
    const double expected_mm=InterpEnergy(candidate,wl);

    std::ostringstream fn;
    fn<<prefix<<"_run"<<run<<".root";
    TFile f(fn.str().c_str(),"READ");
    if (f.IsZombie()) {
      std::cout<<"run"<<run<<" "<<wl<<" nm: missing file : FAIL\n";
      allPass=false;
      continue;
    }

    auto* p=dynamic_cast<TTree*>(f.Get("Photon"));
    auto* rc=dynamic_cast<TTree*>(f.Get("RunConfig"));
    if (!p || !rc || rc->GetEntries()!=1) {
      std::cout<<"run"<<run<<" "<<wl<<" nm: missing Photon/RunConfig : FAIL\n";
      allPass=false;
      continue;
    }

    Int_t runID=0;
    Float_t endL=0;
    Int_t killed=0, oor=0;
    p->SetBranchAddress("RunID",&runID);
    p->SetBranchAddress("EndTrackLength_mm",&endL);
    p->SetBranchAddress("KilledBySafetyCut",&killed);
    p->SetBranchAddress("AcrylicAbsOutOfRange",&oor);

    char schema[128]={0}, revision[256]={0}, coreProfile[128]={0};
    char dataset[512]={0}, status[512]={0}, policy[512]={0};
    Float_t validMin=-1, validMax=-1;
    rc->SetBranchAddress("SchemaVersion",schema);
    rc->SetBranchAddress("ProjectRevision",revision);
    rc->SetBranchAddress("CoreProfile",coreProfile);
    rc->SetBranchAddress("CoreProfileAbsDataset",dataset);
    rc->SetBranchAddress("CoreProfileAbsValidMin_nm",&validMin);
    rc->SetBranchAddress("CoreProfileAbsValidMax_nm",&validMax);
    rc->SetBranchAddress("CoreProfileAbsStatus",status);
    rc->SetBranchAddress("CoreProfileOutOfRangePolicy",policy);
    rc->GetEntry(0);

    Long64_t absorbed=0, cens=0, badSafety=0, badOOR=0;
    double exposure=0.0;
    std::vector<double> observed;
    observed.reserve(p->GetEntries());

    for (Long64_t i=0;i<p->GetEntries();++i) {
      p->GetEntry(i);
      if (killed) ++badSafety;
      if (oor) ++badOOR;
      if (std::isfinite(endL) && endL>=0.0 && endL<censor_mm-1e-4) {
        ++absorbed;
        exposure += endL;
        observed.push_back(endL);
      } else {
        ++cens;
        exposure += censor_mm;
      }
    }

    const double mle=absorbed?exposure/double(absorbed):0.0;
    const double rel=std::abs(mle-expected_mm)/expected_mm;

    const double pAbs=1.0-std::exp(-censor_mm/expected_mm);
    const double expectedAbs=expectedEvents*pAbs;
    const double sigma=std::sqrt(expectedEvents*pAbs*(1.0-pAbs));
    const double z=sigma>1e-12?(double(absorbed)-expectedAbs)/sigma:0.0;

    std::sort(observed.begin(),observed.end());
    double ks=0.0;
    const double n=double(observed.size());
    if (n>0) {
      for (size_t i=0;i<observed.size();++i) {
        const double F=ConditionalExpCDF(observed[i],expected_mm,censor_mm);
        const double lo=double(i)/n;
        const double hi=double(i+1)/n;
        ks=std::max(ks,std::max(std::abs(F-lo),std::abs(hi-F)));
      }
    }
    const double ksLimit=n>0?1.63/std::sqrt(n):0.0;

    const bool metadataOK=
        std::string(schema)=="general-v1.7" &&
        std::string(revision)=="hotfix6_physics_audit5_pmma_stage3c_abs_runtime" &&
        std::string(coreProfile)=="PMMA_SNO_A" &&
        std::string(dataset)=="PMMA_SNO_UVT_SampleA_ABS_TableI_RT_plus_TableII_300to500_Candidate.txt" &&
        std::string(status)=="RuntimeValidationOverride_via_pmmaAbsFile" &&
        std::abs(validMin-300.0f)<1e-4f &&
        std::abs(validMax-500.0f)<1e-4f &&
        std::string(policy).find("Geant4EndpointClamp")!=std::string::npos;

    const bool countOK=p->GetEntries()==expectedEvents;
    const bool mleOK=absorbed>5000 && rel<0.035;
    const bool fracOK=std::abs(z)<5.0;
    const bool shapeOK=absorbed>5000 && ks<ksLimit;
    const bool flagsOK=badSafety==0 && badOOR==0;
    const bool pass=metadataOK && countOK && mleOK && fracOK && shapeOK && flagsOK;
    allPass &= pass;

    std::cout<<std::fixed<<std::setprecision(4)
      <<"run"<<run<<" "<<wl<<" nm"
      <<" expectedL="<<expected_mm<<" mm"
      <<" MLE="<<mle<<" mm"
      <<" rel="<<100.0*rel<<"%"
      <<" abs="<<absorbed<<" cens="<<cens
      <<" zFrac="<<z
      <<" KS="<<ks<<" limit="<<ksLimit
      <<" OOR="<<badOOR<<" safety="<<badSafety
      <<" metadata="<<(metadataOK?"OK":"FAIL")
      <<" : "<<(pass?"PASS":"FAIL")<<"\n";

    if (!metadataOK) {
      std::cout<<"  schema="<<schema
               <<" revision="<<revision
               <<" coreProfile="<<coreProfile
               <<" dataset="<<dataset
               <<" valid=["<<validMin<<","<<validMax<<"]"
               <<" status="<<status
               <<" policy="<<policy<<"\n";
    }
  }

  std::cout<<"PMMA SNO DENSIFIED ABS RUNTIME VALIDATION: "
           <<(allPass?"PASS":"FAIL")<<"\n";
}
