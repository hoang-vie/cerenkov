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
    double wl=0,L=0;
    if (ss>>wl>>L) out.push_back({wl,L});
  }
  return out;
}

double Interp(const std::vector<std::pair<double,double>>& t, double x) {
  if (t.empty()) return std::numeric_limits<double>::quiet_NaN();
  if (x<=t.front().first) return t.front().second;
  if (x>=t.back().first) return t.back().second;
  auto it=std::lower_bound(t.begin(),t.end(),x,
      [](const auto& a,double v){return a.first<v;});
  const auto& b=*it; const auto& a=*(it-1);
  const double f=(x-a.first)/(b.first-a.first);
  return a.second+f*(b.second-a.second);
}

double ConditionalExpCDF(double x,double lambda,double censor) {
  const double den=1.0-std::exp(-censor/lambda);
  if (den<=0) return 0;
  return (1.0-std::exp(-x/lambda))/den;
}

}

void validate_water_rayleigh_pmh_runtime(
    const char* prefix="validation_water_rayleigh_pmh_runtime") {
  const double wavelengths[4]={300.0,400.0,500.0,600.0};
  constexpr double censor_mm=250000.0;
  constexpr long long expectedEvents=20000;
  const auto candidate=
      LoadTable("../data/Water_RAYLEIGH_ZhangHu2009_PMH_21p5C.txt");
  if (candidate.empty()) {
    std::cerr<<"Cannot load candidate Rayleigh table\n";
    return;
  }

  bool allPass=true;
  std::cout<<"\n=== WATER RAYLEIGH PMH RUNTIME VALIDATION ===\n";

  for (int run=0;run<4;++run) {
    const double wl=wavelengths[run];
    const double expected_m=Interp(candidate,wl);
    const double expected_mm=1000.0*expected_m;
    const std::string fn=std::string(prefix)+"_run"+std::to_string(run)+".root";

    TFile f(fn.c_str(),"READ");
    if (f.IsZombie()) {
      std::cerr<<"Cannot open "<<fn<<"\n";
      allPass=false; continue;
    }
    auto* p=dynamic_cast<TTree*>(f.Get("Photon"));
    auto* c=dynamic_cast<TTree*>(f.Get("RunConfig"));
    if (!p || !c || c->GetEntries()!=1) {
      std::cerr<<"Missing Photon/RunConfig in "<<fn<<"\n";
      allPass=false; continue;
    }

    char schema[128]={};
    char revision[256]={};
    char dataset[512]={};
    char status[512]={};
    c->SetBranchAddress("SchemaVersion",schema);
    c->SetBranchAddress("ProjectRevision",revision);
    c->SetBranchAddress("CoreProfileRayleighDataset",dataset);
    c->SetBranchAddress("CoreProfileRayleighStatus",status);
    c->GetEntry(0);

    int nRay=0;
    float first=std::numeric_limits<float>::quiet_NaN();
    p->SetBranchAddress("NRayleigh",&nRay);
    p->SetBranchAddress("FirstRayleighTrackLength_mm",&first);

    long long scat=0,cens=0;
    double exposure=0.0;
    std::vector<double> observed;
    observed.reserve(p->GetEntries());

    for (Long64_t i=0;i<p->GetEntries();++i) {
      p->GetEntry(i);
      if (nRay>0 && std::isfinite(first) && first>=0 && first<censor_mm+1.0) {
        ++scat;
        exposure+=first;
        observed.push_back(first);
      } else {
        ++cens;
        exposure+=censor_mm;
      }
    }

    const double lambdaHat_mm=scat?exposure/double(scat):0.0;
    const double rel=expected_mm>0?std::abs(lambdaHat_mm/expected_mm-1.0):999.0;

    const double pScatter=1.0-std::exp(-censor_mm/expected_mm);
    const double expectedScat=expectedEvents*pScatter;
    const double sigmaScat=std::sqrt(expectedEvents*pScatter*(1.0-pScatter));
    const double z=sigmaScat>0?(double(scat)-expectedScat)/sigmaScat:0.0;

    std::sort(observed.begin(),observed.end());
    double ks=0.0;
    const double n=double(observed.size());
    for (size_t i=0;i<observed.size();++i) {
      const double F=ConditionalExpCDF(observed[i],expected_mm,censor_mm);
      const double lo=double(i)/n;
      const double hi=double(i+1)/n;
      ks=std::max(ks,std::max(std::abs(F-lo),std::abs(hi-F)));
    }
    const double ksLimit=n>0?1.63/std::sqrt(n):0.0; // ~1% one-sample KS threshold

    const bool meta=
        std::string(schema)=="general-v1.7" &&
        std::string(revision)=="hotfix6_physics_audit5_water_stage3e_pmh_default" &&
        std::string(dataset)=="Water_RAYLEIGH_ZhangHu2009_PMH_21p5C.txt" &&
        std::string(status)=="ZhangHu2009_PMH_21p5C_delta0p039_MFPValidated_AngularApproxDocumented";
    const bool countOK=p->GetEntries()==expectedEvents;
    const bool mleOK=scat>1000 && rel<0.04;
    const bool fractionOK=std::abs(z)<5.0;
    const bool shapeOK=scat>1000 && ks<ksLimit;
    const bool pass=meta && countOK && mleOK && fractionOK && shapeOK;
    allPass &= pass;

    std::cout<<std::fixed<<std::setprecision(4)
      <<"run"<<run<<" "<<wl<<" nm"
      <<" expectedL="<<expected_m<<" m"
      <<" MLE="<<lambdaHat_mm/1000.0<<" m"
      <<" rel="<<100.0*rel<<"%"
      <<" scat="<<scat<<" cens="<<cens
      <<" zFrac="<<z
      <<" KS="<<ks<<" limit="<<ksLimit
      <<" metadata="<<(meta?"OK":"FAIL")
      <<" : "<<(pass?"PASS":"FAIL")<<"\n";

    if (!meta) {
      std::cout<<"  SchemaVersion="<<schema
               <<" ProjectRevision="<<revision
               <<" CoreProfileRayleighDataset="<<dataset
               <<" CoreProfileRayleighStatus="<<status<<"\n";
    }
  }

  std::cout<<"WATER RAYLEIGH PMH RUNTIME VALIDATION: "
           <<(allPass?"PASS":"FAIL")<<"\n";
}
