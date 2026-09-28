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

constexpr double HC=1239.8419843320026;

std::vector<std::pair<double,double>>
Load(const char* fn) {
  std::ifstream in(fn);
  std::vector<std::pair<double,double>> v;
  std::string line;

  while (std::getline(in,line)) {
    if (line.empty() || line[0]=='#')
      continue;

    std::istringstream ss(line);
    double x=0,y=0;

    if (ss>>x>>y)
      v.push_back({x,y});
  }

  return v;
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

double Deg(double x) {
  return x*180.0/M_PI;
}

}

void validate_pdms_interface_stage4c_pdms_to_pmma_scan() {

  constexpr Long64_t N=100000;
  constexpr double wl=400.0;

  const double angles[6]={
      0.,30.,45.,46.66,60.,70.
  };

  const auto pdms=
      Load(
        "../data/"
        "OpticalGrease_PDMS_RINDEX_"
        "Roychowdhury2018_Cauchy3Anchor.txt"
      );

  const auto pmma=
      Load(
        "../data/"
        "PMMA_RINDEX_Roychowdhury2020_Cauchy3Anchor.txt"
      );

  const double n1=InterpEnergy(pdms,wl);
  const double n2=InterpEnergy(pmma,wl);

  bool allPass=true;

  std::cout
      <<"\n=== PDMS STAGE 4C "
      <<"PDMS->PMMA FRESNEL/SNELL VALIDATION ===\n";

  for (int run=0;run<12;++run) {

    const int ia=run/2;
    const bool isP=(run%2)==1;

    const double ti=
        angles[ia]*M_PI/180.0;

    const double thetaT=
        std::asin(
          (n1/n2)*std::sin(ti)
        );

    const double ci=std::cos(ti);
    const double ct=std::cos(thetaT);

    const double Rs=
        std::pow(
          (n1*ci-n2*ct)/
          (n1*ci+n2*ct),2);

    const double Rp=
        std::pow(
          (n2*ci-n1*ct)/
          (n2*ci+n1*ct),2);

    const double expectedR=
        isP ? Rp : Rs;

    std::ostringstream fn;

    fn
      <<"validation_pdms_interface_stage4c_"
      <<"pdms_to_pmma_scan_run"
      <<run
      <<".root";

    TFile f(fn.str().c_str(),"READ");

    auto* b=
        dynamic_cast<TTree*>(
            f.Get("BoundaryDebug")
        );

    if (!b) {
      std::cout
          <<"run"<<run
          <<" missing BoundaryDebug : FAIL\n";
      allPass=false;
      continue;
    }

    Int_t step=0;
    Float_t x=0;
    Float_t preX=0,preY=0,preZ=0;
    Float_t postX=0,postY=0,postZ=0;
    char status[128]={0};

    b->SetBranchAddress("StepNumber",&step);
    b->SetBranchAddress("X_mm",&x);

    b->SetBranchAddress("PreDirX",&preX);
    b->SetBranchAddress("PreDirY",&preY);
    b->SetBranchAddress("PreDirZ",&preZ);

    b->SetBranchAddress("PostDirX",&postX);
    b->SetBranchAddress("PostDirY",&postY);
    b->SetBranchAddress("PostDirZ",&postZ);

    b->SetBranchAddress(
        "BoundaryStatusName",status);

    Long64_t first=0;
    Long64_t refl=0;
    Long64_t refr=0;
    Long64_t tir=0;
    Long64_t other=0;

    double maxSnellErr=0.0;
    double maxReflVecErr=0.0;

    for (Long64_t i=0;i<b->GetEntries();++i) {

      b->GetEntry(i);

      const bool firstInterface=
          step==1 &&
          std::abs(double(x)-53.0)<1e-4 &&
          preX<0.0;

      if (!firstInterface)
        continue;

      ++first;

      const std::string s(status);

      if (s=="FresnelReflection") {

        ++refl;

        // Incoming (-cos,+sin,0)
        // reflected (+cos,+sin,0)

        const double ex= std::cos(ti);
        const double ey= std::sin(ti);
        const double ez=0.0;

        const double dr=
            std::sqrt(
              std::pow(postX-ex,2)+
              std::pow(postY-ey,2)+
              std::pow(postZ-ez,2));

        maxReflVecErr=
            std::max(maxReflVecErr,dr);

      } else if (s=="FresnelRefraction") {

        ++refr;

        const double observedTheta=
            std::atan2(
              std::sqrt(
                double(postY)*postY+
                double(postZ)*postZ),
              -double(postX));

        maxSnellErr=
            std::max(
              maxSnellErr,
              std::abs(
                Deg(observedTheta-thetaT)));

      } else if (
          s=="TotalInternalReflection") {

        ++tir;

      } else {

        ++other;
      }
    }

    const double expectedCount=
        N*expectedR;

    const double sigma=
        std::sqrt(
          N*expectedR*(1.0-expectedR));

    double z=0.0;

    if (sigma>1e-12)
      z=
        (double(refl)-expectedCount)/sigma;

    const bool pass=
        first==N &&
        refl+refr==N &&
        tir==0 &&
        other==0 &&
        std::abs(z)<5.0 &&
        maxSnellErr<1e-3 &&
        maxReflVecErr<1e-5;

    allPass &= pass;

    std::cout
      <<std::fixed<<std::setprecision(6)
      <<"run"<<run
      <<" theta="<<angles[ia]
      <<" pol="<<(isP?"p":"s")
      <<" Rexp="<<100.0*expectedR<<"%"
      <<" Rmc="
      <<100.0*double(refl)/double(first)<<"%"
      <<" refl="<<refl
      <<" refr="<<refr
      <<" TIR="<<tir
      <<" z="<<z
      <<" maxSnellErr="<<maxSnellErr
      <<" maxReflVecErr="<<maxReflVecErr
      <<" first="<<first
      <<" other="<<other
      <<" : "<<(pass?"PASS":"FAIL")
      <<"\n";
  }

  std::cout
      <<"PDMS STAGE 4C PDMS->PMMA "
      <<"FRESNEL/SNELL VALIDATION: "
      <<(allPass?"PASS":"FAIL")
      <<"\n";
}
