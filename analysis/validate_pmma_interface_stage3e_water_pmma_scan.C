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

constexpr double kHC = 1239.8419843320026;

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
    e.push_back({kHC/p.first,p.second});

  std::sort(e.begin(),e.end());

  const double E=kHC/wl;

  if (E<=e.front().first)
    return e.front().second;

  if (E>=e.back().first)
    return e.back().second;

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

void validate_pmma_interface_stage3e_water_pmma_scan() {

  constexpr Long64_t N=100000;
  constexpr double wl=400.0;

  const auto water=
      Load("../data/Water_RINDEX_Daimon2007_21p5C.txt");

  const auto pmma=
      Load("../data/PMMA_RINDEX_Roychowdhury2020_Cauchy3Anchor.txt");

  const double n1=InterpEnergy(water,wl);
  const double n2=InterpEnergy(pmma,wl);

  const double angles[5]={0.,30.,50.,60.,70.};

  bool allPass=true;

  std::cout
      << "\n=== PMMA STAGE 3E WATER->PMMA FRESNEL/SNELL VALIDATION ===\n";

  for (int run=0;run<10;++run) {

    const int ia=run/2;
    const bool isP=(run%2)==1;

    const double thetaI=angles[ia]*M_PI/180.0;

    const double sinT=
        (n1/n2)*std::sin(thetaI);

    const double thetaT=
        std::asin(sinT);

    const double ci=std::cos(thetaI);
    const double ct=std::cos(thetaT);

    const double Rs=
        std::pow(
          (n1*ci-n2*ct)/(n1*ci+n2*ct),2);

    const double Rp=
        std::pow(
          (n2*ci-n1*ct)/(n2*ci+n1*ct),2);

    const double expectedR=
        isP ? Rp : Rs;

    std::ostringstream fn;

    fn<<"validation_pmma_interface_stage3e_water_to_pmma_scan_run"
      <<run<<".root";

    TFile f(fn.str().c_str(),"READ");

    auto* b=dynamic_cast<TTree*>(f.Get("BoundaryDebug"));

    if (!b) {
      std::cout<<"run"<<run<<" missing BoundaryDebug : FAIL\n";
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

    b->SetBranchAddress("BoundaryStatusName",status);

    Long64_t first=0;
    Long64_t reflected=0;
    Long64_t refracted=0;
    Long64_t tir=0;
    Long64_t other=0;

    double maxSnellResidualDeg=0.0;
    double maxReflectionVectorResidual=0.0;

    for (Long64_t i=0;i<b->GetEntries();++i) {

      b->GetEntry(i);

      const bool firstInterface =
          step==1 &&
          std::abs(double(x)-50.0)<1e-4 &&
          preX>0.0;

      if (!firstInterface)
        continue;

      ++first;

      const std::string s(status);

      if (s=="FresnelReflection") {

        ++reflected;

        const double ex=-std::cos(thetaI);
        const double ey= std::sin(thetaI);
        const double ez=0.0;

        const double dr=
            std::sqrt(
              std::pow(postX-ex,2)+
              std::pow(postY-ey,2)+
              std::pow(postZ-ez,2));

        maxReflectionVectorResidual=
            std::max(maxReflectionVectorResidual,dr);

      } else if (s=="FresnelRefraction") {

        ++refracted;

        const double observedTheta=
            std::atan2(
              std::sqrt(
                double(postY)*postY+
                double(postZ)*postZ),
              double(postX));

        maxSnellResidualDeg=
            std::max(
              maxSnellResidualDeg,
              std::abs(
                Deg(observedTheta-thetaT)));

      } else if (s=="TotalInternalReflection") {

        ++tir;

      } else {

        ++other;
      }
    }

    const double expectedCount=N*expectedR;

    const double sigma=
        std::sqrt(
          N*expectedR*(1.0-expectedR));

    const double z=
        sigma>1e-12
        ? (double(reflected)-expectedCount)/sigma
        : 0.0;

    const bool countOK=
        first==N;

    const bool statusOK=
        reflected+refracted==N &&
        tir==0 &&
        other==0;

    const bool fresnelOK=
        std::abs(z)<5.0;

    const bool snellOK=
        maxSnellResidualDeg<1e-3;

    const bool reflectionDirOK=
        maxReflectionVectorResidual<1e-5;

    const bool pass=
        countOK &&
        statusOK &&
        fresnelOK &&
        snellOK &&
        reflectionDirOK;

    allPass &= pass;

    std::cout
      <<std::fixed<<std::setprecision(6)
      <<"run"<<run
      <<" theta="<<angles[ia]
      <<" pol="<<(isP?"p":"s")
      <<" Rexp="<<100.0*expectedR<<"%"
      <<" Rmc="<<100.0*double(reflected)/double(first)<<"%"
      <<" z="<<z
      <<" thetaT="<<Deg(thetaT)
      <<" maxSnellErr="<<maxSnellResidualDeg<<" deg"
      <<" maxReflVecErr="<<maxReflectionVectorResidual
      <<" first="<<first
      <<" TIR="<<tir
      <<" other="<<other
      <<" : "<<(pass?"PASS":"FAIL")
      <<"\n";
  }

  std::cout
      <<"PMMA STAGE 3E WATER->PMMA FRESNEL/SNELL VALIDATION: "
      <<(allPass?"PASS":"FAIL")
      <<"\n";
}
