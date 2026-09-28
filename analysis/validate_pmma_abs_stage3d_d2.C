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

std::vector<std::pair<double,double>>
LoadTable(const std::string& path) {
  std::ifstream in(path);
  std::vector<std::pair<double,double>> out;
  std::string line;

  while (std::getline(in,line)) {
    if (line.empty() || line[0]=='#') continue;
    std::istringstream ss(line);
    double wl=0,L=0;
    if (ss>>wl>>L) out.push_back({wl,L});
  }

  std::sort(out.begin(),out.end());
  return out;
}

double InterpEnergy(
    const std::vector<std::pair<double,double>>& t,
    double wl) {

  constexpr double hc=1239.8419843320026;

  std::vector<std::pair<double,double>> e;

  for (const auto& p:t)
    e.push_back({hc/p.first,p.second});

  std::sort(e.begin(),e.end());

  const double x=hc/wl;

  if (x<=e.front().first) return e.front().second;
  if (x>=e.back().first) return e.back().second;

  auto it=std::lower_bound(
      e.begin(),e.end(),std::make_pair(x,-1e300));

  const auto hi=*it;
  const auto lo=*(it-1);

  const double f=(x-lo.first)/(hi.first-lo.first);

  return lo.second+f*(hi.second-lo.second);
}

}

void validate_pmma_abs_stage3d_d2() {

  constexpr double censor=3.0;
  constexpr Long64_t N=50000;

  const std::vector<double> wl={280.,400.,600.};

  struct Policy {
    std::string tag;
    std::string dataset;
  };

  const std::vector<Policy> policies={
    {
      "nominal",
      "PMMA_SNO_UVT_SampleA_ABS_Measured300to500_Model200to900_Candidate.txt"
    },
    {
      "uvopaque",
      "PMMA_SNO_UVT_SampleA_ABS_UVOpaqueBelow300_Model200to900_Candidate.txt"
    },
    {
      "vistransparent",
      "PMMA_SNO_UVT_SampleA_ABS_VisibleTransparentAbove500_Model200to900_Candidate.txt"
    }
  };

  bool allPass=true;

  std::cout
      << "\n=== PMMA STAGE 3D-D2 3-MM RUNTIME VALIDATION ===\n";

  for (const auto& policy:policies) {

    const auto table=LoadTable("../data/"+policy.dataset);

    for (int run=0;run<3;++run) {

      const double expectedL=InterpEnergy(table,wl[run]);
      const double pAbs=1.0-std::exp(-censor/expectedL);

      std::ostringstream fn;
      fn<<"validation_pmma_abs_stage3d_d2_"
        <<policy.tag<<"_run"<<run<<".root";

      TFile f(fn.str().c_str(),"READ");

      auto* p=dynamic_cast<TTree*>(f.Get("Photon"));

      if (!p) {
        std::cout<<policy.tag<<" run"<<run
                 <<" missing Photon : FAIL\n";
        allPass=false;
        continue;
      }

      Float_t endL=0;
      Int_t killed=0,oor=0;

      p->SetBranchAddress("EndTrackLength_mm",&endL);
      p->SetBranchAddress("KilledBySafetyCut",&killed);
      p->SetBranchAddress("AcrylicAbsOutOfRange",&oor);

      Long64_t absorbed=0;
      Long64_t badSafety=0;
      Long64_t oorCount=0;

      for (Long64_t i=0;i<p->GetEntries();++i) {
        p->GetEntry(i);

        if (killed) ++badSafety;
        if (oor) ++oorCount;

        if (std::isfinite(endL) &&
            endL>=0.0 &&
            endL<censor-1e-4)
          ++absorbed;
      }

      const double expectedAbs=N*pAbs;
      const double sigma=
          std::sqrt(N*pAbs*(1.0-pAbs));

      double z=0.0;

      if (sigma>1e-9)
        z=(absorbed-expectedAbs)/sigma;

      const Long64_t expectedOOR =
          (wl[run]<300.0 || wl[run]>500.0) ? N : 0;

      const bool countOK=p->GetEntries()==N;

      const bool fracOK=
          sigma>1e-9
          ? std::abs(z)<5.0
          : absorbed==N;

      const bool oorOK=oorCount==expectedOOR;
      const bool safetyOK=badSafety==0;

      const bool pass=
          countOK && fracOK && oorOK && safetyOK;

      allPass &= pass;

      const double obsSurv=
          1.0-double(absorbed)/double(N);

      std::cout<<std::fixed<<std::setprecision(6)
        <<policy.tag
        <<" "<<wl[run]<<" nm"
        <<" L="<<expectedL<<" mm"
        <<" expectedSurv="<<100.0*(1.0-pAbs)<<"%"
        <<" observedSurv="<<100.0*obsSurv<<"%"
        <<" z="<<z
        <<" OOR="<<oorCount
        <<" safety="<<badSafety
        <<" : "<<(pass?"PASS":"FAIL")
        <<"\n";
    }
  }

  std::cout
      <<"PMMA STAGE 3D-D2 3-MM RUNTIME VALIDATION: "
      <<(allPass?"PASS":"FAIL")
      <<"\n";
}
