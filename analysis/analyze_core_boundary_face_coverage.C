#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {
using Key = std::pair<int,int>;
const char* FaceName(int f){
  static const char* n[]={"+X","-X","+Z","-Z","+Y","-Y"};
  return (f>=0&&f<6)?n[f]:"?";
}
double pct(long long a,long long b){return b?100.0*double(a)/double(b):0.0;}
}

void analyze_core_boundary_face_coverage(
  const char* inputName="study_cry_core_boundary_fullall_1000_run0.root")
{
  TFile f(inputName,"READ");
  auto* cb=(TTree*)f.Get("CoreBoundaryCrossing");
  if(!cb){std::cerr<<"Missing CoreBoundaryCrossing\n";return;}

  int event=0,root=0,face=0;
  cb->SetBranchAddress("EventID",&event);
  cb->SetBranchAddress("RootPrimaryTrackID",&root);
  cb->SetBranchAddress("FaceID",&face);

  std::map<Key,std::array<long long,6>> perPrimary;
  std::array<long long,6> totalByFace{};
  long long total=0,bad=0;

  for(Long64_t i=0;i<cb->GetEntries();++i){
    cb->GetEntry(i);
    if(face<0||face>=6){++bad;continue;}
    ++perPrimary[{event,root}][face];
    ++totalByFace[face];
    ++total;
  }

  int bestSide=0;
  for(int f=1;f<=3;++f)
    if(totalByFace[f]>totalByFace[bestSide]) bestSide=f;

  struct Row{std::string name; std::set<int> faces;};
  std::vector<Row> rows={
    {"bottom_only",{5}},
    {"bottom_plus_best_1_side",{5,bestSide}},
    {"bottom_plus_+X",{5,0}},
    {"bottom_plus_-X",{5,1}},
    {"bottom_plus_+Z",{5,2}},
    {"bottom_plus_-Z",{5,3}},
    {"bottom_plus_X_pair",{5,0,1}},
    {"bottom_plus_Z_pair",{5,2,3}},
    {"all_sides_only",{0,1,2,3}},
    {"bottom_plus_4_sides",{5,0,1,2,3}},
    {"all_6_faces",{0,1,2,3,4,5}}
  };

  std::cout<<std::fixed<<std::setprecision(2);
  std::cout<<"\nCORE-BOUNDARY FACE COVERAGE STUDY\n";
  std::cout<<"Input: "<<inputName<<"\n";
  std::cout<<"Useful root primaries: "<<perPrimary.size()<<"\n";
  std::cout<<"Total firstExit photons: "<<total<<"\n";
  std::cout<<"Bad FaceID: "<<bad<<"\n\n";

  std::cout<<"=== RAW FACE SHARE ===\n";
  for(int f=0;f<6;++f)
    std::cout<<FaceName(f)<<" photons="<<totalByFace[f]
             <<" share="<<pct(totalByFace[f],total)<<"%\n";

  std::cout<<"\nBest single side = "<<FaceName(bestSide)<<"\n\n";

  std::cout<<std::left
           <<std::setw(28)<<"scenario"
           <<std::setw(8)<<"Nface"
           <<std::setw(14)<<"photons"
           <<std::setw(14)<<"photon_cov"
           <<std::setw(14)<<"event_cov"
           <<std::setw(16)<<"dominant_cov"<<"\n";

  for(const auto& r:rows){
    long long capPhot=0,capEvt=0,capDom=0;
    for(const auto& kv:perPrimary){
      const auto& a=kv.second;
      long long cap=0;
      for(int f:r.faces) cap+=a[f];
      if(cap>0) ++capEvt;
      capPhot+=cap;
      long long mx=*std::max_element(a.begin(),a.end());
      bool dom=false;
      for(int f:r.faces) if(a[f]==mx && mx>0) dom=true;
      if(dom) ++capDom;
    }
    std::cout<<std::left
             <<std::setw(28)<<r.name
             <<std::setw(8)<<r.faces.size()
             <<std::setw(14)<<capPhot
             <<std::setw(13)<<pct(capPhot,total)<<"%"
             <<std::setw(13)<<pct(capEvt,perPrimary.size())<<"%"
             <<std::setw(15)<<pct(capDom,perPrimary.size())<<"%"
             <<"\n";
  }

  std::cout<<"\n=== MARGINAL GAIN OVER BOTTOM ONLY ===\n";
  for(int f=0;f<5;++f)
    std::cout<<"add "<<FaceName(f)
             <<" : +"<<totalByFace[f]<<" photons, +"
             <<pct(totalByFace[f],total)<<" percentage points\n";

  std::cout<<"\nNote: geometric boundary coverage only; no sensor area, concentrator acceptance, PDE, or cost applied.\n";
  std::cout<<"\nFACE COVERAGE STRUCTURAL SANITY: "
           <<((bad==0 && total==cb->GetEntries())?"PASS":"FAIL")<<"\n";
}
