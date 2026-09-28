#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <utility>

namespace {
constexpr double PI=3.14159265358979323846;
using Key=std::pair<int,int>;
const char* FN(int f){static const char* n[]={"+X","-X","+Z","-Z","+Y","-Y"}; return (f>=0&&f<6)?n[f]:"?";}
double dg(double r){return r*180.0/PI;}
double cl(double x){return std::max(-1.0,std::min(1.0,x));}
int zb(double t){double e[]={0,15,30,45,60,75,90.000001}; for(int i=0;i<6;++i) if(t>=e[i]&&t<e[i+1]) return i; return -1;}
const char* ZL(int b){static const char* s[]={"0-15","15-30","30-45","45-60","60-75","75-90"}; return (b>=0&&b<6)?s[b]:"?";}
}

void analyze_cry_core_boundary_vs_primary_v2(const char* inputName="study_cry_core_boundary_fullall_1000_run0.root")
{
  TFile f(inputName,"READ");
  auto* pt=(TTree*)f.Get("PrimaryTruth");
  auto* cb=(TTree*)f.Get("CoreBoundaryCrossing");
  if(!pt||!cb){std::cerr<<"Missing required tree\n";return;}

  int e=0,tr=0,pdg=0; float ke=0,dx=0,dy=0,dz=0;
  pt->SetBranchAddress("EventID",&e); pt->SetBranchAddress("TrackID",&tr);
  pt->SetBranchAddress("PDG",&pdg); pt->SetBranchAddress("VertexKE_MeV",&ke);
  pt->SetBranchAddress("DirX",&dx); pt->SetBranchAddress("DirY",&dy); pt->SetBranchAddress("DirZ",&dz);

  struct P{int pdg; double ke,dx,dy,dz,th,ph; int b;};
  std::map<Key,P> Pm; std::map<int,int> mult; std::array<long long,6> nprim{}; long long nonDown=0;
  for(Long64_t i=0;i<pt->GetEntries();++i){
    pt->GetEntry(i); double th=dg(std::acos(cl(-double(dy)))); double ph=dg(std::atan2(double(dz),double(dx))); int b=zb(th);
    Pm[{e,tr}]={pdg,ke,dx,dy,dz,th,ph,b}; ++mult[e]; if(b>=0) ++nprim[b]; if(dy>=0) ++nonDown;
  }

  long long multi=0,maxm=0; for(auto &kv:mult){if(kv.second>1)++multi; maxm=std::max<long long>(maxm,kv.second);}

  int ce=0,root=0,face=0; cb->SetBranchAddress("EventID",&ce); cb->SetBranchAddress("RootPrimaryTrackID",&root); cb->SetBranchAddress("FaceID",&face);
  std::array<std::array<long long,6>,6> by{}; std::array<long long,6> nc{};
  std::map<Key,std::array<long long,6>> per; std::set<Key> useful; long long miss=0,bad=0;
  for(Long64_t i=0;i<cb->GetEntries();++i){
    cb->GetEntry(i); Key k{ce,root}; auto it=Pm.find(k); if(it==Pm.end()){++miss;continue;}
    int b=it->second.b; if(b<0) continue; if(face<0||face>=6){++bad;continue;}
    ++by[b][face]; ++nc[b]; ++per[k][face]; useful.insert(k);
  }

  std::array<std::array<long long,6>,6> dom{}; std::array<long long,6> zero{},ties{};
  for(auto &kv:Pm){
    int b=kv.second.b; if(b<0) continue; auto it=per.find(kv.first); if(it==per.end()){++zero[b];continue;}
    auto a=it->second; auto mx=*std::max_element(a.begin(),a.end()); int win=-1,nw=0;
    for(int ff=0;ff<6;++ff) if(a[ff]==mx){win=ff;++nw;} if(nw==1) ++dom[b][win]; else ++ties[b];
  }

  std::cout<<std::fixed<<std::setprecision(2);
  std::cout<<"\nCRY CORE-BOUNDARY vs ROOT PRIMARY — V2\n";
  std::cout<<"Join: (EventID,RootPrimaryTrackID) -> (EventID,TrackID)\n";
  std::cout<<"Unique EventIDs              : "<<mult.size()<<"\n";
  std::cout<<"PrimaryTruth rows            : "<<pt->GetEntries()<<"\n";
  std::cout<<"Multi-primary events         : "<<multi<<"\n";
  std::cout<<"Max primary multiplicity     : "<<maxm<<"\n";
  std::cout<<"Primaries with firstExit     : "<<useful.size()<<"\n";
  std::cout<<"Crossings missing primary    : "<<miss<<"\n";
  std::cout<<"Non-downward primaries       : "<<nonDown<<"\n";
  std::cout<<"Bad FaceID                   : "<<bad<<"\n\n";

  std::cout<<"=== FACE OCCUPANCY BY PRIMARY ZENITH ===\n";
  std::cout<<std::setw(9)<<"theta"<<std::setw(11)<<"primaries"<<std::setw(12)<<"crossings";
  for(int ff=0;ff<6;++ff) std::cout<<std::setw(10)<<FN(ff); std::cout<<"\n";
  for(int b=0;b<6;++b){
    std::cout<<std::setw(9)<<ZL(b)<<std::setw(11)<<nprim[b]<<std::setw(12)<<nc[b];
    for(int ff=0;ff<6;++ff){double q=nc[b]?100.0*by[b][ff]/double(nc[b]):0; std::cout<<std::setw(9)<<q<<"%";} std::cout<<"\n";
  }

  std::cout<<"\n=== PRIMARY-LEVEL DOMINANT EXIT FACE ===\n";
  std::cout<<std::setw(9)<<"theta"<<std::setw(11)<<"primaries"<<std::setw(9)<<"zero"<<std::setw(9)<<"ties";
  for(int ff=0;ff<6;++ff) std::cout<<std::setw(9)<<FN(ff); std::cout<<"\n";
  for(int b=0;b<6;++b){
    std::cout<<std::setw(9)<<ZL(b)<<std::setw(11)<<nprim[b]<<std::setw(9)<<zero[b]<<std::setw(9)<<ties[b];
    for(int ff=0;ff<6;++ff) std::cout<<std::setw(9)<<dom[b][ff]; std::cout<<"\n";
  }

  std::cout<<"\n=== MULTI-PRIMARY EVENTS ===\n";
  for(auto &kv:mult) if(kv.second>1) std::cout<<"EventID="<<kv.first<<" multiplicity="<<kv.second<<"\n";

  bool pass=(miss==0&&bad==0&&nonDown==0&&mult.size()==1000);
  std::cout<<"\nCRY ROOT-PRIMARY CORRELATION SANITY: "<<(pass?"PASS":"FAIL")<<"\n";
}
