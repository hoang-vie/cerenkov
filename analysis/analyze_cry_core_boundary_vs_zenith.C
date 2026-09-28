#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>

namespace {
constexpr double kPi = 3.14159265358979323846;
const char* FaceName(int f){
  switch(f){case 0:return "+X";case 1:return "-X";case 2:return "+Z";
            case 3:return "-Z";case 4:return "+Y";case 5:return "-Y";
            default:return "?";}
}
double deg(double r){ return r*180.0/kPi; }
double clamp1(double x){ return std::max(-1.0,std::min(1.0,x)); }
int zenithBin(double t){
  const double e[]={0,15,30,45,60,75,90.000001};
  for(int i=0;i<6;++i) if(t>=e[i] && t<e[i+1]) return i;
  return -1;
}
const char* binLabel(int b){
  static const char* x[]={"0-15","15-30","30-45","45-60","60-75","75-90"};
  return (b>=0&&b<6)?x[b]:"?";
}
}

void analyze_cry_core_boundary_vs_zenith(
    const char* inputName="study_cry_core_boundary_100_run0.root")
{
  TFile f(inputName,"READ");
  if(f.IsZombie()){ std::cerr<<"Cannot open "<<inputName<<"\n"; return; }

  auto* pt=dynamic_cast<TTree*>(f.Get("PrimaryTruth"));
  auto* cb=dynamic_cast<TTree*>(f.Get("CoreBoundaryCrossing"));
  if(!pt||!cb){ std::cerr<<"Missing PrimaryTruth or CoreBoundaryCrossing\n"; return; }

  int pEvent=0,pPDG=0;
  float pKE=0,pDX=0,pDY=0,pDZ=0;
  pt->SetBranchAddress("EventID",&pEvent);
  pt->SetBranchAddress("PDG",&pPDG);
  pt->SetBranchAddress("VertexKE_MeV",&pKE);
  pt->SetBranchAddress("DirX",&pDX);
  pt->SetBranchAddress("DirY",&pDY);
  pt->SetBranchAddress("DirZ",&pDZ);

  struct P{int pdg; double ke,dx,dy,dz,theta,phi; int bin;};
  std::map<int,P> prim;
  std::array<long long,6> eventsPerBin{};
  long long nonDownward=0;

  for(Long64_t i=0;i<pt->GetEntries();++i){
    pt->GetEntry(i);
    double theta=deg(std::acos(clamp1(-double(pDY))));   // 0 = straight down (-Y)
    double phi=deg(std::atan2(double(pDZ),double(pDX))); // 0=+X, +90=+Z
    int b=zenithBin(theta);
    prim[pEvent]={pPDG,pKE,pDX,pDY,pDZ,theta,phi,b};
    if(b>=0) ++eventsPerBin[b];
    if(pDY>=0) ++nonDownward;
  }

  int cEvent=0,cFace=0;
  cb->SetBranchAddress("EventID",&cEvent);
  cb->SetBranchAddress("FaceID",&cFace);

  std::array<std::array<long long,6>,6> faceByBin{};
  std::array<long long,6> crossingsPerBin{};
  std::map<int,std::array<long long,6>> eventFaces;
  long long missingPrimary=0,badFace=0;

  for(Long64_t i=0;i<cb->GetEntries();++i){
    cb->GetEntry(i);
    auto it=prim.find(cEvent);
    if(it==prim.end()){ ++missingPrimary; continue; }
    int b=it->second.bin;
    if(b<0||b>=6) continue;
    if(cFace<0||cFace>=6){ ++badFace; continue; }
    ++faceByBin[b][cFace];
    ++crossingsPerBin[b];
    ++eventFaces[cEvent][cFace];
  }

  std::array<std::array<long long,6>,6> dominantByBin{};
  std::array<long long,6> zeroCross{},ties{};

  for(const auto& kv:prim){
    int ev=kv.first,b=kv.second.bin;
    if(b<0||b>=6) continue;
    auto it=eventFaces.find(ev);
    if(it==eventFaces.end()){ ++zeroCross[b]; continue; }
    auto a=it->second;
    long long mx=*std::max_element(a.begin(),a.end());
    if(mx==0){ ++zeroCross[b]; continue; }
    int winner=-1,n=0;
    for(int face=0;face<6;++face) if(a[face]==mx){winner=face;++n;}
    if(n==1) ++dominantByBin[b][winner];
    else ++ties[b];
  }

  std::cout<<std::fixed<<std::setprecision(2);
  std::cout<<"\nCRY CORE-BOUNDARY vs PRIMARY ZENITH\n";
  std::cout<<"Input: "<<inputName<<"\n";
  std::cout<<"theta = acos(-DirY): 0 deg is straight downward (-Y)\n";
  std::cout<<"phi = atan2(DirZ,DirX): 0 deg=+X, +90 deg=+Z\n\n";
  std::cout<<"PrimaryTruth entries              : "<<pt->GetEntries()<<"\n";
  std::cout<<"CoreBoundaryCrossing entries      : "<<cb->GetEntries()<<"\n";
  std::cout<<"Non-downward primaries (DirY>=0) : "<<nonDownward<<"\n";
  std::cout<<"Crossings missing primary         : "<<missingPrimary<<"\n";
  std::cout<<"Bad FaceID                        : "<<badFace<<"\n\n";

  std::cout<<"=== FACE OCCUPANCY BY ZENITH BIN ===\n";
  std::cout<<std::setw(9)<<"theta"<<std::setw(9)<<"events"<<std::setw(12)<<"crossings";
  for(int face=0;face<6;++face) std::cout<<std::setw(10)<<FaceName(face);
  std::cout<<"\n";
  for(int b=0;b<6;++b){
    std::cout<<std::setw(9)<<binLabel(b)<<std::setw(9)<<eventsPerBin[b]
             <<std::setw(12)<<crossingsPerBin[b];
    for(int face=0;face<6;++face){
      double p=crossingsPerBin[b]?100.0*faceByBin[b][face]/double(crossingsPerBin[b]):0.0;
      std::cout<<std::setw(9)<<p<<"%";
    }
    std::cout<<"\n";
  }

  std::cout<<"\n=== EVENT-LEVEL DOMINANT EXIT FACE ===\n";
  std::cout<<std::setw(9)<<"theta"<<std::setw(9)<<"events"<<std::setw(9)<<"zero"
           <<std::setw(9)<<"ties";
  for(int face=0;face<6;++face) std::cout<<std::setw(9)<<FaceName(face);
  std::cout<<"\n";
  for(int b=0;b<6;++b){
    std::cout<<std::setw(9)<<binLabel(b)<<std::setw(9)<<eventsPerBin[b]
             <<std::setw(9)<<zeroCross[b]<<std::setw(9)<<ties[b];
    for(int face=0;face<6;++face) std::cout<<std::setw(9)<<dominantByBin[b][face];
    std::cout<<"\n";
  }

  std::cout<<"\n=== FIRST 15 PRIMARIES ===\n";
  int shown=0;
  for(const auto& kv:prim){
    if(shown++>=15) break;
    const auto& p=kv.second;
    std::cout<<"event="<<kv.first<<" PDG="<<p.pdg<<" KE="<<p.ke<<" MeV"
             <<" dir=("<<p.dx<<","<<p.dy<<","<<p.dz<<")"
             <<" theta="<<p.theta<<" deg phi="<<p.phi<<" deg\n";
  }

  bool pass=(missingPrimary==0 && badFace==0 && nonDownward==0);
  std::cout<<"\nCRY ZENITH CORRELATION SANITY: "<<(pass?"PASS":"FAIL")<<"\n";
}
