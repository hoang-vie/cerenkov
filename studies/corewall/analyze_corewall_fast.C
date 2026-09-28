#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TPad.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace fastcw {

uint64_t Key(int e,int t){
  return (uint64_t(uint32_t(e))<<32) | uint32_t(t);
}

struct TrackInfo{ int parent=0,pdg=0; };

struct Geo{
  double hx=50,hy=50,hz=50;
};

struct Hit{
  int face=-1; double u=0,v=0; bool ambiguous=false,off=false;
};

Hit Classify(double x,double y,double z,const Geo& g,double tol=1e-3){
  Hit h;
  const double rx=std::abs(std::abs(x)-g.hx);
  const double ry=std::abs(std::abs(y)-g.hy);
  const double rz=std::abs(std::abs(z)-g.hz);
  const bool onX=rx<=tol,onY=ry<=tol,onZ=rz<=tol;
  const int n=int(onX)+int(onY)+int(onZ);
  if(n==0){h.off=true;return h;}
  if(n>1){h.ambiguous=true;return h;}
  if(onX){
    if(x>=0){h.face=0;h.u=y;h.v=z;}
    else    {h.face=1;h.u=-y;h.v=z;}
  } else if(onZ){
    if(z>=0){h.face=2;h.u=x;h.v=y;}
    else    {h.face=3;h.u=x;h.v=-y;}
  } else {
    if(y>=0){h.face=4;h.u=x;h.v=-z;}
    else    {h.face=5;h.u=x;h.v=z;}
  }
  return h;
}

const char* FaceName(int f){
  static const char* n[6]={"+X","-X","+Z","-Z","+Y","-Y"};
  return (f>=0&&f<6)?n[f]:"?";
}

struct TimeWindow{
  const char* name; double lo,hi;
};

enum Cat{ALL=0,MU=1,ELECTRON=2,NCAT=3};

struct Set{
  std::unique_ptr<TH2D> h[6];
  long long face[6]={0,0,0,0,0,0};
  long long n=0,amb=0,off=0;
};

void FaceRange(int f,const Geo& g,double& u0,double& u1,double& v0,double& v1){
  if(f==0||f==1){u0=-g.hy;u1=g.hy;v0=-g.hz;v1=g.hz;}
  else if(f==2||f==3){u0=-g.hx;u1=g.hx;v0=-g.hy;v1=g.hy;}
  else {u0=-g.hx;u1=g.hx;v0=-g.hz;v1=g.hz;}
}

void DrawSix(TCanvas& c,Set& s,const TString& title){
  c.Clear(); c.SetFillColor(0);
  const int slot[6]={4,2,1,5,0,3};
  const double xy[6][4]={
    {0.39,0.73,0.61,0.95},{0.39,0.51,0.61,0.73},
    {0.17,0.29,0.39,0.51},{0.39,0.29,0.61,0.51},
    {0.61,0.29,0.83,0.51},{0.39,0.07,0.61,0.29}
  };
  double mx=0;
  for(int f=0;f<6;++f) mx=std::max(mx,s.h[f]->GetMaximum());
  for(int i=0;i<6;++i){
    auto* p=new TPad(Form("p%d",i),"",xy[i][0],xy[i][1],xy[i][2],xy[i][3]);
    p->SetLeftMargin(.18);p->SetRightMargin(.05);p->SetBottomMargin(.18);p->SetTopMargin(.15);
    p->Draw();p->cd();
    auto* h=s.h[slot[i]].get();
    h->SetStats(0);h->SetMinimum(0);if(mx>0)h->SetMaximum(mx);
    h->Draw("COL");
    c.cd();
  }
  TLatex tx;tx.SetNDC();tx.SetTextAlign(22);tx.SetTextSize(.032);
  tx.DrawLatex(.5,.985,title);
}

} // namespace fastcw

void analyze_corewall_fast(
    const char* inputName="validation_corewall_muon_topology_run0.root",
    int eventID=71){

  using namespace fastcw;
  TFile* f=TFile::Open(inputName,"READ");
  if(!f||f->IsZombie()){std::cerr<<"Cannot open "<<inputName<<"\n";return;}

  auto* tPho=dynamic_cast<TTree*>(f->Get("Photon"));
  auto* tTrack=dynamic_cast<TTree*>(f->Get("Track"));
  auto* tCfg=dynamic_cast<TTree*>(f->Get("RunConfig"));
  if(!tPho||!tTrack||!tCfg){std::cerr<<"Missing Photon/Track/RunConfig\n";f->Close();return;}

  Geo g;
  float cx=0,cy=0,cz=0;
  tCfg->SetBranchAddress("CoreX_mm",&cx);
  tCfg->SetBranchAddress("CoreY_mm",&cy);
  tCfg->SetBranchAddress("CoreZ_mm",&cz);
  tCfg->GetEntry(0);
  g.hx=.5*cx;g.hy=.5*cy;g.hz=.5*cz;

  std::unordered_map<uint64_t,TrackInfo> tracks;
  int te=0,tid=0,parent=0,pdg=0;
  tTrack->SetBranchAddress("EventID",&te);
  tTrack->SetBranchAddress("TrackID",&tid);
  tTrack->SetBranchAddress("ParentID",&parent);
  tTrack->SetBranchAddress("PDG",&pdg);
  for(Long64_t i=0;i<tTrack->GetEntries();++i){
    tTrack->GetEntry(i);
    if(te!=eventID) continue;
    tracks[Key(te,tid)]={parent,pdg};
  }

  const std::vector<fastcw::TimeWindow> w={
    {"prompt_e1",0.17000,0.25000},
    {"prompt_e2",0.25000,0.40000},
    {"main_early",0.51197,0.57013},
    {"main_middle",0.57013,0.63725},
    {"main_late",0.63725,0.69379},
    {"mu_band",0.52000,0.69000}
  };
  const char* catName[NCAT]={"all","mu_primary","electron"};

  std::vector<std::vector<Set>> sets;
  sets.resize(w.size());
  for (auto& row : sets) row.resize(NCAT);

  for(size_t iw=0;iw<w.size();++iw){
    for(int c=0;c<NCAT;++c){
      for(int face=0;face<6;++face){
        double u0,u1,v0,v1;FaceRange(face,g,u0,u1,v0,v1);
        sets[iw][c].h[face]=std::make_unique<TH2D>(
          Form("h_%s_%s_f%d",w[iw].name,catName[c],face),
          Form("%s %s %s;Local U (mm);Local V (mm)",w[iw].name,catName[c],FaceName(face)),
          80,u0,u1,80,v0,v1);
        // ROOT histograms are automatically attached to the current TDirectory.
        // We also own them with std::unique_ptr, so detach them to avoid the input
        // TFile/TDirectory deleting the same histogram at shutdown (double free).
        sets[iw][c].h[face]->SetDirectory(nullptr);
      }
    }
  }

  int pe=0,parentPho=0,reached=0;
  float wt=0,wx=0,wy=0,wz=0;
  tPho->SetBranchAddress("EventID",&pe);
  tPho->SetBranchAddress("ParentTrackID",&parentPho);
  tPho->SetBranchAddress("ReachedCoreWall",&reached);
  tPho->SetBranchAddress("FirstWallTime_ns",&wt);
  tPho->SetBranchAddress("FirstWallX_mm",&wx);
  tPho->SetBranchAddress("FirstWallY_mm",&wy);
  tPho->SetBranchAddress("FirstWallZ_mm",&wz);

  long long scanned=0,wall=0;
  for(Long64_t i=0;i<tPho->GetEntries();++i){
    tPho->GetEntry(i);
    if(pe!=eventID) continue;
    ++scanned;
    if(!reached) continue;
    ++wall;

    Cat cat=ALL;
    auto it=tracks.find(Key(pe,parentPho));
    bool isMu=false,isE=false;
    if(it!=tracks.end()){
      isMu=(it->second.parent==0 && std::abs(it->second.pdg)==13);
      isE=(std::abs(it->second.pdg)==11);
    }

    Hit hit=Classify(wx,wy,wz,g);
    for(size_t iw=0;iw<w.size();++iw){
      if(!(wt>=w[iw].lo && wt<w[iw].hi)) continue;

      auto fill=[&](Cat c){
        auto& s=sets[iw][c];
        ++s.n;
        if(hit.off){++s.off;return;}
        if(hit.ambiguous){++s.amb;return;}
        ++s.face[hit.face];
        s.h[hit.face]->Fill(hit.u,hit.v);
      };
      fill(ALL);
      if(isMu) fill(MU);
      if(isE) fill(ELECTRON);
    }
  }

  std::cout<<"\n=== FAST CORE-WALL ANALYSIS ===\n";
  std::cout<<"Input="<<inputName<<" EventID="<<eventID
           <<" PhotonRowsForEvent="<<scanned<<" ReachedWall="<<wall<<"\n\n";
  std::cout<<std::left<<std::setw(14)<<"window"
           <<std::right<<std::setw(8)<<"all"<<std::setw(8)<<"mu"<<std::setw(8)<<"e"
           <<std::setw(10)<<"e/all%"<<std::setw(10)<<"mu/all%"<<"\n";
  for(size_t iw=0;iw<w.size();++iw){
    auto& a=sets[iw][ALL];auto& m=sets[iw][MU];auto& e=sets[iw][ELECTRON];
    const double ef=a.n?100.0*e.n/a.n:0;
    const double mf=a.n?100.0*m.n/a.n:0;
    std::cout<<std::left<<std::setw(14)<<w[iw].name
             <<std::right<<std::setw(8)<<a.n<<std::setw(8)<<m.n<<std::setw(8)<<e.n
             <<std::setw(10)<<std::fixed<<std::setprecision(2)<<ef
             <<std::setw(10)<<mf<<"\n";
  }

  TString base=inputName; auto slash=base.Last('/'); if(slash!=kNPOS) base.Remove(0,slash+1);
  base.ReplaceAll(".root","");
  TString pdf=Form("corewall_fast_%s_evt%d.pdf",base.Data(),eventID);
  TCanvas c("cFast","Fast Core-wall analysis",1400,900);
  c.Print(pdf+"[");
  for(size_t iw=0;iw<w.size();++iw){
    for(int cat=0;cat<NCAT;++cat){
      DrawSix(c,sets[iw][cat],
              Form("Event %d | %s | %s | [%.5f, %.5f) ns",
                   eventID,w[iw].name,catName[cat],w[iw].lo,w[iw].hi));
      c.Print(pdf);
    }
  }
  c.Print(pdf+"]");
  std::cout<<"Wrote "<<pdf<<"\n";
  std::cout<<"FAST CORE-WALL ANALYSIS: DONE\n";
  f->Close();
}