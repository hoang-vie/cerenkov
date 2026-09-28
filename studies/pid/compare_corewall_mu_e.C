#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace pidcw {

struct Event {
  long long nwall=0;
  long long prompt=0;
  long long mainband=0;
  long long face[6]={0,0,0,0,0,0};
  double tsum=0,tsum2=0;
};

int Face(double x,double y,double z,double hx,double hy,double hz,double tol=1e-3){
  const bool X=std::abs(std::abs(x)-hx)<=tol;
  const bool Y=std::abs(std::abs(y)-hy)<=tol;
  const bool Z=std::abs(std::abs(z)-hz)<=tol;
  if(int(X)+int(Y)+int(Z)!=1) return -1;
  if(X) return x>=0?0:1;
  if(Z) return z>=0?2:3;
  return y>=0?4:5;
}

struct Feat{
  int event=0,label=0;
  double Nwall=0,fPrompt=0,fMain=0,fPlusY=0,fMinusY=0,meanT=0,rmsT=0;
};

std::vector<Feat> Load(const char* name,int label){
  TFile* f=TFile::Open(name,"READ");
  if(!f||f->IsZombie()) return {};
  auto* p=dynamic_cast<TTree*>(f->Get("Photon"));
  auto* cfg=dynamic_cast<TTree*>(f->Get("RunConfig"));
  if(!p||!cfg){f->Close();return{};}

  float cx=0,cy=0,cz=0;
  cfg->SetBranchAddress("CoreX_mm",&cx);
  cfg->SetBranchAddress("CoreY_mm",&cy);
  cfg->SetBranchAddress("CoreZ_mm",&cz);
  cfg->GetEntry(0);
  const double hx=.5*cx,hy=.5*cy,hz=.5*cz;

  std::map<int,Event> ev;
  int eid=0,reached=0; float t=0,x=0,y=0,z=0;
  p->SetBranchAddress("EventID",&eid);
  p->SetBranchAddress("ReachedCoreWall",&reached);
  p->SetBranchAddress("FirstWallTime_ns",&t);
  p->SetBranchAddress("FirstWallX_mm",&x);
  p->SetBranchAddress("FirstWallY_mm",&y);
  p->SetBranchAddress("FirstWallZ_mm",&z);

  for(Long64_t i=0;i<p->GetEntries();++i){
    p->GetEntry(i);
    if(!reached) continue;
    auto& e=ev[eid];
    ++e.nwall;
    // Fixed windows are valid only because both controlled samples use the same
    // source position, global t=0 and geometry.
    if(t>=0.17 && t<0.40) ++e.prompt;
    if(t>=0.52 && t<0.69) ++e.mainband;
    int face=Face(x,y,z,hx,hy,hz);
    if(face>=0) ++e.face[face];
    e.tsum+=t; e.tsum2+=double(t)*double(t);
  }

  std::vector<Feat> out;
  for(auto& kv:ev){
    const auto& e=kv.second;
    if(!e.nwall) continue;
    Feat q; q.event=kv.first;q.label=label;q.Nwall=e.nwall;
    q.fPrompt=double(e.prompt)/e.nwall;
    q.fMain=double(e.mainband)/e.nwall;
    q.fPlusY=double(e.face[4])/e.nwall;
    q.fMinusY=double(e.face[5])/e.nwall;
    q.meanT=e.tsum/e.nwall;
    q.rmsT=std::sqrt(std::max(0.0,e.tsum2/e.nwall-q.meanT*q.meanT));
    out.push_back(q);
  }
  f->Close();
  return out;
}

double AUC(const std::vector<double>& mu,const std::vector<double>& el){
  if(mu.empty()||el.empty()) return 0.5;
  long long gt=0,eq=0;
  for(double a:mu) for(double b:el){
    if(a>b) ++gt;
    else if(a==b) ++eq;
  }
  double auc=(gt+0.5*eq)/double(mu.size()*el.size());
  return std::max(auc,1.0-auc); // separation, independent of sign convention
}

void PrintFeature(const char* name,const std::vector<double>& mu,const std::vector<double>& el){
  auto mr=[](const std::vector<double>& v){
    double s=0,s2=0;for(double x:v){s+=x;s2+=x*x;}
    double m=s/v.size(),r=std::sqrt(std::max(0.0,s2/v.size()-m*m));
    return std::pair<double,double>(m,r);
  };
  auto m=mr(mu),e=mr(el);
  std::cout<<std::left<<std::setw(13)<<name
           <<" mu="<<std::setw(11)<<std::fixed<<std::setprecision(5)<<m.first
           <<" +/- "<<std::setw(10)<<m.second
           <<" e="<<std::setw(11)<<e.first
           <<" +/- "<<std::setw(10)<<e.second
           <<" AUCsep="<<std::setprecision(4)<<AUC(mu,el)<<"\n";
}

} // namespace pidcw

void compare_corewall_mu_e(
    const char* muFile="pid_muon_3GeV_run0.root",
    const char* eFile="pid_electron_3GeV_run0.root"){

  using namespace pidcw;
  auto mu=Load(muFile,0), el=Load(eFile,1);
  if(mu.empty()||el.empty()){
    std::cerr<<"Could not load both samples\n";return;
  }

  auto collect=[](const std::vector<Feat>& v,double Feat::*m){
    std::vector<double> o;o.reserve(v.size());
    for(const auto& x:v)o.push_back(x.*m);
    return o;
  };

  std::cout<<"\n=== CONTROLLED CORE-WALL MUON vs ELECTRON ===\n";
  std::cout<<"Definition: same 3 GeV kinetic energy, same centered -Y trajectory/source time.\n";
  std::cout<<"Muon events="<<mu.size()<<" Electron events="<<el.size()<<"\n";
  PrintFeature("Nwall",collect(mu,&Feat::Nwall),collect(el,&Feat::Nwall));
  PrintFeature("fPrompt",collect(mu,&Feat::fPrompt),collect(el,&Feat::fPrompt));
  PrintFeature("fMain",collect(mu,&Feat::fMain),collect(el,&Feat::fMain));
  PrintFeature("fPlusY",collect(mu,&Feat::fPlusY),collect(el,&Feat::fPlusY));
  PrintFeature("fMinusY",collect(mu,&Feat::fMinusY),collect(el,&Feat::fMinusY));
  PrintFeature("meanT",collect(mu,&Feat::meanT),collect(el,&Feat::meanT));
  PrintFeature("rmsT",collect(mu,&Feat::rmsT),collect(el,&Feat::rmsT));

  TH1D hMuN("hMuN","Nwall;N_{wall};Events",60,0,25000);
  TH1D hEN("hEN","Nwall;N_{wall};Events",60,0,25000);
  TH1D hMuR("hMuR","Wall timing RMS;RMS(t_{wall}) (ns);Events",60,0,0.5);
  TH1D hER("hER","Wall timing RMS;RMS(t_{wall}) (ns);Events",60,0,0.5);
  TH1D hMuY("hMuY","-Y face fraction;f_{-Y};Events",60,0,1);
  TH1D hEY("hEY","-Y face fraction;f_{-Y};Events",60,0,1);
  for(const auto& x:mu){hMuN.Fill(x.Nwall);hMuR.Fill(x.rmsT);hMuY.Fill(x.fMinusY);}
  for(const auto& x:el){hEN.Fill(x.Nwall);hER.Fill(x.rmsT);hEY.Fill(x.fMinusY);}

  TCanvas c("cPid","Controlled mu/e comparison",1200,400);c.Divide(3,1);
  auto draw=[&](int pad,TH1D& a,TH1D& b){
    c.cd(pad);a.SetStats(0);b.SetStats(0);a.SetLineWidth(2);b.SetLineWidth(2);
    double mx=std::max(a.GetMaximum(),b.GetMaximum());a.SetMaximum(1.15*mx);
    a.Draw("HIST");b.Draw("HIST SAME");
    auto* leg=new TLegend(.65,.75,.88,.88);leg->AddEntry(&a,"#mu^{-}","l");leg->AddEntry(&b,"e^{-}","l");leg->Draw();
  };
  draw(1,hMuN,hEN);draw(2,hMuR,hER);draw(3,hMuY,hEY);
  c.SaveAs("corewall_pid_mu_e_3GeV.pdf");

  std::cout<<"Wrote corewall_pid_mu_e_3GeV.pdf\n";
  std::cout<<"CONTROLLED MU/E COMPARISON: DONE\n";
}
