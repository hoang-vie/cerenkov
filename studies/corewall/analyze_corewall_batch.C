#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace batchcw {

uint64_t Key(int e,int t){
  return (uint64_t(uint32_t(e))<<32) | uint32_t(t);
}

struct TrackInfo {
  int parent=0;
  int pdg=0;
};

struct EventFeat {
  long long ngen=0;
  long long nwall=0;
  long long nmu=0;
  long long ne=0;
  long long nother=0;

  long long nprompt1=0;
  long long nprompt1_mu=0;
  long long nprompt1_e=0;

  long long nprompt2=0;
  long long nprompt2_mu=0;
  long long nprompt2_e=0;

  long long nmuband=0;
  long long nmuband_mu=0;
  long long nmuband_e=0;

  long long face[6]={0,0,0,0,0,0};

  double tsum=0.0;
  double tsum2=0.0;
};

int FaceFromXYZ(double x,double y,double z,
                double hx,double hy,double hz,double tol=1e-3) {
  const double rx=std::abs(std::abs(x)-hx);
  const double ry=std::abs(std::abs(y)-hy);
  const double rz=std::abs(std::abs(z)-hz);
  const bool onX=rx<=tol,onY=ry<=tol,onZ=rz<=tol;
  const int n=int(onX)+int(onY)+int(onZ);
  if(n!=1) return -1;
  if(onX) return x>=0 ? 0 : 1;
  if(onZ) return z>=0 ? 2 : 3;
  return y>=0 ? 4 : 5;
}

double Frac(long long a,long long b){ return b?double(a)/double(b):0.0; }

} // namespace batchcw

void analyze_corewall_batch(
    const char* inputName="validation_corewall_muon_topology_run0.root") {

  using namespace batchcw;

  TFile* f=TFile::Open(inputName,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"Cannot open "<<inputName<<"\n"; return; }

  auto* tTrack=dynamic_cast<TTree*>(f->Get("Track"));
  auto* tPho=dynamic_cast<TTree*>(f->Get("Photon"));
  auto* tCfg=dynamic_cast<TTree*>(f->Get("RunConfig"));
  if(!tTrack||!tPho||!tCfg){
    std::cerr<<"Missing Track/Photon/RunConfig\n";
    f->Close(); return;
  }

  float cx=0,cy=0,cz=0;
  int boundaryDebug=-1,cerenkovEnabled=-1;
  char revision[256]={};
  tCfg->SetBranchAddress("CoreX_mm",&cx);
  tCfg->SetBranchAddress("CoreY_mm",&cy);
  tCfg->SetBranchAddress("CoreZ_mm",&cz);
  tCfg->SetBranchAddress("BoundaryDebug",&boundaryDebug);
  tCfg->SetBranchAddress("CerenkovEnabled",&cerenkovEnabled);
  if(tCfg->GetBranch("ProjectRevision"))
    tCfg->SetBranchAddress("ProjectRevision",revision);
  tCfg->GetEntry(0);

  const double hx=.5*cx,hy=.5*cy,hz=.5*cz;

  std::unordered_map<uint64_t,TrackInfo> tracks;
  int te=0,tid=0,parent=0,pdg=0;
  tTrack->SetBranchAddress("EventID",&te);
  tTrack->SetBranchAddress("TrackID",&tid);
  tTrack->SetBranchAddress("ParentID",&parent);
  tTrack->SetBranchAddress("PDG",&pdg);
  for(Long64_t i=0;i<tTrack->GetEntries();++i){
    tTrack->GetEntry(i);
    tracks[Key(te,tid)]={parent,pdg};
  }

  std::unordered_map<int,EventFeat> ev;

  int pe=0,parentPho=0,reached=0;
  float wt=0,wx=0,wy=0,wz=0;
  tPho->SetBranchAddress("EventID",&pe);
  tPho->SetBranchAddress("ParentTrackID",&parentPho);
  tPho->SetBranchAddress("ReachedCoreWall",&reached);
  tPho->SetBranchAddress("FirstWallTime_ns",&wt);
  tPho->SetBranchAddress("FirstWallX_mm",&wx);
  tPho->SetBranchAddress("FirstWallY_mm",&wy);
  tPho->SetBranchAddress("FirstWallZ_mm",&wz);

  long long nAmb=0,nOff=0;
  for(Long64_t i=0;i<tPho->GetEntries();++i){
    tPho->GetEntry(i);
    auto& e=ev[pe];
    ++e.ngen;
    if(!reached) continue;
    ++e.nwall;

    bool isMu=false,isE=false;
    auto it=tracks.find(Key(pe,parentPho));
    if(it!=tracks.end()){
      isMu=(it->second.parent==0 && std::abs(it->second.pdg)==13);
      isE=(std::abs(it->second.pdg)==11);
    }
    if(isMu) ++e.nmu;
    else if(isE) ++e.ne;
    else ++e.nother;

    if(wt>=0.17000 && wt<0.25000){
      ++e.nprompt1;
      if(isMu) ++e.nprompt1_mu;
      if(isE) ++e.nprompt1_e;
    }
    if(wt>=0.25000 && wt<0.40000){
      ++e.nprompt2;
      if(isMu) ++e.nprompt2_mu;
      if(isE) ++e.nprompt2_e;
    }
    if(wt>=0.52000 && wt<0.69000){
      ++e.nmuband;
      if(isMu) ++e.nmuband_mu;
      if(isE) ++e.nmuband_e;
    }

    const int face=FaceFromXYZ(wx,wy,wz,hx,hy,hz);
    if(face>=0) ++e.face[face];
    else {
      const double rx=std::abs(std::abs(double(wx))-hx);
      const double ry=std::abs(std::abs(double(wy))-hy);
      const double rz=std::abs(std::abs(double(wz))-hz);
      const int n=int(rx<=1e-3)+int(ry<=1e-3)+int(rz<=1e-3);
      if(n>1) ++nAmb; else ++nOff;
    }

    e.tsum += wt;
    e.tsum2 += double(wt)*double(wt);
  }

  std::vector<int> ids;
  ids.reserve(ev.size());
  for(const auto& kv:ev) ids.push_back(kv.first);
  std::sort(ids.begin(),ids.end());

  TString base=inputName;
  auto slash=base.Last('/'); if(slash!=kNPOS) base.Remove(0,slash+1);
  base.ReplaceAll(".root","");

  TString rootOut=Form("corewall_batch_%s.root",base.Data());
  TString csvOut=Form("corewall_batch_%s.csv",base.Data());
  TString pdfOut=Form("corewall_batch_%s.pdf",base.Data());

  TFile fout(rootOut,"RECREATE");
  TTree ft("CoreWallFeatures","Event-level Core-wall features");

  int eventID=0;
  long long Ngen=0,Nwall=0,Nmu=0,Ne=0;
  double fElectron=0,fMuon=0;
  double fPrompt1=0,fPrompt1Electron=0;
  double fPrompt2=0,fPrompt2Electron=0;
  double fMuBand=0,fMuBandMuon=0;
  double fPlusY=0,fMinusY=0;
  double meanWallTime=0,rmsWallTime=0;

  ft.Branch("EventID",&eventID);
  ft.Branch("Ngen",&Ngen);
  ft.Branch("Nwall",&Nwall);
  ft.Branch("Nmu",&Nmu);
  ft.Branch("Ne",&Ne);
  ft.Branch("fElectron",&fElectron);
  ft.Branch("fMuon",&fMuon);
  ft.Branch("fPrompt1",&fPrompt1);
  ft.Branch("fPrompt1Electron",&fPrompt1Electron);
  ft.Branch("fPrompt2",&fPrompt2);
  ft.Branch("fPrompt2Electron",&fPrompt2Electron);
  ft.Branch("fMuBand",&fMuBand);
  ft.Branch("fMuBandMuon",&fMuBandMuon);
  ft.Branch("fPlusY",&fPlusY);
  ft.Branch("fMinusY",&fMinusY);
  ft.Branch("MeanWallTime_ns",&meanWallTime);
  ft.Branch("RMSWallTime_ns",&rmsWallTime);

  std::ofstream csv(csvOut.Data());
  csv<<"EventID,Ngen,Nwall,Nmu,Ne,fElectron,fMuon,fPrompt1,fPrompt1Electron,"
        "fPrompt2,fPrompt2Electron,fMuBand,fMuBandMuon,fPlusY,fMinusY,"
        "MeanWallTime_ns,RMSWallTime_ns\n";

  TH1D hNwall("hNwall","Core-wall photons/event;N_{wall};Events",60,7000,19000);
  TH1D hFe("hFe","Electron-emitter fraction;N_{e}/N_{wall};Events",50,0,0.6);
  TH1D hFp1("hFp1","Prompt-1 fraction;N(0.17-0.25 ns)/N_{wall};Events",50,0,0.12);
  TH1D hFp1e("hFp1e","Electron purity in prompt-1;N_{e,prompt1}/N_{prompt1};Events",50,0,1.02);
  TH1D hFmu("hFmu","Muon purity in main band;N_{#mu,band}/N_{band};Events",50,0.7,1.02);
  TH1D hTm("hTm","Mean Core-wall time;Mean t_{wall} (ns);Events",50,0.45,0.75);

  for(int id:ids){
    const auto& e=ev[id];
    eventID=id; Ngen=e.ngen; Nwall=e.nwall; Nmu=e.nmu; Ne=e.ne;
    fElectron=Frac(e.ne,e.nwall);
    fMuon=Frac(e.nmu,e.nwall);
    fPrompt1=Frac(e.nprompt1,e.nwall);
    fPrompt1Electron=Frac(e.nprompt1_e,e.nprompt1);
    fPrompt2=Frac(e.nprompt2,e.nwall);
    fPrompt2Electron=Frac(e.nprompt2_e,e.nprompt2);
    fMuBand=Frac(e.nmuband,e.nwall);
    fMuBandMuon=Frac(e.nmuband_mu,e.nmuband);
    fPlusY=Frac(e.face[4],e.nwall);
    fMinusY=Frac(e.face[5],e.nwall);
    meanWallTime=e.nwall?e.tsum/e.nwall:0;
    const double var=e.nwall?e.tsum2/e.nwall-meanWallTime*meanWallTime:0;
    rmsWallTime=std::sqrt(std::max(0.0,var));

    ft.Fill();
    csv<<eventID<<","<<Ngen<<","<<Nwall<<","<<Nmu<<","<<Ne<<","
       <<fElectron<<","<<fMuon<<","<<fPrompt1<<","<<fPrompt1Electron<<","
       <<fPrompt2<<","<<fPrompt2Electron<<","<<fMuBand<<","<<fMuBandMuon<<","
       <<fPlusY<<","<<fMinusY<<","<<meanWallTime<<","<<rmsWallTime<<"\n";

    hNwall.Fill(Nwall);
    hFe.Fill(fElectron);
    hFp1.Fill(fPrompt1);
    if(e.nprompt1>0) hFp1e.Fill(fPrompt1Electron);
    if(e.nmuband>0) hFmu.Fill(fMuBandMuon);
    hTm.Fill(meanWallTime);
  }

  ft.Write();
  hNwall.Write();hFe.Write();hFp1.Write();hFp1e.Write();hFmu.Write();hTm.Write();
  fout.Close();
  csv.close();

  TCanvas c("cBatch","Core-wall batch features",1100,750);
  c.Divide(2,3);
  c.cd(1);hNwall.SetStats(0);hNwall.Draw();
  c.cd(2);hFe.SetStats(0);hFe.Draw();
  c.cd(3);hFp1.SetStats(0);hFp1.Draw();
  c.cd(4);hFp1e.SetStats(0);hFp1e.Draw();
  c.cd(5);hFmu.SetStats(0);hFmu.Draw();
  c.cd(6);hTm.SetStats(0);hTm.Draw();
  c.SaveAs(pdfOut);

  std::cout<<"\n=== CORE-WALL BATCH SUMMARY ===\n";
  std::cout<<"Input: "<<inputName<<"\n";
  std::cout<<"ProjectRevision: "<<(revision[0]?revision:"?")<<"\n";
  std::cout<<"Events="<<ids.size()
           <<" CerenkovEnabled="<<cerenkovEnabled
           <<" BoundaryDebug="<<boundaryDebug<<"\n";
  std::cout<<"Geometry mapping QA: ambiguous="<<nAmb<<" offSurface="<<nOff<<"\n";
  std::cout<<"<Nwall>="<<hNwall.GetMean()<<" RMS="<<hNwall.GetRMS()<<"\n";
  std::cout<<"<electron fraction>="<<hFe.GetMean()<<" RMS="<<hFe.GetRMS()<<"\n";
  std::cout<<"<prompt1 fraction>="<<hFp1.GetMean()<<" RMS="<<hFp1.GetRMS()<<"\n";
  std::cout<<"<prompt1 electron purity>="<<hFp1e.GetMean()<<" RMS="<<hFp1e.GetRMS()<<"\n";
  std::cout<<"<main-band muon purity>="<<hFmu.GetMean()<<" RMS="<<hFmu.GetRMS()<<"\n";
  std::cout<<"<mean wall time>="<<hTm.GetMean()<<" ns RMS(event means)="<<hTm.GetRMS()<<" ns\n";
  std::cout<<"Wrote "<<rootOut<<"\n";
  std::cout<<"Wrote "<<csvOut<<"\n";
  std::cout<<"Wrote "<<pdfOut<<"\n";
  std::cout<<"CORE-WALL BATCH: DONE\n";

  f->Close();
}