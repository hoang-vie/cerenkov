#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLine.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <utility>

void validate_muon_decay(const char* inputName="validation_muon_decay_run0.root") {
  TFile f(inputName,"READ");
  if (f.IsZombie()) { std::cerr << "Cannot open file\n"; return; }
  auto* tr=dynamic_cast<TTree*>(f.Get("Track"));
  if(!tr){std::cerr<<"Missing Track\n";return;}

  int event=0,track=0,pdg=0,parent=0; float properNs=0,vertexKE=0; char endProc[256] = {};
  tr->SetBranchAddress("EventID",&event); tr->SetBranchAddress("TrackID",&track);
  tr->SetBranchAddress("PDG",&pdg); tr->SetBranchAddress("ParentID",&parent);
  tr->SetBranchAddress("VertexKE_MeV",&vertexKE);
  tr->SetBranchAddress("EndProperTime_ns",&properNs); tr->SetBranchAddress("EndProcessName",endProc);

  TH1D hLife("hMuPlusProper","#mu^{+} proper decay time;Proper time (#mus);Decays",160,0,12);
  TH1D hMichel("hMichelPositron","Michel e^{+} kinetic-energy spectrum;Kinetic energy (MeV);e^{+}",160,0,55);
  std::set<std::pair<int,int>> primaryMuPlus;
  for(Long64_t i=0;i<tr->GetEntries();++i){
    tr->GetEntry(i);
    if(parent==0 && pdg==-13){
      primaryMuPlus.insert({event,track});
      if(std::string(endProc).find("Decay")!=std::string::npos && std::isfinite(properNs) && properNs>=0)
        hLife.Fill(properNs/1000.0);
    }
  }
  double maxMichelKE=0;
  for(Long64_t i=0;i<tr->GetEntries();++i){
    tr->GetEntry(i);
    if(pdg==-11 && primaryMuPlus.count({event,parent})){
      hMichel.Fill(vertexKE);
      if(vertexKE>maxMichelKE) maxMichelKE=vertexKE;
    }
  }

  if(hLife.GetEntries()<100){std::cerr<<"Too few primary mu+ decays: "<<hLife.GetEntries()<<"\n";return;}
  TF1 fit("fit","[0]*exp(-x/[1])",0.1,10.0);
  fit.SetParameters(hLife.GetMaximum(),2.197);
  hLife.Fit(&fit,"Q0","",0.1,10.0);
  const double tau=fit.GetParameter(1), err=fit.GetParError(1);
  const double ref=2.1969811;
  const double mMu=105.6583755, mE=0.51099895;
  const double michelTotal=(mMu*mMu+mE*mE)/(2*mMu);
  const double michelKE=michelTotal-mE;
  const bool lifetimePass = err>0 && std::abs(tau-ref) < 5.0*err;
  const bool michelPass = hMichel.GetEntries()>1000 && maxMichelKE>51.5 && maxMichelKE<53.0;

  std::cout<<"Primary mu+ decays="<<hLife.GetEntries()<<"\n"
           <<"Fitted tau="<<tau<<" +/- "<<err<<" us\n"
           <<"Reference free-muon lifetime="<<ref<<" us\n"
           <<"Difference="<<(tau-ref)<<" us";
  if(err>0) std::cout<<"  pull="<<(tau-ref)/err;
  std::cout<<"\nMichel positrons="<<hMichel.GetEntries()
           <<" observed max KE="<<maxMichelKE<<" MeV"
           <<" kinematic endpoint KE~"<<michelKE<<" MeV\n"
           <<"MUON DECAY VALIDATION: "<<((lifetimePass&&michelPass)?"PASS":"FAIL")<<"\n";

  gSystem->mkdir("validation_plots", kTRUE);
  TCanvas c("cMuonDecayValidation","Muon decay validation",1200,500);
  c.Divide(2,1);
  c.cd(1); hLife.Draw(); fit.Draw("SAME");
  c.cd(2); hMichel.Draw(); auto* endpoint=new TLine(michelKE,0,michelKE,hMichel.GetMaximum()*1.05); endpoint->SetLineStyle(2); endpoint->Draw("SAME");
  c.SaveAs("validation_plots/muon_decay_validation.png");
  c.SaveAs("validation_plots/muon_decay_validation.pdf");
  std::cout << "Wrote validation_plots/muon_decay_validation.{png,pdf}\n";
}
