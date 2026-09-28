#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <utility>

void validate_dedx(const char* inputName="validation_dedx_muon_run0.root", int requestedPDG=13, bool combineCharge=true) {
  TFile f(inputName,"READ");
  if (f.IsZombie()) { std::cerr << "Cannot open file\n"; return; }
  auto* tr=dynamic_cast<TTree*>(f.Get("Track"));
  auto* tv=dynamic_cast<TTree*>(f.Get("TrackVolume"));
  if (!tr || !tv) { std::cerr << "Missing Track/TrackVolume\n"; return; }

  int e=0,t=0,pdg=0,parent=0;
  std::map<std::pair<int,int>,std::pair<int,int>> truth;
  tr->SetBranchAddress("EventID",&e); tr->SetBranchAddress("TrackID",&t);
  tr->SetBranchAddress("PDG",&pdg); tr->SetBranchAddress("ParentID",&parent);
  for(Long64_t i=0;i<tr->GetEntries();++i){tr->GetEntry(i);truth[{e,t}]={pdg,parent};}

  int ve=0,vt=0,vol=0; float edep=0,dke=0,len=0;
  tv->SetBranchAddress("EventID",&ve); tv->SetBranchAddress("TrackID",&vt);
  tv->SetBranchAddress("VolumeID",&vol); tv->SetBranchAddress("Edep_MeV",&edep);
  tv->SetBranchAddress("DeltaKE_MeV",&dke); tv->SetBranchAddress("Length_mm",&len);

  TH1::StatOverflows(kTRUE);
  TH1D hPath("hDedxPath","Primary-muon path in Core;Path length (cm);Tracks",160,9.95,10.05);
  TH1D hEdep("hDedxEdep","Local deposited energy per path;E_{dep}/L (MeV/cm);Tracks",160,0,4);
  TH1D hDke("hDedxDke","Primary kinetic-energy loss per path;#DeltaKE/L (MeV/cm);Tracks",240,0,30);

  long long n=0; double sL=0,sL2=0,sE=0,sE2=0,sK=0,sK2=0;
  for(Long64_t i=0;i<tv->GetEntries();++i){
    tv->GetEntry(i); if(vol!=0 || len<=0) continue;
    auto it=truth.find({ve,vt}); if(it==truth.end() || it->second.second!=0) continue;
    const int p=it->second.first;
    if(requestedPDG!=0 && (combineCharge ? std::abs(p)!=std::abs(requestedPDG) : p!=requestedPDG)) continue;
    const double Lcm=len/10.0;
    const double eDx=edep/Lcm;
    const double kDx=dke/Lcm;
    ++n; sL+=Lcm; sL2+=Lcm*Lcm; sE+=eDx; sE2+=eDx*eDx; sK+=kDx; sK2+=kDx*kDx;
    hPath.Fill(Lcm); hEdep.Fill(eDx); hDke.Fill(kDx);
  }
  if(!n){std::cerr<<"No selected primary Core rows\n";return;}
  auto rms=[n](double s,double s2){double m=s/n; return std::sqrt(std::max(0.0,s2/n-m*m));};
  std::cout << "Selected primary Core rows=" << n << "\n"
            << "Mean path=" << sL/n << " cm  RMS=" << rms(sL,sL2) << " cm\n"
            << "Mean local Edep/L=" << sE/n << " MeV/cm  RMS=" << rms(sE,sE2) << "\n"
            << "Mean primary DeltaKE/L=" << sK/n << " MeV/cm  RMS=" << rms(sK,sK2) << "\n"
            << "For stopping-power validation use DeltaKE/L; Edep/L is the local deposited-energy observable.\n";

  gSystem->mkdir("validation_plots", kTRUE);
  TCanvas c("cDedxValidation","dE/dx validation",1500,500);
  c.Divide(3,1);
  c.cd(1); hPath.Draw();
  c.cd(2); hEdep.Draw();
  c.cd(3); gPad->SetLogy(); hDke.SetStats(0); hDke.Draw();
  TLatex dkeText; dkeText.SetNDC(); dkeText.SetTextSize(0.035);
  dkeText.DrawLatex(0.52,0.86,Form("full-sample mean = %.4f MeV/cm",sK/n));
  dkeText.DrawLatex(0.52,0.81,Form("full-sample RMS = %.4f MeV/cm",rms(sK,sK2)));
  dkeText.DrawLatex(0.52,0.76,Form("overflow (>30 MeV/cm) = %.0f",hDke.GetBinContent(hDke.GetNbinsX()+1)));
  c.SaveAs("validation_plots/dedx_validation.png");
  c.SaveAs("validation_plots/dedx_validation.pdf");
  std::cout << "Wrote validation_plots/dedx_validation.{png,pdf}\n";
}
