#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <cstring>
#include <iostream>

void validate_tyvek_reflectivity(const char* file="validation_tyvek_reflectivity_450nm_run0.root") {
  TFile f(file,"READ"); auto* t=dynamic_cast<TTree*>(f.Get("BoundaryDebug")); if(!t){std::cerr<<"Missing BoundaryDebug in "<<file<<"\n";return;}
  int step=0; char status[128]={0}; t->SetBranchAddress("StepNumber",&step); t->SetBranchAddress("BoundaryStatusName",status);
  Long64_t n=0,nRef=0,nAbs=0,nOther=0; for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); if(step!=1) continue; ++n; if(std::strstr(status,"Reflection") || std::strcmp(status,"BackScattering")==0) ++nRef; else if(std::strcmp(status,"Absorption")==0) ++nAbs; else ++nOther; }
  const double R=n?double(nRef)/n:0,expected=0.970,sigma=n?std::sqrt(expected*(1-expected)/n):0; const bool pass=n>45000&&nOther==0&&std::abs(R-expected)<0.005;
  std::cout<<"Tyvek first-hit @450 nm: reflected="<<nRef<<" absorbed="<<nAbs<<" other="<<nOther<<" R_MC="<<R<<" expected="<<expected<<" binomial sigma~"<<sigma<<" : "<<(pass?"PASS":"FAIL")<<"\n";
  std::cout<<"TYVEK REFLECTIVITY VALIDATION: "<<(pass?"PASS":"FAIL")<<"\n";
  gSystem->mkdir("validation_plots", kTRUE); TH1D h("hTyvekR","Tyvek first-hit outcome @450 nm;Outcome;Fraction",3,0.5,3.5); h.GetXaxis()->SetBinLabel(1,"Reflected"); h.GetXaxis()->SetBinLabel(2,"Absorbed"); h.GetXaxis()->SetBinLabel(3,"Other"); h.SetBinContent(1,n?double(nRef)/n:0); h.SetBinContent(2,n?double(nAbs)/n:0); h.SetBinContent(3,n?double(nOther)/n:0); h.SetMaximum(1.05); TCanvas c("cTyvekValidation","Tyvek reflectivity validation",700,500); h.Draw("HIST"); c.SaveAs("validation_plots/tyvek_reflectivity_validation.png"); c.SaveAs("validation_plots/tyvek_reflectivity_validation.pdf"); std::cout<<"Wrote validation_plots/tyvek_reflectivity_validation.{png,pdf}\n";
}
