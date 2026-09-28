#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>

void validate_pde(const char* inputName="validation_pde_450nm_run0.root", double expected=0.400) {
  TFile f(inputName,"READ"); if (f.IsZombie()) { std::cerr << "Cannot open file\n"; return; }
  auto* inc=dynamic_cast<TTree*>(f.Get("SensorIncident")); auto* det=dynamic_cast<TTree*>(f.Get("DetectedSignal")); if (!inc || !det) { std::cerr << "Missing SensorIncident/DetectedSignal\n"; return; }
  const double n=inc->GetEntries(), k=det->GetEntries(); if (n<=0) { std::cerr << "No incident photons\n"; return; }
  const double p=k/n,sigma=std::sqrt(expected*(1.0-expected)/n),pull=sigma>0?(p-expected)/sigma:1e99; const bool pass=n>90000&&std::abs(p-expected)<std::max(5.0*sigma,0.005);
  std::cout << "Nincident=" << n << " Ndetected=" << k << "\n" << "Measured efficiency=" << p << " expected=" << expected << " binomial sigma=" << sigma << " pull=" << pull << "\n" << "DIRECT-SIPM PDE VALIDATION: " << (pass?"PASS":"FAIL") << "\n";
  gSystem->mkdir("validation_plots", kTRUE); double x[2]={1,2},y[2]={p,expected},ex[2]={0,0},ey[2]={sigma,0}; TGraphErrors g(2,x,y,ex,ey); g.SetTitle("S13360-3050CS effective PDE validation @450 nm;Estimate;Detection probability"); g.GetXaxis()->SetLimits(0.5,2.5); g.SetMinimum(std::max(0.0,expected-0.05)); g.SetMaximum(std::min(1.0,expected+0.05)); TCanvas c("cPDEValidation","PDE validation",700,500); g.Draw("AP"); c.SaveAs("validation_plots/pde_450nm_validation.png"); c.SaveAs("validation_plots/pde_450nm_validation.pdf"); std::cout<<"Plot point 1=MC, point 2=expected. Wrote validation_plots/pde_450nm_validation.{png,pdf}\n";
}
