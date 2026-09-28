#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

void validate_optical_bulk(
    const char* absorptionFile="validation_optical_absorption_pmma_run0.root",
    const char* rayleighFile="validation_optical_rayleigh_water_run0.root") {
  bool allPass=true;
  TH1D hAbs("hAbsLength","PMMA absorption track length @300 nm;Track length to absorption (mm);Photons",120,0,180);
  TGraph gRaySurv; gRaySurv.SetName("gRaySurv"); gRaySurv.SetTitle("Water Rayleigh survival @300 nm;Path x (m);P(no Rayleigh before x)");
  double absExpected=24.3902439, rayExpectedM=55.7645;

  {
    TFile f(absorptionFile,"READ");
    auto* t=dynamic_cast<TTree*>(f.Get("Photon"));
    if (!t) { std::cerr<<"Missing Photon in "<<absorptionFile<<"\n"; return; }
    int endReason=0; float length=std::numeric_limits<float>::quiet_NaN();
    t->SetBranchAddress("EndReason",&endReason); t->SetBranchAddress("EndTrackLength_mm",&length);
    Long64_t nAbs=0; double sum=0, sum2=0;
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); if(endReason==2 && std::isfinite(length)){ ++nAbs; sum+=length; sum2+=double(length)*length; hAbs.Fill(length); } }
    const double mean=nAbs?sum/nAbs:0; const double var=nAbs>1?(sum2-nAbs*mean*mean)/(nAbs-1):0; const double sem=nAbs?std::sqrt(std::max(0.0,var)/nAbs):0;
    const bool pass=nAbs>19000 && std::abs(mean-absExpected)/absExpected<0.03; allPass &= pass;
    std::cout<<"PMMA OpAbsorption @300 nm: Nabs="<<nAbs<<" mean free path="<<mean<<" +/- "<<sem<<" mm expected="<<absExpected<<" mm : "<<(pass?"PASS":"FAIL")<<"\n";
  }

  {
    TFile f(rayleighFile,"READ");
    auto* t=dynamic_cast<TTree*>(f.Get("Photon"));
    if (!t) { std::cerr<<"Missing Photon in "<<rayleighFile<<"\n"; return; }
    int nRay=0; float first=std::numeric_limits<float>::quiet_NaN();
    t->SetBranchAddress("NRayleigh",&nRay); t->SetBranchAddress("FirstRayleighTrackLength_mm",&first);
    Long64_t scat=0,cens=0; double exposure=0; constexpr double censor_mm=50000.0;
    std::vector<double> firstOrCensor; firstOrCensor.reserve(t->GetEntries());
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); if(nRay>0 && std::isfinite(first)){ ++scat; exposure+=first; firstOrCensor.push_back(first); } else { ++cens; exposure+=censor_mm; firstOrCensor.push_back(censor_mm); } }
    const double lambdaHat_mm=scat?exposure/scat:0; const double expected_mm=55764.5;
    const bool pass=scat>10000 && std::abs(lambdaHat_mm-expected_mm)/expected_mm<0.03; allPass &= pass;
    std::cout<<"Water OpRayleigh @300 nm: scat="<<scat<<" censored="<<cens<<" lambda_MLE="<<lambdaHat_mm/1000.0<<" m expected="<<expected_mm/1000.0<<" m : "<<(pass?"PASS":"FAIL")<<"\n";
    const int nPts=51;
    for(int j=0;j<nPts;++j){ const double x_mm=censor_mm*j/(nPts-1.0); long long survive=0; for(double v:firstOrCensor) if(v>x_mm || (x_mm>=censor_mm && v>=censor_mm)) ++survive; gRaySurv.SetPoint(j,x_mm/1000.0,double(survive)/firstOrCensor.size()); }
  }

  std::cout<<"OPTICAL BULK PROCESS VALIDATION: "<<(allPass?"PASS":"FAIL")<<"\n";
  gSystem->mkdir("validation_plots", kTRUE);
  TCanvas c("cOpticalBulkValidation","Optical bulk validation",1200,500); c.Divide(2,1);
  c.cd(1); hAbs.Draw(); auto* expPdf=new TF1("expPdf","[0]*exp(-x/[1])",0,180); expPdf->SetParameters(hAbs.GetBinWidth(1)*hAbs.GetEntries()/absExpected,absExpected); expPdf->SetLineStyle(2); expPdf->Draw("SAME");
  c.cd(2); gRaySurv.SetMinimum(0); gRaySurv.SetMaximum(1.05); gRaySurv.Draw("ALP"); auto* expSurv=new TF1("expSurv","exp(-x/[0])",0,50); expSurv->SetParameter(0,rayExpectedM); expSurv->SetLineStyle(2); expSurv->Draw("SAME");
  c.SaveAs("validation_plots/optical_bulk_validation.png"); c.SaveAs("validation_plots/optical_bulk_validation.pdf");
  std::cout << "Wrote validation_plots/optical_bulk_validation.{png,pdf}\n";
}
