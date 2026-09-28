#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TMultiGraph.h>
#include <TSystem.h>
#include <TString.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

void validate_flexible_layout(const char* inputName="validation_flexible_layout_run0.root"){
  TFile f(inputName,"READ"); if(f.IsZombie()){std::cerr<<"Cannot open file\n";return;} auto* t=dynamic_cast<TTree*>(f.Get("ChannelMap")); if(!t){std::cerr<<"Missing ChannelMap\n";return;}
  int channel=0,face=0; float u=0,v=0,nx=0,ny=0,nz=0; t->SetBranchAddress("ChannelID",&channel); t->SetBranchAddress("FaceID",&face); t->SetBranchAddress("LocalU_mm",&u); t->SetBranchAddress("LocalV_mm",&v); t->SetBranchAddress("NormalX",&nx); t->SetBranchAddress("NormalY",&ny); t->SetBranchAddress("NormalZ",&nz);
  const int expected[6]={1,4,3,0,2,2}; const char* labels[6]={"+X","-X","+Z","-Z","+Y","-Y"}; int counts[6]={0}; long long badNormal=0; std::vector<double> us[6],vs[6];
  for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); if(face<0||face>=6){continue;} ++counts[face]; us[face].push_back(u);vs[face].push_back(v); if(std::abs(std::sqrt(nx*nx+ny*ny+nz*nz)-1.0)>1e-5) ++badNormal;}
  bool pass=t->GetEntries()==12&&badNormal==0; for(int i=0;i<6;++i){pass &= counts[i]==expected[i]; std::cout<<labels[i]<<": channels="<<counts[i]<<" expected="<<expected[i]<<"\n";} std::cout<<"bad unit normals="<<badNormal<<"\nFLEXIBLE LAYOUT VALIDATION: "<<(pass?"PASS":"FAIL")<<"\n";
  gSystem->mkdir("validation_plots", kTRUE);
  TCanvas c("cFlexibleLayout","Flexible readout geometry",1200,800);
  c.Divide(3,2);
  for(int i=0;i<6;++i){
    c.cd(i+1);
    if(us[i].empty()){
      auto* frame = gPad->DrawFrame(-1.0,-1.0,1.0,1.0,Form("%s face;Local U (mm);Local V (mm)",labels[i]));
      (void)frame;
      TLatex text; text.SetTextAlign(22); text.SetTextSize(0.05);
      text.DrawLatex(0.0,0.0,"disabled / 0 channels");
      continue;
    }
    TGraph* g=new TGraph(us[i].size());
    for(size_t j=0;j<us[i].size();++j) g->SetPoint(j,us[i][j],vs[i][j]);
    g->SetTitle(Form("%s face;Local U (mm);Local V (mm)",labels[i]));
    g->SetMarkerStyle(20);
    g->Draw("AP");
  }
  c.SaveAs("validation_plots/flexible_layout_validation.png");
  c.SaveAs("validation_plots/flexible_layout_validation.pdf");
  std::cout<<"Wrote validation_plots/flexible_layout_validation.{png,pdf}\n";
}
