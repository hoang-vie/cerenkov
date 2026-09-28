#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

namespace {
bool isReflection(const char* s){ return std::strstr(s,"Reflection") || std::strcmp(s,"BackScattering")==0; }
struct Result { Long64_t n=0,nRef=0,nTrans=0,nTIR=0,nOther=0; double outAngleSum=0; };
Result readFirstBoundary(const std::string& fn){
  TFile f(fn.c_str(),"READ"); auto* t=dynamic_cast<TTree*>(f.Get("BoundaryDebug")); Result r;
  if(!t){ std::cerr<<"Missing BoundaryDebug in "<<fn<<"\n"; return r; }
  int step=0; float px=0,py=0,pz=0; char status[128]={0};
  t->SetBranchAddress("StepNumber",&step); t->SetBranchAddress("PostDirX",&px); t->SetBranchAddress("PostDirY",&py); t->SetBranchAddress("PostDirZ",&pz); t->SetBranchAddress("BoundaryStatusName",status);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); if(step!=1) continue; ++r.n; if(std::strcmp(status,"TotalInternalReflection")==0){++r.nTIR;++r.nRef;} else if(isReflection(status)) ++r.nRef; else if(std::strcmp(status,"FresnelRefraction")==0 || std::strcmp(status,"Transmission")==0){ ++r.nTrans; r.outAngleSum += std::atan2(std::abs((double)py),std::abs((double)px))*180.0/M_PI; } else ++r.nOther; }
  return r;
}
double fresnelS(double n1,double n2,double ti){ const double st=n1/n2*std::sin(ti); const double tt=std::asin(st); const double rs=(n1*std::cos(ti)-n2*std::cos(tt))/(n1*std::cos(ti)+n2*std::cos(tt)); return rs*rs; }
double fresnelP(double n1,double n2,double ti){ const double st=n1/n2*std::sin(ti); const double tt=std::asin(st); const double rp=(n2*std::cos(ti)-n1*std::cos(tt))/(n2*std::cos(ti)+n1*std::cos(tt)); return rp*rp; }
}

void validate_fresnel_tir(const char* prefix="validation_fresnel_tir_pmma_air"){
  const double n1=1.492,n2=1.0003,ti30=30.0*M_PI/180.0; const double Rs=fresnelS(n1,n2,ti30), Rp=fresnelP(n1,n2,ti30); const double thetaT=std::asin(n1/n2*std::sin(ti30))*180.0/M_PI; const double critical=std::asin(n2/n1)*180.0/M_PI;
  bool allPass=true; double rMC[3]={0},rRef[3]={Rs,Rp,1.0},thetaMC[2]={0,0};
  for(int run=0;run<3;++run){ const auto r=readFirstBoundary(std::string(prefix)+"_run"+std::to_string(run)+".root"); const double refl=r.n?double(r.nRef)/r.n:0; rMC[run]=refl; bool pass=false;
    if(run==0){ thetaMC[0]=r.nTrans?r.outAngleSum/r.nTrans:0; pass=r.n>25000&&r.nOther==0&&std::abs(refl-Rs)<0.01&&r.nTrans>0&&std::abs(thetaMC[0]-thetaT)<0.2; std::cout<<"30deg s: R_MC="<<refl<<" R_Fresnel="<<Rs<<" theta_t_MC="<<thetaMC[0]<<" theta_t_Snell="<<thetaT; }
    else if(run==1){ thetaMC[1]=r.nTrans?r.outAngleSum/r.nTrans:0; pass=r.n>25000&&r.nOther==0&&std::abs(refl-Rp)<0.01&&r.nTrans>0&&std::abs(thetaMC[1]-thetaT)<0.2; std::cout<<"30deg p: R_MC="<<refl<<" R_Fresnel="<<Rp<<" theta_t_MC="<<thetaMC[1]<<" theta_t_Snell="<<thetaT; }
    else { pass=r.n>9000&&r.nTIR==r.n&&r.nOther==0; std::cout<<"50deg TIR: TIR="<<r.nTIR<<"/"<<r.n<<" theta_critical="<<critical<<" deg"; }
    allPass&=pass; std::cout<<" : "<<(pass?"PASS":"FAIL")<<"\n"; }
  std::cout<<"FRESNEL/SNELL/TIR VALIDATION: "<<(allPass?"PASS":"FAIL")<<"\n";
  gSystem->mkdir("validation_plots", kTRUE);
  TH1D hMC("hFresnelMC","Reflection/TIR probability;Case;Probability",3,0.5,3.5), hRef("hFresnelRef","",3,0.5,3.5);
  const char* labels[3]={"30deg s","30deg p","50deg TIR"}; for(int i=0;i<3;++i){hMC.GetXaxis()->SetBinLabel(i+1,labels[i]); hRef.GetXaxis()->SetBinLabel(i+1,labels[i]); hMC.SetBinContent(i+1,rMC[i]); hRef.SetBinContent(i+1,rRef[i]);}
  TH1D hTheta("hSnellTheta","Transmitted angle at 30deg incidence;Polarization;#theta_{t} (deg)",2,0.5,2.5), hThetaRef("hSnellThetaRef","",2,0.5,2.5); hTheta.GetXaxis()->SetBinLabel(1,"s"); hTheta.GetXaxis()->SetBinLabel(2,"p"); hThetaRef.GetXaxis()->SetBinLabel(1,"s"); hThetaRef.GetXaxis()->SetBinLabel(2,"p"); for(int i=0;i<2;++i){hTheta.SetBinContent(i+1,thetaMC[i]);hThetaRef.SetBinContent(i+1,thetaT);}
  TCanvas c("cFresnelValidation","Fresnel/Snell/TIR validation",1200,500); c.Divide(2,1); c.cd(1); hMC.SetMinimum(0); hMC.SetMaximum(1.1); hMC.Draw("HIST"); hRef.SetLineStyle(2); hRef.Draw("HIST SAME"); c.cd(2); hTheta.SetMinimum(0); hTheta.SetMaximum(90); hTheta.Draw("HIST"); hThetaRef.SetLineStyle(2); hThetaRef.Draw("HIST SAME"); c.SaveAs("validation_plots/fresnel_tir_validation.png"); c.SaveAs("validation_plots/fresnel_tir_validation.pdf");
  std::cout << "Wrote validation_plots/fresnel_tir_validation.{png,pdf}\n";
}
