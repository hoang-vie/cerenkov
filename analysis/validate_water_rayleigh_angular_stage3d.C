#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace {

double Clamp(double x) { return std::max(-1.0,std::min(1.0,x)); }

double CdfCosThetaG4(double u) {
  u=Clamp(u);
  return 0.5 + 3.0*u/8.0 + u*u*u/8.0;
}

double CdfPolAxisG4(double q) {
  q=Clamp(q);
  return 0.5 + 3.0*q/4.0 - q*q*q/4.0;
}

template<class F>
double KSDistance(std::vector<double> x, F cdf) {
  if(x.empty()) return 1.0;
  std::sort(x.begin(),x.end());
  const double n=double(x.size());
  double d=0.0;
  for(std::size_t i=0;i<x.size();++i){
    const double f=cdf(x[i]);
    const double lo=double(i)/n;
    const double hi=double(i+1)/n;
    d=std::max(d,std::max(std::abs(f-lo),std::abs(hi-f)));
  }
  return d;
}

}

void validate_water_rayleigh_angular_stage3d(
    const char* inputName="validation_water_rayleigh_angular_stage3d_run0.root") {
  TFile* f=TFile::Open(inputName,"READ");
  if(!f || f->IsZombie()){
    std::cerr<<"Cannot open "<<inputName<<"\n";
    return;
  }

  auto* dbg=dynamic_cast<TTree*>(f->Get("BoundaryDebug"));
  auto* pho=dynamic_cast<TTree*>(f->Get("Photon"));
  auto* cfg=dynamic_cast<TTree*>(f->Get("RunConfig"));
  if(!dbg || !pho || !cfg){
    std::cerr<<"Missing BoundaryDebug, Photon, or RunConfig\n";
    f->Close(); return;
  }

  char schema[128]={};
  char revision[256]={};
  char dataset[256]={};
  char rayStatus[256]={};
  cfg->SetBranchAddress("SchemaVersion",schema);
  cfg->SetBranchAddress("ProjectRevision",revision);
  cfg->SetBranchAddress("CoreProfileRayleighDataset",dataset);
  cfg->SetBranchAddress("CoreProfileRayleighStatus",rayStatus);
  cfg->GetEntry(0);

  int nRay=0;
  pho->SetBranchAddress("NRayleigh",&nRay);
  Long64_t photonsWithRayleigh=0;
  for(Long64_t i=0;i<pho->GetEntries();++i){
    pho->GetEntry(i);
    if(nRay>0) ++photonsWithRayleigh;
  }

  int status=0;
  float preDx=0,preDy=0,preDz=0,postDx=0,postDy=0,postDz=0;
  float prePx=0,prePy=0,prePz=0,postPx=0,postPy=0,postPz=0;
  char statusName[128]={};
  dbg->SetBranchAddress("BoundaryStatus",&status);
  dbg->SetBranchAddress("BoundaryStatusName",statusName);
  dbg->SetBranchAddress("PreDirX",&preDx); dbg->SetBranchAddress("PreDirY",&preDy); dbg->SetBranchAddress("PreDirZ",&preDz);
  dbg->SetBranchAddress("PostDirX",&postDx); dbg->SetBranchAddress("PostDirY",&postDy); dbg->SetBranchAddress("PostDirZ",&postDz);
  dbg->SetBranchAddress("PrePolX",&prePx); dbg->SetBranchAddress("PrePolY",&prePy); dbg->SetBranchAddress("PrePolZ",&prePz);
  dbg->SetBranchAddress("PostPolX",&postPx); dbg->SetBranchAddress("PostPolY",&postPy); dbg->SetBranchAddress("PostPolZ",&postPz);

  std::vector<double> uvals,qvals;
  uvals.reserve(40000); qvals.reserve(40000);
  TH1D hU("hRayleighCosTheta","First Rayleigh scatter;cos(#theta);rows",80,-1,1);
  TH1D hQ("hRayleighPolAxis","First Rayleigh scatter;post direction #upoint pre polarization;rows",80,-1,1);

  double sumU=0,sumU2=0,sumQ=0,sumQ2=0;
  double maxInputDirResidual=0,maxInputPolResidual=0;
  double maxPostTransverse=0,maxProjectionResidual=0;
  Long64_t tagged=0,badSentinel=0,badProjection=0;

  for(Long64_t i=0;i<dbg->GetEntries();++i){
    dbg->GetEntry(i);
    if(std::string(statusName)!="FirstRayleighScatter") continue;
    ++tagged;
    if(status!=-1001) ++badSentinel;

    const double pdi=std::sqrt(double(preDx)*preDx+double(preDy)*preDy+double(preDz)*preDz);
    const double pdo=std::sqrt(double(postDx)*postDx+double(postDy)*postDy+double(postDz)*postDz);
    const double ppi=std::sqrt(double(prePx)*prePx+double(prePy)*prePy+double(prePz)*prePz);
    const double ppo=std::sqrt(double(postPx)*postPx+double(postPy)*postPy+double(postPz)*postPz);
    if(pdi<=0 || pdo<=0 || ppi<=0 || ppo<=0){ ++badProjection; continue; }

    const double idx=preDx/pdi, idy=preDy/pdi, idz=preDz/pdi;
    const double odx=postDx/pdo, ody=postDy/pdo, odz=postDz/pdo;
    const double ipx=prePx/ppi, ipy=prePy/ppi, ipz=prePz/ppi;
    const double opx=postPx/ppo, opy=postPy/ppo, opz=postPz/ppo;

    maxInputDirResidual=std::max(maxInputDirResidual,std::sqrt(idx*idx+idy*idy+(idz-1)*(idz-1)));
    maxInputPolResidual=std::max(maxInputPolResidual,std::sqrt((ipx-1)*(ipx-1)+ipy*ipy+ipz*ipz));

    const double u=Clamp(idx*odx+idy*ody+idz*odz);
    const double q=Clamp(ipx*odx+ipy*ody+ipz*odz);
    uvals.push_back(u); qvals.push_back(q);
    hU.Fill(u); hQ.Fill(q);
    sumU+=u; sumU2+=u*u; sumQ+=q; sumQ2+=q*q;

    const double transverse=std::abs(odx*opx+ody*opy+odz*opz);
    maxPostTransverse=std::max(maxPostTransverse,transverse);

    // G4OpRayleigh constructs post polarization from the projection of the
    // pre polarization onto the plane perpendicular to post momentum, with
    // an arbitrary +/- sign. Check |dot| against that projected unit vector.
    double ex=ipx-q*odx, ey=ipy-q*ody, ez=ipz-q*odz;
    const double em=std::sqrt(ex*ex+ey*ey+ez*ez);
    if(em>1e-8){
      ex/=em; ey/=em; ez/=em;
      const double align=std::abs(ex*opx+ey*opy+ez*opz);
      maxProjectionResidual=std::max(maxProjectionResidual,std::abs(1.0-align));
    }
  }

  const double n=double(uvals.size());
  const double meanU=n?sumU/n:0;
  const double meanU2=n?sumU2/n:0;
  const double meanQ=n?sumQ/n:0;
  const double meanQ2=n?sumQ2/n:0;
  const double ksU=KSDistance(uvals,CdfCosThetaG4);
  const double ksQ=KSDistance(qvals,CdfPolAxisG4);
  const double ksLimit=n>0?2.0/std::sqrt(n):0;

  const bool metadataOK=
      std::string(schema)=="general-v1.7" &&
      std::string(revision)=="hotfix6_physics_audit5_water_stage3e_pmh_default" &&
      std::string(dataset)=="Water_RAYLEIGH_ZhangHu2009_PMH_21p5C.txt" &&
      std::string(rayStatus)=="ZhangHu2009_PMH_21p5C_delta0p039_MFPValidated_AngularApproxDocumented";
  const bool countOK=tagged==photonsWithRayleigh && tagged>20000;
  const bool shapeOK=ksU<ksLimit && ksQ<ksLimit && std::abs(meanU)<0.015 &&
                     std::abs(meanQ)<0.015 && std::abs(meanU2-0.4)<0.01 &&
                     std::abs(meanQ2-0.2)<0.01;
  const bool polOK=badProjection==0 && maxInputDirResidual<1e-5 &&
                   maxInputPolResidual<1e-5 && maxPostTransverse<1e-5 &&
                   maxProjectionResidual<1e-5;
  const bool pass=metadataOK && countOK && badSentinel==0 && shapeOK && polOK;

  std::cout<<"\n=== WATER RAYLEIGH ANGULAR STAGE-3D VALIDATION ===\n";
  std::cout<<"tagged first scatters="<<tagged
           <<" photons NRayleigh>0="<<photonsWithRayleigh
           <<" count="<<(countOK?"OK":"FAIL")<<"\n";
  std::cout<<"cos(theta): mean="<<meanU<<" <u^2>="<<meanU2
           <<" expected <u^2>=0.4 KS="<<ksU<<" limit="<<ksLimit<<"\n";
  std::cout<<"q=postDir.prePol: mean="<<meanQ<<" <q^2>="<<meanQ2
           <<" expected <q^2>=0.2 KS="<<ksQ<<" limit="<<ksLimit<<"\n";
  std::cout<<"polarization: max|postDir.postPol|="<<maxPostTransverse
           <<" max projection residual="<<maxProjectionResidual
           <<" inputDirResidual="<<maxInputDirResidual
           <<" inputPolResidual="<<maxInputPolResidual<<"\n";
  std::cout<<"metadata="<<(metadataOK?"OK":"FAIL")
           <<" sentinelBad="<<badSentinel<<"\n";

  gSystem->mkdir("validation_plots",kTRUE);
  TCanvas c("cRayleighAngular3D","Rayleigh angular Stage 3D",1200,500);
  c.Divide(2,1);
  c.cd(1); hU.Draw("E");
  TF1 fu("fu","[0]*(3./8.)*(1+x*x)",-1,1);
  fu.SetParameter(0,tagged*hU.GetBinWidth(1)); fu.Draw("SAME");
  c.cd(2); hQ.Draw("E");
  TF1 fq("fq","[0]*(3./4.)*(1-x*x)",-1,1);
  fq.SetParameter(0,tagged*hQ.GetBinWidth(1)); fq.Draw("SAME");
  c.SaveAs("validation_plots/water_rayleigh_angular_stage3d.png");
  c.SaveAs("validation_plots/water_rayleigh_angular_stage3d.pdf");

  std::cout<<"WATER RAYLEIGH GEANT4 ANGULAR LAW: "<<(pass?"PASS":"FAIL")<<"\n";
  std::cout<<"WATER DELTA=0.039 PHYSICAL PHASE FUNCTION: NOT_EXPLICITLY_MODELED\n";
  f->Close();
}
