#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

void validate_photon_birth_truth(
    const char* inputName="validation_photon_birth_truth_run0.root") {
  TFile* f = TFile::Open(inputName, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << inputName << "\n";
    return;
  }

  auto* p = dynamic_cast<TTree*>(f->Get("Photon"));
  auto* cfg = dynamic_cast<TTree*>(f->Get("RunConfig"));
  if (!p || !cfg) {
    std::cerr << "Missing Photon or RunConfig\n";
    f->Close();
    return;
  }

  int boundaryDebug=-1;
  char revision[256] = {};
  cfg->SetBranchAddress("BoundaryDebug",&boundaryDebug);
  if (cfg->GetBranch("ProjectRevision"))
    cfg->SetBranchAddress("ProjectRevision",revision);
  cfg->GetEntry(0);

  int event=0, photon=0, reached=0;
  float bt=0,bx=0,by=0,bz=0,bdx=0,bdy=0,bdz=0;
  float wt=0,wx=0,wy=0,wz=0;
  p->SetBranchAddress("EventID",&event);
  p->SetBranchAddress("PhotonID",&photon);
  p->SetBranchAddress("BirthTime_ns",&bt);
  p->SetBranchAddress("BirthX_mm",&bx);
  p->SetBranchAddress("BirthY_mm",&by);
  p->SetBranchAddress("BirthZ_mm",&bz);
  p->SetBranchAddress("BirthDirX",&bdx);
  p->SetBranchAddress("BirthDirY",&bdy);
  p->SetBranchAddress("BirthDirZ",&bdz);
  p->SetBranchAddress("ReachedCoreWall",&reached);
  p->SetBranchAddress("FirstWallTime_ns",&wt);
  p->SetBranchAddress("FirstWallX_mm",&wx);
  p->SetBranchAddress("FirstWallY_mm",&wy);
  p->SetBranchAddress("FirstWallZ_mm",&wz);

  TH1D hDt("hDt","First-wall flight time;t_{wall}-t_{birth} (ns);Photons",100,0,0.5);
  TH1D hWallX("hWallX","First wall X;X_{wall} (mm);Photons",100,49.9,50.1);

  long long n=0,nReached=0,nBadTime=0,nBadBirthPos=0,nBadBirthDir=0,nBadWall=0;
  double maxBirthTimeResidual=0,maxBirthPos=0,maxBirthDirResidual=0,maxWallResidual=0;
  double minDt=std::numeric_limits<double>::infinity();
  double maxDt=-std::numeric_limits<double>::infinity();

  for (Long64_t i=0;i<p->GetEntries();++i) {
    p->GetEntry(i);
    ++n;
    if (reached) ++nReached;

    maxBirthTimeResidual=std::max(maxBirthTimeResidual,std::abs(double(bt)-7.0));
    maxBirthPos=std::max(maxBirthPos,std::max({std::abs(double(bx)),
                                               std::abs(double(by)),
                                               std::abs(double(bz))}));
    maxBirthDirResidual=std::max(maxBirthDirResidual,std::max({
        std::abs(double(bdx)-1.0),std::abs(double(bdy)),std::abs(double(bdz))}));

    if (std::abs(double(bt)-7.0)>1e-5) ++nBadTime;
    if (std::max({std::abs(double(bx)),std::abs(double(by)),std::abs(double(bz))})>1e-4)
      ++nBadBirthPos;
    if (std::max({std::abs(double(bdx)-1.0),std::abs(double(bdy)),std::abs(double(bdz))})>1e-6)
      ++nBadBirthDir;

    if (reached) {
      const double wallResidual=std::max({std::abs(double(wx)-50.0),
                                          std::abs(double(wy)),
                                          std::abs(double(wz))});
      maxWallResidual=std::max(maxWallResidual,wallResidual);
      if (wallResidual>1e-3) ++nBadWall;

      const double dt=double(wt)-double(bt);
      minDt=std::min(minDt,dt);
      maxDt=std::max(maxDt,dt);
      hDt.Fill(dt);
      hWallX.Fill(wx);
    }
  }

  // 50 mm in vacuum already takes ~0.1668 ns. Water must be slower than vacuum.
  // This broad lower bound catches the old implementation, where BirthTime was
  // recorded at the first-step endpoint and dt was ~0 for this one-step ray.
  const bool flightTimePhysical =
      nReached==n && hDt.GetEntries()==n && minDt>0.16 && maxDt<1.0;

  const bool pass =
      n==100 && nReached==100 && boundaryDebug==0 &&
      nBadTime==0 && nBadBirthPos==0 && nBadBirthDir==0 && nBadWall==0 &&
      flightTimePhysical;

  std::cout << "\n=== PHOTON BIRTH TRUTH VALIDATION ===\n";
  std::cout << "Input: " << inputName << "\n";
  std::cout << "ProjectRevision: " << (revision[0] ? revision : "?") << "\n";
  std::cout << "BoundaryDebug=" << boundaryDebug << "\n";
  std::cout << "Photon rows=" << n << " reachedWall=" << nReached << "\n";
  std::cout << "max |BirthTime-7ns| = " << maxBirthTimeResidual << " ns\n";
  std::cout << "max |BirthXYZ| = " << maxBirthPos << " mm\n";
  std::cout << "max BirthDir residual = " << maxBirthDirResidual << "\n";
  std::cout << "max wall residual from (+X,50mm) = " << maxWallResidual << " mm\n";
  std::cout << "wall-birth dt: min=" << minDt
            << " mean=" << hDt.GetMean() << " max=" << maxDt << " ns\n";
  std::cout << "PHOTON BIRTH TRUTH VALIDATION: " << (pass?"PASS":"FAIL") << "\n";

  gSystem->mkdir("validation_plots",true);
  TCanvas c("c","Photon birth truth validation",1200,500);
  c.Divide(2,1);
  c.cd(1); hDt.SetStats(0); hDt.Draw();
  c.cd(2); hWallX.SetStats(0); hWallX.Draw();
  c.SaveAs("validation_plots/photon_birth_truth_validation.png");
  c.SaveAs("validation_plots/photon_birth_truth_validation.pdf");

  f->Close();
}