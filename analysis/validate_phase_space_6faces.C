#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <string>

void validate_phase_space_6faces(const char* prefix="validation_phase_space_6faces") {
  const int expectedChannel[6] = {0,1,2,3,4,5}; const int oppositeChannel[6] = {1,0,3,2,5,4}; const char* faceLabel[6]={"+X","-X","+Z","-Z","+Y","-Y"}; bool allPass = true;
  TH1D hPrimary("hPhasePrimary","6-face destination counts;Injected face;Photons",6,0.5,6.5), hOpp("hPhaseOpposite","",6,0.5,6.5), hMaxZ("hPhaseMaxZ","Maximum entrance-plane residual;Injected face;max |z_{entrance}| (mm)",6,0.5,6.5), hMaxAngle("hPhaseMaxAngle","Maximum local incident-angle residual;Injected face;max |#theta| (deg)",6,0.5,6.5);
  for(int b=1;b<=6;++b){ hPrimary.GetXaxis()->SetBinLabel(b,faceLabel[b-1]); hOpp.GetXaxis()->SetBinLabel(b,faceLabel[b-1]); hMaxZ.GetXaxis()->SetBinLabel(b,faceLabel[b-1]); hMaxAngle.GetXaxis()->SetBinLabel(b,faceLabel[b-1]); }
  for (int run=0; run<6; ++run) {
    const std::string fn = std::string(prefix) + "_run" + std::to_string(run) + ".root"; TFile f(fn.c_str(), "READ"); auto* t = dynamic_cast<TTree*>(f.Get("SensorIncident")); if (!t) { std::cerr << "run" << run << ": missing SensorIncident in " << fn << "\n"; allPass = false; continue; }
    int channel=-1, forward=0; float ex=0,ey=0,ez=0,dx=0,dy=0,dz=0,angle=0; t->SetBranchAddress("ChannelID",&channel); t->SetBranchAddress("EntranceX_mm",&ex); t->SetBranchAddress("EntranceY_mm",&ey); t->SetBranchAddress("EntranceZ_mm",&ez); t->SetBranchAddress("LocalDirX",&dx); t->SetBranchAddress("LocalDirY",&dy); t->SetBranchAddress("LocalDirZ",&dz); t->SetBranchAddress("IncidentAngle_deg",&angle); t->SetBranchAddress("Forward",&forward);
    Long64_t badGeom=0,nPrimaryFace=0,nOppositeFace=0,nOtherFace=0; double maxAbsZ=0,maxTransDir=0,maxAngle=0; const auto n=t->GetEntries();
    for (Long64_t i=0;i<n;++i) { t->GetEntry(i); maxAbsZ=std::max(maxAbsZ,std::abs((double)ez)); maxTransDir=std::max(maxTransDir,std::hypot((double)dx,(double)dy)); maxAngle=std::max(maxAngle,std::abs((double)angle)); if (channel == expectedChannel[run]) ++nPrimaryFace; else if (channel == oppositeChannel[run]) ++nOppositeFace; else ++nOtherFace; const bool geomBad = forward != 1 || std::abs(ez) > 1e-3 || dz < 0.999999 || std::hypot(dx,dy) > 1e-5 || std::abs(angle) > 1e-3; if (geomBad) { if (badGeom < 5) std::cout << "  bad geometry row " << i << ": ch=" << channel << " forward=" << forward << " entranceZ=" << ez << " dir=(" << dx << "," << dy << "," << dz << ") angle=" << angle << " deg\n"; ++badGeom; } }
    const bool pass = n > 900 && nPrimaryFace > 900 && nOtherFace == 0 && badGeom == 0; allPass &= pass; hPrimary.SetBinContent(run+1,nPrimaryFace); hOpp.SetBinContent(run+1,nOppositeFace); hMaxZ.SetBinContent(run+1,maxAbsZ); hMaxAngle.SetBinContent(run+1,maxAngle);
    std::cout << "run" << run << " expectedChannel=" << expectedChannel[run] << " entries=" << n << " primary=" << nPrimaryFace << " opposite(Fresnel)=" << nOppositeFace << " other=" << nOtherFace << " badGeom=" << badGeom << " max|EntranceZ|=" << maxAbsZ << " mm max transverse local-dir=" << maxTransDir << " max angle=" << maxAngle << " deg : " << (pass?"PASS":"FAIL") << "\n";
  }
  std::cout << "PHASE-SPACE 6-FACE VALIDATION: " << (allPass?"PASS":"FAIL") << "\n";
  gSystem->mkdir("validation_plots", kTRUE); TCanvas c("cPhaseSpaceValidation","6-face phase-space validation",1200,800); c.Divide(2,2); c.cd(1); hPrimary.Draw("HIST"); hOpp.SetLineStyle(2); hOpp.Draw("HIST SAME"); c.cd(2); hMaxZ.Draw("HIST"); c.cd(3); hMaxAngle.Draw("HIST"); c.cd(4); TH1D hOppFrac("hPhaseOppFrac","Opposite-face Fresnel return;Injected face;Fraction",6,0.5,6.5); for(int b=1;b<=6;++b){hOppFrac.GetXaxis()->SetBinLabel(b,faceLabel[b-1]); const double tot=hPrimary.GetBinContent(b)+hOpp.GetBinContent(b); hOppFrac.SetBinContent(b,tot>0?hOpp.GetBinContent(b)/tot:0);} hOppFrac.Draw("HIST"); c.SaveAs("validation_plots/phase_space_6faces_validation.png"); c.SaveAs("validation_plots/phase_space_6faces_validation.pdf"); std::cout<<"Wrote validation_plots/phase_space_6faces_validation.{png,pdf}\n";
}
