#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <limits>
#include <utility>
#include <vector>
#include <fstream>
#include <sstream>

namespace {
std::vector<std::pair<double,double>> LoadRindexTable(const char* path){
  std::ifstream in(path); std::vector<std::pair<double,double>> out; std::string line;
  while(std::getline(in,line)){
    if(line.empty() || line[0]=='#') continue;
    std::istringstream ss(line); double wl=0,n=0; if(ss>>wl>>n) out.push_back({wl,n});
  }
  return out;
}
double Interp(const std::vector<std::pair<double,double>>& tab,double x){
  if(tab.empty()) return std::numeric_limits<double>::quiet_NaN();
  if(x<=tab.front().first) return tab.front().second;
  if(x>=tab.back().first) return tab.back().second;
  auto it=std::lower_bound(tab.begin(),tab.end(),x,[](const auto& a,double v){return a.first<v;});
  const auto& b=*it; const auto& a=*(it-1); const double f=(x-a.first)/(b.first-a.first); return a.second+f*(b.second-a.second);
}
}

void validate_frank_tamm(const char* inputName="validation_frank_tamm_water_run0.root") {
  TFile f(inputName,"READ");
  if (f.IsZombie()) { std::cerr << "Cannot open file\n"; return; }
  auto* tr=dynamic_cast<TTree*>(f.Get("Track"));
  auto* tv=dynamic_cast<TTree*>(f.Get("TrackVolume"));
  auto* ph=dynamic_cast<TTree*>(f.Get("Photon"));
  if (!tr || !tv || !ph) { std::cerr << "Missing Track/TrackVolume/Photon\n"; return; }

  int e=0,t=0,p=0,par=0; float vertexKE=0;
  std::map<std::pair<int,int>,int> pdg,parent; std::map<std::pair<int,int>,double> vke;
  tr->SetBranchAddress("EventID",&e); tr->SetBranchAddress("TrackID",&t);
  tr->SetBranchAddress("PDG",&p); tr->SetBranchAddress("ParentID",&par); tr->SetBranchAddress("VertexKE_MeV",&vertexKE);
  for (Long64_t i=0;i<tr->GetEntries();++i) { tr->GetEntry(i); pdg[{e,t}]=p; parent[{e,t}]=par; vke[{e,t}]=vertexKE; }

  std::map<std::pair<int,int>,double> length;
  int vol=0; float L=0;
  tv->SetBranchAddress("EventID",&e); tv->SetBranchAddress("TrackID",&t);
  tv->SetBranchAddress("VolumeID",&vol); tv->SetBranchAddress("Length_mm",&L);
  for (Long64_t i=0;i<tv->GetEntries();++i) {
    tv->GetEntry(i); if (vol==0) length[{e,t}] += L;
  }

  std::map<std::pair<int,int>,long long> npho;
  int pe=0,ppar=0,origin=0;
  float angle=0,wl=0,emitterBeta=std::numeric_limits<float>::quiet_NaN();
  ph->SetBranchAddress("EventID",&pe);
  ph->SetBranchAddress("ParentTrackID",&ppar);
  ph->SetBranchAddress("OriginVolumeID",&origin);
  ph->SetBranchAddress("Wavelength_nm",&wl);
  ph->SetBranchAddress("EmissionAngle_deg",&angle);
  const bool hasEmitterBeta = ph->GetBranch("EmitterBeta") != nullptr;
  if (hasEmitterBeta) ph->SetBranchAddress("EmitterBeta",&emitterBeta);
  TH1D hWl("hFrankTammWl","Generated Cherenkov spectrum;Wavelength (nm);Photons",140,200,900);
  TH2D hAngleWl("hFrankTammAngleWl","Primary-#mu Cherenkov angle vs wavelength;Wavelength (nm);Emission angle (deg)",140,200,900,120,35,55);
  TH1D hAngleResidual("hFrankTammResidual","Primary-#mu Cherenkov-angle residual;#theta_{MC}-#theta_{theory} (deg);Photons",200,-0.10,0.10);
  const auto rindex = LoadRindexTable("../data/Water_RINDEX_Daimon2007_21p5C.txt");
  for (Long64_t i=0;i<ph->GetEntries();++i) {
    ph->GetEntry(i);
    if (origin==0 && std::isfinite(angle)) {
      npho[{pe,ppar}]++;
      hWl.Fill(wl);
      const auto key=std::make_pair(pe,ppar);
      auto ip=pdg.find(key), ipa=parent.find(key);
      if(ip!=pdg.end() && ipa!=parent.end() && ipa->second==0 && std::abs(ip->second)==13){
        hAngleWl.Fill(wl,angle);
        const double nr=Interp(rindex,wl);
        if(hasEmitterBeta && std::isfinite(emitterBeta) && std::isfinite(nr) && emitterBeta*nr>1.0){
          const double theory=std::acos(1.0/(double(emitterBeta)*nr))*180.0/M_PI;
          hAngleResidual.Fill(angle-theory);
        }
      }
    }
  }

  TH1D hYield("hFrankTammYield","Cherenkov yield per primary muon;N_{#gamma}/L (photons/cm);Tracks",180,650,1000);
  double sum=0, sum2=0; long long n=0;
  for (const auto& kv : length) {
    const auto key=kv.first;
    auto ip=pdg.find(key), ipa=parent.find(key);
    if (ip==pdg.end() || ipa==parent.end()) continue;
    if (ipa->second!=0 || std::abs(ip->second)!=13 || kv.second<=0) continue;
    const double y = npho[key] / (kv.second/10.0);
    sum += y; sum2 += y*y; ++n; hYield.Fill(y);
  }
  if (!n) { std::cerr << "No primary-muon Core path rows\n"; return; }
  const double mean=sum/n;
  const double rms=std::sqrt(std::max(0.0,sum2/n-mean*mean));
  const double sem=rms/std::sqrt(double(n));
  const double expected=817.0;
  const bool pass=std::abs(mean-expected)/expected < 0.03;
  std::cout << "Tracks=" << n << " mean N_gamma/L=" << mean
            << " photons/cm, RMS=" << rms << ", SEM=" << sem << "\n"
            << "Independent Frank-Tamm expectation~" << expected << " photons/cm\n"
            << "FRANK-TAMM YIELD VALIDATION: " << (pass?"PASS":"FAIL") << "\n";

  gSystem->mkdir("validation_plots", kTRUE);
  TCanvas c("cFrankTammValidation","Frank-Tamm validation",1500,900);
  c.Divide(2,2);
  c.cd(1); hYield.Draw(); auto* yieldLine=new TLine(expected,0,expected,hYield.GetMaximum()*1.05); yieldLine->SetLineStyle(2); yieldLine->Draw("SAME");
  c.cd(2); hWl.Draw();
  c.cd(3); hAngleWl.Draw("COLZ");
  c.cd(4); hAngleResidual.Draw(); auto* zeroLine=new TLine(0,0,0,hAngleResidual.GetMaximum()*1.05); zeroLine->SetLineStyle(2); zeroLine->Draw("SAME");
  if (hasEmitterBeta) {
    std::cout << "Primary-muon exact angle residual using Photon.EmitterBeta: mean="
              << hAngleResidual.GetMean() << " deg, RMS=" << hAngleResidual.GetRMS()
              << " deg, entries=" << hAngleResidual.GetEntries() << "\n";
  } else {
    std::cout << "WARNING: Photon.EmitterBeta is absent; exact angle residual not evaluated.\n";
  }
  c.SaveAs("validation_plots/frank_tamm_validation.png");
  c.SaveAs("validation_plots/frank_tamm_validation.pdf");
  std::cout << "Wrote validation_plots/frank_tamm_validation.{png,pdf}\n";
}
