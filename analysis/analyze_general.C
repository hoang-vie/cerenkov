#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TString.h>
#include <TTree.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLine.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

uint64_t Key(int eventID, int trackID) {
  return (uint64_t(uint32_t(eventID)) << 32) | uint32_t(trackID);
}

std::vector<std::pair<double,double>> LoadRindexTable(const std::string& path) {
  std::ifstream in(path);
  std::vector<std::pair<double,double>> out;
  std::string line;
  while (std::getline(in,line)) {
    if (line.empty() || line[0]=='#') continue;
    std::istringstream ss(line);
    double wl=0,n=0;
    if (ss>>wl>>n) out.push_back({wl,n});
  }
  return out;
}

double InterpRindex(const std::vector<std::pair<double,double>>& tab, double wl) {
  if (tab.empty()) return std::numeric_limits<double>::quiet_NaN();
  if (wl<=tab.front().first) return tab.front().second;
  if (wl>=tab.back().first) return tab.back().second;
  auto it=std::lower_bound(tab.begin(),tab.end(),wl,
      [](const auto& a,double v){ return a.first < v; });
  const auto& b=*it; const auto& a=*(it-1);
  const double f=(wl-a.first)/(b.first-a.first);
  return a.second+f*(b.second-a.second);
}

struct TrackInfo {
  int parent = 0;
  int root = 0;
  int pdg = 0;
  float vertexKE = 0.f;
  float dirX = 0.f, dirY = 0.f, dirZ = 0.f;
  float endProperNs = 0.f;
  std::string endProcess;
  std::string creatorProcess;
};

bool PdgMatch(int pdg, int requested, bool combineCharge) {
  if (requested == 0) return true;
  return combineCharge ? std::abs(pdg) == std::abs(requested) : pdg == requested;
}

TString ParticleLabel(int pdg, bool combineCharge) {
  if (pdg == 0) return "all particles";
  if (combineCharge) {
    if (std::abs(pdg)==11) return "e^{#pm}";
    if (std::abs(pdg)==13) return "#mu^{#pm}";
    if (std::abs(pdg)==211) return "#pi^{#pm}";
    if (std::abs(pdg)==321) return "K^{#pm}";
    if (std::abs(pdg)==2212) return "p/#bar{p}";
    return Form("|PDG|=%d", std::abs(pdg));
  }
  if (pdg==11) return "e^{-}";
  if (pdg==-11) return "e^{+}";
  if (pdg==13) return "#mu^{-}";
  if (pdg==-13) return "#mu^{+}";
  if (pdg==211) return "#pi^{+}";
  if (pdg==-211) return "#pi^{-}";
  if (pdg==321) return "K^{+}";
  if (pdg==-321) return "K^{-}";
  if (pdg==2212) return "p";
  if (pdg==-2212) return "#bar{p}";
  return Form("PDG=%d", pdg);
}

bool SelectPhotonSource(const std::unordered_map<uint64_t,TrackInfo>& tracks,
                        int eventID, int parentTrackID, int sourceTrackID,
                        int selectedPDG, bool combineCharge,
                        const std::string& sourceMode) {
  if (selectedPDG == 0) return true;
  if (sourceMode == "ancestor") {
    auto sourceIt = tracks.find(Key(eventID, sourceTrackID));
    if (sourceIt == tracks.end()) return false;
    auto rootIt = tracks.find(Key(eventID, sourceIt->second.root));
    return rootIt != tracks.end() &&
           PdgMatch(rootIt->second.pdg, selectedPDG, combineCharge);
  }

  auto parentIt = tracks.find(Key(eventID, parentTrackID));
  if (parentIt == tracks.end()) return false; // optical parent => e.g. WLS, not direct emission
  if (sourceMode == "emitter") {
    return PdgMatch(parentIt->second.pdg, selectedPDG, combineCharge);
  }
  return parentIt->second.parent == 0 &&
         PdgMatch(parentIt->second.pdg, selectedPDG, combineCharge);
}

TString SourceLabel(const TString& particle, const std::string& mode) {
  if (particle == "all particles") return "all optical photons";
  if (mode == "ancestor") return particle + " primary-initiated cascade";
  if (mode == "emitter") return particle + " emitter";
  return "primary " + particle + " emitted";
}

} // namespace

void analyze_general(const char* inputName,
                     int selectedPDG = 13,
                     bool combineCharge = true,
                     const char* sourceModeC = "primaryEmitter") {
  const std::string sourceMode(sourceModeC ? sourceModeC : "primaryEmitter");
  if (sourceMode != "primaryEmitter" && sourceMode != "emitter" && sourceMode != "ancestor") {
    std::cerr << "sourceMode must be primaryEmitter, emitter, or ancestor\n";
    return;
  }

  TFile* f = TFile::Open(inputName, "READ");
  if (!f || f->IsZombie()) { std::cerr << "Cannot open " << inputName << "\n"; return; }

  auto* tTrack = dynamic_cast<TTree*>(f->Get("Track"));
  auto* tTV = dynamic_cast<TTree*>(f->Get("TrackVolume"));
  auto* tPho = dynamic_cast<TTree*>(f->Get("Photon"));
  auto* tInc = dynamic_cast<TTree*>(f->Get("SensorIncident"));
  auto* tDet = dynamic_cast<TTree*>(f->Get("DetectedSignal"));
  auto* tMap = dynamic_cast<TTree*>(f->Get("ChannelMap"));
  auto* tTrig = dynamic_cast<TTree*>(f->Get("TriggerEvent"));
  auto* tPrimary = dynamic_cast<TTree*>(f->Get("PrimaryTruth"));
  auto* tCfg = dynamic_cast<TTree*>(f->Get("RunConfig"));
  if (!tTrack || !tTV || !tPho || !tInc || !tDet || !tMap || !tTrig) {
    std::cerr << "Missing one or more required general-schema ntuples.\n";
    return;
  }

  TString particle = ParticleLabel(selectedPDG, combineCharge);
  TString sourceLabel = SourceLabel(particle, sourceMode);
  TString base(inputName);
  Ssiz_t slash = base.Last('/'); if (slash != kNPOS) base.Remove(0, slash+1);
  base.ReplaceAll(".root", "");
  TString pdf = Form("report_%s_PDG%d_%s.pdf", base.Data(), selectedPDG, sourceMode.c_str());

  std::string coreMaterial;
  std::string coreRindexPath;
  if (tCfg && tCfg->GetEntries()>0) {
    char coreMat[128] = {};
    char waterRi[256] = {};
    char acrylicRi[256] = {};
    if (tCfg->GetBranch("CoreMaterial")) tCfg->SetBranchAddress("CoreMaterial",coreMat);
    if (tCfg->GetBranch("WaterRindexDataset")) tCfg->SetBranchAddress("WaterRindexDataset",waterRi);
    if (tCfg->GetBranch("AcrylicRindexDataset")) tCfg->SetBranchAddress("AcrylicRindexDataset",acrylicRi);
    tCfg->GetEntry(0);
    coreMaterial = coreMat;
    if (coreMaterial=="Water" && waterRi[0]) coreRindexPath = std::string("../data/")+waterRi;
    else if (coreMaterial=="Acrylic" && acrylicRi[0]) coreRindexPath = std::string("../data/")+acrylicRi;
  }
  const auto coreRindex = coreRindexPath.empty() ? std::vector<std::pair<double,double>>{} : LoadRindexTable(coreRindexPath);

  std::unordered_map<uint64_t,TrackInfo> tracks;
  int te=0, tid=0, parent=0, root=0, pdg=0;
  float vke=0, dx=0,dy=0,dz=0, proper=0;
  char endProc[256] = {};
  char creatorProc[256] = {};
  tTrack->SetBranchAddress("EventID",&te);
  tTrack->SetBranchAddress("TrackID",&tid);
  tTrack->SetBranchAddress("ParentID",&parent);
  tTrack->SetBranchAddress("RootPrimaryTrackID",&root);
  tTrack->SetBranchAddress("PDG",&pdg);
  tTrack->SetBranchAddress("VertexKE_MeV",&vke);
  tTrack->SetBranchAddress("DirX",&dx); tTrack->SetBranchAddress("DirY",&dy); tTrack->SetBranchAddress("DirZ",&dz);
  tTrack->SetBranchAddress("EndProperTime_ns",&proper);
  tTrack->SetBranchAddress("EndProcessName",endProc);
  tTrack->SetBranchAddress("CreatorProcessName",creatorProc);

  TH1D* hPrimaryE = new TH1D("hPrimaryE","Selected primary kinetic energy;Kinetic energy (MeV);Tracks",160,0,8000);
  TH1D* hZenith = new TH1D("hZenith","Selected primary zenith cosine;cos(#theta_{zenith});Tracks",100,-1,1);
  TH1D* hLifetime = new TH1D("hLifetime","Muon decay proper time;Proper time (#mus);Decays",120,0,12);
  TH1D* hMichel = new TH1D("hMichel","Electron/positron kinetic energy from selected-muon decay;Kinetic energy (MeV);Tracks",110,0,55);

  for (Long64_t i=0;i<tTrack->GetEntries();++i) {
    tTrack->GetEntry(i);
    TrackInfo info;
    info.parent=parent; info.root=root; info.pdg=pdg; info.vertexKE=vke;
    info.dirX=dx; info.dirY=dy; info.dirZ=dz; info.endProperNs=proper;
    info.endProcess = endProc;
    info.creatorProcess = creatorProc;
    tracks[Key(te,tid)] = info;

    if (parent==0 && PdgMatch(pdg,selectedPDG,combineCharge) &&
        std::abs(pdg)==13 && std::string(endProc).find("Decay")!=std::string::npos) {
      hLifetime->Fill(proper/1000.0);
    }
  }

  if (tPrimary) {
    int pre=0, prID=0, pPDG=0;
    float pKE=0, pDx=0, pDy=0, pDz=0;
    tPrimary->SetBranchAddress("EventID",&pre);
    tPrimary->SetBranchAddress("TrackID",&prID);
    tPrimary->SetBranchAddress("PDG",&pPDG);
    tPrimary->SetBranchAddress("VertexKE_MeV",&pKE);
    tPrimary->SetBranchAddress("DirX",&pDx);
    tPrimary->SetBranchAddress("DirY",&pDy);
    tPrimary->SetBranchAddress("DirZ",&pDz);
    for (Long64_t i=0;i<tPrimary->GetEntries();++i) {
      tPrimary->GetEntry(i);
      if (!PdgMatch(pPDG,selectedPDG,combineCharge)) continue;
      hPrimaryE->Fill(pKE);
      const double mag=std::sqrt(pDx*pDx+pDy*pDy+pDz*pDz);
      if (mag>0) hZenith->Fill(-pDy/mag);
    }
  } else {
    for (const auto& kv : tracks) {
      const auto& tr = kv.second;
      if (tr.parent!=0 || !PdgMatch(tr.pdg,selectedPDG,combineCharge)) continue;
      hPrimaryE->Fill(tr.vertexKE);
      const double mag=std::sqrt(tr.dirX*tr.dirX+tr.dirY*tr.dirY+tr.dirZ*tr.dirZ);
      if (mag>0) hZenith->Fill(-tr.dirY/mag);
    }
  }

  for (const auto& kv : tracks) {
    const auto& tr = kv.second;
    if (std::abs(tr.pdg)!=11 || tr.parent==0 || tr.creatorProcess.find("Decay")==std::string::npos) continue;
    const uint32_t eid = uint32_t(kv.first >> 32);
    auto pit = tracks.find(Key(int(eid),tr.parent));
    if (pit!=tracks.end() && std::abs(pit->second.pdg)==13 &&
        PdgMatch(pit->second.pdg,selectedPDG,combineCharge)) hMichel->Fill(tr.vertexKE);
  }

  TH1D* hCorePath = new TH1D("hCorePath",Form("Primary path length in Core -- primary %s;Path length (cm);Tracks",particle.Data()),160,0,20);
  TH1D* hLocalEdepDx = new TH1D("hLocalEdepDx",Form("Local deposited energy per path in Core -- primary %s;E_{dep}/L (MeV/cm);Track-volume rows",particle.Data()),160,0,8);
  TH1D* hStoppingPower = new TH1D("hStoppingPower",Form("Primary kinetic-energy loss per path in Core -- primary %s;#DeltaKE/L (MeV/cm);Track-volume rows",particle.Data()),240,0,30);
  TH2D* hStoppingPowerE = new TH2D("hStoppingPowerE",Form("Kinetic energy vs kinetic-energy loss per path -- primary %s;Vertex KE (MeV);#DeltaKE/L (MeV/cm)",particle.Data()),120,0,6000,180,0,30);
  int ve=0, vt=0, vv=0; float ved=0,vdelta=0,vlen=0;
  tTV->SetBranchAddress("EventID",&ve); tTV->SetBranchAddress("TrackID",&vt);
  tTV->SetBranchAddress("VolumeID",&vv); tTV->SetBranchAddress("Edep_MeV",&ved);
  tTV->SetBranchAddress("DeltaKE_MeV",&vdelta); tTV->SetBranchAddress("Length_mm",&vlen);
  for (Long64_t i=0;i<tTV->GetEntries();++i) {
    tTV->GetEntry(i); if (vv!=0 || vlen<=0) continue;
    auto it=tracks.find(Key(ve,vt)); if (it==tracks.end()) continue;
    if (it->second.parent!=0 || !PdgMatch(it->second.pdg,selectedPDG,combineCharge)) continue;
    const double localEdepDx=ved/(vlen/10.0);
    const double stoppingPower=vdelta/(vlen/10.0);
    hCorePath->Fill(vlen/10.0);
    hLocalEdepDx->Fill(localEdepDx);
    hStoppingPower->Fill(stoppingPower); hStoppingPowerE->Fill(it->second.vertexKE,stoppingPower);
  }

  int pe=0,pid=0,parentPho=0,src=0,origin=0,reached=0,extracted=0,nref=0,ntir=0,ntyv=0,nray=0;
  float wl=0,btime=0,angle=0,emitterBeta=std::numeric_limits<float>::quiet_NaN(),wallT=0,wx=0,wy=0,wz=0;
  tPho->SetBranchAddress("EventID",&pe); tPho->SetBranchAddress("PhotonID",&pid);
  tPho->SetBranchAddress("ParentTrackID",&parentPho); tPho->SetBranchAddress("SourceTrackID",&src); tPho->SetBranchAddress("OriginVolumeID",&origin);
  tPho->SetBranchAddress("Wavelength_nm",&wl); tPho->SetBranchAddress("BirthTime_ns",&btime);
  tPho->SetBranchAddress("EmissionAngle_deg",&angle);
  const bool hasEmitterBeta = tPho->GetBranch("EmitterBeta") != nullptr;
  if (hasEmitterBeta) tPho->SetBranchAddress("EmitterBeta",&emitterBeta);
  tPho->SetBranchAddress("ReachedCoreWall",&reached);
  tPho->SetBranchAddress("FirstWallTime_ns",&wallT); tPho->SetBranchAddress("FirstWallX_mm",&wx);
  tPho->SetBranchAddress("FirstWallY_mm",&wy); tPho->SetBranchAddress("FirstWallZ_mm",&wz);
  tPho->SetBranchAddress("ExtractedToReadout",&extracted); tPho->SetBranchAddress("NBoundaryReflections",&nref);
  tPho->SetBranchAddress("NTIR",&ntir); tPho->SetBranchAddress("NTyvekReflections",&ntyv); tPho->SetBranchAddress("NRayleigh",&nray);

  TH1D* hGenWl = new TH1D("hGenWl",Form("Optical photons -- %s;Wavelength (nm);Photons",sourceLabel.Data()),120,200,900);
  TH1D* hCkAngle = new TH1D("hCkAngle",Form("Cherenkov angle -- %s;Emission angle (deg);Photons",sourceLabel.Data()),100,0,90);
  TH2D* hAngleWl = new TH2D("hAngleWl",Form("Cherenkov angle vs wavelength -- %s;Wavelength (nm);Emission angle (deg)",sourceLabel.Data()),120,200,900,100,0,90);
  TH1D* hEmitterBeta = new TH1D("hEmitterBeta",Form("Cerenkov emitter #beta -- %s;Emitter #beta;Photons",sourceLabel.Data()),120,0.95,1.0005);
  TH1D* hAngleResidual = new TH1D("hAngleResidual",Form("Core Cerenkov-angle residual -- %s;#theta_{MC}-#theta_{theory} (deg);Photons",sourceLabel.Data()),240,-0.12,0.12);
  TH1D* hNRef = new TH1D("hNRef","Boundary reflections per selected photon;N reflections;Photons",80,-.5,79.5);
  TH1D* hNTIR = new TH1D("hNTIR","TIR count per selected photon;N TIR;Photons",60,-.5,59.5);
  TH1D* hNTyv = new TH1D("hNTyv","Tyvek reflections per selected photon;N Tyvek reflections;Photons",60,-.5,59.5);
  TH1D* hNRay = new TH1D("hNRay","Rayleigh scatters per selected photon;N Rayleigh;Photons",10,-.5,9.5);
  TH2D* hWallXZ = new TH2D("hWallXZ","First Core-boundary position (projection);X (mm);Z (mm)",120,-100,100,120,-100,100);

  std::unordered_set<uint64_t> selectedPhotons;
  std::unordered_map<uint64_t,float> photonWl;
  long long nGenerated=0,nWall=0,nExtract=0;
  long long nCoreGenerated=0,nCoreWall=0,nCoreExtract=0;
  std::unordered_set<uint64_t> selectedCorePhotons;
  for (Long64_t i=0;i<tPho->GetEntries();++i) {
    tPho->GetEntry(i);
    if (!SelectPhotonSource(tracks,pe,parentPho,src,selectedPDG,combineCharge,sourceMode)) continue;
    const auto k=Key(pe,pid); selectedPhotons.insert(k); photonWl[k]=wl;
    ++nGenerated; if(reached) ++nWall; if(extracted) ++nExtract;
    if(origin==0){ ++nCoreGenerated; selectedCorePhotons.insert(k); if(reached) ++nCoreWall; if(extracted) ++nCoreExtract; }
    hGenWl->Fill(wl);
    if (std::isfinite(angle)) {
      hCkAngle->Fill(angle); hAngleWl->Fill(wl,angle);
      if (hasEmitterBeta && std::isfinite(emitterBeta)) {
        hEmitterBeta->Fill(emitterBeta);
        if (origin==0 && !coreRindex.empty()) {
          const double nr=InterpRindex(coreRindex,wl);
          if (std::isfinite(nr) && emitterBeta*nr>1.0) {
            const double theory=std::acos(1.0/(double(emitterBeta)*nr))*180.0/M_PI;
            hAngleResidual->Fill(angle-theory);
          }
        }
      }
    }
    hNRef->Fill(nref); hNTIR->Fill(ntir); hNTyv->Fill(ntyv); hNRay->Fill(nray);
    if (reached) hWallXZ->Fill(wx,wz);
  }

  int ie=0,ipid=0,ich=0; float iwl=0,it=0,ix=0,iy=0,ia=0;
  tInc->SetBranchAddress("EventID",&ie); tInc->SetBranchAddress("PhotonID",&ipid); tInc->SetBranchAddress("ChannelID",&ich);
  tInc->SetBranchAddress("Wavelength_nm",&iwl); tInc->SetBranchAddress("Time_ns",&it);
  tInc->SetBranchAddress("LocalX_mm",&ix); tInc->SetBranchAddress("LocalY_mm",&iy); tInc->SetBranchAddress("IncidentAngle_deg",&ia);
  TH1D* hIncWl = new TH1D("hIncWl","Readout incident-attempt spectrum;Wavelength (nm);Incident attempts",120,200,900);
  TH1D* hIncUniqueWl = new TH1D("hIncUniqueWl","Photons reaching readout at least once;Wavelength (nm);Unique photons",120,200,900);
  TH1D* hIncAngle = new TH1D("hIncAngle","Readout incident angle;Incident angle (deg);Photons",100,0,90);
  TH2D* hIncWlAngle = new TH2D("hIncWlAngle","Wavelength vs readout incident angle;Incident angle (deg);Wavelength (nm)",100,0,90,120,200,900);
  TH2D* hApertureXY = new TH2D("hApertureXY","Readout entrance local phase space position;Local X (mm);Local Y (mm)",120,-30,30,120,-30,30);
  std::unordered_set<uint64_t> incidentPhotons;
  std::unordered_set<uint64_t> coreIncidentPhotons;
  std::unordered_map<int,double> selectedIncByChannel;
  for(Long64_t i=0;i<tInc->GetEntries();++i){
    tInc->GetEntry(i); auto k=Key(ie,ipid); if(!selectedPhotons.count(k))continue;
    const bool firstIncident = incidentPhotons.insert(k).second;
    if(selectedCorePhotons.count(k)) coreIncidentPhotons.insert(k);
    selectedIncByChannel[ich] += 1.0;
    hIncWl->Fill(iwl);
    if(firstIncident) hIncUniqueWl->Fill(iwl);
    hIncAngle->Fill(ia); hIncWlAngle->Fill(ia,iwl); hApertureXY->Fill(ix,iy);
  }

  // Trigger reference times are loaded before detected-photon timing so the timing plot
  // is relative to the selected top-scintillator trigger pulse, not an arbitrary absolute clock.
  int ge=0,trig=0; float etop=0,ebot=0,ttop=0,tbot=0;
  tTrig->SetBranchAddress("EventID",&ge); tTrig->SetBranchAddress("TriggerPassed",&trig);
  tTrig->SetBranchAddress("EdepTop_MeV",&etop); tTrig->SetBranchAddress("EdepBottom_MeV",&ebot);
  tTrig->SetBranchAddress("TimeTop_ns",&ttop); tTrig->SetBranchAddress("TimeBottom_ns",&tbot);
  std::unordered_map<int,float> topTriggerTime;
  for(Long64_t i=0;i<tTrig->GetEntries();++i){
    tTrig->GetEntry(i);
    if(ttop>=0) topTriggerTime[ge]=ttop;
  }

  int de=0,dpid=0,dch=0; float dt=0;
  tDet->SetBranchAddress("EventID",&de); tDet->SetBranchAddress("PhotonID",&dpid); tDet->SetBranchAddress("ChannelID",&dch); tDet->SetBranchAddress("Time_ns",&dt);
  TH1D* hDetWl = new TH1D("hDetWl","Detected spectrum;Wavelength (nm);Detected photons",120,200,900);
  TH1D* hDetTime = new TH1D("hDetTime","Detected photon delay relative to top trigger;t_{det}-t_{top} (ns);Detected photons",240,-5,55);
  std::unordered_set<uint64_t> detectedPhotons;
  std::unordered_set<uint64_t> coreDetectedPhotons;
  std::unordered_map<int,double> selectedDetByChannel;
  for(Long64_t i=0;i<tDet->GetEntries();++i){
    tDet->GetEntry(i);
    auto k=Key(de,dpid); if(!selectedPhotons.count(k))continue;
    detectedPhotons.insert(k); if(selectedCorePhotons.count(k)) coreDetectedPhotons.insert(k); selectedDetByChannel[dch] += 1.0;
    auto w=photonWl.find(k); if(w!=photonWl.end())hDetWl->Fill(w->second);
    auto tt=topTriggerTime.find(de); if(tt!=topTriggerTime.end()) hDetTime->Fill(dt-tt->second);
  }

  TH1D* hTransport = (TH1D*)hGenWl->Clone("hTransportEfficiency"); hTransport->Reset();
  hTransport->SetTitle("Transport reach probability vs wavelength;Wavelength (nm);N_{unique reach}/N_{generated}");
  TH1D* hSensor = (TH1D*)hGenWl->Clone("hSensorEfficiency"); hSensor->Reset();
  hSensor->SetTitle("Detection probability after reaching readout vs wavelength;Wavelength (nm);N_{detected}/N_{unique reach}");
  TH1D* hTotal = (TH1D*)hGenWl->Clone("hTotalEfficiency"); hTotal->Reset();
  hTotal->SetTitle("Total detection efficiency vs wavelength;Wavelength (nm);N_{detected}/N_{generated}");
  for(int b=1;b<=hGenWl->GetNbinsX();++b){
    const double g=hGenWl->GetBinContent(b);
    const double reach=hIncUniqueWl->GetBinContent(b);
    const double d=hDetWl->GetBinContent(b);
    if(g>0){hTransport->SetBinContent(b,reach/g);hTotal->SetBinContent(b,d/g);}
    if(reach>0)hSensor->SetBinContent(b,d/reach);
  }

  int maxChannel=-1, mapCh=0; tMap->SetBranchAddress("ChannelID",&mapCh);
  for(Long64_t i=0;i<tMap->GetEntries();++i){tMap->GetEntry(i);maxChannel=std::max(maxChannel,mapCh);}
  TH1D* hOccDet = new TH1D("hOccDet",Form("Detected occupancy by channel -- %s;Channel ID;Detected photons",sourceLabel.Data()),std::max(1,maxChannel+1),-.5,std::max(.5,maxChannel+.5));
  TH1D* hOccInc = new TH1D("hOccInc",Form("Incident occupancy by channel -- %s;Channel ID;Incident photons",sourceLabel.Data()),std::max(1,maxChannel+1),-.5,std::max(.5,maxChannel+.5));
  for(const auto& kv : selectedIncByChannel) hOccInc->SetBinContent(kv.first+1,kv.second);
  for(const auto& kv : selectedDetByChannel) hOccDet->SetBinContent(kv.first+1,kv.second);

  TH2D* hTrigE = new TH2D("hTrigE","Total event energy deposited in trigger scintillators;Top total Edep (MeV);Bottom total Edep (MeV)",120,0,20,120,0,20);
  TH1D* hTrigDt = new TH1D("hTrigDt","Top-bottom scintillator time difference;t_{bottom}-t_{top} (ns);Events",200,-100,100);
  long long nAll=0,nTrig=0;
  for(Long64_t i=0;i<tTrig->GetEntries();++i){tTrig->GetEntry(i);++nAll;if(trig)++nTrig;hTrigE->Fill(etop,ebot);if(ttop>=0&&tbot>=0)hTrigDt->Fill(tbot-ttop);}

  TH1D* hFunnel = new TH1D("hFunnel",Form("Core-origin transport funnel -- %s;Stage;Photons",sourceLabel.Data()),5,0,5);
  hFunnel->GetXaxis()->SetBinLabel(1,"Generated in Core");hFunnel->GetXaxis()->SetBinLabel(2,"Reach Core boundary");hFunnel->GetXaxis()->SetBinLabel(3,"Extracted to readout");hFunnel->GetXaxis()->SetBinLabel(4,"Readout incident");hFunnel->GetXaxis()->SetBinLabel(5,"Detected");
  hFunnel->SetBinContent(1,nCoreGenerated);hFunnel->SetBinContent(2,nCoreWall);hFunnel->SetBinContent(3,nCoreExtract);hFunnel->SetBinContent(4,coreIncidentPhotons.size());hFunnel->SetBinContent(5,coreDetectedPhotons.size());

  TCanvas c("c","General Cherenkov analysis",1200,700); c.Print(pdf+"[");
  c.Clear(); TLatex tx; tx.SetNDC();tx.SetTextSize(.035);tx.DrawLatex(.08,.90,Form("General analysis: %s",sourceLabel.Data()));tx.DrawLatex(.08,.83,Form("Input: %s",inputName));tx.DrawLatex(.08,.76,Form("Generated selected photons: %lld",nGenerated));tx.DrawLatex(.08,.69,Form("Trigger rows: %lld, passed: %lld (%.3f%%)",nAll,nTrig,nAll?100.0*nTrig/nAll:0));tx.DrawLatex(.08,.62,Form("Readout-incident selected photons: %zu",incidentPhotons.size()));tx.DrawLatex(.08,.55,Form("Detected selected photons: %zu",detectedPhotons.size()));c.Print(pdf);
  c.Clear();hPrimaryE->Draw();c.Print(pdf); c.Clear();hZenith->Draw();c.Print(pdf);
  c.Clear();hCorePath->Draw();c.Print(pdf);
  c.Clear();hLocalEdepDx->Draw();c.Print(pdf);
  c.Clear();c.SetLogy(1);hStoppingPower->Draw();c.Print(pdf);c.SetLogy(0);
  c.Clear();hStoppingPowerE->Draw("COLZ");c.Print(pdf);
  c.Clear();hGenWl->Draw();c.Print(pdf); c.Clear();hCkAngle->Draw();c.Print(pdf); c.Clear();hAngleWl->Draw("COLZ");c.Print(pdf);
  if (hasEmitterBeta) { c.Clear();hEmitterBeta->Draw();c.Print(pdf); }
  if (hAngleResidual->GetEntries()>0) { c.Clear();hAngleResidual->Draw(); TLine zero(0,0,0,hAngleResidual->GetMaximum()*1.05); zero.SetLineStyle(2); zero.Draw("SAME"); c.Print(pdf); }
  c.Clear();hFunnel->Draw("HIST");c.Print(pdf);
  c.Clear();hIncWl->Draw();c.Print(pdf); c.Clear();hIncUniqueWl->Draw();c.Print(pdf); c.Clear();hIncAngle->Draw();c.Print(pdf); c.Clear();hIncWlAngle->Draw("COLZ");c.Print(pdf); c.Clear();hApertureXY->Draw("COLZ");c.Print(pdf);
  c.Clear();hDetWl->Draw();c.Print(pdf); c.Clear();hDetTime->Draw();c.Print(pdf);
  c.Clear();hTransport->SetMaximum(1.05);hTransport->Draw();c.Print(pdf); c.Clear();hSensor->SetMaximum(1.05);hSensor->Draw();c.Print(pdf); c.Clear();hTotal->SetMaximum(1.05);hTotal->Draw();c.Print(pdf);
  c.Clear();hOccInc->Draw("HIST");c.Print(pdf); c.Clear();hOccDet->Draw("HIST");c.Print(pdf);
  c.Clear();hNRef->Draw();c.Print(pdf);c.Clear();hNTIR->Draw();c.Print(pdf);c.Clear();hNTyv->Draw();c.Print(pdf);c.Clear();hNRay->Draw();c.Print(pdf);c.Clear();hWallXZ->Draw("COLZ");c.Print(pdf);
  c.Clear();hTrigE->Draw("COLZ");c.Print(pdf);c.Clear();hTrigDt->Draw();c.Print(pdf);
  if (std::abs(selectedPDG)==13) {
    c.Clear();
    if (selectedPDG==-13 && !combineCharge)
      hLifetime->SetTitle("#mu^{+} decay proper time (free-decay validation);Proper time (#mus);Decays");
    else
      hLifetime->SetTitle("Selected-muon decay proper time (no free-lifetime fit for mixed/negative muons);Proper time (#mus);Decays");
    hLifetime->Draw();
    // Only auto-fit the free lifetime for an explicit mu+ selection.  mu- in matter can
    // disappear by nuclear capture, so a mixed-charge or mu- sample must not be fit
    // to the free-muon lifetime as if it were the same observable.
    if(selectedPDG==-13 && !combineCharge && hLifetime->GetEntries()>100){
      TF1 expo("muLife","[0]*exp(-x/[1])",0.2,10);
      expo.SetParameters(hLifetime->GetMaximum(),2.2);
      hLifetime->Fit(&expo,"Q","",0.2,10); expo.Draw("SAME");
    }
    c.Print(pdf);
    c.Clear();hMichel->Draw();c.Print(pdf);
  }
  c.Print(pdf+"]");
  std::cout << "Integrated analysis summary:\n"
            << "  Core path mean = " << hCorePath->GetMean() << " cm\n"
            << "  Local Edep/L mean = " << hLocalEdepDx->GetMean() << " MeV/cm\n"
            << "  DeltaKE/L histogram in-range mean = " << hStoppingPower->GetMean() << " MeV/cm\n";
  if (hAngleResidual->GetEntries()>0)
    std::cout << "  Core Cherenkov angle residual = " << hAngleResidual->GetMean()
              << " deg RMS=" << hAngleResidual->GetRMS()
              << " entries=" << hAngleResidual->GetEntries() << "\n";
  std::cout << "Wrote " << pdf << "\n";
  f->Close();
}
