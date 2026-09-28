#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <algorithm>
#include <cmath>

EventAction::EventAction(RunAction* runAction)
  : fRunAction(runAction) {}

void EventAction::BeginOfEventAction(const G4Event*) {
  fTracks.clear();
  fTrackVolumes.clear();
  fPhotons.clear();
  fRootPrimaryByAnyTrack.clear();
  fSourceTrackByAnyTrack.clear();
  fSensorIncidents.clear();
  fDetections.clear();
  fBoundaryDebug.clear();
  fCoreBoundaryCrossings.clear();
  fCoreBoundaryRecorded.clear();
  fChannels.clear();
  fScintTopDeposits.clear();
  fScintBottomDeposits.clear();

  fEdepTopMeV = 0.0;
  fEdepBottomMeV = 0.0;
  fTriggerTopTimeNs = -1.0;
  fTriggerBottomTimeNs = -1.0;
  fTriggerPassed = false;
}

G4int EventAction::ResolveRootPrimary(G4int trackID) const {
  auto it = fRootPrimaryByAnyTrack.find(trackID);
  return it == fRootPrimaryByAnyTrack.end() ? 0 : it->second;
}

G4int EventAction::ResolveSourceTrack(G4int trackID) const {
  auto it = fSourceTrackByAnyTrack.find(trackID);
  return it == fSourceTrackByAnyTrack.end() ? 0 : it->second;
}

TrackRecord& EventAction::RegisterTrack(const G4Track* track, G4int birthVolumeID) {
  const G4int id = track->GetTrackID();
  auto found = fTracks.find(id);
  if (found != fTracks.end()) return found->second;

  TrackRecord r;
  r.trackID = id;
  r.parentID = track->GetParentID();
  r.rootPrimaryTrackID = (r.parentID == 0) ? id : ResolveRootPrimary(r.parentID);

  const auto* def = track->GetDefinition();
  r.pdg = def ? def->GetPDGEncoding() : 0;
  r.chargeE = def ? static_cast<G4float>(def->GetPDGCharge()/eplus) : 0.0f;

  const auto* creator = track->GetCreatorProcess();
  if (creator) {
    r.creatorProcessType = creator->GetProcessType();
    r.creatorProcessSubType = creator->GetProcessSubType();
    r.creatorProcessName = creator->GetProcessName();
  }

  r.birthVolumeID = birthVolumeID;
  r.vertexKE_MeV = static_cast<G4float>(track->GetVertexKineticEnergy()/MeV);
  r.vertexTime_ns = static_cast<G4float>(track->GetGlobalTime()/ns);
  const auto vpos = track->GetVertexPosition();
  const auto vdir = track->GetVertexMomentumDirection();
  r.vertexX_mm = static_cast<G4float>(vpos.x()/mm);
  r.vertexY_mm = static_cast<G4float>(vpos.y()/mm);
  r.vertexZ_mm = static_cast<G4float>(vpos.z()/mm);
  r.dirX = static_cast<G4float>(vdir.x());
  r.dirY = static_cast<G4float>(vdir.y());
  r.dirZ = static_cast<G4float>(vdir.z());

  fRootPrimaryByAnyTrack[id] = r.rootPrimaryTrackID;
  fSourceTrackByAnyTrack[id] = id;

  auto inserted = fTracks.emplace(id, r);
  return inserted.first->second;
}

void EventAction::UpdateTrackEnd(const G4Track* track, G4int endVolumeID,
                                 const G4VProcess* endProcess) {
  auto it = fTracks.find(track->GetTrackID());
  if (it == fTracks.end()) return;
  auto& r = it->second;
  r.endTime_ns = static_cast<G4float>(track->GetGlobalTime()/ns);
  r.endProperTime_ns = static_cast<G4float>(track->GetProperTime()/ns);
  const auto pos = track->GetPosition();
  r.endX_mm = static_cast<G4float>(pos.x()/mm);
  r.endY_mm = static_cast<G4float>(pos.y()/mm);
  r.endZ_mm = static_cast<G4float>(pos.z()/mm);
  r.endVolumeID = endVolumeID;
  r.endTrackStatus = static_cast<G4int>(track->GetTrackStatus());
  if (endProcess) {
    r.endProcessType = endProcess->GetProcessType();
    r.endProcessSubType = endProcess->GetProcessSubType();
    r.endProcessName = endProcess->GetProcessName();
  }
}

void EventAction::AddTrackVolume(G4int trackID, G4int volumeID,
                                 G4double edep, G4double deltaKE, G4double length) {
  if (volumeID < 0) return;
  const auto key = std::make_pair(trackID, volumeID);
  auto& r = fTrackVolumes[key];
  r.trackID = trackID;
  r.volumeID = volumeID;
  r.edep_MeV += static_cast<G4float>(edep/MeV);
  r.deltaKE_MeV += static_cast<G4float>(deltaKE/MeV);
  r.length_mm += static_cast<G4float>(length/mm);
  r.nSteps += 1;
}

void EventAction::SetCoreEntry(G4int trackID, const G4ThreeVector& pos) {
  auto it = fTracks.find(trackID);
  if (it == fTracks.end() || it->second.coreEntrySet) return;
  it->second.coreEntrySet = 1;
  it->second.coreEntryX_mm = static_cast<G4float>(pos.x()/mm);
  it->second.coreEntryY_mm = static_cast<G4float>(pos.y()/mm);
  it->second.coreEntryZ_mm = static_cast<G4float>(pos.z()/mm);
}

void EventAction::SetCoreExit(G4int trackID, const G4ThreeVector& pos) {
  auto it = fTracks.find(trackID);
  if (it == fTracks.end()) return;
  it->second.coreExitSet = 1;
  it->second.coreExitX_mm = static_cast<G4float>(pos.x()/mm);
  it->second.coreExitY_mm = static_cast<G4float>(pos.y()/mm);
  it->second.coreExitZ_mm = static_cast<G4float>(pos.z()/mm);
}

PhotonRecord& EventAction::RegisterPhoton(const G4Track* track, G4int originVolumeID,
                                          G4double wavelengthNm,
                                          G4double emissionAngleDeg,
                                          G4double emitterBeta,
                                          G4double birthGlobalTime,
                                          const G4ThreeVector& birthPosition,
                                          const G4ThreeVector& birthDirection) {
  const G4int id = track->GetTrackID();
  auto found = fPhotons.find(id);
  if (found != fPhotons.end()) return found->second;

  PhotonRecord p;
  p.photonID = id;
  p.parentTrackID = track->GetParentID();
  p.sourceTrackID = ResolveSourceTrack(p.parentTrackID);
  p.originVolumeID = originVolumeID;
  p.wavelength_nm = static_cast<G4float>(wavelengthNm);
  // RegisterPhoton is called from optical step #1.  At UserSteppingAction time
  // the G4Track already represents the transported step state, so using
  // track->GetGlobalTime()/GetPosition() here would record the END of the first
  // optical step, not the photon creation state.  The caller therefore passes
  // the first-step PRE point explicitly.
  p.birthTime_ns = static_cast<G4float>(birthGlobalTime/ns);
  p.emissionAngle_deg = static_cast<G4float>(emissionAngleDeg);
  p.emitterBeta = static_cast<G4float>(emitterBeta);

  const auto* creator = track->GetCreatorProcess();
  if (creator) {
    p.creatorProcessType = creator->GetProcessType();
    p.creatorProcessSubType = creator->GetProcessSubType();
  }

  const auto* det = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  if (det && det->GetPhotonBirthDetail() == "full") {
    p.birthX_mm = static_cast<G4float>(birthPosition.x()/mm);
    p.birthY_mm = static_cast<G4float>(birthPosition.y()/mm);
    p.birthZ_mm = static_cast<G4float>(birthPosition.z()/mm);
    p.birthDirX = static_cast<G4float>(birthDirection.x());
    p.birthDirY = static_cast<G4float>(birthDirection.y());
    p.birthDirZ = static_cast<G4float>(birthDirection.z());
  }

  const G4int root = ResolveRootPrimary(p.parentTrackID);
  fRootPrimaryByAnyTrack[id] = root;
  fSourceTrackByAnyTrack[id] = p.sourceTrackID;

  auto inserted = fPhotons.emplace(id, p);
  return inserted.first->second;
}

PhotonRecord* EventAction::GetPhoton(G4int photonID) {
  auto it = fPhotons.find(photonID);
  return it == fPhotons.end() ? nullptr : &it->second;
}

void EventAction::AddSensorIncident(const SensorIncidentRecord& hit) {
  fSensorIncidents.push_back(hit);
  auto& c = fChannels[hit.channelID];
  c.nIncident += 1;
  if (std::isnan(c.firstIncidentTime_ns) || hit.time_ns < c.firstIncidentTime_ns) {
    c.firstIncidentTime_ns = hit.time_ns;
  }
}

void EventAction::AddDetection(const DetectionRecord& hit) {
  fDetections.push_back(hit);
  auto& c = fChannels[hit.channelID];
  c.nDetected += 1;
  if (std::isnan(c.firstDetectedTime_ns) || hit.time_ns < c.firstDetectedTime_ns) {
    c.firstDetectedTime_ns = hit.time_ns;
  }
}

void EventAction::AddBoundaryDebug(const BoundaryDebugRecord& hit) {
  fBoundaryDebug.push_back(hit);
}

G4bool EventAction::AddCoreBoundaryCrossing(
    const CoreBoundaryCrossingRecord& hit) {
  if (hit.photonID <= 0) return false;
  if (!fCoreBoundaryRecorded.insert(hit.photonID).second) return false;
  fCoreBoundaryCrossings.push_back(hit);
  return true;
}

void EventAction::AddScintTop(G4int trackID, G4double edepMeV, G4double timeNs) {
  fEdepTopMeV += edepMeV;
  fScintTopDeposits.push_back({trackID,edepMeV,timeNs});
}

void EventAction::AddScintBottom(G4int trackID, G4double edepMeV, G4double timeNs) {
  fEdepBottomMeV += edepMeV;
  fScintBottomDeposits.push_back({trackID,edepMeV,timeNs});
}

void EventAction::EndOfEventAction(const G4Event* event) {
  const auto* det = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  const G4double threshold = det ? det->GetScintThreshold()/MeV : 0.0;
  const G4double coincidence = det ? det->GetCoincidenceWindow()/ns : 100.0;
  const G4double integration = det ? det->GetScintIntegrationWindow()/ns : 100.0;

  struct Pulse { G4double t=0.0; G4double e=0.0; };
  auto buildPulses = [integration](std::vector<ScintDeposit> deps) {
    std::vector<Pulse> pulses;
    if (deps.empty()) return pulses;
    std::sort(deps.begin(), deps.end(),
              [](const ScintDeposit& a, const ScintDeposit& b){ return a.timeNs < b.timeNs; });
    Pulse cur{deps.front().timeNs,deps.front().edepMeV};
    for (std::size_t i=1; i<deps.size(); ++i) {
      if (deps[i].timeNs - cur.t <= integration) cur.e += deps[i].edepMeV;
      else {
        pulses.push_back(cur);
        cur = {deps[i].timeNs,deps[i].edepMeV};
      }
    }
    pulses.push_back(cur);
    return pulses;
  };

  const auto top = buildPulses(fScintTopDeposits);
  const auto bottom = buildPulses(fScintBottomDeposits);
  G4double bestTriggerTime = 1.0e300;
  G4double bestDt = 1.0e300;
  for (const auto& t : top) {
    if (t.e < threshold) continue;
    for (const auto& b : bottom) {
      if (b.e < threshold) continue;
      const G4double dt = std::abs(t.t-b.t);
      if (dt > coincidence) continue;
      const G4double triggerTime = std::max(t.t,b.t);
      if (triggerTime < bestTriggerTime ||
          (triggerTime == bestTriggerTime && dt < bestDt)) {
        bestTriggerTime = triggerTime;
        bestDt = dt;
        fTriggerPassed = true;
        fTriggerTopTimeNs = t.t;
        fTriggerBottomTimeNs = b.t;
      }
    }
  }

  const G4String outputMode = det ? det->GetOutputEventMode() : "thinAll";
  const G4bool details = (outputMode == "fullAll") || fTriggerPassed;
  const G4bool writeThin = (outputMode != "triggered") || fTriggerPassed;

  if (fRunAction) fRunAction->CountEvent(fTriggerPassed, details);

  const G4int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
  const G4int eventID = event->GetEventID();
  auto* am = G4AnalysisManager::Instance();

  G4int nPrimary = 0;
  for (const auto& kv : fTracks) if (kv.second.parentID == 0) ++nPrimary;

  G4float coreEdep = 0.0f;
  for (const auto& kv : fTrackVolumes) {
    if (kv.second.volumeID == 0) coreEdep += kv.second.edep_MeV;
  }

  if (writeThin) {
    for (const auto& d : fScintTopDeposits) {
      am->FillNtupleIColumn(kScintDeposit,0,runID);
      am->FillNtupleIColumn(kScintDeposit,1,eventID);
      am->FillNtupleIColumn(kScintDeposit,2,0);
      am->FillNtupleIColumn(kScintDeposit,3,d.trackID);
      am->FillNtupleFColumn(kScintDeposit,4,static_cast<G4float>(d.timeNs));
      am->FillNtupleFColumn(kScintDeposit,5,static_cast<G4float>(d.edepMeV));
      am->AddNtupleRow(kScintDeposit);
    }
    for (const auto& d : fScintBottomDeposits) {
      am->FillNtupleIColumn(kScintDeposit,0,runID);
      am->FillNtupleIColumn(kScintDeposit,1,eventID);
      am->FillNtupleIColumn(kScintDeposit,2,1);
      am->FillNtupleIColumn(kScintDeposit,3,d.trackID);
      am->FillNtupleFColumn(kScintDeposit,4,static_cast<G4float>(d.timeNs));
      am->FillNtupleFColumn(kScintDeposit,5,static_cast<G4float>(d.edepMeV));
      am->AddNtupleRow(kScintDeposit);
    }

    for (const auto& kv : fTracks) {
      const auto& r = kv.second;
      if (r.parentID != 0) continue;
      am->FillNtupleIColumn(kPrimaryTruth,0,runID);
      am->FillNtupleIColumn(kPrimaryTruth,1,eventID);
      am->FillNtupleIColumn(kPrimaryTruth,2,r.trackID);
      am->FillNtupleIColumn(kPrimaryTruth,3,r.pdg);
      am->FillNtupleFColumn(kPrimaryTruth,4,r.vertexKE_MeV);
      am->FillNtupleFColumn(kPrimaryTruth,5,r.vertexTime_ns);
      am->FillNtupleFColumn(kPrimaryTruth,6,r.vertexX_mm);
      am->FillNtupleFColumn(kPrimaryTruth,7,r.vertexY_mm);
      am->FillNtupleFColumn(kPrimaryTruth,8,r.vertexZ_mm);
      am->FillNtupleFColumn(kPrimaryTruth,9,r.dirX);
      am->FillNtupleFColumn(kPrimaryTruth,10,r.dirY);
      am->FillNtupleFColumn(kPrimaryTruth,11,r.dirZ);
      am->AddNtupleRow(kPrimaryTruth);
    }

    am->FillNtupleIColumn(kTriggerEvent,0,runID);
    am->FillNtupleIColumn(kTriggerEvent,1,eventID);
    am->FillNtupleIColumn(kTriggerEvent,2,fTriggerPassed?1:0);
    am->FillNtupleFColumn(kTriggerEvent,3,static_cast<G4float>(fEdepTopMeV));
    am->FillNtupleFColumn(kTriggerEvent,4,static_cast<G4float>(fEdepBottomMeV));
    am->FillNtupleFColumn(kTriggerEvent,5,static_cast<G4float>(fTriggerTopTimeNs));
    am->FillNtupleFColumn(kTriggerEvent,6,static_cast<G4float>(fTriggerBottomTimeNs));
    am->FillNtupleIColumn(kTriggerEvent,7,nPrimary);
    am->FillNtupleIColumn(kTriggerEvent,8,static_cast<G4int>(fTracks.size()));
    am->FillNtupleIColumn(kTriggerEvent,9,static_cast<G4int>(fPhotons.size()));
    am->FillNtupleIColumn(kTriggerEvent,10,static_cast<G4int>(fSensorIncidents.size()));
    am->FillNtupleIColumn(kTriggerEvent,11,static_cast<G4int>(fDetections.size()));
    am->FillNtupleFColumn(kTriggerEvent,12,coreEdep);
    am->AddNtupleRow(kTriggerEvent);
  }

  if (!details) return;

  for (const auto& kv : fTracks) {
    const auto& r = kv.second;
    am->FillNtupleIColumn(kTrack,0,runID);
    am->FillNtupleIColumn(kTrack,1,eventID);
    am->FillNtupleIColumn(kTrack,2,r.trackID);
    am->FillNtupleIColumn(kTrack,3,r.parentID);
    am->FillNtupleIColumn(kTrack,4,r.rootPrimaryTrackID);
    am->FillNtupleIColumn(kTrack,5,r.pdg);
    am->FillNtupleFColumn(kTrack,6,r.chargeE);
    am->FillNtupleIColumn(kTrack,7,r.creatorProcessType);
    am->FillNtupleIColumn(kTrack,8,r.creatorProcessSubType);
    am->FillNtupleSColumn(kTrack,9,r.creatorProcessName);
    am->FillNtupleIColumn(kTrack,10,r.birthVolumeID);
    am->FillNtupleFColumn(kTrack,11,r.vertexKE_MeV);
    am->FillNtupleFColumn(kTrack,12,r.vertexTime_ns);
    am->FillNtupleFColumn(kTrack,13,r.vertexX_mm);
    am->FillNtupleFColumn(kTrack,14,r.vertexY_mm);
    am->FillNtupleFColumn(kTrack,15,r.vertexZ_mm);
    am->FillNtupleFColumn(kTrack,16,r.dirX);
    am->FillNtupleFColumn(kTrack,17,r.dirY);
    am->FillNtupleFColumn(kTrack,18,r.dirZ);
    am->FillNtupleIColumn(kTrack,19,r.coreEntrySet);
    am->FillNtupleFColumn(kTrack,20,r.coreEntryX_mm);
    am->FillNtupleFColumn(kTrack,21,r.coreEntryY_mm);
    am->FillNtupleFColumn(kTrack,22,r.coreEntryZ_mm);
    am->FillNtupleIColumn(kTrack,23,r.coreExitSet);
    am->FillNtupleFColumn(kTrack,24,r.coreExitX_mm);
    am->FillNtupleFColumn(kTrack,25,r.coreExitY_mm);
    am->FillNtupleFColumn(kTrack,26,r.coreExitZ_mm);
    am->FillNtupleFColumn(kTrack,27,r.endTime_ns);
    am->FillNtupleFColumn(kTrack,28,r.endProperTime_ns);
    am->FillNtupleFColumn(kTrack,29,r.endX_mm);
    am->FillNtupleFColumn(kTrack,30,r.endY_mm);
    am->FillNtupleFColumn(kTrack,31,r.endZ_mm);
    am->FillNtupleIColumn(kTrack,32,r.endVolumeID);
    am->FillNtupleIColumn(kTrack,33,r.endProcessType);
    am->FillNtupleIColumn(kTrack,34,r.endProcessSubType);
    am->FillNtupleSColumn(kTrack,35,r.endProcessName);
    am->FillNtupleIColumn(kTrack,36,r.endTrackStatus);
    am->AddNtupleRow(kTrack);
  }

  for (const auto& kv : fTrackVolumes) {
    const auto& r = kv.second;
    am->FillNtupleIColumn(kTrackVolume,0,runID);
    am->FillNtupleIColumn(kTrackVolume,1,eventID);
    am->FillNtupleIColumn(kTrackVolume,2,r.trackID);
    am->FillNtupleIColumn(kTrackVolume,3,r.volumeID);
    am->FillNtupleFColumn(kTrackVolume,4,r.edep_MeV);
    am->FillNtupleFColumn(kTrackVolume,5,r.deltaKE_MeV);
    am->FillNtupleFColumn(kTrackVolume,6,r.length_mm);
    am->FillNtupleIColumn(kTrackVolume,7,r.nSteps);
    am->AddNtupleRow(kTrackVolume);
  }

  for (const auto& kv : fPhotons) {
    const auto& p = kv.second;
    am->FillNtupleIColumn(kPhoton,0,runID);
    am->FillNtupleIColumn(kPhoton,1,eventID);
    am->FillNtupleIColumn(kPhoton,2,p.photonID);
    am->FillNtupleIColumn(kPhoton,3,p.parentTrackID);
    am->FillNtupleIColumn(kPhoton,4,p.sourceTrackID);
    am->FillNtupleIColumn(kPhoton,5,p.creatorProcessType);
    am->FillNtupleIColumn(kPhoton,6,p.creatorProcessSubType);
    am->FillNtupleIColumn(kPhoton,7,p.originVolumeID);
    am->FillNtupleFColumn(kPhoton,8,p.wavelength_nm);
    am->FillNtupleFColumn(kPhoton,9,p.birthTime_ns);
    am->FillNtupleFColumn(kPhoton,10,p.emissionAngle_deg);
    am->FillNtupleFColumn(kPhoton,11,p.birthX_mm);
    am->FillNtupleFColumn(kPhoton,12,p.birthY_mm);
    am->FillNtupleFColumn(kPhoton,13,p.birthZ_mm);
    am->FillNtupleFColumn(kPhoton,14,p.birthDirX);
    am->FillNtupleFColumn(kPhoton,15,p.birthDirY);
    am->FillNtupleFColumn(kPhoton,16,p.birthDirZ);
    am->FillNtupleIColumn(kPhoton,17,p.reachedCoreWall);
    am->FillNtupleFColumn(kPhoton,18,p.firstWallTime_ns);
    am->FillNtupleFColumn(kPhoton,19,p.firstWallX_mm);
    am->FillNtupleFColumn(kPhoton,20,p.firstWallY_mm);
    am->FillNtupleFColumn(kPhoton,21,p.firstWallZ_mm);
    am->FillNtupleIColumn(kPhoton,22,p.extractedToReadout);
    am->FillNtupleFColumn(kPhoton,23,p.extractionTime_ns);
    am->FillNtupleIColumn(kPhoton,24,p.nBoundaryReflections);
    am->FillNtupleIColumn(kPhoton,25,p.nTIR);
    am->FillNtupleIColumn(kPhoton,26,p.nTyvekReflections);
    am->FillNtupleIColumn(kPhoton,27,p.nRayleigh);
    am->FillNtupleIColumn(kPhoton,28,p.endReason);
    am->FillNtupleIColumn(kPhoton,29,p.endVolumeID);
    am->FillNtupleFColumn(kPhoton,30,p.endTime_ns);
    am->FillNtupleIColumn(kPhoton,31,p.killedBySafetyCut);
    am->FillNtupleFColumn(kPhoton,32,p.firstRayleighTrackLength_mm);
    am->FillNtupleFColumn(kPhoton,33,p.endTrackLength_mm);
    am->FillNtupleIColumn(kPhoton,34,p.acrylicAbsOutOfRange);
    am->FillNtupleFColumn(kPhoton,35,p.emitterBeta);
    am->AddNtupleRow(kPhoton);
  }

  for (const auto& h : fSensorIncidents) {
    am->FillNtupleIColumn(kSensorIncident,0,runID);
    am->FillNtupleIColumn(kSensorIncident,1,eventID);
    am->FillNtupleIColumn(kSensorIncident,2,h.photonID);
    am->FillNtupleIColumn(kSensorIncident,3,h.channelID);
    am->FillNtupleFColumn(kSensorIncident,4,h.wavelength_nm);
    am->FillNtupleFColumn(kSensorIncident,5,h.time_ns);
    am->FillNtupleFColumn(kSensorIncident,6,h.localX_mm);
    am->FillNtupleFColumn(kSensorIncident,7,h.localY_mm);
    am->FillNtupleFColumn(kSensorIncident,8,h.localZ_mm);
    am->FillNtupleFColumn(kSensorIncident,9,h.localDirX);
    am->FillNtupleFColumn(kSensorIncident,10,h.localDirY);
    am->FillNtupleFColumn(kSensorIncident,11,h.localDirZ);
    am->FillNtupleFColumn(kSensorIncident,12,h.localPolX);
    am->FillNtupleFColumn(kSensorIncident,13,h.localPolY);
    am->FillNtupleFColumn(kSensorIncident,14,h.localPolZ);
    am->FillNtupleFColumn(kSensorIncident,15,h.incidentAngle_deg);
    am->FillNtupleIColumn(kSensorIncident,16,h.nBoundaryReflections);
    am->FillNtupleIColumn(kSensorIncident,17,h.nTIR);
    am->FillNtupleIColumn(kSensorIncident,18,h.nTyvekReflections);
    am->FillNtupleIColumn(kSensorIncident,19,h.nRayleigh);
    am->FillNtupleFColumn(kSensorIncident,20,h.entranceX_mm);
    am->FillNtupleFColumn(kSensorIncident,21,h.entranceY_mm);
    am->FillNtupleFColumn(kSensorIncident,22,h.entranceZ_mm);
    am->FillNtupleIColumn(kSensorIncident,23,h.forward);
    am->AddNtupleRow(kSensorIncident);
  }

  for (const auto& d : fDetections) {
    am->FillNtupleIColumn(kDetectedSignal,0,runID);
    am->FillNtupleIColumn(kDetectedSignal,1,eventID);
    am->FillNtupleIColumn(kDetectedSignal,2,d.photonID);
    am->FillNtupleIColumn(kDetectedSignal,3,d.channelID);
    am->FillNtupleFColumn(kDetectedSignal,4,d.time_ns);
    am->AddNtupleRow(kDetectedSignal);
  }

  for (const auto& h : fCoreBoundaryCrossings) {
    am->FillNtupleIColumn(kCoreBoundaryCrossing,0,runID);
    am->FillNtupleIColumn(kCoreBoundaryCrossing,1,eventID);
    am->FillNtupleIColumn(kCoreBoundaryCrossing,2,h.photonID);
    am->FillNtupleIColumn(kCoreBoundaryCrossing,3,h.parentTrackID);
    am->FillNtupleIColumn(kCoreBoundaryCrossing,4,h.sourceTrackID);
    am->FillNtupleIColumn(kCoreBoundaryCrossing,5,h.rootPrimaryTrackID);
    am->FillNtupleIColumn(kCoreBoundaryCrossing,6,h.originVolumeID);
    am->FillNtupleIColumn(kCoreBoundaryCrossing,7,h.faceID);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,8,h.x_mm);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,9,h.y_mm);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,10,h.z_mm);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,11,h.localU_mm);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,12,h.localV_mm);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,13,h.dirX);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,14,h.dirY);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,15,h.dirZ);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,16,h.cosThetaNormal);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,17,h.wavelength_nm);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,18,h.time_ns);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,19,h.polX);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,20,h.polY);
    am->FillNtupleFColumn(kCoreBoundaryCrossing,21,h.polZ);
    am->AddNtupleRow(kCoreBoundaryCrossing);
  }

  for (const auto& kv : fChannels) {
    const auto& c = kv.second;
    am->FillNtupleIColumn(kChannelEvent,0,runID);
    am->FillNtupleIColumn(kChannelEvent,1,eventID);
    am->FillNtupleIColumn(kChannelEvent,2,kv.first);
    am->FillNtupleIColumn(kChannelEvent,3,c.nIncident);
    am->FillNtupleIColumn(kChannelEvent,4,c.nDetected);
    am->FillNtupleFColumn(kChannelEvent,5,c.firstIncidentTime_ns);
    am->FillNtupleFColumn(kChannelEvent,6,c.firstDetectedTime_ns);
    am->AddNtupleRow(kChannelEvent);
  }

  if (det && det->GetBoundaryDebug()) {
    for (const auto& b : fBoundaryDebug) {
      am->FillNtupleIColumn(kBoundaryDebug,0,runID);
      am->FillNtupleIColumn(kBoundaryDebug,1,eventID);
      am->FillNtupleIColumn(kBoundaryDebug,2,b.photonID);
      am->FillNtupleIColumn(kBoundaryDebug,3,b.stepNumber);
      am->FillNtupleIColumn(kBoundaryDebug,4,b.preVolumeID);
      am->FillNtupleIColumn(kBoundaryDebug,5,b.postVolumeID);
      am->FillNtupleIColumn(kBoundaryDebug,6,b.boundaryStatus);
      am->FillNtupleFColumn(kBoundaryDebug,7,b.time_ns);
      am->FillNtupleFColumn(kBoundaryDebug,8,b.x_mm);
      am->FillNtupleFColumn(kBoundaryDebug,9,b.y_mm);
      am->FillNtupleFColumn(kBoundaryDebug,10,b.z_mm);
      am->FillNtupleFColumn(kBoundaryDebug,11,b.preDirX);
      am->FillNtupleFColumn(kBoundaryDebug,12,b.preDirY);
      am->FillNtupleFColumn(kBoundaryDebug,13,b.preDirZ);
      am->FillNtupleFColumn(kBoundaryDebug,14,b.postDirX);
      am->FillNtupleFColumn(kBoundaryDebug,15,b.postDirY);
      am->FillNtupleFColumn(kBoundaryDebug,16,b.postDirZ);
      am->FillNtupleFColumn(kBoundaryDebug,17,b.prePolX);
      am->FillNtupleFColumn(kBoundaryDebug,18,b.prePolY);
      am->FillNtupleFColumn(kBoundaryDebug,19,b.prePolZ);
      am->FillNtupleFColumn(kBoundaryDebug,20,b.postPolX);
      am->FillNtupleFColumn(kBoundaryDebug,21,b.postPolY);
      am->FillNtupleFColumn(kBoundaryDebug,22,b.postPolZ);
      am->FillNtupleSColumn(kBoundaryDebug,23,b.boundaryStatusName);
      am->AddNtupleRow(kBoundaryDebug);
    }
  }
}
