#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4LogicalVolume.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

TrackingAction::TrackingAction(EventAction* eventAction)
  : fEventAction(eventAction) {}

const DetectorConstruction* TrackingAction::GetDetector() const {
  auto* rm = G4RunManager::GetRunManager();
  return rm ? static_cast<const DetectorConstruction*>(
                  rm->GetUserDetectorConstruction())
            : nullptr;
}

void TrackingAction::PreUserTrackingAction(const G4Track* track) {
  if (!track || !fEventAction) return;
  if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) return;

  const auto* detector = GetDetector();
  const auto* pv = track->GetVolume();
  const auto* lv = pv ? pv->GetLogicalVolume() : nullptr;
  const G4int volumeID = detector ? detector->GetVolumeID(lv) : -1;
  fEventAction->RegisterTrack(track, volumeID);
}

void TrackingAction::PostUserTrackingAction(const G4Track* track) {
  if (!track || !fEventAction) return;
  if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) return;

  const auto* detector = GetDetector();
  const auto* pv = track->GetVolume();
  const auto* lv = pv ? pv->GetLogicalVolume() : nullptr;
  G4int volumeID = detector ? detector->GetVolumeID(lv) : -1;

  const G4VProcess* endProcess = nullptr;
  const auto* step = track->GetStep();
  if (step) {
    const auto* post = step->GetPostStepPoint();
    if (post) {
      endProcess = post->GetProcessDefinedStep();
      const auto* postPV = post->GetPhysicalVolume();
      const auto* postLV = postPV ? postPV->GetLogicalVolume() : nullptr;
      volumeID = detector ? detector->GetVolumeID(postLV) : -1;
    }
  }

  fEventAction->UpdateTrackEnd(track, volumeID, endProcess);
}
