#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4VUserTrackInformation.hh"
#include "globals.hh"

#include <cstdint>
#include <limits>

class DetectorConstruction;
class EventAction;
class G4LogicalVolume;
class G4OpBoundaryProcess;
class G4Step;

class PhotonInfo : public G4VUserTrackInformation {
public:
  PhotonInfo() = default;
  ~PhotonInfo() override = default;

  G4double emissionAngleDeg = std::numeric_limits<G4double>::quiet_NaN();
  // Exact beta convention used by G4Cerenkov for photon-angle sampling in the
  // parent step: 0.5*(preStepBeta + postStepBeta). NaN for non-Cerenkov photons.
  G4double emitterBeta = std::numeric_limits<G4double>::quiet_NaN();
};

class SteppingAction : public G4UserSteppingAction {
public:
  explicit SteppingAction(EventAction* eventAction);
  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step*) override;

private:
  void RefreshGeometryHandles();
  G4OpBoundaryProcess* GetBoundaryProcess();

  EventAction* fEventAction = nullptr;
  const DetectorConstruction* fDetector = nullptr;
  std::uint64_t fGeometryVersion = std::numeric_limits<std::uint64_t>::max();

  G4LogicalVolume* fCoreLV = nullptr;
  G4LogicalVolume* fAcrylicLV = nullptr;
  G4LogicalVolume* fScintTLV = nullptr;
  G4LogicalVolume* fScintBLV = nullptr;
  G4LogicalVolume* fTyvekLV = nullptr;
  G4LogicalVolume* fCouplingLV = nullptr;
  G4LogicalVolume* fSensorLV = nullptr;
  G4double fAcrylicAbsMinNm = 0.0;
  G4double fAcrylicAbsMaxNm = 0.0;
  G4OpBoundaryProcess* fBoundaryProcess = nullptr;
};

#endif
