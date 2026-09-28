#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

class G4Run;

enum NtupleID {
  kTriggerEvent = 0,
  kScintDeposit = 1,
  kTrack = 2,
  kTrackVolume = 3,
  kPhoton = 4,
  kSensorIncident = 5,
  kDetectedSignal = 6,
  kChannelEvent = 7,
  kRunConfig = 8,
  kChannelMap = 9,
  kVolumeMap = 10,
  kRunStatistics = 11,
  kBoundaryDebug = 12,
  kReadoutPlaneConfig = 13,
  kPrimaryTruth = 14,
  kCoreBoundaryCrossing = 15
};

class RunAction : public G4UserRunAction {
public:
  RunAction();
  ~RunAction() override;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

  void CountEvent(G4bool triggered, G4bool detailsWritten);

private:
  G4Accumulable<G4int> fGeneratedEvents {0};
  G4Accumulable<G4int> fTriggeredEvents {0};
  G4Accumulable<G4int> fDetailedEvents {0};
  G4Accumulable<G4double> fCRYSimulatedTime_s {0.0};
  G4double fCRYTimeStart_s = 0.0;
};

#endif
