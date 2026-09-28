#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class DetectorConstruction;
class EventAction;
class G4Track;

// Registers generic non-optical track truth at tracking start/end.  Keeping
// registration here (rather than waiting for step 1) makes Track semantics
// exact even for a track that is created but never takes a transport step.
class TrackingAction : public G4UserTrackingAction {
public:
  explicit TrackingAction(EventAction* eventAction);
  ~TrackingAction() override = default;

  void PreUserTrackingAction(const G4Track* track) override;
  void PostUserTrackingAction(const G4Track* track) override;

private:
  const DetectorConstruction* GetDetector() const;
  EventAction* fEventAction = nullptr;
};

#endif
