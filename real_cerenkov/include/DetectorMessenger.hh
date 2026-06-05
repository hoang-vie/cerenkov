#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;

class DetectorMessenger : public G4UImessenger {
public:
  DetectorMessenger(DetectorConstruction* Det);
  virtual ~DetectorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

private:
  DetectorConstruction* fDetectorConstruction;
  G4UIcmdWithAString* fCouplingFileCmd;
  G4UIcmdWithAString* fTyvekPosCmd;
  G4UIcmdWith3VectorAndUnit* fScintSizeCmd;
  G4UIcmdWith3VectorAndUnit* fSipmSizeCmd;
};

#endif