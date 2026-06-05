#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIcmdWithAString;

class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction* );
    virtual ~DetectorMessenger();
    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    DetectorConstruction* fDetectorConstruction;
    G4UIcmdWithAString* fCouplingFileCmd;
};
#endif