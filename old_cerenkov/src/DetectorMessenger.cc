#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIcmdWithAString.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
 : G4UImessenger(), fDetectorConstruction(Det)
{
  fCouplingFileCmd = new G4UIcmdWithAString("/detector/setCouplingFile", this);
  fCouplingFileCmd->SetGuidance("Nhap duong dan file txt chua RINDEX cua lop vo SD2.");
  fCouplingFileCmd->SetParameterName("filepath", false);
  fCouplingFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger() {
  delete fCouplingFileCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == fCouplingFileCmd) {
    fDetectorConstruction->SetCouplingFile(newValue);
  }
}