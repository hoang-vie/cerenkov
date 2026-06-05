#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
 : G4UImessenger(), fDetectorConstruction(Det)
{
  fCouplingFileCmd = new G4UIcmdWithAString("/detector/setCouplingFile", this);
  fCouplingFileCmd->SetGuidance("Set RINDEX file for Coupling layer.");
  fCouplingFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fTyvekPosCmd = new G4UIcmdWithAString("/detector/tyvekPosition", this);
  fTyvekPosCmd->SetGuidance("Select Tyvek setup: Inside or Outside");
  fTyvekPosCmd->AvailableForStates(G4State_PreInit);

  fScintSizeCmd = new G4UIcmdWith3VectorAndUnit("/detector/scintSize", this);
  fScintSizeCmd->SetGuidance("Set the full dimensions of Scintillator (X Y Z).");
  fScintSizeCmd->SetParameterName("sizeX", "sizeY", "sizeZ", false);
  fScintSizeCmd->SetDefaultUnit("mm");
  fScintSizeCmd->AvailableForStates(G4State_PreInit);

  fSipmSizeCmd = new G4UIcmdWith3VectorAndUnit("/detector/sipmSize", this);
  fSipmSizeCmd->SetGuidance("Set the full dimensions of SiPM (X Y Z).");
  fSipmSizeCmd->SetParameterName("sizeX", "sizeY", "sizeZ", false);
  fSipmSizeCmd->SetDefaultUnit("mm");
  fSipmSizeCmd->AvailableForStates(G4State_PreInit);
}

DetectorMessenger::~DetectorMessenger() {
  delete fCouplingFileCmd;
  delete fTyvekPosCmd;
  delete fScintSizeCmd;
  delete fSipmSizeCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == fCouplingFileCmd) fDetectorConstruction->SetCouplingFile(newValue);
  else if (command == fTyvekPosCmd) fDetectorConstruction->SetTyvekPosition(newValue);
  else if (command == fScintSizeCmd) fDetectorConstruction->SetScintSize(fScintSizeCmd->GetNew3VectorValue(newValue));
  else if (command == fSipmSizeCmd) fDetectorConstruction->SetSipmSize(fSipmSizeCmd->GetNew3VectorValue(newValue));
}