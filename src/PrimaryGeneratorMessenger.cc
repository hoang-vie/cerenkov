#include "PrimaryGeneratorMessenger.hh"

#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"
#include "PrimaryGeneratorAction.hh"
#include <iomanip>
#include "G4SystemOfUnits.hh"
PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* action)
  : fAction(action) {
  fCRYDir = new G4UIdirectory("/CRY/");
  fCRYDir->SetGuidance("CRY configuration.");

  fFileCmd = new G4UIcmdWithAString("/CRY/file", this);
  fFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fInputCmd = new G4UIcmdWithAString("/CRY/input", this);
  fInputCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fUpdateCmd = new G4UIcmdWithoutParameter("/CRY/update", this);
  fUpdateCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fTimeModeCmd = new G4UIcmdWithAString("/CRY/timeMode", this);
  fTimeModeCmd->SetCandidates("relative absolute zero");
  fTimeModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGeneratorDir = new G4UIdirectory("/generator/");
  fGeneratorDir->SetGuidance("General primary generator settings.");

  fModeCmd = new G4UIcmdWithAString("/generator/mode", this);
  fModeCmd->SetCandidates("CRY Gun");
  fModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fParticleCmd = new G4UIcmdWithAString("/generator/particle", this);
  fParticleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEnergyCmd = new G4UIcmdWithADoubleAndUnit("/generator/energy", this);
  fEnergyCmd->SetUnitCategory("Energy");
  fEnergyCmd->SetDefaultUnit("GeV");
  fEnergyCmd->SetParameterName("gunEnergy", false);
  fEnergyCmd->SetRange("gunEnergy>=0");
  fEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/generator/position", this);
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->SetDefaultUnit("mm");
  fPositionCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDirectionCmd = new G4UIcmdWith3Vector("/generator/direction", this);
  fDirectionCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fPolarizationCmd = new G4UIcmdWith3Vector("/generator/polarization", this);
  fPolarizationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGunTimeCmd = new G4UIcmdWithADoubleAndUnit("/generator/time", this);
  fGunTimeCmd->SetUnitCategory("Time");
  fGunTimeCmd->SetDefaultUnit("ns");
  fGunTimeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fUseAutoYCmd = new G4UIcmdWithABool("/generator/useAutoY", this);
  fUseAutoYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() {
  delete fUseAutoYCmd;
  delete fGunTimeCmd;
  delete fPolarizationCmd;
  delete fDirectionCmd;
  delete fPositionCmd;
  delete fEnergyCmd;
  delete fParticleCmd;
  delete fModeCmd;
  delete fGeneratorDir;
  delete fTimeModeCmd;
  delete fUpdateCmd;
  delete fInputCmd;
  delete fFileCmd;
  delete fCRYDir;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String value) {
  if (command == fInputCmd) {
    fAction->InputCRY();
    fCRYInput += value;
    fCRYInput += " ";
  } else if (command == fUpdateCmd) {
    fAction->UpdateCRY(&fCRYInput);
    fCRYInput.clear();
  } else if (command == fFileCmd) {
    fAction->CRYFromFile(value);
  } else if (command == fTimeModeCmd) {
    fAction->SetTimeMode(value);
  } else if (command == fModeCmd) {
    fAction->SetGeneratorMode(value);
  } else if (command == fParticleCmd) {
    fAction->SetGunParticle(value);
  } else if (command == fEnergyCmd) {
    const G4double parsed = fEnergyCmd->GetNewDoubleValue(value);
    fAction->SetGunEnergy(parsed);
  } else if (command == fPositionCmd) {
    fAction->SetGunPosition(fPositionCmd->GetNew3VectorValue(value));
  } else if (command == fDirectionCmd) {
    fAction->SetGunDirection(fDirectionCmd->GetNew3VectorValue(value));
  } else if (command == fPolarizationCmd) {
    fAction->SetGunPolarization(fPolarizationCmd->GetNew3VectorValue(value));
  } else if (command == fGunTimeCmd) {
    fAction->SetGunTime(fGunTimeCmd->GetNewDoubleValue(value));
  } else if (command == fUseAutoYCmd) {
    fAction->SetGunUseAutoY(fUseAutoYCmd->GetNewBoolValue(value));
  }
}
