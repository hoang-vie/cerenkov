#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <string>

class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIdirectory;
class PrimaryGeneratorAction;

class PrimaryGeneratorMessenger : public G4UImessenger {
public:
  explicit PrimaryGeneratorMessenger(PrimaryGeneratorAction* action);
  ~PrimaryGeneratorMessenger() override;
  void SetNewValue(G4UIcommand*, G4String) override;

private:
  PrimaryGeneratorAction* fAction = nullptr;
  G4UIdirectory* fCRYDir = nullptr;
  G4UIdirectory* fGeneratorDir = nullptr;

  G4UIcmdWithAString* fFileCmd = nullptr;
  G4UIcmdWithAString* fInputCmd = nullptr;
  G4UIcmdWithoutParameter* fUpdateCmd = nullptr;
  G4UIcmdWithAString* fTimeModeCmd = nullptr;

  G4UIcmdWithAString* fModeCmd = nullptr;
  G4UIcmdWithAString* fParticleCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fEnergyCmd = nullptr;
  G4UIcmdWith3VectorAndUnit* fPositionCmd = nullptr;
  G4UIcmdWith3Vector* fDirectionCmd = nullptr;
  G4UIcmdWith3Vector* fPolarizationCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fGunTimeCmd = nullptr;
  G4UIcmdWithABool* fUseAutoYCmd = nullptr;

  std::string fCRYInput;
};

#endif
