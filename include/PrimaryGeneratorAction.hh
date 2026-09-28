#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <string>
#include <vector>

class CRYGenerator;
class CRYParticle;
class CRYSetup;
class G4Event;
class G4ParticleGun;
class G4ParticleTable;
class PrimaryGeneratorMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  explicit PrimaryGeneratorAction(const char* inputfile);
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event*) override;

  void InputCRY();
  void UpdateCRY(std::string* input);
  void CRYFromFile(G4String fileName);
  void SetTimeMode(const G4String& mode);

  void SetGeneratorMode(const G4String& mode);
  void SetGunParticle(const G4String& particleName);
  void SetGunEnergy(G4double energy) { fGunEnergy = energy; }
  void SetGunPosition(const G4ThreeVector& pos) { fGunPosition = pos; fGunUseAutoY = false; }
  void SetGunDirection(const G4ThreeVector& dir);
  void SetGunPolarization(const G4ThreeVector& pol);
  void SetGunTime(G4double time) { fGunTime = time; }
  void SetGunUseAutoY(G4bool value) { fGunUseAutoY = value; }

  G4String GetGeneratorMode() const { return fGeneratorMode; }
  G4String GetTimeMode() const { return fTimeMode; }
  G4String GetCRYConfigString() const { return fCryConfigString; }
  G4String GetGunParticleName() const { return fGunParticleName; }
  G4double GetGunEnergy() const { return fGunEnergy; }
  G4ThreeVector GetGunPosition() const { return fGunPosition; }
  G4ThreeVector GetGunDirection() const { return fGunDirection; }
  G4ThreeVector GetGunPolarization() const { return fGunPolarization; }
  G4double GetGunTime() const { return fGunTime; }
  G4bool GetGunUseAutoY() const { return fGunUseAutoY; }
  G4double GetCRYTimeSimulated() const;
  const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

private:
  void ConfigureCRY(const std::string& setupString);
  G4String CRYDataDirectory() const;

  G4ParticleGun* fParticleGun = nullptr;
  std::vector<CRYParticle*>* fParticles = nullptr;
  G4ParticleTable* fParticleTable = nullptr;
  CRYGenerator* fGenerator = nullptr;
  CRYSetup* fCrySetup = nullptr;
  PrimaryGeneratorMessenger* fMessenger = nullptr;
  G4int fInputState = -1;

  G4String fGeneratorMode = "CRY";
  G4String fTimeMode = "relative";
  G4String fCryConfigString;

  G4String fGunParticleName = "mu-";
  G4double fGunEnergy = 3.0*GeV;
  G4ThreeVector fGunPosition = G4ThreeVector();
  G4ThreeVector fGunDirection = G4ThreeVector(0.,-1.,0.);
  G4ThreeVector fGunPolarization = G4ThreeVector(1.,0.,0.);
  G4double fGunTime = 0.0;
  G4bool fGunUseAutoY = true;
};

#endif
