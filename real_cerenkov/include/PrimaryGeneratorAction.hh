
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYSetup.h"
#include "CRYUtils.h"
#include "G4DataVector.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "RNGWrapper.hh"
#include "Randomize.hh"
#include "globals.hh"
#include <vector>

class G4ParticleGun;
class G4Event;
class G4Box;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued
/// in front of the phantom across 80% of the (X,Y) phantom size.

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction(const char *inputfile);
  virtual ~PrimaryGeneratorAction();

  // method from the base class
  virtual void GeneratePrimaries(G4Event *);
  void InputCRY();
  void UpdateCRY(std::string *MessInput);
  void CRYFromFile(G4String newValue);

  // method to access particle gun
  const G4ParticleGun *GetParticleGun() const { return fParticleGun; }

private:
  G4ParticleGun *fParticleGun; // pointer a to G4 gun class
  G4Box *fEnvelopeBox;
  std::vector<CRYParticle *> *vect; // vector of generated particles
  G4ParticleTable *particleTable;
  CRYGenerator *gen;
  G4int InputState;
  PrimaryGeneratorMessenger *gunMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
