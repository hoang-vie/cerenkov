#include "PrimaryGeneratorAction.hh"

#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYSetup.h"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "RNGWrapper.hh"
#include "Randomize.hh"

#include <algorithm>
#include <fstream>
#include <sstream>

#ifndef CRY_DATA_DIR
#define CRY_DATA_DIR "./CRY-1.7-prefix/src/CRY-1.7/data"
#endif

PrimaryGeneratorAction::PrimaryGeneratorAction(const char* inputfile)
  : fParticleGun(new G4ParticleGun(1)),
    fParticles(new std::vector<CRYParticle*>),
    fParticleTable(G4ParticleTable::GetParticleTable()) {
  fMessenger = new PrimaryGeneratorMessenger(this);

  if (inputfile && std::string(inputfile).size() > 0) {
    CRYFromFile(inputfile);
  }
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  if (fParticles) {
    for (auto* p : *fParticles) delete p;
  }
  delete fMessenger;
  delete fParticleGun;
  delete fParticles;
  delete fGenerator;
  delete fCrySetup;
}

G4String PrimaryGeneratorAction::CRYDataDirectory() const {
  return CRY_DATA_DIR;
}


G4double PrimaryGeneratorAction::GetCRYTimeSimulated() const {
  return fGenerator ? fGenerator->timeSimulated() : 0.0;
}

void PrimaryGeneratorAction::ConfigureCRY(const std::string& setupString) {
  auto* setup = new CRYSetup(setupString, CRYDataDirectory().c_str());
  RNGWrapper<CLHEP::HepRandomEngine>::set(
      CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
  setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);

  delete fGenerator;
  delete fCrySetup;
  fCrySetup = setup;
  fGenerator = new CRYGenerator(setup);
  fCryConfigString = setupString;
  fInputState = 0;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  if (fGeneratorMode == "Gun") {
    auto* particle = fParticleTable->FindParticle(fGunParticleName);
    if (!particle) {
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "UnknownGunParticle",
                  RunMustBeAborted, "Requested particle is not in G4ParticleTable.");
      return;
    }

    G4ThreeVector position = fGunPosition;
    if (fGunUseAutoY) {
      const auto* det = static_cast<const DetectorConstruction*>(
          G4RunManager::GetRunManager()->GetUserDetectorConstruction());
      if (det) position.setY(det->GetGunPosY());
    }
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleEnergy(fGunEnergy);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->SetParticleMomentumDirection(fGunDirection);
    if (particle == G4OpticalPhoton::OpticalPhotonDefinition()) {
      G4ThreeVector pol = fGunPolarization - fGunDirection * fGunPolarization.dot(fGunDirection);
      if (pol.mag2() <= 0.) {
        pol = fGunDirection.orthogonal();
      }
      fParticleGun->SetParticlePolarization(pol.unit());
    }
    fParticleGun->SetParticleTime(fGunTime);
    fParticleGun->GeneratePrimaryVertex(event);
    return;
  }

  if (fInputState != 0 || !fGenerator) {
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "CRYNotInitialized",
                RunMustBeAborted,
                "CRY mode requires /CRY/input ... followed by /CRY/update, or /CRY/file.");
    return;
  }

  for (auto* p : *fParticles) delete p;
  fParticles->clear();
  fGenerator->genEvent(fParticles);

  G4double tMin = 0.0;
  if (!fParticles->empty()) {
    tMin = (*fParticles)[0]->t();
    for (auto* p : *fParticles) tMin = std::min(tMin, p->t());
  }

  const auto* det = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  const G4double gunY = det ? det->GetGunPosY() : 120.0*cm;

  for (auto* p : *fParticles) {
    auto* definition = fParticleTable->FindParticle(p->PDGid());
    if (!definition) continue;

    fParticleGun->SetParticleDefinition(definition);
    fParticleGun->SetParticleEnergy(p->ke()*MeV);
    fParticleGun->SetParticlePosition(G4ThreeVector(p->x()*m, gunY, p->y()*m));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(p->u(), p->w(), p->v()));

    G4double time = 0.0;
    if (fTimeMode == "absolute") time = p->t()*s;
    else if (fTimeMode == "relative") time = (p->t()-tMin)*s;
    fParticleGun->SetParticleTime(time);
    fParticleGun->GeneratePrimaryVertex(event);
  }

  for (auto* p : *fParticles) delete p;
  fParticles->clear();
}

void PrimaryGeneratorAction::InputCRY() {
  fInputState = 1;
}

void PrimaryGeneratorAction::UpdateCRY(std::string* input) {
  if (!input) return;
  ConfigureCRY(*input);
}

void PrimaryGeneratorAction::CRYFromFile(G4String fileName) {
  std::ifstream in(fileName);
  if (!in) {
    fInputState = -1;
    G4Exception("PrimaryGeneratorAction::CRYFromFile", "CRYFileOpenFailed",
                JustWarning, ("Cannot open CRY input file: " + fileName).c_str());
    return;
  }
  std::ostringstream ss;
  ss << in.rdbuf();
  ConfigureCRY(ss.str());
}

void PrimaryGeneratorAction::SetTimeMode(const G4String& mode) {
  if (mode == "relative" || mode == "absolute" || mode == "zero") {
    fTimeMode = mode;
    return;
  }
  G4Exception("PrimaryGeneratorAction::SetTimeMode", "InvalidTimeMode", JustWarning,
              "Supported modes: relative absolute zero.");
}

void PrimaryGeneratorAction::SetGeneratorMode(const G4String& mode) {
  if (mode == "CRY" || mode == "cry") {
    fGeneratorMode = "CRY";
  } else if (mode == "Gun" || mode == "gun") {
    fGeneratorMode = "Gun";
  } else {
    G4Exception("PrimaryGeneratorAction::SetGeneratorMode", "InvalidGeneratorMode",
                JustWarning, "Supported modes: CRY Gun.");
  }
}

void PrimaryGeneratorAction::SetGunParticle(const G4String& particleName) {
  fGunParticleName = particleName;
}

void PrimaryGeneratorAction::SetGunDirection(const G4ThreeVector& dir) {
  if (dir.mag2() <= 0.) {
    G4Exception("PrimaryGeneratorAction::SetGunDirection", "ZeroDirection",
                JustWarning, "Gun direction must be non-zero.");
    return;
  }
  fGunDirection = dir.unit();
}

void PrimaryGeneratorAction::SetGunPolarization(const G4ThreeVector& pol) {
  if (pol.mag2() <= 0.) {
    G4Exception("PrimaryGeneratorAction::SetGunPolarization", "ZeroPolarization",
                JustWarning, "Gun polarization must be non-zero.");
    return;
  }
  fGunPolarization = pol.unit();
}
