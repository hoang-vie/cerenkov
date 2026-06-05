#include "PrimaryGeneratorAction.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include "CRYSetup.h"
#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYUtils.h"


/*
static void GenThetaPhi(double& theta, double& phi) {
  double cosTheta;
  do { cosTheta = G4UniformRand(); } while (G4UniformRand() > cosTheta * cosTheta);
  theta = std::acos(cosTheta);  
  phi   = G4UniformRand() * 2.0 * CLHEP::pi;
}
*/
PrimaryGeneratorAction::PrimaryGeneratorAction(const char* inputfile)
  : G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr), fEnvelopeBox(nullptr)
{
  fParticleGun = new G4ParticleGun();
  std::ifstream inputFile;
  inputFile.open(inputfile, std::ios::in);
  char buffer[1000];

  if (inputFile.fail()) {
    InputState = -1;
  } else {
    std::string setupString("");
    while (!inputFile.getline(buffer, 1000).eof()) {
      setupString.append(buffer);
      setupString.append(" ");
    }
    CRYSetup* setup = new CRYSetup(setupString, "../data");
    gen = new CRYGenerator(setup);
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState = 0;
  }
  vect = new std::vector<CRYParticle*>;
  particleTable = G4ParticleTable::GetParticleTable();
  gunMessenger = new PrimaryGeneratorMessenger(this);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() { 
  delete fParticleGun; 
  delete vect;
  delete gen;
  delete gunMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  if (InputState != 0) {
    G4Exception("PrimaryGeneratorAction", "1", RunMustBeAborted, "CRY not initialized");
  }

  vect->clear();
  gen->genEvent(vect);

  for (unsigned j = 0; j < vect->size(); j++) {
    G4double particleEnergy = (*vect)[j]->ke() * MeV;
    G4double particlePosX = (*vect)[j]->x() * m;
    G4double particlePosY =  10.*cm;            
    G4double particlePosZ = (*vect)[j]->y() * m;
    G4double particleTime = (*vect)[j]->t();
    /*
    double theta, phi;
    GenThetaPhi(theta, phi);

    G4double ux =  std::sin(theta) * std::cos(phi);
    G4double uy = -std::cos(theta);   
    G4double uz =  std::sin(theta) * std::sin(phi); 
    */

    G4double ux = (*vect)[j]->u();
    G4double uy = (*vect)[j]->v();
    G4double uz = (*vect)[j]->w();

    fParticleGun->SetParticleDefinition(particleTable->FindParticle((*vect)[j]->PDGid()));
    fParticleGun->SetParticleEnergy(particleEnergy);
    fParticleGun->SetParticlePosition(G4ThreeVector(particlePosX, particlePosY, particlePosZ));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));
    fParticleGun->SetParticleTime(particleTime);

    fParticleGun->GeneratePrimaryVertex(anEvent);
    delete (*vect)[j];
  }
}

void PrimaryGeneratorAction::InputCRY() { InputState = 1; }

void PrimaryGeneratorAction::UpdateCRY(std::string* MessInput) {
  CRYSetup* setup = new CRYSetup(*MessInput, "./CRY-1.7-prefix/src/CRY-1.7/data");
  gen = new CRYGenerator(setup);
  RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
  setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
  InputState = 0;
}

void PrimaryGeneratorAction::CRYFromFile(G4String newValue) {
  std::ifstream inputFile;
  inputFile.open(newValue, std::ios::in);
  char buffer[1000];
  if (inputFile.fail()) {
    InputState = -1;
  } else {
    std::string setupString("");
    while (!inputFile.getline(buffer, 1000).eof()) {
      setupString.append(buffer); setupString.append(" ");
    }
    CRYSetup* setup = new CRYSetup(setupString, "../data");
    gen = new CRYGenerator(setup);
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState = 0;
  }
}