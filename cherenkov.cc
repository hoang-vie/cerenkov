#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "SimulationMetadata.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include <chrono>
#include <cstdlib>
#include <random>

namespace {
long MakeSeed(int argc, char** argv) {
  if (argc >= 3) {
    const long requested = std::strtol(argv[2], nullptr, 10);
    if (requested > 0) return requested;
  }
  const auto ticks = static_cast<unsigned long long>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
  std::random_device rd;
  unsigned long long mixed = ticks ^
      (static_cast<unsigned long long>(rd()) << 32) ^ rd();
  long seed = static_cast<long>(mixed % 900000000ULL) + 1L;
  return seed;
}
}

int main(int argc, char** argv) {
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc, argv);

  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  const long initialSeed = MakeSeed(argc, argv);
  CherenkovMetadata::InitialSeed.store(initialSeed);
  CLHEP::HepRandom::setTheSeed(initialSeed);

#ifdef G4MULTITHREADED
  auto* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(2);
#else
  auto* runManager = new G4RunManager;
#endif

  runManager->SetUserInitialization(new DetectorConstruction());

  auto* physicsList = new FTFP_BERT;
  physicsList->SetVerboseLevel(0);
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  physicsList->RegisterPhysics(new G4OpticalPhysics());
  physicsList->SetDefaultCutValue(0.1*mm);

  auto* optical = G4OpticalParameters::Instance();
  optical->SetProcessActivation("Cerenkov", true);
  optical->SetProcessActivation("Scintillation", false);
  optical->SetProcessActivation("OpAbsorption", true);
  optical->SetProcessActivation("OpRayleigh", true);
  optical->SetProcessActivation("OpMieHG", false);
  optical->SetProcessActivation("OpBoundary", true);
  optical->SetProcessActivation("OpWLS", false);
  optical->SetProcessActivation("OpWLS2", false);

  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new ActionInitialization());

  G4VisManager* visManager = nullptr;
  if (ui) {
    visManager = new G4VisExecutive;
    visManager->Initialize();
  }

  auto* uiManager = G4UImanager::GetUIpointer();
  if (ui) {
    uiManager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  } else if (argc >= 2) {
    uiManager->ApplyCommand("/control/execute " + G4String(argv[1]));
  }

  delete visManager;
  delete runManager;
  return 0;
}
