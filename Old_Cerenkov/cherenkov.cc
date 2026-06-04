/// cherenkov.cc  —  Main program cho mô phỏng Cherenkov trong bể nước

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4OpticalParameters.hh" 
#include "G4OpticalPhysics.hh"        
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "Randomize.hh"
#include <ctime>

int main(int argc, char** argv) {

  // ── Interactive / batch mode ──────────────────────────────────────────────
  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  // ── Random engine ─────────────────────────────────────────────────────────
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  time_t timeStart;
  time(&timeStart);
  CLHEP::HepRandom::setTheSeed((unsigned long)timeStart);

  // ── Run manager ───────────────────────────────────────────────────────────
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(2);
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // ── Detector ──────────────────────────────────────────────────────────────
  runManager->SetUserInitialization(new DetectorConstruction());

  // ── Physics list ──────────────────────────────────────────────────────────
  auto physicsList = new FTFP_BERT;
  physicsList->SetVerboseLevel(0);
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

  // Đăng ký G4OpticalPhysics (Chỉ 1 lần duy nhất)
  auto opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);
  
  // Cấu hình tham số quang học
  G4OpticalParameters* optParams = G4OpticalParameters::Instance();
  optParams->SetProcessActivation("Cerenkov",      true);  
  optParams->SetProcessActivation("Scintillation", false); 
  optParams->SetProcessActivation("OpAbsorption",  true);  
  optParams->SetProcessActivation("OpRayleigh",    true);  
  optParams->SetProcessActivation("OpMieHG",       false); 
  optParams->SetProcessActivation("OpBoundary",    true);  
  optParams->SetProcessActivation("OpWLS",         false); 

  // Nạp Physics List vào RunManager
  runManager->SetUserInitialization(physicsList);

  // ── User actions ──────────────────────────────────────────────────────────
  runManager->SetUserInitialization(new ActionInitialization());

  // ── Visualization ─────────────────────────────────────────────────────────
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // ── Run ───────────────────────────────────────────────────────────────────
  if (!ui) {
    // Batch mode
    UImanager->ApplyCommand("/control/execute " + G4String(argv[1]));
  } else {
    // Interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
  return 0;
}