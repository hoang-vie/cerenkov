#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

RunAction::RunAction() : G4UserRunAction()
{
  auto anaMan = G4AnalysisManager::Instance();
  anaMan->SetVerboseLevel(0);
  anaMan->SetNtupleMerging(true);
  anaMan->SetDefaultFileType("root");

  // TREE 0: Truth
  anaMan->CreateNtuple("Truth", "Primary Cosmic Ray");
  anaMan->CreateNtupleIColumn("EventID");          // 0
  anaMan->CreateNtupleIColumn("TrackID");          // 1
  anaMan->CreateNtupleIColumn("PDG_Code");         // 2 [MỚI] 
  anaMan->CreateNtupleDColumn("Initial_E_MeV");    // 3 [MỚI]
  anaMan->CreateNtupleDColumn("Edep_Water_MeV");   // 4 [MỚI]
  anaMan->CreateNtupleDColumn("Len_Water_mm");     // 5 [MỚI]
  anaMan->CreateNtupleDColumn("Edep_Acrylic_MeV"); // 6 [MỚI]
  anaMan->CreateNtupleDColumn("Len_Acrylic_mm");   // 7 [MỚI]
  anaMan->CreateNtupleIColumn("OptPhotonsGen");    // 8
  anaMan->CreateNtupleDColumn("DirX");             // 9
  anaMan->CreateNtupleDColumn("DirY");             // 10
  anaMan->CreateNtupleDColumn("DirZ");             // 11
  anaMan->FinishNtuple(0);

  // TREE 1: Secondary 
  anaMan->CreateNtuple("Secondary", "Delta Electrons");
  anaMan->CreateNtupleIColumn("EventID");          
  anaMan->CreateNtupleIColumn("TrackID");          
  anaMan->CreateNtupleIColumn("ParentID");         
  anaMan->CreateNtupleIColumn("PDG_Code");         
  anaMan->CreateNtupleDColumn("Energy_MeV");       
  anaMan->CreateNtupleSColumn("Process");          
  anaMan->FinishNtuple(1);

  // TREE 2: InnerData 
  anaMan->CreateNtuple("InnerData", "Photons hitting Water-Acrylic boundary");
  anaMan->CreateNtupleIColumn("EventID");          // 0
  anaMan->CreateNtupleIColumn("TrackID");          // 1
  anaMan->CreateNtupleIColumn("ParentID");         // 2
  anaMan->CreateNtupleDColumn("Energy_eV");        // 3
  anaMan->CreateNtupleDColumn("Wavelength_nm");    // 4
  anaMan->CreateNtupleDColumn("HitX_mm");          // 5
  anaMan->CreateNtupleDColumn("HitY_mm");          // 6
  anaMan->CreateNtupleDColumn("HitZ_mm");          // 7
  anaMan->CreateNtupleDColumn("DirX");             // 8
  anaMan->CreateNtupleDColumn("DirY");             // 9
  anaMan->CreateNtupleDColumn("DirZ");             // 10
  anaMan->CreateNtupleDColumn("Theta_deg");        // 11
  anaMan->CreateNtupleDColumn("Phi_deg");          // 12
  anaMan->CreateNtupleDColumn("AngleToPlane_deg"); // 13
  anaMan->CreateNtupleDColumn("Time_ns");          // 14
  anaMan->CreateNtupleIColumn("SurfaceID");        // 15
  anaMan->CreateNtupleDColumn("NormalX");          // 16
  anaMan->CreateNtupleDColumn("NormalY");          // 17
  anaMan->CreateNtupleDColumn("NormalZ");          // 18
  anaMan->CreateNtupleDColumn("PolX");             // 19
  anaMan->CreateNtupleDColumn("PolY");             // 20
  anaMan->CreateNtupleDColumn("PolZ");             // 21
  anaMan->CreateNtupleIColumn("OriginVol");        // 22 [MỚI] 0=Water, 1=Acrylic
  anaMan->FinishNtuple(2);

  // TREE 3: OuterData
  anaMan->CreateNtuple("OuterData", "Photons exiting Acrylic to Coupling");
  anaMan->CreateNtupleIColumn("EventID");          // 0
  anaMan->CreateNtupleIColumn("TrackID");          // 1
  anaMan->CreateNtupleIColumn("ParentID");         // 2
  anaMan->CreateNtupleDColumn("Energy_eV");        // 3
  anaMan->CreateNtupleDColumn("Wavelength_nm");    // 4
  anaMan->CreateNtupleDColumn("HitX_mm");          // 5
  anaMan->CreateNtupleDColumn("HitY_mm");          // 6
  anaMan->CreateNtupleDColumn("HitZ_mm");          // 7
  anaMan->CreateNtupleDColumn("DirX");             // 8
  anaMan->CreateNtupleDColumn("DirY");             // 9
  anaMan->CreateNtupleDColumn("DirZ");             // 10
  anaMan->CreateNtupleDColumn("Theta_deg");        // 11
  anaMan->CreateNtupleDColumn("Phi_deg");          // 12
  anaMan->CreateNtupleDColumn("AngleToPlane_deg"); // 13
  anaMan->CreateNtupleDColumn("Time_ns");          // 14
  anaMan->CreateNtupleIColumn("SurfaceID");        // 15
  anaMan->CreateNtupleDColumn("NormalX");          // 16
  anaMan->CreateNtupleDColumn("NormalY");          // 17
  anaMan->CreateNtupleDColumn("NormalZ");          // 18
  anaMan->CreateNtupleDColumn("PolX");             // 19
  anaMan->CreateNtupleDColumn("PolY");             // 20
  anaMan->CreateNtupleDColumn("PolZ");             // 21
  anaMan->CreateNtupleIColumn("OriginVol");        // 22 [MỚI]
  anaMan->FinishNtuple(3);
  
  // TREE 4: OpticalTrace
  anaMan->CreateNtuple("OpticalTrace", "Trace directions");
  anaMan->CreateNtupleIColumn("EventID");          // 0
  anaMan->CreateNtupleIColumn("TrackID");          // 1
  anaMan->CreateNtupleIColumn("BoundaryID");       // 2
  anaMan->CreateNtupleIColumn("IsReflected");      // 3
  anaMan->CreateNtupleDColumn("HitX_mm");          // 4
  anaMan->CreateNtupleDColumn("HitY_mm");          // 5
  anaMan->CreateNtupleDColumn("HitZ_mm");          // 6
  anaMan->CreateNtupleDColumn("InDirX");           // 7
  anaMan->CreateNtupleDColumn("InDirY");           // 8
  anaMan->CreateNtupleDColumn("InDirZ");           // 9
  anaMan->CreateNtupleDColumn("OutDirX");          // 10
  anaMan->CreateNtupleDColumn("OutDirY");          // 11
  anaMan->CreateNtupleDColumn("OutDirZ");          // 12
  anaMan->CreateNtupleIColumn("OriginVol");        // 13 [MỚI]
  anaMan->FinishNtuple(4);

  // TREE 5: EventSummary (Bổ sung Link với Primary)
  anaMan->CreateNtuple("EventSummary", "Counters");
  anaMan->CreateNtupleIColumn("EventID");          // 0
  anaMan->CreateNtupleIColumn("Primary_PDG");      // 1 [MỚI]
  anaMan->CreateNtupleDColumn("Primary_E_MeV");    // 2 [MỚI]
  anaMan->CreateNtupleDColumn("Edep_Water_MeV");   // 3 [MỚI]
  anaMan->CreateNtupleIColumn("DeltaRays_Count");  // 4 [MỚI]
  anaMan->CreateNtupleIColumn("GenWater");         // 5
  anaMan->CreateNtupleIColumn("GenAcrylic");       // 6
  anaMan->CreateNtupleIColumn("GenCoupling");      // 7 [MỚI]
  anaMan->CreateNtupleIColumn("W2A_Unique");       // 8
  anaMan->CreateNtupleIColumn("W2A_Total");        // 9
  anaMan->CreateNtupleIColumn("A2W_Total");        // 10 
  anaMan->CreateNtupleIColumn("A2C_Unique");       // 11
  anaMan->CreateNtupleIColumn("A2C_Total");        // 12 
  anaMan->CreateNtupleIColumn("C2A_Total");        // 13 
  anaMan->CreateNtupleIColumn("C2Sensor_Unique");  // 14 
  anaMan->CreateNtupleIColumn("Detected_Total");   // 15 
  anaMan->FinishNtuple(5);

  // TREE 6: DetectedSignal 
  anaMan->CreateNtuple("DetectedSignal", "Photons absorbed by SiPM");
  anaMan->CreateNtupleIColumn("EventID");          // 0
  anaMan->CreateNtupleIColumn("TrackID");          // 1
  anaMan->CreateNtupleDColumn("Wavelength_nm");    // 2
  anaMan->CreateNtupleDColumn("Time_ns");          // 3
  anaMan->CreateNtupleDColumn("HitX_mm");          // 4
  anaMan->CreateNtupleDColumn("HitY_mm");          // 5
  anaMan->CreateNtupleDColumn("HitZ_mm");          // 6
  anaMan->CreateNtupleIColumn("OriginVol");        // 7 [MỚI]
  anaMan->CreateNtupleIColumn("BounceCount");      // 8 [MỚI]
  anaMan->FinishNtuple(6);
}

RunAction::~RunAction() { delete G4AnalysisManager::Instance(); }
void RunAction::BeginOfRunAction(const G4Run*) {
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  G4AnalysisManager::Instance()->OpenFile("cherenkov_water_data"); 
}
void RunAction::EndOfRunAction(const G4Run* run) {
  if (run->GetNumberOfEvent() == 0) return;
  G4AnalysisManager::Instance()->Write();
  G4AnalysisManager::Instance()->CloseFile();
}