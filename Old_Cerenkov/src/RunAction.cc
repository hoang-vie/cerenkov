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

  // ==============================================================
  // TREE 0: Truth (Lưu thông tin hạt sơ cấp)
  // ==============================================================
  anaMan->CreateNtuple("Truth", "Primary Cosmic Ray Info");
  anaMan->CreateNtupleIColumn("EventID");          
  anaMan->CreateNtupleSColumn("Particle");         
  anaMan->CreateNtupleDColumn("Energy_GeV");       
  anaMan->CreateNtupleDColumn("StartX_mm");        
  anaMan->CreateNtupleDColumn("StartY_mm");        
  anaMan->CreateNtupleDColumn("StartZ_mm");        
  anaMan->CreateNtupleDColumn("GlobalTime_ns");    
  anaMan->CreateNtupleDColumn("Theta_deg");        
  anaMan->CreateNtupleDColumn("Phi_deg");          
  anaMan->FinishNtuple(0);

  // ==============================================================
  // TREE 1: Secondary (Bộ lọc electron/positron thứ cấp)
  // ==============================================================
  anaMan->CreateNtuple("Secondary", "Delta Electrons & Cherenkov Candidates");
  anaMan->CreateNtupleIColumn("EventID");          
  anaMan->CreateNtupleIColumn("TrackID");          
  anaMan->CreateNtupleIColumn("ParentID");         
  anaMan->CreateNtupleIColumn("PDG_Code");         
  anaMan->CreateNtupleDColumn("Energy_MeV");       
  anaMan->CreateNtupleSColumn("Process");          
  anaMan->FinishNtuple(1);

  // ==============================================================
  // TREE 2: InnerData (Nước -> Acrylic)
  // ==============================================================
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
  anaMan->CreateNtupleDColumn("AngleToPlane_deg"); // 13: Góc TỚI (Incidence)
  anaMan->CreateNtupleDColumn("Time_ns");          // 14
  anaMan->CreateNtupleIColumn("SurfaceID");        // 15
  anaMan->CreateNtupleDColumn("NormalX");          // 16: Thêm Vector Pháp tuyến
  anaMan->CreateNtupleDColumn("NormalY");          // 17
  anaMan->CreateNtupleDColumn("NormalZ");          // 18
  anaMan->CreateNtupleDColumn("PolX");             // 19: Thêm Vector Phân cực
  anaMan->CreateNtupleDColumn("PolY");             // 20
  anaMan->CreateNtupleDColumn("PolZ");             // 21
  anaMan->FinishNtuple(2);

  // ==============================================================
  // TREE 3: OuterData (Acrylic -> Air)
  // ==============================================================
  anaMan->CreateNtuple("OuterData", "Photons exiting Acrylic to Air");
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
  anaMan->CreateNtupleDColumn("AngleToPlane_deg"); // 13: Góc KHÚC XẠ (Refraction)
  anaMan->CreateNtupleDColumn("Time_ns");          // 14
  anaMan->CreateNtupleIColumn("SurfaceID");        // 15
  anaMan->CreateNtupleDColumn("NormalX");          // 16
  anaMan->CreateNtupleDColumn("NormalY");          // 17
  anaMan->CreateNtupleDColumn("NormalZ");          // 18
  anaMan->CreateNtupleDColumn("PolX");             // 19
  anaMan->CreateNtupleDColumn("PolY");             // 20
  anaMan->CreateNtupleDColumn("PolZ");             // 21
  anaMan->FinishNtuple(3);
  
  // ==============================================================
  // TREE 4: OpticalTrace (Theo dõi tương tác ranh giới chi tiết)
  // ==============================================================
  anaMan->CreateNtuple("OpticalTrace", "Trace directions for Snell/Fresnel laws");
  anaMan->CreateNtupleIColumn("EventID");          
  anaMan->CreateNtupleIColumn("TrackID");          
  anaMan->CreateNtupleIColumn("BoundaryID");       
  anaMan->CreateNtupleIColumn("IsReflected");      
  anaMan->CreateNtupleDColumn("HitX_mm");          
  anaMan->CreateNtupleDColumn("HitY_mm");          
  anaMan->CreateNtupleDColumn("HitZ_mm");          
  anaMan->CreateNtupleDColumn("InDirX");           
  anaMan->CreateNtupleDColumn("InDirY");           
  anaMan->CreateNtupleDColumn("InDirZ");           
  anaMan->CreateNtupleDColumn("OutDirX");          
  anaMan->CreateNtupleDColumn("OutDirY");          
  anaMan->CreateNtupleDColumn("OutDirZ");          
  anaMan->FinishNtuple(4);

  // ==============================================================
  // TREE 5: EventSummary (Hệ thống biến đếm Counters)
  // ==============================================================
  anaMan->CreateNtuple("EventSummary", "Granular counters per event");
  anaMan->CreateNtupleIColumn("EventID");          // 0
  anaMan->CreateNtupleIColumn("N_Gen");            // 1: Tổng sinh ra trong nước
  anaMan->CreateNtupleIColumn("N_W2A");            // 2: Tổng từ Nước đâm vào Acrylic
  anaMan->CreateNtupleIColumn("N_A2W");            // 3: Tổng từ Acrylic dội ngược lại Nước
  anaMan->CreateNtupleIColumn("N_A2Air");          // 4: Tổng lọt ra ngoài Không khí thành công
  anaMan->FinishNtuple(5);
}

RunAction::~RunAction() { delete G4AnalysisManager::Instance(); }

void RunAction::BeginOfRunAction(const G4Run*) {
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  auto anaMan = G4AnalysisManager::Instance();
  anaMan->OpenFile("cherenkov_water_data"); 
}

void RunAction::EndOfRunAction(const G4Run* run) {
  if (run->GetNumberOfEvent() == 0) return;
  auto anaMan = G4AnalysisManager::Instance();
  anaMan->Write();
  anaMan->CloseFile();
}