#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include <string>

RunAction::RunAction() : G4UserRunAction()
{
  auto anaMan = G4AnalysisManager::Instance();
  anaMan->SetVerboseLevel(0);
  anaMan->SetNtupleMerging(true);
  anaMan->SetDefaultFileType("root");

  // ĐÃ BỎ QUA TREE 0 (TRUTH) VÌ ĐÃ GỘP VÀO EVENT SUMMARY (Tree 3 bên dưới)

  // TREE 0: Secondary (Delta Electrons) - Dùng để xem nhiễu điện tử
  anaMan->CreateNtuple("Secondary", "Delta Electrons");
  anaMan->CreateNtupleIColumn("EventID");          
  anaMan->CreateNtupleIColumn("TrackID");          
  anaMan->CreateNtupleIColumn("ParentID");         
  anaMan->CreateNtupleIColumn("PDG_Code");         
  anaMan->CreateNtupleDColumn("Energy_MeV");       
  anaMan->CreateNtupleSColumn("Process");          
  anaMan->FinishNtuple(0);

  // TREE 1: InnerData (W2A - Lọt qua lỗ Tyvek)
  anaMan->CreateNtuple("InnerData", "Photons hitting Water-Acrylic boundary");
  anaMan->CreateNtupleIColumn("EventID");          
  anaMan->CreateNtupleIColumn("TrackID");          
  anaMan->CreateNtupleDColumn("Wavelength_nm");    
  anaMan->CreateNtupleDColumn("HitX_mm");          
  anaMan->CreateNtupleDColumn("HitY_mm");          
  anaMan->CreateNtupleDColumn("HitZ_mm");          
  anaMan->CreateNtupleDColumn("Theta_deg");        
  anaMan->CreateNtupleDColumn("Phi_deg");          
  anaMan->CreateNtupleDColumn("Angle_In_deg");     // Góc tới thực tế (Khảo sát Fresnel)
  anaMan->CreateNtupleDColumn("Time_ns");          
  anaMan->CreateNtupleIColumn("OriginVol");         
  anaMan->CreateNtupleIColumn("Wall_ID");          // 0=+X, 1=-X, 2=+Z, 3=-Z (Thay cho SurfaceID)
  anaMan->CreateNtupleIColumn("SiPM_ID");          // Kênh cửa sổ (0-15), -1 nếu trượt
  anaMan->CreateNtupleIColumn("Accumulated_Bounces"); // Số lần đập Tyvek trước khi lọt
  anaMan->FinishNtuple(1);

  // TREE 2: OuterData (A2C - Thoát khỏi vách Nhựa)
  anaMan->CreateNtuple("OuterData", "Photons exiting Acrylic to Coupling");
  anaMan->CreateNtupleIColumn("EventID");          
  anaMan->CreateNtupleIColumn("TrackID");          
  anaMan->CreateNtupleDColumn("Wavelength_nm");    
  anaMan->CreateNtupleDColumn("HitX_mm");          
  anaMan->CreateNtupleDColumn("HitY_mm");          
  anaMan->CreateNtupleDColumn("HitZ_mm");          
  anaMan->CreateNtupleDColumn("Theta_deg");        
  anaMan->CreateNtupleDColumn("Phi_deg");          
  anaMan->CreateNtupleDColumn("Angle_In_deg");     
  anaMan->CreateNtupleDColumn("Time_ns");          
  anaMan->CreateNtupleIColumn("OriginVol");         
  anaMan->CreateNtupleIColumn("Wall_ID");          
  anaMan->CreateNtupleIColumn("SiPM_ID");          
  anaMan->CreateNtupleIColumn("Accumulated_Bounces"); 
  anaMan->FinishNtuple(2);

  // TREE 3: EventSummary (Cây dữ liệu Trung tâm)
  anaMan->CreateNtuple("EventSummary", "Counters and Primary Truth");
  anaMan->CreateNtupleIColumn("EventID");          
  anaMan->CreateNtupleIColumn("Primary_PDG");      
  anaMan->CreateNtupleDColumn("Primary_E_MeV");    
  anaMan->CreateNtupleDColumn("Primary_Theta_deg");// Góc thiên đỉnh của Muon
  anaMan->CreateNtupleDColumn("Primary_Phi_deg");  // Góc phương vị của Muon
  anaMan->CreateNtupleDColumn("Edep_Water_MeV");   
  anaMan->CreateNtupleIColumn("GenWater");         
  anaMan->CreateNtupleIColumn("GenAcrylic");       
  anaMan->CreateNtupleIColumn("W2A_Total");        // Đo rò rỉ bẫy quang (So sánh với A2C)
  anaMan->CreateNtupleIColumn("A2C_Total");        
  anaMan->CreateNtupleIColumn("Total_TIR_Count");  // Đếm tổng TIR (Waveguide indicator)
  anaMan->CreateNtupleIColumn("Total_Refl_Count"); // Đếm tổng đập nhả màng phản xạ Tyvek
  anaMan->CreateNtupleIColumn("Detected_Total");   
  anaMan->CreateNtupleIColumn("IsCoincidence");    // Trigger Hodoscope
  // Tạo 16 cột biểu diễn y hệt mạch DAQ 16 kênh
  for(int i=0; i<16; i++) {
      anaMan->CreateNtupleIColumn("Hit_SiPM_" + std::to_string(i));
  }
  anaMan->FinishNtuple(3);

  // TREE 4: DetectedSignal (Danh sách photon thu được)
  anaMan->CreateNtuple("DetectedSignal", "Photons absorbed by SiPM");
  anaMan->CreateNtupleIColumn("EventID");          
  anaMan->CreateNtupleIColumn("TrackID");          
  anaMan->CreateNtupleDColumn("Wavelength_nm");    
  anaMan->CreateNtupleDColumn("Time_ns");          
  anaMan->CreateNtupleDColumn("HitX_mm");          // Giữ theo yêu cầu của bạn
  anaMan->CreateNtupleDColumn("HitY_mm");          
  anaMan->CreateNtupleDColumn("HitZ_mm");          
  anaMan->CreateNtupleIColumn("OriginVol");        
  anaMan->CreateNtupleIColumn("BounceCount");      // Trả lại đúng ý nghĩa (số lần đập nhả)
  anaMan->CreateNtupleIColumn("SiPM_ID");          // Kênh nào đã bắt được hạt (0-15)
  anaMan->CreateNtupleDColumn("TrackLength_mm");   // Quãng đường sống sót của photon
  anaMan->FinishNtuple(4);
}

RunAction::~RunAction() { delete G4AnalysisManager::Instance(); }

void RunAction::BeginOfRunAction(const G4Run*) {
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  G4AnalysisManager::Instance()->OpenFile(); 
}
void RunAction::EndOfRunAction(const G4Run* run) {
  if (run->GetNumberOfEvent() == 0) return;
  G4AnalysisManager::Instance()->Write();
  G4AnalysisManager::Instance()->CloseFile();
}