#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(RunAction* runAction) : G4UserEventAction(), fRunAction(runAction) {}
EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*) {
    fPrimaryMap.clear(); 
    fBounceMap.clear(); 
    fGenWaterSet.clear(); 
    fGenAcrylicSet.clear(); 
    fGenCouplingSet.clear(); 
    
    fDeltaRayCount = 0;
    fHitScintTop = false; 
    fHitScintBottom = false;

    fW2A_Unique.clear(); fW2A_Total = 0; fA2W_Total = 0;
    fA2C_Unique.clear(); fA2C_Total = 0; fC2A_Total = 0;
    fC2Sensor_Unique.clear(); fC2Sensor_Total = 0; 
    fSensorDetect_Unique.clear();

    // Reset các biến đếm mới
    fTotalTIR = 0;
    fTotalReflection = 0;
    fSiPMHits.assign(16, 0); // Khởi tạo mảng 16 kênh SiPM với giá trị 0
}

void EventAction::EndOfEventAction(const G4Event* evt) {
    G4int evtID = evt->GetEventID();
    auto anaMan = G4AnalysisManager::Instance();
  
    // GHI TREE 0: Delta Electrons (Hạt thứ cấp) được ghi trực tiếp bên SteppingAction

    // 1. TRÍCH XUẤT THÔNG TIN HẠT MUON SƠ CẤP
    G4int priPdg = 0; 
    G4double priE = 0.0, edepW = 0.0;
    G4double priTheta = 0.0, priPhi = 0.0;
    
    // Quét tìm hạt Muon (+-13) trước tiên
    for (const auto& [trk, data] : fPrimaryMap) {
        if (std::abs(data.pdg) == 13) {
            priPdg = data.pdg; 
            priE = data.initialEnergy; 
            edepW = data.edepWater;
            priTheta = data.dir.theta() / deg;
            priPhi = data.dir.phi() / deg;
            break; 
        }
    }
    // Dự phòng an toàn: Nếu CRY không sinh Muon ở event này, lấy hạt sơ cấp đầu tiên
    if (priPdg == 0 && !fPrimaryMap.empty()) {
        auto firstPrimary = fPrimaryMap.begin()->second;
        priPdg = firstPrimary.pdg; 
        priE = firstPrimary.initialEnergy; 
        edepW = firstPrimary.edepWater;
        priTheta = firstPrimary.dir.theta() / deg;
        priPhi = firstPrimary.dir.phi() / deg;
    }

    // 2. GHI TREE 3: EventSummary (TTree số 3 theo đúng RunAction.cc)
    anaMan->FillNtupleIColumn(3, 0, evtID);
    anaMan->FillNtupleIColumn(3, 1, priPdg);
    anaMan->FillNtupleDColumn(3, 2, priE / MeV);
    anaMan->FillNtupleDColumn(3, 3, priTheta);
    anaMan->FillNtupleDColumn(3, 4, priPhi);
    anaMan->FillNtupleDColumn(3, 5, edepW / MeV);
    anaMan->FillNtupleIColumn(3, 6, fGenWaterSet.size());
    anaMan->FillNtupleIColumn(3, 7, fGenAcrylicSet.size());
    anaMan->FillNtupleIColumn(3, 8, fW2A_Total);
    anaMan->FillNtupleIColumn(3, 9, fA2C_Total);
    anaMan->FillNtupleIColumn(3, 10, fTotalTIR);
    anaMan->FillNtupleIColumn(3, 11, fTotalReflection);
    anaMan->FillNtupleIColumn(3, 12, fSensorDetect_Unique.size());
    anaMan->FillNtupleIColumn(3, 13, IsCoincidence() ? 1 : 0);
    
    // Điền 16 kênh SiPM vào các cột từ 14 đến 29
    for (int i = 0; i < 16; i++) {
        anaMan->FillNtupleIColumn(3, 14 + i, fSiPMHits[i]);
    }
    anaMan->AddNtupleRow(3);

    // IN LOG THEO DÕI
    if (evtID % 100 == 0) {
        G4cout << ">>> Event " << evtID 
           << " | Coincidence: " << (IsCoincidence() ? "YES" : "NO")
           << " | TIR Bounces: " << fTotalTIR
           << " | Detected: " << fSensorDetect_Unique.size() << G4endl;
    }
}