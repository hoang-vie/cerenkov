#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(RunAction* runAction)
  : G4UserEventAction(), fRunAction(runAction) {}

EventAction::~EventAction() {}


void EventAction::BeginOfEventAction(const G4Event*) {
    fPrimaryMap.clear();
    fBounceMap.clear(); 
    fGenWaterSet.clear(); 
    fGenAcrylicSet.clear();
    fGenCouplingSet.clear(); 
    fDeltaRayCount = 0;

    fW2A_Unique.clear(); fW2A_Total = 0;
    fA2W_Total = 0;
    fA2C_Unique.clear(); fA2C_Total = 0;
    fC2A_Total = 0;
    fC2Sensor_Unique.clear(); fC2Sensor_Total = 0;
    fSensorDetect_Unique.clear();
}

void EventAction::EndOfEventAction(const G4Event* evt) {
    G4int evtID = evt->GetEventID();
    auto anaMan = G4AnalysisManager::Instance();
  
    // Ghi TREE 0: Bổ sung Vector hướng bay (DirX, DirY, DirZ)
    for (auto const& [trkID, data] : fPrimaryMap) {
      anaMan->FillNtupleIColumn(0, 0, evtID);
      anaMan->FillNtupleIColumn(0, 1, trkID);
      anaMan->FillNtupleIColumn(0, 2, data.pdg);
      anaMan->FillNtupleDColumn(0, 3, data.initialEnergy / MeV);
      anaMan->FillNtupleDColumn(0, 4, data.edepWater / MeV);
      anaMan->FillNtupleDColumn(0, 5, data.lenWater / mm);
      anaMan->FillNtupleDColumn(0, 6, data.edepAcrylic / MeV);
      anaMan->FillNtupleDColumn(0, 7, data.lenAcrylic / mm);
      anaMan->FillNtupleIColumn(0, 8, data.nPhotonsGen);
      anaMan->FillNtupleDColumn(0, 9, data.dir.x());
      anaMan->FillNtupleDColumn(0, 10, data.dir.y());
      anaMan->FillNtupleDColumn(0, 11, data.dir.z());
      anaMan->AddNtupleRow(0);
    }

    // TRÍCH XUẤT THÔNG TIN MUON CHO TREE 5 (Lấy TrackID = 1 làm chuẩn)
    G4int priPdg = 0; G4double priE = 0.0, edepW = 0.0;
    if (fPrimaryMap.find(1) != fPrimaryMap.end()) {
        priPdg = fPrimaryMap[1].pdg;
        priE = fPrimaryMap[1].initialEnergy;
        edepW = fPrimaryMap[1].edepWater;
    }

    // GHI TREE 5
    anaMan->FillNtupleIColumn(5, 0, evtID);
    anaMan->FillNtupleIColumn(5, 1, priPdg);
    anaMan->FillNtupleDColumn(5, 2, priE / MeV);
    anaMan->FillNtupleDColumn(5, 3, edepW / MeV);
    anaMan->FillNtupleIColumn(5, 4, fDeltaRayCount);
    anaMan->FillNtupleIColumn(5, 5, fGenWaterSet.size());
    anaMan->FillNtupleIColumn(5, 6, fGenAcrylicSet.size());
    anaMan->FillNtupleIColumn(5, 7, fGenCouplingSet.size());
    anaMan->FillNtupleIColumn(5, 8, fW2A_Unique.size());
    anaMan->FillNtupleIColumn(5, 9, fW2A_Total);
    anaMan->FillNtupleIColumn(5, 10, fA2C_Unique.size());
    anaMan->FillNtupleIColumn(5, 11, fA2C_Total);
    anaMan->FillNtupleIColumn(5, 12, fC2Sensor_Unique.size());
    anaMan->FillNtupleIColumn(5, 13, fSensorDetect_Unique.size());
    anaMan->AddNtupleRow(5);

    if (evtID % 100 == 0) {
        G4cout << ">>> Event " << evtID 
           << " | Gen(W/A): " << fGenWaterSet.size() << "/" << fGenAcrylicSet.size()
           << " | SD1(U/T): " << fW2A_Unique.size() << "/" << fW2A_Total
           << " | SD2(U/T): " << fA2C_Unique.size() << "/" << fA2C_Total
           << " | Detected: " << fSensorDetect_Unique.size() << G4endl;
    }
}