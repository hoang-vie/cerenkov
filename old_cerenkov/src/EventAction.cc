#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

EventAction::EventAction(RunAction* runAction)
  : G4UserEventAction(), fRunAction(runAction) {}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*) {
  // Xóa sạch lịch sử TrackID khi bắt đầu Event mới
  fGenSet.clear();
  fW2ASet.clear();
  fA2WSet.clear();
  fA2AirSet.clear();
}

void EventAction::EndOfEventAction(const G4Event* evt) {
  G4int evtID = evt->GetEventID();
  auto anaMan = G4AnalysisManager::Instance();
  
  // Ghi tổng kết số liệu đếm chuẩn xác vào TREE 5
  anaMan->FillNtupleIColumn(5, 0, evtID);
  anaMan->FillNtupleIColumn(5, 1, fGenSet.size());
  anaMan->FillNtupleIColumn(5, 2, fW2ASet.size());
  anaMan->FillNtupleIColumn(5, 3, fA2WSet.size());
  anaMan->FillNtupleIColumn(5, 4, fA2AirSet.size());
  anaMan->AddNtupleRow(5);

  // In log ra màn hình
  if (evtID % 100 == 0) {
    G4cout << ">>> Event " << evtID 
           << " | Gen Photons: " << fGenSet.size()
           << " | Hit SD1 (W2A): " << fW2ASet.size()
           << " | Reflected (A2W): " << fA2WSet.size()
           << " | Escaped SD2 (A2Air): " << fA2AirSet.size() << G4endl;
  }
}