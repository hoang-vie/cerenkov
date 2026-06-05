#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <set> // THÊM THƯ VIỆN NÀY

class RunAction;

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    // Cập nhật các hàm đếm nhận đầu vào là TrackID
    void AddGen(G4int trackID)   { fGenSet.insert(trackID); }
    void AddW2A(G4int trackID)   { fW2ASet.insert(trackID); }
    void AddA2W(G4int trackID)   { fA2WSet.insert(trackID); }
    void AddA2Air(G4int trackID) { fA2AirSet.insert(trackID); }

  private:
    RunAction* fRunAction;
    
    // Thay thế các biến đếm nguyên thủy bằng std::set
    std::set<G4int> fGenSet;
    std::set<G4int> fW2ASet;
    std::set<G4int> fA2WSet;
    std::set<G4int> fA2AirSet;
};

#endif