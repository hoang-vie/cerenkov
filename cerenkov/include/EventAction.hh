#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <set>
#include <map>

class RunAction;

// Cấu trúc lưu lịch sử của một hạt sơ cấp (Tia vũ trụ) khi bay qua nước
struct PrimaryData {
    G4int pdg = 0;
    G4double initialEnergy = 0.0;
    G4ThreeVector dir;
    G4double edepWater = 0.0;
    G4double lenWater = 0.0;
    G4double edepAcrylic = 0.0;
    G4double lenAcrylic = 0.0;
    G4int nPhotonsGen = 0;
};

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    // Đếm photon sinh ra
    void AddGenWater(G4int trkID)   { fGenWaterSet.insert(trkID); }
    void AddGenAcrylic(G4int trkID) { fGenAcrylicSet.insert(trkID); }
    void AddGenCoupling(G4int id) { fGenCouplingSet.insert(id); }

    // Trích xuất nguồn gốc (0 = Water, 1 = Acrylic, 2 = Gel/Air;  -1 = Khác)
    G4int GetOriginVol(G4int trkID) {
        if (fGenWaterSet.count(trkID)) return 0;
        if (fGenAcrylicSet.count(trkID)) return 1;
        if (fGenCouplingSet.count(trkID)) return 2; 
        return -1; 
    }

    // Hệ thống đếm đập nhả (Bounce Counter)
    void AddBounce(G4int trkID) { fBounceMap[trkID]++; }
    G4int GetBounce(G4int trkID) { return fBounceMap[trkID]; }

    // Đếm tương tác ranh giới (Unique & Total)
    void AddW2A(G4int trkID)   { fW2A_Unique.insert(trkID); fW2A_Total++; }
    void AddA2W(G4int trkID)   { fA2W_Total++; }
    void AddA2C(G4int trkID)   { fA2C_Unique.insert(trkID); fA2C_Total++; }
    void AddC2A(G4int trkID)   { fC2A_Total++; }
    void AddC2Sensor(G4int trkID) { fC2Sensor_Unique.insert(trkID); fC2Sensor_Total++; }
    void AddSensorDetect(G4int trkID) { fSensorDetect_Unique.insert(trkID); }

    // ── [ĐÃ SỬA & BỔ SUNG] Các hàm xử lý tracking Hạt sơ cấp & Thứ cấp ──
    void SetPrimaryInfo(G4int trkID, G4int pdg, G4double E, G4ThreeVector dir) {
        fPrimaryMap[trkID].pdg = pdg;
        fPrimaryMap[trkID].initialEnergy = E;
        fPrimaryMap[trkID].dir = dir;
    }
    
    void AddPrimaryStepWater(G4int trkID, G4double edep, G4double len) {
        fPrimaryMap[trkID].edepWater += edep;
        fPrimaryMap[trkID].lenWater += len;
    }
    
    void AddPrimaryStepAcrylic(G4int trkID, G4double edep, G4double len) {
        fPrimaryMap[trkID].edepAcrylic += edep;
        fPrimaryMap[trkID].lenAcrylic += len;
    }
    
    void AddPhotonToPrimary(G4int parentID) {
        if (fPrimaryMap.find(parentID) != fPrimaryMap.end()) {
            fPrimaryMap[parentID].nPhotonsGen++;
        }
    }
    
    void AddDeltaRay() { fDeltaRayCount++; }

  private:
    RunAction* fRunAction;
    
    std::map<G4int, PrimaryData> fPrimaryMap;
    std::map<G4int, G4int> fBounceMap; 
    G4int fDeltaRayCount = 0;

    std::set<G4int> fGenWaterSet;
    std::set<G4int> fGenAcrylicSet;
    std::set<G4int> fGenCouplingSet;

    std::set<G4int> fW2A_Unique; G4int fW2A_Total;
    G4int fA2W_Total;
    std::set<G4int> fA2C_Unique; G4int fA2C_Total;
    G4int fC2A_Total;
    std::set<G4int> fC2Sensor_Unique; G4int fC2Sensor_Total;
    std::set<G4int> fSensorDetect_Unique;
};

#endif