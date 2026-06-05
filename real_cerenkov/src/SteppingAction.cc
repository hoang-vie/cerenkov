#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"

// =========================================================================
// HÀM NỘI SUY: Tìm ID của Cửa sổ SiPM dựa trên Tọa độ (Khắc phục lỗi CopyNumber)
// =========================================================================
G4int GetSiPMIDFromPosition(G4ThreeVector pos) {
    G4double tol = 12.55 * mm; 
    std::vector<G4double> grid = {-25.0 * mm, 25.0 * mm};
    G4int copyNo = 0;
    
    for (int face = 0; face < 4; face++) {
        for (G4double yPos : grid) {
            for (G4double hPos : grid) {
                bool match = false;
                if (face == 0 && pos.x() > 49.0*mm) {         // Mặt +X
                    if (std::abs(pos.y() - yPos) <= tol && std::abs(pos.z() - hPos) <= tol) match = true;
                } else if (face == 1 && pos.x() < -49.0*mm) { // Mặt -X
                    if (std::abs(pos.y() - yPos) <= tol && std::abs(pos.z() - hPos) <= tol) match = true;
                } else if (face == 2 && pos.z() > 49.0*mm) {  // Mặt +Z
                    if (std::abs(pos.y() - yPos) <= tol && std::abs(pos.x() - hPos) <= tol) match = true;
                } else if (face == 3 && pos.z() < -49.0*mm) { // Mặt -Z
                    if (std::abs(pos.y() - yPos) <= tol && std::abs(pos.x() - hPos) <= tol) match = true;
                }
                if (match) return copyNo;
                copyNo++;
            }
        }
    }
    return -1; // Trượt ra ngoài vùng cửa sổ
}

SteppingAction::SteppingAction(EventAction* eventAction) : G4UserSteppingAction(), fEventAction(eventAction) {}
SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
  auto track     = step->GetTrack();
  auto particle  = track->GetParticleDefinition();
  auto prePoint  = step->GetPreStepPoint();
  auto postPoint = step->GetPostStepPoint();
  auto anaMan    = G4AnalysisManager::Instance();
  G4int eventID  = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  G4int trkID    = track->GetTrackID();

  G4String preVolName = prePoint->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();
  G4String postVolName = postPoint->GetTouchable() && postPoint->GetTouchable()->GetVolume() ? postPoint->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName() : "Unknown";

  // 1. TRACKING PRIMARY & HODOSCOPE TRIGGER
  if (track->GetParentID() == 0) {
      if (track->GetCurrentStepNumber() == 1) fEventAction->SetPrimaryInfo(trkID, particle->GetPDGEncoding(), prePoint->GetKineticEnergy(), track->GetMomentumDirection());
      G4double edep = step->GetTotalEnergyDeposit();
      G4double stepLen = step->GetStepLength();
      if (edep > 0 || stepLen > 0) {
          if (preVolName == "WaterBox") fEventAction->AddPrimaryStepWater(trkID, edep, stepLen);
          else if (preVolName == "AcrylicBox") fEventAction->AddPrimaryStepAcrylic(trkID, edep, stepLen);
          else if (preVolName == "ScintLogic") {
              // Discriminator Threshold: 0.5 MeV
              if (edep > 0.5 * MeV) {
                  if (prePoint->GetTouchable()->GetVolume()->GetName() == "ScintTop") fEventAction->SetHitScintTop();
                  if (prePoint->GetTouchable()->GetVolume()->GetName() == "ScintBottom") fEventAction->SetHitScintBottom();
              }
          }
      }
  }

  // 2. PHÂN LOẠI NGUỒN CHERENKOV
  if (particle == G4OpticalPhoton::OpticalPhotonDefinition() && track->GetCurrentStepNumber() == 1) {
      if (preVolName == "WaterBox") {
          fEventAction->AddGenWater(trkID);
          if (track->GetCreatorProcess() && track->GetCreatorProcess()->GetProcessName() == "Cerenkov") fEventAction->AddPhotonToPrimary(track->GetParentID());
      } else if (preVolName == "AcrylicBox") { fEventAction->AddGenAcrylic(trkID); }
      else if (preVolName == "CouplingBox") { fEventAction->AddGenCoupling(trkID); }
  }

  // 3. SECONDARY DELTA ELECTRONS (Ghi vào Tree 0)
  if (track->GetParentID() > 0 && track->GetCurrentStepNumber() == 1) {
      G4double eKin = track->GetKineticEnergy();
      if (particle->GetPDGCharge() != 0.0 && eKin > 0.0) { 
          fEventAction->AddDeltaRay(); 
          G4String processName = track->GetCreatorProcess() ? track->GetCreatorProcess()->GetProcessName() : "Unknown";
          anaMan->FillNtupleIColumn(0, 0, eventID); 
          anaMan->FillNtupleIColumn(0, 1, track->GetTrackID()); 
          anaMan->FillNtupleIColumn(0, 2, track->GetParentID());
          anaMan->FillNtupleIColumn(0, 3, particle->GetPDGEncoding()); 
          anaMan->FillNtupleDColumn(0, 4, eKin / MeV); 
          anaMan->FillNtupleSColumn(0, 5, processName);
          anaMan->AddNtupleRow(0);
      }
  }

  // 4. QUANG HỌC TẠI RANH GIỚI
  if (particle == G4OpticalPhoton::OpticalPhotonDefinition() && postPoint->GetStepStatus() == fGeomBoundary) {
    G4OpBoundaryProcessStatus boundaryStatus = Undefined;
    G4ProcessManager* pm = particle->GetProcessManager();
    if (pm) {
        G4ProcessVector* pv = pm->GetPostStepProcessVector(typeDoIt);
        for (G4int i = 0; i < pv->entries(); i++) {
            if ((*pv)[i]->GetProcessName() == "OpBoundary") { boundaryStatus = ((G4OpBoundaryProcess*)(*pv)[i])->GetStatus(); break; }
        }
    }

    // Đếm Bounce (Phản xạ toàn phần & Phản xạ Tyvek)
    if (boundaryStatus == TotalInternalReflection) {
        fEventAction->AddTIR();
        fEventAction->AddBounce(trkID);
    } else if (boundaryStatus == FresnelReflection || boundaryStatus == SpikeReflection || 
               boundaryStatus == LambertianReflection || boundaryStatus == LobeReflection) {
        fEventAction->AddReflection();
        fEventAction->AddBounce(trkID);
    }

    G4int boundary_id = 0;     // 1=W2A, 2=A2C
    G4int treeID = -1;         // Tree 1 hoặc 2
    G4int sipm_id = -1; 
    auto pos = postPoint->GetPosition();
    G4ThreeVector dir_in = prePoint->GetMomentumDirection();
    
    // Tính Góc tới và Wall_ID
    G4int wall_id = -1;
    G4double angle_in_deg = -1.0;
    if (std::abs(pos.x()) > 49.0*mm) {
        wall_id = (pos.x() > 0) ? 0 : 1; 
        angle_in_deg = std::acos(std::abs(dir_in.x())) / deg; // Vuông góc mặt X
    } else if (std::abs(pos.z()) > 49.0*mm) {
        wall_id = (pos.z() > 0) ? 2 : 3;
        angle_in_deg = std::acos(std::abs(dir_in.z())) / deg; // Vuông góc mặt Z
    }

    // Bắt sự kiện xuyên biên
    if (boundaryStatus == Transmission || boundaryStatus == FresnelRefraction) {
        if (preVolName == "WaterBox" && postVolName == "AcrylicBox") {
            fEventAction->AddW2A(trkID);
            boundary_id = 1; treeID = 1; // Tree 1
            sipm_id = GetSiPMIDFromPosition(pos); 
        } 
        else if (preVolName == "AcrylicBox" && postVolName == "CouplingBox") {
            fEventAction->AddA2C(trkID);
            boundary_id = 2; treeID = 2; // Tree 2
            sipm_id = postPoint->GetTouchable()->GetCopyNumber(); 
        }
    }

    G4int origin = fEventAction->GetOriginVol(trkID); 

    // GHI NHẬN TREE 1 (W2A) và TREE 2 (A2C)
    if (boundary_id > 0) {
      anaMan->FillNtupleIColumn(treeID, 0, eventID); 
      anaMan->FillNtupleIColumn(treeID, 1, trkID); 
      anaMan->FillNtupleDColumn(treeID, 2, (1239.84193 * eV * nm) / prePoint->GetTotalEnergy() / nm);
      anaMan->FillNtupleDColumn(treeID, 3, pos.x() / mm); 
      anaMan->FillNtupleDColumn(treeID, 4, pos.y() / mm); 
      anaMan->FillNtupleDColumn(treeID, 5, pos.z() / mm);
      anaMan->FillNtupleDColumn(treeID, 6, dir_in.theta() / deg); 
      anaMan->FillNtupleDColumn(treeID, 7, dir_in.phi() / deg);
      anaMan->FillNtupleDColumn(treeID, 8, angle_in_deg); 
      anaMan->FillNtupleDColumn(treeID, 9, postPoint->GetGlobalTime() / ns);
      anaMan->FillNtupleIColumn(treeID, 10, origin); 
      anaMan->FillNtupleIColumn(treeID, 11, wall_id);
      anaMan->FillNtupleIColumn(treeID, 12, sipm_id);
      anaMan->FillNtupleIColumn(treeID, 13, fEventAction->GetBounce(trkID)); // Bounce tích lũy
      anaMan->AddNtupleRow(treeID);
    }

    // 5. GHI NHẬN TÍN HIỆU TẠI SENSOR (TREE 4)
    if (boundaryStatus == Detection) {
        G4int final_sipm_id = postPoint->GetTouchable()->GetCopyNumber();
        fEventAction->AddSensorDetect(trkID, final_sipm_id);
        
        anaMan->FillNtupleIColumn(4, 0, eventID); 
        anaMan->FillNtupleIColumn(4, 1, trkID);
        anaMan->FillNtupleDColumn(4, 2, (1239.84193 * eV * nm) / prePoint->GetTotalEnergy() / nm);
        anaMan->FillNtupleDColumn(4, 3, postPoint->GetGlobalTime() / ns);
        anaMan->FillNtupleDColumn(4, 4, pos.x() / mm); 
        anaMan->FillNtupleDColumn(4, 5, pos.y() / mm); 
        anaMan->FillNtupleDColumn(4, 6, pos.z() / mm);
        anaMan->FillNtupleIColumn(4, 7, origin); 
        anaMan->FillNtupleIColumn(4, 8, fEventAction->GetBounce(trkID)); 
        anaMan->FillNtupleIColumn(4, 9, final_sipm_id); 
        anaMan->FillNtupleDColumn(4, 10, track->GetTrackLength() / mm);
        anaMan->AddNtupleRow(4);
        
        track->SetTrackStatus(fStopAndKill);
    }
  }
}