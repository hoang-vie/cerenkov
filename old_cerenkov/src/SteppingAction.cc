#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4VProcess.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
  : G4UserSteppingAction(), fEventAction(eventAction) {}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  auto track     = step->GetTrack();
  auto particle  = track->GetParticleDefinition();
  auto prePoint  = step->GetPreStepPoint();
  auto postPoint = step->GetPostStepPoint();
  auto anaMan    = G4AnalysisManager::Instance();
  G4int eventID  = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  // ── 1. ĐẾM PHOTON CHERENKOV KHI VỪA SINH RA ──
  if (particle == G4OpticalPhoton::OpticalPhotonDefinition() && track->GetCurrentStepNumber() == 1) {
      fEventAction->AddGen(track->GetTrackID()); 
  }

  // ── 2. BỘ LỌC HẠT THỨ CẤP CHUYÊN BIỆT (ELECTRON / POSITRON) ──
  if (track->GetParentID() > 0 && track->GetCurrentStepNumber() == 1) {
      G4int pdgCode = particle->GetPDGEncoding();
      G4double eKin = track->GetKineticEnergy();
      if ((std::abs(pdgCode) == 11) && eKin > 0.26 * MeV) {
          G4String processName = track->GetCreatorProcess() ? track->GetCreatorProcess()->GetProcessName() : "Unknown";
          anaMan->FillNtupleIColumn(1, 0, eventID);
          anaMan->FillNtupleIColumn(1, 1, track->GetTrackID());
          anaMan->FillNtupleIColumn(1, 2, track->GetParentID());
          anaMan->FillNtupleIColumn(1, 3, pdgCode);
          anaMan->FillNtupleDColumn(1, 4, eKin / MeV);
          anaMan->FillNtupleSColumn(1, 5, processName);
          anaMan->AddNtupleRow(1);
      }
  }

  // ── 3. XỬ LÝ QUANG HỌC TẠI CÁC RANH GIỚI HÌNH HỌC ──
  if (particle == G4OpticalPhoton::OpticalPhotonDefinition() && postPoint->GetStepStatus() == fGeomBoundary)
  {
    G4String preVolName = prePoint->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();
    G4String postVolName = postPoint->GetTouchable() && postPoint->GetTouchable()->GetVolume() ? 
                           postPoint->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName() : "Unknown";

    G4int trackID = track->GetTrackID();
    
    if (preVolName == "WaterBox" && postVolName == "AcrylicBox") {
        fEventAction->AddW2A(trackID);
    } else if (preVolName == "AcrylicBox" && postVolName == "WaterBox") {
        fEventAction->AddA2W(trackID);
    } else if (preVolName == "AcrylicBox" && postVolName == "CouplingBox") {  // Đã sửa thành CouplingBox
        fEventAction->AddA2Air(trackID); // Vẫn dùng biến tên là A2Air trong EventAction để tiện cho ROOT macro
    }

    G4double energy = prePoint->GetTotalEnergy();
    G4double wavelength = (1239.84193 * eV * nm) / energy; 
    
    G4OpBoundaryProcessStatus boundaryStatus = Undefined;
    G4ProcessManager* pm = particle->GetProcessManager();
    if (pm) {
        G4ProcessVector* pv = pm->GetPostStepProcessVector(typeDoIt);
        for (G4int i = 0; i < pv->entries(); i++) {
            G4VProcess* p = (*pv)[i];
            if (p->GetProcessName() == "OpBoundary") {
                boundaryStatus = ((G4OpBoundaryProcess*)p)->GetStatus();
                break;
            }
        }
    }

    auto pos = postPoint->GetPosition();
    G4ThreeVector dir_in = prePoint->GetMomentumDirection();
    G4ThreeVector dir_out = track->GetMomentumDirection(); 
    G4ThreeVector pol = track->GetPolarization(); 
    
    G4double theta = dir_in.theta();
    G4double phi   = dir_in.phi();

    G4Navigator* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    G4bool validNormal;
    G4ThreeVector localNormal = navigator->GetLocalExitNormal(&validNormal);
    G4ThreeVector globalNormal = validNormal ? 
        prePoint->GetTouchableHandle()->GetHistory()->GetTopTransform().Inverse().TransformAxis(localNormal) : G4ThreeVector();
    
    G4double angle_in = validNormal ? std::acos(std::abs(dir_in.dot(globalNormal))) : 0.;
    G4double angle_out = validNormal ? std::acos(std::abs(dir_out.dot(globalNormal))) : 0.;

    G4int surfaceID = 0;
    if (validNormal) {
        G4ThreeVector absNormal(std::abs(globalNormal.x()), std::abs(globalNormal.y()), std::abs(globalNormal.z()));
        if (absNormal.y() > 0.9) surfaceID = (globalNormal.y() > 0) ? 1 : 2;
        else if (absNormal.z() > 0.9) surfaceID = (globalNormal.z() > 0) ? 3 : 4;
        else if (absNormal.x() > 0.9) surfaceID = (globalNormal.x() > 0) ? 6 : 5;
    }

    // =========================================================================
    // XÁC ĐỊNH BỀ MẶT (BOUNDARY) ĐỂ GHI DỮ LIỆU
    // =========================================================================
    G4int boundary_id = 0;
    if ((preVolName == "WaterBox" && postVolName == "AcrylicBox") || 
        (preVolName == "AcrylicBox" && postVolName == "WaterBox")) {
        boundary_id = 1; // SD1
    } 
    else if ((preVolName == "AcrylicBox" && postVolName == "CouplingBox") || 
             (preVolName == "CouplingBox" && postVolName == "AcrylicBox")) {
        boundary_id = 2; // SD2 (Bây giờ là giao giữa Acrylic và lớp Coupling)
    }

    G4int is_reflected = (boundaryStatus == TotalInternalReflection || boundaryStatus == FresnelReflection) ? 1 : 0;

    // GHI TREE 4 (OpticalTrace)
    if (boundary_id > 0) {
        anaMan->FillNtupleIColumn(4, 0, eventID);
        anaMan->FillNtupleIColumn(4, 1, track->GetTrackID());
        anaMan->FillNtupleIColumn(4, 2, boundary_id);
        anaMan->FillNtupleIColumn(4, 3, is_reflected);
        anaMan->FillNtupleDColumn(4, 4, pos.x() / mm);
        anaMan->FillNtupleDColumn(4, 5, pos.y() / mm);
        anaMan->FillNtupleDColumn(4, 6, pos.z() / mm);
        anaMan->FillNtupleDColumn(4, 7, dir_in.x());
        anaMan->FillNtupleDColumn(4, 8, dir_in.y());
        anaMan->FillNtupleDColumn(4, 9, dir_in.z());
        anaMan->FillNtupleDColumn(4, 10, dir_out.x());
        anaMan->FillNtupleDColumn(4, 11, dir_out.y());
        anaMan->FillNtupleDColumn(4, 12, dir_out.z());
        anaMan->AddNtupleRow(4);
    }

    // GHI TREE 2 (Qua SD1 thành công)
    if (boundary_id == 1 && is_reflected == 0 && preVolName == "WaterBox") {
      anaMan->FillNtupleIColumn(2, 0, eventID);
      anaMan->FillNtupleIColumn(2, 1, track->GetTrackID());
      anaMan->FillNtupleIColumn(2, 2, track->GetParentID());
      anaMan->FillNtupleDColumn(2, 3, energy / eV);
      anaMan->FillNtupleDColumn(2, 4, wavelength / nm); 
      anaMan->FillNtupleDColumn(2, 5, pos.x() / mm);
      anaMan->FillNtupleDColumn(2, 6, pos.y() / mm);
      anaMan->FillNtupleDColumn(2, 7, pos.z() / mm);
      anaMan->FillNtupleDColumn(2, 8, dir_in.x());
      anaMan->FillNtupleDColumn(2, 9, dir_in.y());
      anaMan->FillNtupleDColumn(2, 10, dir_in.z());
      anaMan->FillNtupleDColumn(2, 11, theta / deg);
      anaMan->FillNtupleDColumn(2, 12, phi / deg);
      anaMan->FillNtupleDColumn(2, 13, angle_in / deg); 
      anaMan->FillNtupleDColumn(2, 14, postPoint->GetGlobalTime() / ns);
      anaMan->FillNtupleIColumn(2, 15, surfaceID);
      anaMan->FillNtupleDColumn(2, 16, globalNormal.x());
      anaMan->FillNtupleDColumn(2, 17, globalNormal.y());
      anaMan->FillNtupleDColumn(2, 18, globalNormal.z());
      anaMan->FillNtupleDColumn(2, 19, pol.x());
      anaMan->FillNtupleDColumn(2, 20, pol.y());
      anaMan->FillNtupleDColumn(2, 21, pol.z());
      anaMan->AddNtupleRow(2);
    }
    // GHI TREE 3 (Qua SD2 thành công)
    else if (boundary_id == 2 && is_reflected == 0 && postVolName == "CouplingBox") {
      anaMan->FillNtupleIColumn(3, 0, eventID);
      anaMan->FillNtupleIColumn(3, 1, track->GetTrackID());
      anaMan->FillNtupleIColumn(3, 2, track->GetParentID());
      anaMan->FillNtupleDColumn(3, 3, energy / eV);
      anaMan->FillNtupleDColumn(3, 4, wavelength / nm);
      anaMan->FillNtupleDColumn(3, 5, pos.x() / mm);
      anaMan->FillNtupleDColumn(3, 6, pos.y() / mm);
      anaMan->FillNtupleDColumn(3, 7, pos.z() / mm);
      anaMan->FillNtupleDColumn(3, 8, dir_in.x());
      anaMan->FillNtupleDColumn(3, 9, dir_in.y());
      anaMan->FillNtupleDColumn(3, 10, dir_in.z());
      anaMan->FillNtupleDColumn(3, 11, theta / deg);
      anaMan->FillNtupleDColumn(3, 12, phi / deg);
      anaMan->FillNtupleDColumn(3, 13, angle_out / deg); 
      anaMan->FillNtupleDColumn(3, 14, postPoint->GetGlobalTime() / ns);
      anaMan->FillNtupleIColumn(3, 15, surfaceID);
      anaMan->FillNtupleDColumn(3, 16, globalNormal.x());
      anaMan->FillNtupleDColumn(3, 17, globalNormal.y());
      anaMan->FillNtupleDColumn(3, 18, globalNormal.z());
      anaMan->FillNtupleDColumn(3, 19, pol.x());
      anaMan->FillNtupleDColumn(3, 20, pol.y());
      anaMan->FillNtupleDColumn(3, 21, pol.z());
      anaMan->AddNtupleRow(3);

      // KILL PHOTON NGAY LẬP TỨC ĐỂ KHÔNG TÍNH SD3 (Coupling -> Không khí)
      track->SetTrackStatus(fStopAndKill); 
    }
  }
}