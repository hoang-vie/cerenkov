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
#include "G4LogicalVolumeStore.hh"
#include "G4Box.hh"

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
  G4int trkID    = track->GetTrackID();

  G4String preVolName = prePoint->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();
  G4String postVolName = postPoint->GetTouchable() && postPoint->GetTouchable()->GetVolume() ? 
                         postPoint->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName() : "Unknown";

  // ── [1] TRACKING HẠT SƠ CẤP (CỘNG DỒN NĂNG LƯỢNG) ──
  if (track->GetParentID() == 0) {
      if (track->GetCurrentStepNumber() == 1) {
          fEventAction->SetPrimaryInfo(trkID, particle->GetPDGEncoding(), prePoint->GetKineticEnergy(), track->GetMomentumDirection());
      }
      
      G4double edep = step->GetTotalEnergyDeposit();
      G4double stepLen = step->GetStepLength();
      
      if (edep > 0 || stepLen > 0) {
          if (preVolName == "WaterBox") {
              fEventAction->AddPrimaryStepWater(trkID, edep, stepLen);
          } else if (preVolName == "AcrylicBox") {
              fEventAction->AddPrimaryStepAcrylic(trkID, edep, stepLen);
          }
      }
  }

  // ── [2] PHÂN LOẠI VÀ ĐẾM NGUỒN CHERENKOV ──
  if (particle == G4OpticalPhoton::OpticalPhotonDefinition() && track->GetCurrentStepNumber() == 1) {
      if (preVolName == "WaterBox") {
          fEventAction->AddGenWater(trkID);
          if (track->GetCreatorProcess() && track->GetCreatorProcess()->GetProcessName() == "Cerenkov") {
              fEventAction->AddPhotonToPrimary(track->GetParentID());
          }
      } else if (preVolName == "AcrylicBox") {
          fEventAction->AddGenAcrylic(trkID);
      } else if (preVolName == "CouplingBox") {
          fEventAction->AddGenCoupling(trkID); 
          track->SetTrackStatus(fStopAndKill); // [QUAN TRỌNG] Tiêu diệt nhiễu từ Gel/Air
      }
  }

  // ── [3] BỘ LỌC HẠT THỨ CẤP (Bắt mọi hạt mang điện) ──
  if (track->GetParentID() > 0 && track->GetCurrentStepNumber() == 1) {
      G4double eKin = track->GetKineticEnergy();
      
      if (particle->GetPDGCharge() != 0.0 && eKin > 0.0) { 
          G4String processName = track->GetCreatorProcess() ? track->GetCreatorProcess()->GetProcessName() : "Unknown";
          anaMan->FillNtupleIColumn(1, 0, eventID);
          anaMan->FillNtupleIColumn(1, 1, track->GetTrackID());
          anaMan->FillNtupleIColumn(1, 2, track->GetParentID());
          anaMan->FillNtupleIColumn(1, 3, particle->GetPDGEncoding());
          anaMan->FillNtupleDColumn(1, 4, eKin / MeV);
          anaMan->FillNtupleSColumn(1, 5, processName);
          anaMan->AddNtupleRow(1);
      }
  }

  // ── [4] XỬ LÝ QUANG HỌC TẠI RANH GIỚI ──
  if (particle == G4OpticalPhoton::OpticalPhotonDefinition() && postPoint->GetStepStatus() == fGeomBoundary)
  {
    if (preVolName == "WaterBox" && postVolName == "AcrylicBox") fEventAction->AddW2A(trkID);
    else if (preVolName == "AcrylicBox" && postVolName == "WaterBox") fEventAction->AddA2W(trkID);
    else if (preVolName == "AcrylicBox" && postVolName == "CouplingBox") fEventAction->AddA2C(trkID);
    else if (preVolName == "CouplingBox" && postVolName == "AcrylicBox") fEventAction->AddC2A(trkID);
    else if (preVolName == "CouplingBox" && postVolName == "SensorBox") fEventAction->AddC2Sensor(trkID);

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

    if (boundaryStatus == TotalInternalReflection || 
        boundaryStatus == FresnelReflection ||
        boundaryStatus == SpikeReflection) {
        fEventAction->AddBounce(trkID);
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
    
    G4double dot_in  = std::min(1.0, std::abs(dir_in.dot(globalNormal)));
    G4double dot_out = std::min(1.0, std::abs(dir_out.dot(globalNormal)));
    G4double angle_in  = validNormal ? std::acos(dot_in) : 0.;
    G4double angle_out = validNormal ? std::acos(dot_out) : 0.;

    G4int surfaceID = 0;
    if (validNormal) {
        G4double ax = std::abs(globalNormal.x());
        G4double ay = std::abs(globalNormal.y());
        G4double az = std::abs(globalNormal.z());
        
        if (ay >= ax && ay >= az) surfaceID = (globalNormal.y() > 0) ? 1 : 2;      
        else if (az >= ax && az >= ay) surfaceID = (globalNormal.z() > 0) ? 3 : 4; 
        else surfaceID = (globalNormal.x() > 0) ? 6 : 5;                           
    }

    G4int boundary_id = 0;
    static G4double wHL = -1.0, aHL = -1.0;
    if (wHL < 0) {
        auto wLV = G4LogicalVolumeStore::GetInstance()->GetVolume("WaterBox");
        auto aLV = G4LogicalVolumeStore::GetInstance()->GetVolume("AcrylicBox");
        if (wLV) wHL = ((G4Box*)wLV->GetSolid())->GetXHalfLength();
        if (aLV) aHL = ((G4Box*)aLV->GetSolid())->GetXHalfLength();
    }
    if (wHL > 0 && aHL > 0) {
        G4double rMax = std::max({std::abs(pos.x()), std::abs(pos.y()), std::abs(pos.z())});
        if (std::abs(rMax - wHL) < 0.5) boundary_id = 1;      // SD1
        else if (std::abs(rMax - aHL) < 0.5) boundary_id = 2; // SD2
    }

    G4int is_reflected = (boundaryStatus == TotalInternalReflection || boundaryStatus == FresnelReflection) ? 1 : 0;
    G4int origin = fEventAction->GetOriginVol(trkID); // [MỚI] Lấy thẻ xuất xứ

    if (boundary_id > 0) {
        anaMan->FillNtupleIColumn(4, 0, eventID);
        anaMan->FillNtupleIColumn(4, 1, trkID);
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
        anaMan->FillNtupleIColumn(4, 13, origin); // [MỚI]
        anaMan->AddNtupleRow(4);
    }

    if (boundary_id == 1 && is_reflected == 0 && preVolName == "WaterBox") {
      anaMan->FillNtupleIColumn(2, 0, eventID);
      anaMan->FillNtupleIColumn(2, 1, trkID);
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
      anaMan->FillNtupleIColumn(2, 22, origin); // [MỚI]
      anaMan->AddNtupleRow(2);
    }
    else if (boundary_id == 2 && is_reflected == 0 && preVolName == "AcrylicBox") {
      anaMan->FillNtupleIColumn(3, 0, eventID);
      anaMan->FillNtupleIColumn(3, 1, trkID);
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
      anaMan->FillNtupleIColumn(3, 22, origin); // [MỚI]
      anaMan->AddNtupleRow(3);
    }

    // ── [5] KHI PHOTON ĐƯỢC SIPM (SENSORBOX) HẤP THỤ THÀNH CÔNG ──
    if (boundaryStatus == Detection) {
        fEventAction->AddSensorDetect(trkID);
        
        anaMan->FillNtupleIColumn(6, 0, eventID);
        anaMan->FillNtupleIColumn(6, 1, trkID);
        anaMan->FillNtupleDColumn(6, 2, wavelength / nm);
        anaMan->FillNtupleDColumn(6, 3, postPoint->GetGlobalTime() / ns);
        anaMan->FillNtupleDColumn(6, 4, pos.x() / mm);
        anaMan->FillNtupleDColumn(6, 5, pos.y() / mm);
        anaMan->FillNtupleDColumn(6, 6, pos.z() / mm);
        anaMan->FillNtupleIColumn(6, 7, origin); // [MỚI]
        anaMan->FillNtupleIColumn(6, 8, fEventAction->GetBounce(trkID)); // [MỚI]
        anaMan->AddNtupleRow(6);

        track->SetTrackStatus(fStopAndKill);
    }
  }
}