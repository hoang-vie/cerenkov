#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4AffineTransform.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

bool IsReflectionStatus(G4OpBoundaryProcessStatus s) {
  switch (s) {
    case FresnelReflection:
    case TotalInternalReflection:
    case LambertianReflection:
    case LobeReflection:
    case SpikeReflection:
    case BackScattering:
      return true;
    default:
      return false;
  }
}

bool IsTransmissionStatus(G4OpBoundaryProcessStatus s) {
  // Only statuses that actually represent passage through the boundary count as
  // optical extraction.  In particular, diagnostic/error statuses such as
  // Undefined, NoRINDEX or StepTooSmall must never be promoted to "extracted".
  return s == FresnelRefraction || s == Transmission;
}

const char* BoundaryStatusName(G4OpBoundaryProcessStatus s) {
  switch (s) {
    case Undefined: return "Undefined";
    case FresnelRefraction: return "FresnelRefraction";
    case FresnelReflection: return "FresnelReflection";
    case TotalInternalReflection: return "TotalInternalReflection";
    case LambertianReflection: return "LambertianReflection";
    case LobeReflection: return "LobeReflection";
    case SpikeReflection: return "SpikeReflection";
    case BackScattering: return "BackScattering";
    case Absorption: return "Absorption";
    case Detection: return "Detection";
    case Transmission: return "Transmission";
    default: return "Other";
  }
}

G4double WavelengthNm(const G4StepPoint* point) {
  if (!point || point->GetTotalEnergy() <= 0.) {
    return std::numeric_limits<G4double>::quiet_NaN();
  }
  const G4double hc = h_Planck * c_light;
  return (hc / point->GetTotalEnergy()) / nm;
}

} // namespace

SteppingAction::SteppingAction(EventAction* eventAction)
  : fEventAction(eventAction) {
  fDetector = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
}

void SteppingAction::RefreshGeometryHandles() {
  if (!fDetector) {
    fDetector = static_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  }
  if (!fDetector) return;

  const auto version = fDetector->GetGeometryVersion();
  if (version == fGeometryVersion) return;

  fCoreLV = fDetector->GetCoreLV();
  fAcrylicLV = fDetector->GetAcrylicLV();
  fScintTLV = fDetector->GetScintTopLV();
  fScintBLV = fDetector->GetScintBottomLV();
  fTyvekLV = fDetector->GetTyvekLV();
  fCouplingLV = fDetector->GetCouplingLV();
  fSensorLV = fDetector->GetSensorLV();
  const auto acrylicAbsSupportRange =
      fDetector->GetDatasetWavelengthRange(fDetector->GetAcrylicAbsDataset());

  // AcrylicAbsOutOfRange means outside the empirical/measured PMMA ABS
  // validity domain, not outside the numerical transport-table support.
  const auto* pmmaProfile = fDetector->GetMaterialProfile("PMMA_SNO_A");

  if (pmmaProfile &&
      pmmaProfile->absMeasuredMinNm > 0.0 &&
      pmmaProfile->absMeasuredMaxNm > pmmaProfile->absMeasuredMinNm) {
    fAcrylicAbsMinNm = pmmaProfile->absMeasuredMinNm;
    fAcrylicAbsMaxNm = pmmaProfile->absMeasuredMaxNm;
  } else {
    // Backward-compatible fallback.
    fAcrylicAbsMinNm = acrylicAbsSupportRange.first;
    fAcrylicAbsMaxNm = acrylicAbsSupportRange.second;
  }
  fGeometryVersion = version;
}

G4OpBoundaryProcess* SteppingAction::GetBoundaryProcess() {
  if (fBoundaryProcess) return fBoundaryProcess;
  auto* pm = G4OpticalPhoton::OpticalPhotonDefinition()->GetProcessManager();
  if (!pm) return nullptr;
  auto* processes = pm->GetProcessList();
  for (G4int i = 0; i < pm->GetProcessListLength(); ++i) {
    auto* candidate = dynamic_cast<G4OpBoundaryProcess*>((*processes)[i]);
    if (candidate) {
      fBoundaryProcess = candidate;
      break;
    }
  }
  return fBoundaryProcess;
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
  if (!step || !fEventAction) return;
  RefreshGeometryHandles();
  if (!fDetector) return;

  auto* track = step->GetTrack();
  auto* pre = step->GetPreStepPoint();
  auto* post = step->GetPostStepPoint();
  if (!track || !pre || !post) return;

  auto* prePV = pre->GetPhysicalVolume();
  auto* postPV = post->GetPhysicalVolume();
  auto* preLV = prePV ? prePV->GetLogicalVolume() : nullptr;
  auto* postLV = postPV ? postPV->GetLogicalVolume() : nullptr;
  const G4int preVolumeID = fDetector->GetVolumeID(preLV);
  const G4int postVolumeID = fDetector->GetVolumeID(postLV);

  const bool isOptical =
      track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition();

  if (!isOptical) {
    // Non-optical truth is registered in TrackingAction::PreUserTrackingAction(),
    // including tracks that never take a step.  Do not duplicate registration here.
    const G4double deltaKE = pre->GetKineticEnergy() - post->GetKineticEnergy();
    fEventAction->AddTrackVolume(track->GetTrackID(), preVolumeID,
                                 step->GetTotalEnergyDeposit(),
                                 deltaKE,
                                 step->GetStepLength());

    if (preLV == fScintTLV && step->GetTotalEnergyDeposit() > 0.) {
      fEventAction->AddScintTop(track->GetTrackID(),
                                step->GetTotalEnergyDeposit()/MeV,
                                pre->GetGlobalTime()/ns);
    }
    if (preLV == fScintBLV && step->GetTotalEnergyDeposit() > 0.) {
      fEventAction->AddScintBottom(track->GetTrackID(),
                                   step->GetTotalEnergyDeposit()/MeV,
                                   pre->GetGlobalTime()/ns);
    }

    if (preLV != fCoreLV && postLV == fCoreLV) {
      fEventAction->SetCoreEntry(track->GetTrackID(), post->GetPosition());
    }
    if (preLV == fCoreLV && postLV != fCoreLV) {
      fEventAction->SetCoreExit(track->GetTrackID(), post->GetPosition());
    }

    const auto* secondaries = step->GetSecondaryInCurrentStep();
    if (secondaries) {
      G4ThreeVector referenceDir = step->GetDeltaPosition();
      if (referenceDir.mag2() > 0.) referenceDir = referenceDir.unit();
      else referenceDir = pre->GetMomentumDirection();

      for (const auto* secondary : *secondaries) {
        if (!secondary || secondary->GetDefinition() !=
                              G4OpticalPhoton::OpticalPhotonDefinition()) continue;
        auto* mutableSecondary = const_cast<G4Track*>(secondary);
        auto* info = dynamic_cast<PhotonInfo*>(mutableSecondary->GetUserInformation());
        if (!info) {
          info = new PhotonInfo();
          mutableSecondary->SetUserInformation(info);
        }
        const auto* creator = secondary->GetCreatorProcess();
        if (creator && creator->GetProcessName() == "Cerenkov") {
          info->emissionAngleDeg =
              secondary->GetMomentumDirection().angle(referenceDir)/deg;
          // G4Cerenkov::PostStepDoIt uses beta=(beta_pre+beta_post)/2 when
          // sampling cos(theta)=1/(beta*n(E)). Store that same step-level beta
          // so the recorded angle can be validated photon-by-photon without
          // reconstructing beta approximately from the track vertex.
          info->emitterBeta = 0.5*(pre->GetBeta() + post->GetBeta());
        }
      }
    }

    if (track->GetTrackStatus() != fAlive || !postPV) {
      fEventAction->UpdateTrackEnd(track, postVolumeID,
                                   post->GetProcessDefinedStep());
    }
    return;
  }

  auto* info = dynamic_cast<PhotonInfo*>(track->GetUserInformation());
  if (!info) {
    info = new PhotonInfo();
    track->SetUserInformation(info);
  }

  const G4double wavelengthNm = WavelengthNm(pre);
  if (track->GetCurrentStepNumber() == 1) {
    // PRE point of optical step #1 is the exact photon creation state.
    // Do not use the G4Track state inside RegisterPhoton: UserSteppingAction is
    // invoked after the step has been transported, so the track state is already
    // at the first-step endpoint.
    fEventAction->RegisterPhoton(track, preVolumeID, wavelengthNm,
                                 info->emissionAngleDeg, info->emitterBeta,
                                 pre->GetGlobalTime(),
                                 pre->GetPosition(),
                                 pre->GetMomentumDirection());
  }
  auto* photon = fEventAction->GetPhoton(track->GetTrackID());
  if (!photon) return;

  const bool pmmaStep =
      preLV == fAcrylicLV || postLV == fAcrylicLV ||
      (fDetector->GetCoreMaterialName() == "Acrylic" &&
       (preLV == fCoreLV || postLV == fCoreLV));

  // Numerical-only tolerance for wavelength reconstructed from photon energy.
  // Does not extend the measured optical-data validity range.
  constexpr G4double kOpticalRangeNumericalTolNm = 1e-6;

  if (pmmaStep && std::isfinite(wavelengthNm) &&
      (wavelengthNm < fAcrylicAbsMinNm - kOpticalRangeNumericalTolNm ||
      wavelengthNm > fAcrylicAbsMaxNm + kOpticalRangeNumericalTolNm)) {
    photon->acrylicAbsOutOfRange = 1;
  }

  const auto* stepProcess = post->GetProcessDefinedStep();
  if (stepProcess && stepProcess->GetProcessName() == "OpRayleigh") {
    if (photon->nRayleigh == 0) {
      photon->firstRayleighTrackLength_mm =
          static_cast<G4float>(track->GetTrackLength()/mm);

      // Stage 3D validation instrumentation: reuse the existing debug-only
      // BoundaryDebug row shape to capture the PRE/POST direction and
      // polarization of the FIRST Rayleigh scatter.  The sentinel status/name
      // keeps these rows unambiguous and the production schema unchanged.
      // No RNG calls are added here.
      if (fDetector->GetBoundaryDebug()) {
        BoundaryDebugRecord debug;
        debug.photonID = track->GetTrackID();
        debug.stepNumber = track->GetCurrentStepNumber();
        debug.preVolumeID = preVolumeID;
        debug.postVolumeID = postVolumeID;
        debug.boundaryStatus = -1001;
        debug.boundaryStatusName = "FirstRayleighScatter";
        debug.time_ns = static_cast<G4float>(post->GetGlobalTime()/ns);
        debug.x_mm = static_cast<G4float>(post->GetPosition().x()/mm);
        debug.y_mm = static_cast<G4float>(post->GetPosition().y()/mm);
        debug.z_mm = static_cast<G4float>(post->GetPosition().z()/mm);
        const auto preDir = pre->GetMomentumDirection();
        const auto postDir = post->GetMomentumDirection();
        const auto prePol = pre->GetPolarization();
        const auto postPol = post->GetPolarization();
        debug.preDirX = static_cast<G4float>(preDir.x());
        debug.preDirY = static_cast<G4float>(preDir.y());
        debug.preDirZ = static_cast<G4float>(preDir.z());
        debug.postDirX = static_cast<G4float>(postDir.x());
        debug.postDirY = static_cast<G4float>(postDir.y());
        debug.postDirZ = static_cast<G4float>(postDir.z());
        debug.prePolX = static_cast<G4float>(prePol.x());
        debug.prePolY = static_cast<G4float>(prePol.y());
        debug.prePolZ = static_cast<G4float>(prePol.z());
        debug.postPolX = static_cast<G4float>(postPol.x());
        debug.postPolY = static_cast<G4float>(postPol.y());
        debug.postPolZ = static_cast<G4float>(postPol.z());
        fEventAction->AddBoundaryDebug(debug);
      }
    }
    photon->nRayleigh += 1;
  }

  const bool atBoundary = post->GetStepStatus() == fGeomBoundary;
  G4OpBoundaryProcessStatus boundaryStatus = Undefined;
  if (atBoundary) {
    auto* boundary = GetBoundaryProcess();
    if (boundary) boundaryStatus = boundary->GetStatus();

    if (fDetector->GetBoundaryDebug()) {
      BoundaryDebugRecord debug;
      debug.photonID = track->GetTrackID();
      debug.stepNumber = track->GetCurrentStepNumber();
      debug.preVolumeID = preVolumeID;
      debug.postVolumeID = postVolumeID;
      debug.boundaryStatus = static_cast<G4int>(boundaryStatus);
      debug.boundaryStatusName = BoundaryStatusName(boundaryStatus);
      debug.time_ns = static_cast<G4float>(post->GetGlobalTime()/ns);
      debug.x_mm = static_cast<G4float>(post->GetPosition().x()/mm);
      debug.y_mm = static_cast<G4float>(post->GetPosition().y()/mm);
      debug.z_mm = static_cast<G4float>(post->GetPosition().z()/mm);
      const auto preDir = pre->GetMomentumDirection();
      const auto postDir = post->GetMomentumDirection();
      const auto prePol = pre->GetPolarization();
      const auto postPol = post->GetPolarization();
      debug.preDirX = static_cast<G4float>(preDir.x());
      debug.preDirY = static_cast<G4float>(preDir.y());
      debug.preDirZ = static_cast<G4float>(preDir.z());
      debug.postDirX = static_cast<G4float>(postDir.x());
      debug.postDirY = static_cast<G4float>(postDir.y());
      debug.postDirZ = static_cast<G4float>(postDir.z());
      debug.prePolX = static_cast<G4float>(prePol.x());
      debug.prePolY = static_cast<G4float>(prePol.y());
      debug.prePolZ = static_cast<G4float>(prePol.z());
      debug.postPolX = static_cast<G4float>(postPol.x());
      debug.postPolY = static_cast<G4float>(postPol.y());
      debug.postPolZ = static_cast<G4float>(postPol.z());
      fEventAction->AddBoundaryDebug(debug);
    }

    if (IsReflectionStatus(boundaryStatus)) {
      photon->nBoundaryReflections += 1;
      if (boundaryStatus == TotalInternalReflection) photon->nTIR += 1;
      if (preLV == fTyvekLV || postLV == fTyvekLV) {
        photon->nTyvekReflections += 1;
      }
    }

    if (preLV == fCoreLV && !photon->reachedCoreWall) {
      photon->reachedCoreWall = 1;
      photon->firstWallTime_ns = static_cast<G4float>(post->GetGlobalTime()/ns);
      photon->firstWallX_mm = static_cast<G4float>(post->GetPosition().x()/mm);
      photon->firstWallY_mm = static_cast<G4float>(post->GetPosition().y()/mm);
      photon->firstWallZ_mm = static_cast<G4float>(post->GetPosition().z()/mm);
    }

    // CoreBoundaryCrossing is intentionally stricter than Photon.FirstWall*:
    // it records only the first successful optical transmission from Core to
    // a non-Core volume. Reflections/TIR can set ReachedCoreWall but must not
    // create a CoreBoundaryCrossing row.
    if (fDetector->GetCoreBoundaryMode() == "firstExit" &&
        preLV == fCoreLV && postLV != fCoreLV &&
        IsTransmissionStatus(boundaryStatus)) {
      CoreBoundaryCrossingRecord crossing;
      crossing.photonID = track->GetTrackID();
      crossing.parentTrackID = photon->parentTrackID;
      crossing.sourceTrackID = photon->sourceTrackID;
      crossing.rootPrimaryTrackID =
          fEventAction->ResolveRootPrimary(track->GetTrackID());
      crossing.originVolumeID = photon->originVolumeID;

      const auto position = post->GetPosition();
      crossing.faceID = fDetector->ResolveCoreFace(position);
      const auto uv = fDetector->GetFaceLocalUV(crossing.faceID, position);
      crossing.x_mm = static_cast<G4float>(position.x()/mm);
      crossing.y_mm = static_cast<G4float>(position.y()/mm);
      crossing.z_mm = static_cast<G4float>(position.z()/mm);
      crossing.localU_mm = static_cast<G4float>(uv.first/mm);
      crossing.localV_mm = static_cast<G4float>(uv.second/mm);

      const auto postDir = post->GetMomentumDirection();
      crossing.dirX = static_cast<G4float>(postDir.x());
      crossing.dirY = static_cast<G4float>(postDir.y());
      crossing.dirZ = static_cast<G4float>(postDir.z());
      crossing.cosThetaNormal = static_cast<G4float>(std::clamp(
          postDir.dot(fDetector->GetFaceNormal(crossing.faceID)), -1.0, 1.0));
      crossing.wavelength_nm = static_cast<G4float>(WavelengthNm(post));
      crossing.time_ns = static_cast<G4float>(post->GetGlobalTime()/ns);

      if (fDetector->GetCoreBoundaryDetail() == "phaseSpace") {
        const auto postPol = post->GetPolarization();
        crossing.polX = static_cast<G4float>(postPol.x());
        crossing.polY = static_cast<G4float>(postPol.y());
        crossing.polZ = static_cast<G4float>(postPol.z());
      }

      fEventAction->AddCoreBoundaryCrossing(crossing);
    }

    const bool candidateExtraction =
        ((preLV == fAcrylicLV && postLV == fCouplingLV) ||
         (!fAcrylicLV && preLV == fCoreLV && postLV == fCouplingLV));
    if (candidateExtraction && IsTransmissionStatus(boundaryStatus) &&
        !photon->extractedToReadout) {
      photon->extractedToReadout = 1;
      photon->extractionTime_ns = static_cast<G4float>(post->GetGlobalTime()/ns);
    }

    if (preLV == fCouplingLV && postLV == fSensorLV && postPV) {
      SensorIncidentRecord incident;
      incident.photonID = track->GetTrackID();
      incident.channelID = postPV->GetCopyNo();
      incident.wavelength_nm = static_cast<G4float>(wavelengthNm);
      incident.time_ns = static_cast<G4float>(post->GetGlobalTime()/ns);

      const auto* touchable = post->GetTouchable();
      if (touchable && touchable->GetHistory()) {
        const G4AffineTransform& transform =
            touchable->GetHistory()->GetTopTransform();
        const auto localPos = transform.TransformPoint(post->GetPosition());
        const auto localDir = transform.TransformAxis(pre->GetMomentumDirection());
        const auto localPol = transform.TransformAxis(pre->GetPolarization());
        incident.localX_mm = static_cast<G4float>(localPos.x()/mm);
        incident.localY_mm = static_cast<G4float>(localPos.y()/mm);
        incident.localZ_mm = static_cast<G4float>(localPos.z()/mm);
        incident.entranceX_mm = incident.localX_mm;
        incident.entranceY_mm = incident.localY_mm;
        incident.entranceZ_mm = static_cast<G4float>(
            (localPos.z() + 0.5*fDetector->GetReadoutSize().z())/mm);
        // Local +z points outwards on all six readout faces.  Do not take
        // abs(localDir.z()): a wrong face transform must remain visible.
        const G4double signedCosTheta = std::clamp(localDir.z(), -1.0, 1.0);
        incident.forward = signedCosTheta > 0.0 ? 1 : 0;
        incident.incidentAngle_deg =
            static_cast<G4float>(std::acos(signedCosTheta)/deg);

        G4String detail = fDetector->GetIncidentDetail();
        if (detail == "auto") {
          detail = (fDetector->GetSensorMode() == "ConcentratorInput")
                       ? "phaseSpace" : "minimal";
        }
        if (detail == "phaseSpace") {
          incident.localDirX = static_cast<G4float>(localDir.x());
          incident.localDirY = static_cast<G4float>(localDir.y());
          incident.localDirZ = static_cast<G4float>(localDir.z());
          incident.localPolX = static_cast<G4float>(localPol.x());
          incident.localPolY = static_cast<G4float>(localPol.y());
          incident.localPolZ = static_cast<G4float>(localPol.z());
        } else {
          const auto nan = std::numeric_limits<G4float>::quiet_NaN();
          incident.localDirX = nan;
          incident.localDirY = nan;
          incident.localDirZ = nan;
          incident.localPolX = nan;
          incident.localPolY = nan;
          incident.localPolZ = nan;
        }
      }
      incident.nBoundaryReflections = photon->nBoundaryReflections;
      incident.nTIR = photon->nTIR;
      incident.nTyvekReflections = photon->nTyvekReflections;
      incident.nRayleigh = photon->nRayleigh;
      fEventAction->AddSensorIncident(incident);
    }

    if (boundaryStatus == Detection &&
        preLV == fCouplingLV && postLV == fSensorLV && postPV) {
      DetectionRecord detection;
      detection.photonID = track->GetTrackID();
      detection.channelID = postPV->GetCopyNo();
      detection.time_ns = static_cast<G4float>(post->GetGlobalTime()/ns);
      fEventAction->AddDetection(detection);
      photon->endReason = 1;
      photon->endVolumeID = postVolumeID;
      photon->endTime_ns = static_cast<G4float>(post->GetGlobalTime()/ns);
      photon->endTrackLength_mm = static_cast<G4float>(track->GetTrackLength()/mm);
      track->SetTrackStatus(fStopAndKill);
      return;
    }
  }

  const G4int maxBounces = fDetector->GetMaxOpticalBounces();
  if (maxBounces > 0 && photon->nBoundaryReflections >= maxBounces) {
    photon->killedBySafetyCut = 1;
    photon->endReason = 4;
    photon->endVolumeID = postVolumeID;
    photon->endTime_ns = static_cast<G4float>(post->GetGlobalTime()/ns);
    photon->endTrackLength_mm = static_cast<G4float>(track->GetTrackLength()/mm);
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  if (track->GetTrackStatus() != fAlive || !postPV) {
    photon->endVolumeID = postVolumeID;
    photon->endTime_ns = static_cast<G4float>(post->GetGlobalTime()/ns);
    photon->endTrackLength_mm = static_cast<G4float>(track->GetTrackLength()/mm);
    if (!postPV) photon->endReason = 3;
    else if (stepProcess && stepProcess->GetProcessName() == "OpAbsorption") {
      photon->endReason = 2;
    } else if (atBoundary && boundaryStatus == Absorption) {
      photon->endReason = 5;
    } else if (photon->endReason == 0) {
      photon->endReason = 6;
    }
  }
}
