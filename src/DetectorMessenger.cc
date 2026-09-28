#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4StateManager.hh"
#include "G4Threading.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIdirectory.hh"

#include <algorithm>
#include <cctype>
#include <sstream>

namespace {
G4bool ParseBoolToken(G4String token, G4bool fallback) {
  std::transform(token.begin(), token.end(), token.begin(),
                 [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  if (token == "1" || token == "true" || token == "on" || token == "yes") return true;
  if (token == "0" || token == "false" || token == "off" || token == "no") return false;
  return fallback;
}
}

DetectorMessenger::DetectorMessenger(DetectorConstruction* detector)
  : fDetector(detector) {
  fDetectorDir = new G4UIdirectory("/detector/");
  fReadoutDir = new G4UIdirectory("/readout/");
  fTriggerDir = new G4UIdirectory("/trigger/");
  fOutputDir = new G4UIdirectory("/output/");

  fTyvekPositionCmd = new G4UIcmdWithAString("/detector/tyvekPosition", this);
  fTyvekPositionCmd->SetCandidates("Inside Outside");

  fCoreMaterialCmd = new G4UIcmdWithAString("/detector/coreMaterial", this);
  fCoreMaterialCmd->SetCandidates("Water Acrylic PMMA");
  fCoreMaterialCmd->SetGuidance(
      "Convenience material alias. Prefer /detector/coreProfile for reproducible curated profiles.");

  fCoreProfileCmd = new G4UIcmdWithAString("/detector/coreProfile", this);
  fCoreProfileCmd->SetCandidates("Water_Default PMMA_SNO_A");
  fCoreProfileCmd->SetGuidance(
      "Canonical curated bulk-material profile for the Core.");

  fCoreSizeCmd = new G4UIcmdWith3VectorAndUnit("/detector/coreSize", this);
  fCoreSizeCmd->SetDefaultUnit("mm");

  fScintSizeCmd = new G4UIcmdWith3VectorAndUnit("/detector/scintSize", this);
  fScintSizeCmd->SetDefaultUnit("mm");

  fReadoutSizeCmd = new G4UIcmdWith3VectorAndUnit("/detector/readoutSize", this);
  fReadoutSizeCmd->SetDefaultUnit("mm");

  fLegacySipmSizeCmd = new G4UIcmdWith3VectorAndUnit("/detector/sipmSize", this);
  fLegacySipmSizeCmd->SetDefaultUnit("mm");

  fLegacySipmPitchCmd = new G4UIcmdWithADoubleAndUnit("/detector/sipmPitch", this);
  fLegacySipmPitchCmd->SetDefaultUnit("mm");

  fHoleSizeCmd = new G4UIcmdWith3VectorAndUnit("/detector/tyvekHoleSize", this);
  fHoleSizeCmd->SetDefaultUnit("mm");
  fHoleSizeCmd->SetGuidance(
      "Legacy 3-vector interface: X/Y define the Tyvek readout opening; Z is ignored by geometry and is not stored in RunConfig.");

  fCouplingSizeCmd = new G4UIcmdWith3VectorAndUnit("/detector/couplingSize", this);
  fCouplingSizeCmd->SetDefaultUnit("mm");

  fOuterShellThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/outerShellThickness", this);
  fOuterShellThicknessCmd->SetDefaultUnit("mm");
  fOuterShellThicknessCmd->SetParameterName("outerShellThickness", false);
  fOuterShellThicknessCmd->SetRange("outerShellThickness>=0");
  fOuterShellThicknessCmd->SetGuidance("Thickness of the optional outer shell.");

  fUseOuterShellCmd = new G4UIcmdWithABool("/detector/useOuterShell", this);
  fUseOuterShellCmd->SetGuidance("Enable/disable the optional outer shell.");

  fOuterShellProfileCmd = new G4UIcmdWithAString("/detector/outerShellProfile", this);
  fOuterShellProfileCmd->SetCandidates("PMMA_SNO_A");
  fOuterShellProfileCmd->SetGuidance(
      "Curated bulk-material profile for OuterShell. R2C enables PMMA_SNO_A only.");

  // Legacy command aliases retained for existing macros and pre/post-R1 regression.
  fLegacyAcrylicThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/setAcrylicThickness", this);
  fLegacyAcrylicThicknessCmd->SetDefaultUnit("mm");
  fLegacyAcrylicThicknessCmd->SetGuidance("Legacy alias of /detector/outerShellThickness.");

  fLegacyUseAcrylicShellCmd = new G4UIcmdWithABool("/detector/useAcrylicShell", this);
  fLegacyUseAcrylicShellCmd->SetGuidance("Legacy alias of /detector/useOuterShell.");

  fTyvekThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/setTyvekThickness", this);
  fTyvekThicknessCmd->SetDefaultUnit("mm");

  fLegacySiPMOffsetCmd = new G4UIcmdWithADoubleAndUnit("/detector/setSiPMOffsetY", this);
  fLegacySiPMOffsetCmd->SetDefaultUnit("mm");

  fScintGapCmd = new G4UIcmdWithADoubleAndUnit("/detector/scintGap", this);
  fScintGapCmd->SetDefaultUnit("mm");

  fSourceDistanceCmd = new G4UIcmdWithADoubleAndUnit("/detector/sourceDistance", this);
  fSourceDistanceCmd->SetDefaultUnit("mm");

  fCouplingProfileCmd = new G4UIcmdWithAString("/detector/couplingProfile", this);
  fCouplingProfileCmd->SetCandidates("PDMS_Generic");
  fCouplingProfileCmd->SetGuidance(
      "Curated bulk-material profile for optical coupling. R2D enables PDMS_Generic only.");

  fCouplingFileCmd = new G4UIcmdWithAString("/detector/couplingFile", this);
  fCouplingFileCmd->SetGuidance(
      "Legacy expert override of the active coupling profile RINDEX dataset.");

  fWaterRayleighFileCmd =
      new G4UIcmdWithAString("/detector/waterRayleighFile", this);
  fWaterRayleighFileCmd->SetGuidance(
      "Validation/expert override of Water_Default RAYLEIGH using a data/ basename. "
      "Production Water_Default is unchanged unless this command is explicitly used.");

  fPmmaAbsFileCmd =
      new G4UIcmdWithAString("/detector/pmmaAbsFile", this);
  fPmmaAbsFileCmd->SetGuidance(
      "Validation/expert override of PMMA_SNO_A ABSLENGTH using a data/ basename. "
      "Production PMMA_SNO_A is unchanged unless this command is explicitly used.");

  fSensorModeCmd = new G4UIcmdWithAString("/detector/sensorMode", this);
  fSensorModeCmd->SetCandidates("DirectSiPM SiPM ConcentratorInput");

  fUniformGridCmd = new G4UIcmdWithAString("/readout/uniformGrid", this);
  fUniformGridCmd->SetGuidance("nU nV pitchU_mm pitchV_mm");

  fPlaneGridCmd = new G4UIcmdWithAString("/readout/planeGrid", this);
  fPlaneGridCmd->SetGuidance("face enabled nU nV pitchU_mm pitchV_mm offsetU_mm offsetV_mm");

  fPlaneAutoCmd = new G4UIcmdWithAString("/readout/planeAuto", this);
  fPlaneAutoCmd->SetGuidance("face enabled gapU_mm gapV_mm marginU_mm marginV_mm");

  fScintThresholdCmd = new G4UIcmdWithADoubleAndUnit("/trigger/scintThreshold", this);
  fScintThresholdCmd->SetDefaultUnit("MeV");
  fScintThresholdCmd->SetParameterName("scintThreshold", false);
  fScintThresholdCmd->SetRange("scintThreshold>=0");

  fCoincidenceWindowCmd = new G4UIcmdWithADoubleAndUnit("/trigger/coincidenceWindow", this);
  fCoincidenceWindowCmd->SetDefaultUnit("ns");
  fCoincidenceWindowCmd->SetParameterName("coincidenceWindow", false);
  fCoincidenceWindowCmd->SetRange("coincidenceWindow>=0");

  fIntegrationWindowCmd = new G4UIcmdWithADoubleAndUnit("/trigger/integrationWindow", this);
  fIntegrationWindowCmd->SetDefaultUnit("ns");
  fIntegrationWindowCmd->SetParameterName("integrationWindow", false);
  fIntegrationWindowCmd->SetRange("integrationWindow>=0");

  fEventModeCmd = new G4UIcmdWithAString("/output/eventMode", this);
  fEventModeCmd->SetCandidates("triggered thinAll fullAll");

  fIncidentDetailCmd = new G4UIcmdWithAString("/output/incidentDetail", this);
  fIncidentDetailCmd->SetCandidates("auto minimal phaseSpace");

  fPhotonBirthDetailCmd = new G4UIcmdWithAString("/output/photonBirthDetail", this);
  fPhotonBirthDetailCmd->SetCandidates("basic full");

  fCoreBoundaryModeCmd = new G4UIcmdWithAString("/output/coreBoundaryMode", this);
  fCoreBoundaryModeCmd->SetCandidates("off firstExit");
  fCoreBoundaryModeCmd->SetGuidance(
      "Core-boundary optical recorder: off or firstExit (first successful Core transmission only).");

  fCoreBoundaryDetailCmd = new G4UIcmdWithAString("/output/coreBoundaryDetail", this);
  fCoreBoundaryDetailCmd->SetCandidates("basic phaseSpace");
  fCoreBoundaryDetailCmd->SetGuidance(
      "Core-boundary detail: basic stores position/direction/time/wavelength; phaseSpace also stores polarization.");

  fBoundaryDebugCmd = new G4UIcmdWithABool("/output/boundaryDebug", this);

  fMaxBouncesCmd = new G4UIcmdWithAnInteger("/output/maxOpticalBounces", this);
  fMaxBouncesCmd->SetParameterName("maxOpticalBounces", false);
  fMaxBouncesCmd->SetRange("maxOpticalBounces>=0");

  fFileNameCmd = new G4UIcmdWithAString("/output/fileName", this);

  fCompressionCmd = new G4UIcmdWithAnInteger("/output/compressionLevel", this);
  fCompressionCmd->SetParameterName("compressionLevel", false);
  fCompressionCmd->SetRange("compressionLevel>=0 && compressionLevel<=9");

  for (auto* cmd : {
      static_cast<G4UIcommand*>(fTyvekPositionCmd),
      static_cast<G4UIcommand*>(fCoreMaterialCmd),
      static_cast<G4UIcommand*>(fCoreProfileCmd),
      static_cast<G4UIcommand*>(fCoreSizeCmd),
      static_cast<G4UIcommand*>(fScintSizeCmd),
      static_cast<G4UIcommand*>(fReadoutSizeCmd),
      static_cast<G4UIcommand*>(fLegacySipmSizeCmd),
      static_cast<G4UIcommand*>(fLegacySipmPitchCmd),
      static_cast<G4UIcommand*>(fHoleSizeCmd),
      static_cast<G4UIcommand*>(fCouplingSizeCmd),
      static_cast<G4UIcommand*>(fOuterShellThicknessCmd),
      static_cast<G4UIcommand*>(fUseOuterShellCmd),
      static_cast<G4UIcommand*>(fOuterShellProfileCmd),
      static_cast<G4UIcommand*>(fLegacyAcrylicThicknessCmd),
      static_cast<G4UIcommand*>(fLegacyUseAcrylicShellCmd),
      static_cast<G4UIcommand*>(fTyvekThicknessCmd),
      static_cast<G4UIcommand*>(fLegacySiPMOffsetCmd),
      static_cast<G4UIcommand*>(fScintGapCmd),
      static_cast<G4UIcommand*>(fSourceDistanceCmd),
      static_cast<G4UIcommand*>(fCouplingProfileCmd),
      static_cast<G4UIcommand*>(fSensorModeCmd),
      static_cast<G4UIcommand*>(fUniformGridCmd),
      static_cast<G4UIcommand*>(fPlaneGridCmd),
      static_cast<G4UIcommand*>(fPlaneAutoCmd)}) {
    cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  }

  fCouplingFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fWaterRayleighFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fPmmaAbsFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fScintThresholdCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fCoincidenceWindowCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fIntegrationWindowCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fEventModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fIncidentDetailCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fPhotonBirthDetailCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fCoreBoundaryModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fCoreBoundaryDetailCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fBoundaryDebugCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fMaxBouncesCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fCompressionCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger() {
  delete fTyvekPositionCmd;
  delete fCoreMaterialCmd;
  delete fCoreProfileCmd;
  delete fCoreSizeCmd;
  delete fScintSizeCmd;
  delete fReadoutSizeCmd;
  delete fLegacySipmSizeCmd;
  delete fLegacySipmPitchCmd;
  delete fHoleSizeCmd;
  delete fCouplingSizeCmd;
  delete fOuterShellThicknessCmd;
  delete fUseOuterShellCmd;
  delete fOuterShellProfileCmd;
  delete fLegacyAcrylicThicknessCmd;
  delete fLegacyUseAcrylicShellCmd;
  delete fTyvekThicknessCmd;
  delete fLegacySiPMOffsetCmd;
  delete fScintGapCmd;
  delete fSourceDistanceCmd;
  delete fCouplingProfileCmd;
  delete fCouplingFileCmd;
  delete fWaterRayleighFileCmd;
  delete fPmmaAbsFileCmd;
  delete fSensorModeCmd;
  delete fUniformGridCmd;
  delete fPlaneGridCmd;
  delete fPlaneAutoCmd;
  delete fScintThresholdCmd;
  delete fCoincidenceWindowCmd;
  delete fIntegrationWindowCmd;
  delete fEventModeCmd;
  delete fIncidentDetailCmd;
  delete fPhotonBirthDetailCmd;
  delete fCoreBoundaryModeCmd;
  delete fCoreBoundaryDetailCmd;
  delete fBoundaryDebugCmd;
  delete fMaxBouncesCmd;
  delete fFileNameCmd;
  delete fCompressionCmd;
  delete fDetectorDir;
  delete fReadoutDir;
  delete fTriggerDir;
  delete fOutputDir;
}

G4bool DetectorMessenger::GeometryChangeAllowed() const {
  const auto state = G4StateManager::GetStateManager()->GetCurrentState();
  if (state == G4State_PreInit) return true;
  if (state == G4State_Idle && G4Threading::IsMultithreadedApplication()) {
    G4Exception("DetectorMessenger", "DynamicGeometryMT", JustWarning,
      "Geometry changes at Idle are disabled in multi-threaded mode. "
      "Set geometry/readout commands before /run/initialize or start a new process.");
    return false;
  }
  return state == G4State_Idle;
}

void DetectorMessenger::GeometryChanged() const {
  if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_Idle) {
    // Geant4 <= 11.3 does not clear these pointer tables in ReinitializeGeometry().
    // Cleaning explicitly is also harmless on newer releases because the tables are
    // simply empty when ReinitializeGeometry() performs its own cleanup.
    G4LogicalBorderSurface::CleanSurfaceTable();
    G4LogicalSkinSurface::CleanSurfaceTable();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String value) {
  const G4bool geometryCommand =
      command == fTyvekPositionCmd || command == fCoreMaterialCmd ||
      command == fCoreProfileCmd || command == fCoreSizeCmd || command == fScintSizeCmd ||
      command == fReadoutSizeCmd || command == fLegacySipmSizeCmd || command == fLegacySipmPitchCmd ||
      command == fHoleSizeCmd || command == fCouplingSizeCmd ||
      command == fOuterShellThicknessCmd || command == fUseOuterShellCmd ||
      command == fOuterShellProfileCmd ||
      command == fLegacyAcrylicThicknessCmd || command == fLegacyUseAcrylicShellCmd ||
      command == fTyvekThicknessCmd || command == fLegacySiPMOffsetCmd ||
      command == fScintGapCmd || command == fSourceDistanceCmd ||
      command == fCouplingProfileCmd || command == fSensorModeCmd || command == fUniformGridCmd ||
      command == fPlaneGridCmd || command == fPlaneAutoCmd;

  if (geometryCommand && !GeometryChangeAllowed()) return;

  // Material-property mutation after initialization is allowed only in a
  // sequential Idle session.  In MT production, shared material tables and
  // worker physics state make between-run mutation easy to misinterpret; use
  // a new process/configuration instead (same policy as geometry scans).
  if ((command == fCouplingFileCmd || command == fWaterRayleighFileCmd ||
       command == fPmmaAbsFileCmd) &&
      G4StateManager::GetStateManager()->GetCurrentState() == G4State_Idle &&
      G4Threading::IsMultithreadedApplication()) {
    G4Exception("DetectorMessenger", "DynamicMaterialMT", JustWarning,
      "Changing optical material tables at Idle is disabled in multi-threaded mode. "
      "Set the override before /run/initialize or start a new process.");
    return;
  }

  if (command == fTyvekPositionCmd) fDetector->SetTyvekPosition(value);
  else if (command == fCoreMaterialCmd) fDetector->SetCoreMaterial(value);
  else if (command == fCoreProfileCmd) fDetector->SetCoreProfile(value);
  else if (command == fCoreSizeCmd) fDetector->SetCoreSize(fCoreSizeCmd->GetNew3VectorValue(value));
  else if (command == fScintSizeCmd) fDetector->SetScintSize(fScintSizeCmd->GetNew3VectorValue(value));
  else if (command == fReadoutSizeCmd) fDetector->SetReadoutSize(fReadoutSizeCmd->GetNew3VectorValue(value));
  else if (command == fLegacySipmSizeCmd) fDetector->SetReadoutSize(fLegacySipmSizeCmd->GetNew3VectorValue(value));
  else if (command == fLegacySipmPitchCmd) fDetector->SetSipmPitch(fLegacySipmPitchCmd->GetNewDoubleValue(value));
  else if (command == fHoleSizeCmd) fDetector->SetTyvekHoleSize(fHoleSizeCmd->GetNew3VectorValue(value));
  else if (command == fCouplingSizeCmd) fDetector->SetCouplingSize(fCouplingSizeCmd->GetNew3VectorValue(value));
  else if (command == fOuterShellThicknessCmd) fDetector->SetOuterShellThickness(fOuterShellThicknessCmd->GetNewDoubleValue(value));
  else if (command == fUseOuterShellCmd) fDetector->SetUseOuterShell(fUseOuterShellCmd->GetNewBoolValue(value));
  else if (command == fOuterShellProfileCmd) fDetector->SetOuterShellProfile(value);
  else if (command == fLegacyAcrylicThicknessCmd) fDetector->SetAcrylicThickness(fLegacyAcrylicThicknessCmd->GetNewDoubleValue(value));
  else if (command == fLegacyUseAcrylicShellCmd) fDetector->SetUseAcrylicShell(fLegacyUseAcrylicShellCmd->GetNewBoolValue(value));
  else if (command == fTyvekThicknessCmd) fDetector->SetTyvekThickness(fTyvekThicknessCmd->GetNewDoubleValue(value));
  else if (command == fLegacySiPMOffsetCmd) fDetector->SetSiPMOffsetY(fLegacySiPMOffsetCmd->GetNewDoubleValue(value));
  else if (command == fScintGapCmd) fDetector->SetScintGap(fScintGapCmd->GetNewDoubleValue(value));
  else if (command == fSourceDistanceCmd) fDetector->SetSourceDistance(fSourceDistanceCmd->GetNewDoubleValue(value));
  else if (command == fCouplingProfileCmd) fDetector->SetCouplingProfile(value);
  else if (command == fWaterRayleighFileCmd) fDetector->SetWaterRayleighFile(value);
  else if (command == fPmmaAbsFileCmd) fDetector->SetPmmaAbsFile(value);
  else if (command == fSensorModeCmd) fDetector->SetSensorMode(value);
  else if (command == fUniformGridCmd) {
    std::istringstream iss(value);
    G4int nU=0, nV=0; G4double pU=0, pV=0;
    if (iss >> nU >> nV >> pU >> pV) fDetector->SetUniformGrid(nU,nV,pU*mm,pV*mm);
    else G4Exception("DetectorMessenger", "BadUniformGrid", JustWarning,
                     "Use: /readout/uniformGrid nU nV pitchU_mm pitchV_mm");
  }
  else if (command == fPlaneGridCmd) {
    std::istringstream iss(value);
    G4String faceToken, enabledToken;
    G4int nU=0, nV=0; G4double pU=0,pV=0,oU=0,oV=0;
    if (iss >> faceToken >> enabledToken >> nU >> nV >> pU >> pV >> oU >> oV) {
      const G4int face = fDetector->ParseFace(faceToken);
      if (face >= 0) fDetector->SetPlaneGrid(face, ParseBoolToken(enabledToken,true),
                                             nU,nV,pU*mm,pV*mm,oU*mm,oV*mm);
      else G4Exception("DetectorMessenger","BadFace",JustWarning,"Face must be +X -X +Z -Z +Y or -Y.");
    } else {
      G4Exception("DetectorMessenger","BadPlaneGrid",JustWarning,
                  "Use: /readout/planeGrid face enabled nU nV pitchU_mm pitchV_mm offsetU_mm offsetV_mm");
    }
  }
  else if (command == fPlaneAutoCmd) {
    std::istringstream iss(value);
    G4String faceToken, enabledToken;
    G4double gU=0,gV=0,mU=0,mV=0;
    if (iss >> faceToken >> enabledToken >> gU >> gV >> mU >> mV) {
      const G4int face = fDetector->ParseFace(faceToken);
      if (face >= 0) fDetector->SetPlaneAutoFill(face, ParseBoolToken(enabledToken,true),
                                                 gU*mm,gV*mm,mU*mm,mV*mm);
      else G4Exception("DetectorMessenger","BadFace",JustWarning,"Face must be +X -X +Z -Z +Y or -Y.");
    } else {
      G4Exception("DetectorMessenger","BadPlaneAuto",JustWarning,
                  "Use: /readout/planeAuto face enabled gapU_mm gapV_mm marginU_mm marginV_mm");
    }
  }
  else if (command == fCouplingFileCmd) fDetector->SetCouplingFile(value);
  else if (command == fScintThresholdCmd) fDetector->SetScintThreshold(fScintThresholdCmd->GetNewDoubleValue(value));
  else if (command == fCoincidenceWindowCmd) fDetector->SetCoincidenceWindow(fCoincidenceWindowCmd->GetNewDoubleValue(value));
  else if (command == fIntegrationWindowCmd) fDetector->SetScintIntegrationWindow(fIntegrationWindowCmd->GetNewDoubleValue(value));
  else if (command == fEventModeCmd) fDetector->SetOutputEventMode(value);
  else if (command == fIncidentDetailCmd) fDetector->SetIncidentDetail(value);
  else if (command == fPhotonBirthDetailCmd) fDetector->SetPhotonBirthDetail(value);
  else if (command == fCoreBoundaryModeCmd) fDetector->SetCoreBoundaryMode(value);
  else if (command == fCoreBoundaryDetailCmd) fDetector->SetCoreBoundaryDetail(value);
  else if (command == fBoundaryDebugCmd) fDetector->SetBoundaryDebug(fBoundaryDebugCmd->GetNewBoolValue(value));
  else if (command == fMaxBouncesCmd) fDetector->SetMaxOpticalBounces(fMaxBouncesCmd->GetNewIntValue(value));
  else if (command == fFileNameCmd) fDetector->SetOutputFileName(value);
  else if (command == fCompressionCmd) fDetector->SetCompressionLevel(fCompressionCmd->GetNewIntValue(value));

  if (geometryCommand) GeometryChanged();
}
