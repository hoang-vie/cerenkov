#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3VectorAndUnit;

class DetectorMessenger : public G4UImessenger {
public:
  explicit DetectorMessenger(DetectorConstruction* detector);
  ~DetectorMessenger() override;
  void SetNewValue(G4UIcommand* command, G4String value) override;

private:
  G4bool GeometryChangeAllowed() const;
  void GeometryChanged() const;

  DetectorConstruction* fDetector = nullptr;

  G4UIdirectory* fDetectorDir = nullptr;
  G4UIdirectory* fReadoutDir = nullptr;
  G4UIdirectory* fTriggerDir = nullptr;
  G4UIdirectory* fOutputDir = nullptr;

  G4UIcmdWithAString* fTyvekPositionCmd = nullptr;
  G4UIcmdWithAString* fCoreMaterialCmd = nullptr;
  G4UIcmdWithAString* fCoreProfileCmd = nullptr;
  G4UIcmdWith3VectorAndUnit* fCoreSizeCmd = nullptr;
  G4UIcmdWith3VectorAndUnit* fScintSizeCmd = nullptr;
  G4UIcmdWith3VectorAndUnit* fReadoutSizeCmd = nullptr;
  G4UIcmdWith3VectorAndUnit* fLegacySipmSizeCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fLegacySipmPitchCmd = nullptr;
  G4UIcmdWith3VectorAndUnit* fHoleSizeCmd = nullptr;
  G4UIcmdWith3VectorAndUnit* fCouplingSizeCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fOuterShellThicknessCmd = nullptr;
  G4UIcmdWithABool* fUseOuterShellCmd = nullptr;
  G4UIcmdWithAString* fOuterShellProfileCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fLegacyAcrylicThicknessCmd = nullptr;
  G4UIcmdWithABool* fLegacyUseAcrylicShellCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fTyvekThicknessCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fLegacySiPMOffsetCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fScintGapCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fSourceDistanceCmd = nullptr;
  G4UIcmdWithAString* fCouplingProfileCmd = nullptr;
  G4UIcmdWithAString* fCouplingFileCmd = nullptr;
  G4UIcmdWithAString* fWaterRayleighFileCmd = nullptr;
  G4UIcmdWithAString* fPmmaAbsFileCmd = nullptr;
  G4UIcmdWithAString* fSensorModeCmd = nullptr;

  G4UIcmdWithAString* fUniformGridCmd = nullptr;
  G4UIcmdWithAString* fPlaneGridCmd = nullptr;
  G4UIcmdWithAString* fPlaneAutoCmd = nullptr;

  G4UIcmdWithADoubleAndUnit* fScintThresholdCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fCoincidenceWindowCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fIntegrationWindowCmd = nullptr;

  G4UIcmdWithAString* fEventModeCmd = nullptr;
  G4UIcmdWithAString* fIncidentDetailCmd = nullptr;
  G4UIcmdWithAString* fPhotonBirthDetailCmd = nullptr;
  G4UIcmdWithAString* fCoreBoundaryModeCmd = nullptr;
  G4UIcmdWithAString* fCoreBoundaryDetailCmd = nullptr;
  G4UIcmdWithABool* fBoundaryDebugCmd = nullptr;
  G4UIcmdWithAnInteger* fMaxBouncesCmd = nullptr;
  G4UIcmdWithAString* fFileNameCmd = nullptr;
  G4UIcmdWithAnInteger* fCompressionCmd = nullptr;
};

#endif
