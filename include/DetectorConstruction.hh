#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "MaterialProfile.hh"

#include <array>
#include <cstdint>
#include <map>
#include <utility>
#include <unordered_map>
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

enum ReadoutFace {
  kPlusX  = 0,
  kMinusX = 1,
  kPlusZ  = 2,
  kMinusZ = 3,
  kPlusY  = 4,
  kMinusY = 5
};

struct ReadoutPlaneConfig {
  G4bool enabled = true;
  G4bool autoFill = false;
  G4int nU = 2;
  G4int nV = 2;
  G4double pitchU = 50.0*mm;
  G4double pitchV = 50.0*mm;
  G4double offsetU = 0.0;
  G4double offsetV = 0.0;
  G4double gapU = 0.0;
  G4double gapV = 0.0;
  G4double marginU = 0.0;
  G4double marginV = 0.0;
};


struct VolumeDescriptor {
  G4int volumeID = -1;
  G4String name;
  G4String materialName;
};

struct ChannelDescriptor {
  G4int channelID = -1;
  G4int faceID = -1;
  G4int gridU = -1;
  G4int gridV = -1;
  G4double localU = 0.0;
  G4double localV = 0.0;
  G4ThreeVector center = G4ThreeVector();          // readout-volume center (legacy)
  G4ThreeVector entranceCenter = G4ThreeVector();  // physical/virtual entrance plane
  G4ThreeVector normal = G4ThreeVector();
  G4double sizeU = 0.0;
  G4double sizeV = 0.0;
};

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  void SetCouplingFile(G4String filepath);
  // Validation/expert override of Water_Default RAYLEIGH only.  This does
  // not change the profile identity or production default dataset.
  void SetWaterRayleighFile(G4String filename);
  // Validation/expert override of PMMA_SNO_A ABSLENGTH only. This leaves the
  // production profile identity and default dataset unchanged unless invoked.
  void SetPmmaAbsFile(G4String filename);
  void SetCouplingProfile(G4String profileName);
  void SetTyvekPosition(const G4String& pos) { fTyvekConfig = pos; }
  void SetCoreSize(const G4ThreeVector& size) { fCoreHalf = size / 2.0; }
  void SetCoreMaterial(G4String name);
  void SetCoreProfile(G4String profileName);
  void SetScintSize(const G4ThreeVector& size) { fScintSize = size; }
  void SetReadoutSize(const G4ThreeVector& size) { fReadoutSize = size; }
  void SetSipmSize(const G4ThreeVector& size) { SetReadoutSize(size); }
  void SetTyvekHoleSize(const G4ThreeVector& size) { fTyvekHoleSize = size; }
  void SetCouplingSize(const G4ThreeVector& size) { fCouplingSize = size; }
  // Generic outer-shell geometry controls.  The geometry role is independent
  // of the bulk-material profile; R2C routes the current PMMA shell through
  // the curated MaterialProfile registry.
  void SetOuterShellThickness(G4double value) { fOuterShellThickness = value; }
  void SetUseOuterShell(G4bool value) { fUseOuterShell = value; }
  void SetOuterShellProfile(G4String profileName);

  // Legacy aliases retained so existing macros remain regression-compatible.
  void SetAcrylicThickness(G4double value) { SetOuterShellThickness(value); }
  void SetUseAcrylicShell(G4bool value) { SetUseOuterShell(value); }
  void SetTyvekThickness(G4double value) { fTyvekThickness = value; }
  void SetScintGap(G4double value) { fScintGap = value; }
  void SetSourceDistance(G4double value) { fSourceDist = value; }

  void SetSensorMode(G4String mode);
  void SetScintThreshold(G4double value) { fScintThreshold = value; }
  void SetCoincidenceWindow(G4double value) { fCoincidenceWindow = value; }
  void SetScintIntegrationWindow(G4double value) { fScintIntegrationWindow = value; }

  void SetOutputEventMode(G4String mode);
  void SetIncidentDetail(G4String mode);
  void SetPhotonBirthDetail(G4String mode);
  void SetCoreBoundaryMode(G4String mode);
  void SetCoreBoundaryDetail(G4String mode);
  void SetBoundaryDebug(G4bool value) { fBoundaryDebug = value; }
  void SetMaxOpticalBounces(G4int value) { fMaxOpticalBounces = value < 0 ? 0 : value; }
  void SetOutputFileName(const G4String& value) { fOutputFileName = value; }
  void SetCompressionLevel(G4int value);

  void SetUniformGrid(G4int nU, G4int nV, G4double pitchU, G4double pitchV);
  void SetPlaneGrid(G4int face, G4bool enabled, G4int nU, G4int nV,
                    G4double pitchU, G4double pitchV,
                    G4double offsetU, G4double offsetV);
  void SetPlaneAutoFill(G4int face, G4bool enabled,
                        G4double gapU, G4double gapV,
                        G4double marginU, G4double marginV);
  G4int ParseFace(const G4String& faceName) const;
  G4String FaceName(G4int face) const;
  G4int ResolveCoreFace(const G4ThreeVector& boundaryPosition) const;
  std::pair<G4double,G4double> GetFaceLocalUV(
      G4int face, const G4ThreeVector& worldPosition) const;
  G4ThreeVector GetFaceNormal(G4int face) const { return FaceNormal(face); }

  void SetSipmPitch(G4double halfSpacing) {
    SetUniformGrid(2, 2, 2.0*halfSpacing, 2.0*halfSpacing);
  }
  void SetSiPMOffsetY(G4double offset) {
    for (auto& p : fPlanes) p.offsetV = offset;
  }

  G4double GetGunPosY() const;
  G4String GetCoreMaterialName() const { return fCoreMaterialName; }
  G4String GetCoreProfileName() const { return fCoreProfileName; }
  G4String GetMaterialProfileRegistryVersion() const { return "curated-bulk-v2_absMeasuredRange"; }
  const MaterialProfile* GetMaterialProfile(const G4String& profileName) const;
  G4bool GetUseOuterShell() const { return fUseOuterShell; }
  G4String GetOuterShellProfileName() const { return fOuterShellProfileName; }
  G4String GetCouplingProfileName() const { return fCouplingProfileName; }
  G4bool GetUseAcrylicShell() const { return GetUseOuterShell(); } // legacy alias
  G4ThreeVector GetCoreFullSize() const { return fCoreHalf * 2.0; }
  G4ThreeVector GetScintSize() const { return fScintSize; }
  G4ThreeVector GetReadoutSize() const { return fReadoutSize; }
  G4ThreeVector GetSipmSize() const { return fReadoutSize; }
  G4double GetSipmThickness() const { return fReadoutSize.z(); }
  G4ThreeVector GetTyvekHoleSize() const { return fTyvekHoleSize; }
  G4ThreeVector GetCouplingSize() const { return fCouplingSize; }
  G4double GetOuterShellThickness() const { return fOuterShellThickness; }
  G4double GetAcrylicThickness() const { return GetOuterShellThickness(); } // legacy alias
  G4double GetTyvekThickness() const { return fTyvekThickness; }
  G4double GetScintGap() const { return fScintGap; }
  G4double GetSourceDistance() const { return fSourceDist; }
  G4String GetTyvekConfig() const { return fTyvekConfig; }
  G4String GetSensorMode() const { return fSensorMode; }
  G4double GetScintThreshold() const { return fScintThreshold; }
  G4double GetCoincidenceWindow() const { return fCoincidenceWindow; }
  G4double GetScintIntegrationWindow() const { return fScintIntegrationWindow; }

  G4String GetOutputEventMode() const { return fOutputEventMode; }
  G4String GetIncidentDetail() const { return fIncidentDetail; }
  G4String GetPhotonBirthDetail() const { return fPhotonBirthDetail; }
  G4String GetCoreBoundaryMode() const { return fCoreBoundaryMode; }
  G4String GetCoreBoundaryDetail() const { return fCoreBoundaryDetail; }
  G4bool GetBoundaryDebug() const { return fBoundaryDebug; }
  G4int GetMaxOpticalBounces() const { return fMaxOpticalBounces; }
  G4String GetOutputFileName() const { return fOutputFileName; }
  G4int GetCompressionLevel() const { return fCompressionLevel; }

  const std::array<ReadoutPlaneConfig,6>& GetPlaneConfigs() const { return fPlanes; }
  const std::vector<ChannelDescriptor>& GetChannelMap() const { return fChannels; }
  const std::vector<VolumeDescriptor>& GetVolumeMap() const { return fVolumeDescriptors; }
  G4int GetChannelCount() const { return static_cast<G4int>(fChannels.size()); }
  const ChannelDescriptor* GetChannelDescriptor(G4int id) const;

  std::uint64_t GetGeometryVersion() const { return fGeometryVersion; }
  G4int GetVolumeID(const G4LogicalVolume* lv) const;

  G4LogicalVolume* GetCoreLV() const { return fLogicCore; }
  G4LogicalVolume* GetOuterShellLV() const { return fLogicOuterShell; }
  G4LogicalVolume* GetAcrylicLV() const { return GetOuterShellLV(); } // legacy alias
  G4LogicalVolume* GetScintTopLV() const { return fLogicScintTop; }
  G4LogicalVolume* GetScintBottomLV() const { return fLogicScintBottom; }
  G4LogicalVolume* GetTyvekLV() const { return fLogicTyvek; }
  G4LogicalVolume* GetCouplingLV() const { return fLogicCoupling; }
  G4LogicalVolume* GetSensorLV() const { return fLogicSensor; }
  G4LogicalVolume* GetWorldLV() const { return fLogicWorld; }

  G4String GetAirRindexDataset() const { return fAirRindexDataset; }
  G4String GetWaterRindexDataset() const;
  G4String GetWaterAbsDataset() const;
  G4String GetWaterRayleighDataset() const;
  G4String GetAcrylicRindexDataset() const;
  G4String GetAcrylicAbsDataset() const;
  G4String GetCouplingRindexDataset() const;
  G4String GetCouplingAbsDataset() const;
  G4String GetSensorPdeDataset() const { return fSensorPdeDataset; }
  G4String GetTyvekDataset() const { return fTyvekDataset; }
  std::pair<G4double,G4double> GetDatasetWavelengthRange(const G4String& dataset) const;
  G4String GetOpticalDataFingerprint() const;

private:
  void DefineMaterials();
  void RegisterMaterialProfile(const MaterialProfile& profile);
  MaterialProfile* FindMaterialProfile(const G4String& profileName);
  void ApplyMaterialProfileOptics(MaterialProfile& profile);
  void ReadOpticalData(const G4String& filename,
                       std::vector<G4double>& x,
                       std::vector<G4double>& y,
                       G4bool unitInterval = false);
  G4double ComputeScintTopY() const;
  G4String DataPath(const G4String& fileName) const;
  std::vector<std::pair<G4double,G4double>> BuildPlaneCoordinates(
      G4int face, G4double hx, G4double hy, G4double hz) const;
  G4RotationMatrix* BuildFaceRotation(G4int face);
  G4ThreeVector FaceNormal(G4int face) const;
  void RegisterVolume(G4LogicalVolume* lv, G4int id, const G4String& name);
  void ValidateConfiguration() const;

  DetectorMessenger* fMessenger = nullptr;
  G4bool fCheckOverlaps = true;

  G4Material* fCoreMat = nullptr;
  G4Material* fWater = nullptr;
  G4Material* fAcrylic = nullptr;
  G4Material* fOuterShellMat = nullptr;  // resolved from fOuterShellProfileName
  G4Material* fAir = nullptr;
  G4Material* fCouplingMat = nullptr;
  G4Material* fSensorMat = nullptr;
  G4Material* fScintMat = nullptr;
  G4Material* fDummyTyvekMat = nullptr;

  // Curated bulk-material registry. Core, OuterShell and Coupling resolve
  // through profiles; role-specific geometry remains separate from material data.
  std::map<G4String,MaterialProfile> fMaterialProfiles;

  G4LogicalVolume* fLogicWorld = nullptr;
  G4LogicalVolume* fLogicCore = nullptr;
  G4LogicalVolume* fLogicOuterShell = nullptr;
  G4LogicalVolume* fLogicScintTop = nullptr;
  G4LogicalVolume* fLogicScintBottom = nullptr;
  G4LogicalVolume* fLogicTyvek = nullptr;
  G4LogicalVolume* fLogicCoupling = nullptr;
  G4LogicalVolume* fLogicSensor = nullptr;

  std::vector<G4RotationMatrix*> fRotations;
  std::vector<ChannelDescriptor> fChannels;
  std::unordered_map<const G4LogicalVolume*,G4int> fVolumeIds;
  std::vector<VolumeDescriptor> fVolumeDescriptors;
  std::array<ReadoutPlaneConfig,6> fPlanes;
  std::uint64_t fGeometryVersion = 0;

  G4String fTyvekConfig = "Inside";
  G4String fCoreMaterialName = "Water";          // legacy human-readable label
  G4String fCoreProfileName = "Water_Default";  // canonical material profile
  G4String fOuterShellProfileName = "PMMA_SNO_A";
  G4String fCouplingProfileName = "PDMS_Generic";
  G4ThreeVector fCoreHalf;
  G4ThreeVector fScintSize;
  G4ThreeVector fReadoutSize;
  G4ThreeVector fTyvekHoleSize;
  G4ThreeVector fCouplingSize;
  G4double fOuterShellThickness = 3.0*mm;
  G4bool fUseOuterShell = true;
  G4double fTyvekThickness = 0.1*mm;
  G4double fScintGap = 15.0*mm;
  G4double fSourceDist = 20.0*mm;

  G4String fSensorMode = "DirectSiPM";
  G4double fScintThreshold = 0.0*MeV;
  G4double fCoincidenceWindow = 100.0*ns;
  G4double fScintIntegrationWindow = 100.0*ns;

  G4String fOutputEventMode = "thinAll";
  G4String fIncidentDetail = "auto";
  G4String fPhotonBirthDetail = "basic";
  G4String fCoreBoundaryMode = "off";
  G4String fCoreBoundaryDetail = "basic";
  G4bool fBoundaryDebug = false;
  G4int fMaxOpticalBounces = 0;
  G4String fOutputFileName = "output_cerenkov";
  G4int fCompressionLevel = 6;

  G4String fAirRindexDataset = "Air_RINDEX.txt";
  G4String fSensorPdeDataset = "SiPM_PDE.txt";
  G4String fTyvekDataset = "Tyvek_REF.txt";
};

#endif
