#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4Exception.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MultiUnion.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string>

#ifndef OPTICAL_DATA_DIR
#define OPTICAL_DATA_DIR "../data"
#endif

DetectorConstruction::DetectorConstruction() {
  fCoreHalf = G4ThreeVector(50.0*mm, 50.0*mm, 50.0*mm);
  fScintSize = G4ThreeVector(110.0*mm, 10.0*mm, 110.0*mm);
  fReadoutSize = G4ThreeVector(25.0*mm, 25.0*mm, 1.0*mm);
  fTyvekHoleSize = G4ThreeVector(25.0*mm, 25.0*mm, 0.1*mm);
  fCouplingSize = G4ThreeVector(25.0*mm, 25.0*mm, 1.0*mm);

  for (auto& p : fPlanes) {
    p.enabled = true;
    p.autoFill = false;
    p.nU = 2;
    p.nV = 2;
    p.pitchU = 50.0*mm;
    p.pitchV = 50.0*mm;
  }

  fMessenger = new DetectorMessenger(this);
  DefineMaterials();
}

DetectorConstruction::~DetectorConstruction() {
  delete fMessenger;
  for (auto* rot : fRotations) delete rot;
}

G4String DetectorConstruction::DataPath(const G4String& fileName) const {
  return G4String(OPTICAL_DATA_DIR) + "/" + fileName;
}

void DetectorConstruction::ReadOpticalData(const G4String& filename,
                                           std::vector<G4double>& photonEnergies,
                                           std::vector<G4double>& properties,
                                           G4bool unitInterval) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    G4Exception("DetectorConstruction::ReadOpticalData", "FileOpenError",
                FatalException, ("Cannot open optical data file: " + filename).c_str());
  }

  struct DataPoint { G4double wavelengthNm; G4double energy; G4double value; };
  std::vector<DataPoint> data;
  std::set<G4double> wavelengths;
  G4String line;
  G4int lineNo = 0;
  while (std::getline(file, line)) {
    ++lineNo;
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == G4String::npos || line[first] == '#') continue;

    std::istringstream iss(line);
    G4double wavelength = 0.0;
    G4double value = 0.0;
    if (!(iss >> wavelength >> value)) {
      G4Exception("DetectorConstruction::ReadOpticalData", "MalformedData",
                  FatalException,
                  ("Malformed optical-data row in " + filename +
                   " at line " + std::to_string(lineNo)).c_str());
    }
    if (!std::isfinite(wavelength) || wavelength <= 0.0 ||
        !std::isfinite(value)) {
      G4Exception("DetectorConstruction::ReadOpticalData", "NonPhysicalData",
                  FatalException,
                  ("Non-finite/non-positive wavelength or non-finite value in " +
                   filename + " at line " + std::to_string(lineNo)).c_str());
    }
    if (unitInterval) {
      if (value < 0.0 || value > 1.0) {
        G4Exception("DetectorConstruction::ReadOpticalData", "ProbabilityOutOfRange",
                    FatalException,
                    ("Probability/reflectivity outside [0,1] in " + filename +
                     " at line " + std::to_string(lineNo)).c_str());
      }
    } else if (value <= 0.0) {
      G4Exception("DetectorConstruction::ReadOpticalData", "NonPositiveProperty",
                  FatalException,
                  ("Optical property must be > 0 in " + filename +
                   " at line " + std::to_string(lineNo)).c_str());
    }
    if (!wavelengths.insert(wavelength).second) {
      G4Exception("DetectorConstruction::ReadOpticalData", "DuplicateWavelength",
                  FatalException,
                  ("Duplicate wavelength in " + filename +
                   " at line " + std::to_string(lineNo)).c_str());
    }
    data.push_back({wavelength,
                    (1239.84193 * eV * nm) / (wavelength * nm),
                    value});
  }

  if (data.size() < 2) {
    G4Exception("DetectorConstruction::ReadOpticalData", "InsufficientData",
                FatalException, ("Need at least two optical data points in: " + filename).c_str());
  }

  std::sort(data.begin(), data.end(),
            [](const DataPoint& a, const DataPoint& b) { return a.energy < b.energy; });

  photonEnergies.clear();
  properties.clear();
  photonEnergies.reserve(data.size());
  properties.reserve(data.size());
  for (const auto& point : data) {
    photonEnergies.push_back(point.energy);
    properties.push_back(point.value);
  }
}

std::pair<G4double,G4double>
DetectorConstruction::GetDatasetWavelengthRange(const G4String& dataset) const {
  G4String path = dataset;
  std::ifstream direct(path);
  if (!direct.is_open()) {
    direct.clear();
    path = DataPath(dataset);
    direct.open(path);
  }
  if (!direct.is_open()) {
    G4Exception("DetectorConstruction::GetDatasetWavelengthRange", "FileOpenError",
                FatalException, ("Cannot open optical data file: " + dataset).c_str());
  }

  G4double minNm = std::numeric_limits<G4double>::infinity();
  G4double maxNm = -std::numeric_limits<G4double>::infinity();
  G4String line;
  while (std::getline(direct, line)) {
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == G4String::npos || line[first] == '#') continue;
    std::istringstream iss(line);
    G4double wavelength = 0.0, value = 0.0;
    if (!(iss >> wavelength >> value) || !std::isfinite(wavelength) ||
        wavelength <= 0.0) continue;
    minNm = std::min(minNm, wavelength);
    maxNm = std::max(maxNm, wavelength);
  }
  if (!std::isfinite(minNm) || !std::isfinite(maxNm)) {
    G4Exception("DetectorConstruction::GetDatasetWavelengthRange", "NoValidData",
                FatalException, ("No valid wavelengths in: " + dataset).c_str());
  }
  return {minNm, maxNm};
}


G4String DetectorConstruction::GetOpticalDataFingerprint() const {
  // FNV-1a-64 over dataset label + exact file bytes, in a fixed order.  This
  // is a provenance fingerprint (not a cryptographic security primitive).
  // It makes a ROOT RunConfig distinguish two runs even if a data file was
  // later edited without changing its filename.
  const std::array<G4String,10> datasets = {{
      fAirRindexDataset, GetWaterRindexDataset(), GetWaterAbsDataset(),
      GetWaterRayleighDataset(), GetAcrylicRindexDataset(), GetAcrylicAbsDataset(),
      GetCouplingRindexDataset(), GetCouplingAbsDataset(), fSensorPdeDataset,
      fTyvekDataset}};

  std::uint64_t hash = 14695981039346656037ULL;
  constexpr std::uint64_t prime = 1099511628211ULL;
  auto mixByte = [&](unsigned char c) {
    hash ^= static_cast<std::uint64_t>(c);
    hash *= prime;
  };

  for (const auto& dataset : datasets) {
    for (const char c : dataset) mixByte(static_cast<unsigned char>(c));
    mixByte(0xff);

    G4String path = dataset;
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) {
      in.clear();
      path = DataPath(dataset);
      in.open(path, std::ios::binary);
    }
    if (!in.is_open()) {
      G4Exception("DetectorConstruction::GetOpticalDataFingerprint",
                  "FileOpenError", FatalException,
                  ("Cannot open optical data file: " + dataset).c_str());
    }
    char c = 0;
    while (in.get(c)) mixByte(static_cast<unsigned char>(c));
    mixByte(0x00);
  }

  std::ostringstream os;
  os << std::hex << std::setfill('0') << std::setw(16) << hash;
  return os.str();
}

const MaterialProfile*
DetectorConstruction::GetMaterialProfile(const G4String& profileName) const {
  auto it = fMaterialProfiles.find(profileName);
  if (it != fMaterialProfiles.end()) return &it->second;

  G4String norm = profileName;
  std::transform(norm.begin(), norm.end(), norm.begin(),
                 [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  for (const auto& kv : fMaterialProfiles) {
    G4String candidate = kv.first;
    std::transform(candidate.begin(), candidate.end(), candidate.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    if (candidate == norm) return &kv.second;
  }
  return nullptr;
}

MaterialProfile*
DetectorConstruction::FindMaterialProfile(const G4String& profileName) {
  const auto* profile =
      static_cast<const DetectorConstruction*>(this)->GetMaterialProfile(profileName);
  return const_cast<MaterialProfile*>(profile);
}

void DetectorConstruction::RegisterMaterialProfile(const MaterialProfile& profile) {
  if (profile.profileName.empty() || !profile.material) {
    G4Exception("DetectorConstruction::RegisterMaterialProfile",
                "InvalidMaterialProfile", FatalException,
                "Material profiles require a non-empty profileName and G4Material pointer.");
  }
  if (fMaterialProfiles.find(profile.profileName) != fMaterialProfiles.end()) {
    G4Exception("DetectorConstruction::RegisterMaterialProfile",
                "DuplicateMaterialProfile", FatalException,
                ("Duplicate material profile: " + profile.profileName).c_str());
  }
  fMaterialProfiles.emplace(profile.profileName, profile);
}

void DetectorConstruction::ApplyMaterialProfileOptics(MaterialProfile& profile) {
  auto* mpt = new G4MaterialPropertiesTable();

  if (!profile.rindexDataset.empty()) {
    std::vector<G4double> e, n;
    ReadOpticalData(DataPath(profile.rindexDataset), e, n);
    mpt->AddProperty("RINDEX", e, n);
  }

  if (!profile.absDataset.empty()) {
    std::vector<G4double> e, l;
    ReadOpticalData(DataPath(profile.absDataset), e, l);
    for (auto& x : l) x *= profile.absLengthUnit;
    mpt->AddProperty("ABSLENGTH", e, l);
  }

  if (!profile.rayleighDataset.empty()) {
    std::vector<G4double> e, l;
    ReadOpticalData(DataPath(profile.rayleighDataset), e, l);
    for (auto& x : l) x *= profile.rayleighLengthUnit;
    mpt->AddProperty("RAYLEIGH", e, l);
  }

  profile.material->SetMaterialPropertiesTable(mpt);
}

G4String DetectorConstruction::GetWaterRindexDataset() const {
  const auto* p = GetMaterialProfile("Water_Default");
  return p ? p->rindexDataset : "";
}

G4String DetectorConstruction::GetWaterAbsDataset() const {
  const auto* p = GetMaterialProfile("Water_Default");
  return p ? p->absDataset : "";
}

G4String DetectorConstruction::GetWaterRayleighDataset() const {
  const auto* p = GetMaterialProfile("Water_Default");
  return p ? p->rayleighDataset : "";
}

G4String DetectorConstruction::GetAcrylicRindexDataset() const {
  const auto* p = GetMaterialProfile("PMMA_SNO_A");
  return p ? p->rindexDataset : "";
}

G4String DetectorConstruction::GetAcrylicAbsDataset() const {
  const auto* p = GetMaterialProfile("PMMA_SNO_A");
  return p ? p->absDataset : "";
}

G4String DetectorConstruction::GetCouplingRindexDataset() const {
  const auto* p = GetMaterialProfile(fCouplingProfileName);
  return p ? p->rindexDataset : "";
}

G4String DetectorConstruction::GetCouplingAbsDataset() const {
  const auto* p = GetMaterialProfile(fCouplingProfileName);
  return p ? p->absDataset : "";
}

void DetectorConstruction::DefineMaterials() {
  auto* nist = G4NistManager::Instance();

  fAir = nist->FindOrBuildMaterial("G4_AIR");
  fWater = nist->FindOrBuildMaterial("G4_WATER");
  fAcrylic = nist->FindOrBuildMaterial("G4_PLEXIGLASS");
  // Generic silicone/PDMS optical coupling baseline.  Do not alias the
  // coupling to fused silica: doing so corrupts material provenance and the
  // charged-particle interaction model if a secondary ever traverses it.
  if (auto* existing = G4Material::GetMaterial("OpticalGrease_PDMS_Generic", false)) {
    fCouplingMat = existing;
  } else {
    auto* H  = nist->FindOrBuildElement("H");
    auto* C  = nist->FindOrBuildElement("C");
    auto* O  = nist->FindOrBuildElement("O");
    auto* Si = nist->FindOrBuildElement("Si");
    fCouplingMat = new G4Material("OpticalGrease_PDMS_Generic",
                                  0.97*g/cm3, 4);
    // PDMS repeat unit: [-Si(CH3)2-O-] = C2 H6 O Si.
    fCouplingMat->AddElement(C,  2);
    fCouplingMat->AddElement(H,  6);
    fCouplingMat->AddElement(O,  1);
    fCouplingMat->AddElement(Si, 1);
  }
  fSensorMat = nist->FindOrBuildMaterial("G4_Si");
  fScintMat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  {
    auto* mpt = new G4MaterialPropertiesTable();
    std::vector<G4double> e, v;
    ReadOpticalData(DataPath(fAirRindexDataset), e, v);
    mpt->AddProperty("RINDEX", e, v);
    fAir->SetMaterialPropertiesTable(mpt);
  }

  // Register curated bulk-material profiles. Water_Default Rayleigh was
  // promoted in Stage 3E from the historical seawater-like baseline to the
  // validated Zhang-Hu PMH pure-water molecular-scattering table; the other
  // profile datasets and unit conventions are unchanged.
  fMaterialProfiles.clear();

  MaterialProfile water;
  water.profileName = "Water_Default";
  water.materialLabel = "Water";
  water.material = fWater;
  water.rindexDataset = "Water_RINDEX_Daimon2007_21p5C.txt";
  water.absDataset = "Water_ABS_Fewell2019_Pure23C.txt";
  water.rayleighDataset = "Water_RAYLEIGH_ZhangHu2009_PMH_21p5C.txt";
  water.absLengthUnit = m;
  water.rayleighLengthUnit = m;
  water.rindexStatus = "Daimon2007_21p5C_SellmeierSampling";
  water.absStatus = "Fewell2019_PureWater";
  water.rayleighStatus = "ZhangHu2009_PMH_21p5C_delta0p039_MFPValidated_AngularApproxDocumented";
  water.outOfRangePolicy = "Geant4EndpointClamp";
  RegisterMaterialProfile(water);

  MaterialProfile pmma;
  pmma.profileName = "PMMA_SNO_A";
  pmma.materialLabel = "Acrylic";
  pmma.material = fAcrylic;
  pmma.rindexDataset = "PMMA_RINDEX_Roychowdhury2020_Cauchy3Anchor.txt";
  pmma.absDataset = "PMMA_SNO_UVT_SampleA_ABS_Measured300to500_Model200to900.txt";
  pmma.absLengthUnit = mm;

  // Direct empirical validity of the SNO Sample-A absorption data.
  // Numerical transport support may extend beyond this interval using
  // explicitly-labelled policy models.
  pmma.absMeasuredMinNm = 300.0;
  pmma.absMeasuredMaxNm = 500.0;
  pmma.rindexStatus = "Roychowdhury2020_Cauchy3Anchor_Approximation";
  pmma.absStatus = "SNO_UVT_SampleA_Measured300to500_plus_explicit_transport_model_200to900";
  pmma.rayleighStatus = "NotConfigured";
  pmma.outOfRangePolicy =
      "Measured300to500;"
      "UV250to300_AnchoredExp300to320;"
      "UV200to250_Hold250;"
      "Visible500to900_Hold500;"
      "OutsideSupport_Geant4EndpointClamp";
  RegisterMaterialProfile(pmma);

  MaterialProfile pdms;
  pdms.profileName = "PDMS_Generic";
  pdms.materialLabel = "OpticalGrease_PDMS_Generic";
  pdms.material = fCouplingMat;
  pdms.rindexDataset = "OpticalGrease_PDMS_RINDEX_Roychowdhury2018_Cauchy3Anchor.txt";
  pdms.absDataset = "OpticalGrease_PDMS_ABS_Generic10m.txt";
  pdms.absLengthUnit = m;
  pdms.rindexStatus = "Roychowdhury2018_Cauchy3Anchor_Approximation";
  pdms.absStatus = "Generic10m_Placeholder_NotVendorMeasurement";
  pdms.rayleighStatus = "NotConfigured";
  pdms.outOfRangePolicy = "Geant4EndpointClamp";
  RegisterMaterialProfile(pdms);

  ApplyMaterialProfileOptics(*FindMaterialProfile("Water_Default"));
  ApplyMaterialProfileOptics(*FindMaterialProfile("PMMA_SNO_A"));
  ApplyMaterialProfileOptics(*FindMaterialProfile("PDMS_Generic"));

  // OuterShell resolves through the curated profile registry.  R2E changes
  // metadata only; PMMA_SNO_A remains the sole validated shell profile here.
  SetOuterShellProfile(fOuterShellProfileName);

  // Coupling also resolves through the curated registry.  R2E changes metadata
  // only; PDMS_Generic remains the sole enabled generic coupling profile here.
  SetCouplingProfile(fCouplingProfileName);

  fDummyTyvekMat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");

  // R2B: Core now resolves through the curated profile registry.  The legacy
  // /detector/coreMaterial command is retained and maps onto these profiles.
  SetCoreProfile(fCoreProfileName);
}

void DetectorConstruction::SetCoreProfile(G4String profileName) {
  auto* profile = FindMaterialProfile(profileName);
  if (!profile) {
    G4Exception("DetectorConstruction::SetCoreProfile", "InvalidCoreProfile",
                JustWarning,
                "Supported core profiles are Water_Default and PMMA_SNO_A.");
    return;
  }
  fCoreProfileName = profile->profileName;
  fCoreMaterialName = profile->materialLabel;
  fCoreMat = profile->material;
}

void DetectorConstruction::SetCoreMaterial(G4String name) {
  G4String norm = name;
  std::transform(norm.begin(), norm.end(), norm.begin(),
                 [](unsigned char c){ return static_cast<char>(std::tolower(c)); });

  if (norm == "water") {
    SetCoreProfile("Water_Default");
  } else if (norm == "acrylic" || norm == "pmma") {
    SetCoreProfile("PMMA_SNO_A");
  } else {
    G4Exception("DetectorConstruction::SetCoreMaterial", "InvalidCoreMaterial",
                JustWarning,
                "Supported core materials are Water and Acrylic/PMMA.");
  }
}

void DetectorConstruction::SetOuterShellProfile(G4String profileName) {
  auto* profile = FindMaterialProfile(profileName);
  if (!profile) {
    G4Exception("DetectorConstruction::SetOuterShellProfile",
                "InvalidOuterShellProfile", JustWarning,
                "Unknown outer-shell material profile.");
    return;
  }

  // R2E adds profile-aware metadata while deliberately keeping role eligibility
  // unchanged.  Expanding OuterShell to additional curated profiles is a
  // separate physics/regression step rather than part of this schema refactor.
  if (profile->profileName != "PMMA_SNO_A") {
    G4Exception("DetectorConstruction::SetOuterShellProfile",
                "OuterShellProfileNotYetEnabled", JustWarning,
                "R2E keeps PMMA_SNO_A as the only enabled OuterShell profile. "
                "Additional curated shell profiles require a separate validation step.");
    return;
  }

  fOuterShellProfileName = profile->profileName;
  fOuterShellMat = profile->material;
}

void DetectorConstruction::SetCouplingProfile(G4String profileName) {
  auto* profile = FindMaterialProfile(profileName);
  if (!profile) {
    G4Exception("DetectorConstruction::SetCouplingProfile",
                "InvalidCouplingProfile", JustWarning,
                "Unknown coupling material profile.");
    return;
  }

  // R2E keeps the already-declared generic PDMS model as the only coupling
  // choice.  Additional coupling products require their own measured datasets
  // and a separate role-validation step rather than arbitrary material strings.
  if (profile->profileName != "PDMS_Generic") {
    G4Exception("DetectorConstruction::SetCouplingProfile",
                "CouplingProfileNotYetEnabled", JustWarning,
                "R2E keeps PDMS_Generic as the only enabled Coupling profile.");
    return;
  }

  fCouplingProfileName = profile->profileName;
  fCouplingMat = profile->material;
}

void DetectorConstruction::SetSensorMode(G4String mode) {
  G4String norm = mode;
  std::transform(norm.begin(), norm.end(), norm.begin(),
                 [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  if (norm == "sipm" || norm == "directsipm" || norm == "direct") {
    fSensorMode = "DirectSiPM";
  } else if (norm == "concentratorinput" || norm == "concentrator" || norm == "aperture") {
    fSensorMode = "ConcentratorInput";
  } else {
    G4Exception("DetectorConstruction::SetSensorMode", "InvalidSensorMode",
                JustWarning, "Use DirectSiPM or ConcentratorInput.");
  }
}

void DetectorConstruction::SetOutputEventMode(G4String mode) {
  G4String norm = mode;
  std::transform(norm.begin(), norm.end(), norm.begin(),
                 [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  if (norm == "triggered" || norm == "triggeredonly") fOutputEventMode = "triggered";
  else if (norm == "thinall") fOutputEventMode = "thinAll";
  else if (norm == "fullall") fOutputEventMode = "fullAll";
  else G4Exception("DetectorConstruction::SetOutputEventMode", "InvalidOutputMode",
                   JustWarning, "Use triggered, thinAll, or fullAll.");
}

void DetectorConstruction::SetIncidentDetail(G4String mode) {
  G4String norm = mode;
  std::transform(norm.begin(), norm.end(), norm.begin(),
                 [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  if (norm == "auto") fIncidentDetail = "auto";
  else if (norm == "minimal") fIncidentDetail = "minimal";
  else if (norm == "phasespace" || norm == "full") fIncidentDetail = "phaseSpace";
  else G4Exception("DetectorConstruction::SetIncidentDetail", "InvalidIncidentDetail",
                   JustWarning, "Use auto, minimal, or phaseSpace.");
}

void DetectorConstruction::SetPhotonBirthDetail(G4String mode) {
  G4String norm = mode;
  std::transform(norm.begin(), norm.end(), norm.begin(),
                 [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  if (norm == "basic") fPhotonBirthDetail = "basic";
  else if (norm == "full") fPhotonBirthDetail = "full";
  else G4Exception("DetectorConstruction::SetPhotonBirthDetail", "InvalidPhotonBirthDetail",
                   JustWarning, "Use basic or full.");
}

void DetectorConstruction::SetCoreBoundaryMode(G4String mode) {
  G4String norm = mode;
  std::transform(norm.begin(), norm.end(), norm.begin(),
                 [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  if (norm == "off" || norm == "none") fCoreBoundaryMode = "off";
  else if (norm == "firstexit" || norm == "first") fCoreBoundaryMode = "firstExit";
  else G4Exception("DetectorConstruction::SetCoreBoundaryMode", "InvalidCoreBoundaryMode",
                   JustWarning, "Use off or firstExit.");
}

void DetectorConstruction::SetCoreBoundaryDetail(G4String mode) {
  G4String norm = mode;
  std::transform(norm.begin(), norm.end(), norm.begin(),
                 [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  if (norm == "basic" || norm == "minimal") fCoreBoundaryDetail = "basic";
  else if (norm == "phasespace" || norm == "full") fCoreBoundaryDetail = "phaseSpace";
  else G4Exception("DetectorConstruction::SetCoreBoundaryDetail", "InvalidCoreBoundaryDetail",
                   JustWarning, "Use basic or phaseSpace.");
}

void DetectorConstruction::SetCompressionLevel(G4int value) {
  fCompressionLevel = std::max(0, std::min(9, value));
}

void DetectorConstruction::SetUniformGrid(G4int nU, G4int nV,
                                          G4double pitchU, G4double pitchV) {
  nU = std::max(0, nU);
  nV = std::max(0, nV);
  for (auto& p : fPlanes) {
    p.enabled = (nU > 0 && nV > 0);
    p.autoFill = false;
    p.nU = nU;
    p.nV = nV;
    p.pitchU = pitchU;
    p.pitchV = pitchV;
    p.offsetU = 0.0;
    p.offsetV = 0.0;
  }
}

void DetectorConstruction::SetPlaneGrid(G4int face, G4bool enabled,
                                        G4int nU, G4int nV,
                                        G4double pitchU, G4double pitchV,
                                        G4double offsetU, G4double offsetV) {
  if (face < 0 || face >= 6) return;
  auto& p = fPlanes[face];
  p.enabled = enabled && nU > 0 && nV > 0;
  p.autoFill = false;
  p.nU = std::max(0, nU);
  p.nV = std::max(0, nV);
  p.pitchU = pitchU;
  p.pitchV = pitchV;
  p.offsetU = offsetU;
  p.offsetV = offsetV;
}

void DetectorConstruction::SetPlaneAutoFill(G4int face, G4bool enabled,
                                            G4double gapU, G4double gapV,
                                            G4double marginU, G4double marginV) {
  if (face < 0 || face >= 6) return;
  auto& p = fPlanes[face];
  p.enabled = enabled;
  p.autoFill = true;
  // In AutoFill mode the realized counts/pitches are geometry-derived.
  // Keep configured fixed-grid fields neutral; ChannelMap is authoritative.
  p.nU = 0;
  p.nV = 0;
  p.pitchU = 0.0;
  p.pitchV = 0.0;
  p.gapU = std::max(0.0, gapU);
  p.gapV = std::max(0.0, gapV);
  p.marginU = std::max(0.0, marginU);
  p.marginV = std::max(0.0, marginV);
  p.offsetU = 0.0;
  p.offsetV = 0.0;
}

G4int DetectorConstruction::ParseFace(const G4String& faceName) const {
  G4String s = faceName;
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c){ return static_cast<char>(std::toupper(c)); });
  if (s == "+X" || s == "PX" || s == "PLUSX") return kPlusX;
  if (s == "-X" || s == "NX" || s == "MINUSX") return kMinusX;
  if (s == "+Z" || s == "PZ" || s == "PLUSZ") return kPlusZ;
  if (s == "-Z" || s == "NZ" || s == "MINUSZ") return kMinusZ;
  if (s == "+Y" || s == "PY" || s == "PLUSY") return kPlusY;
  if (s == "-Y" || s == "NY" || s == "MINUSY") return kMinusY;
  return -1;
}

G4String DetectorConstruction::FaceName(G4int face) const {
  static const char* names[6] = {"+X", "-X", "+Z", "-Z", "+Y", "-Y"};
  return (face >= 0 && face < 6) ? names[face] : "Unknown";
}

G4int DetectorConstruction::ResolveCoreFace(
    const G4ThreeVector& boundaryPosition) const {
  const G4double hx = fCoreHalf.x();
  const G4double hy = fCoreHalf.y();
  const G4double hz = fCoreHalf.z();

  struct Candidate { G4double residual; G4int face; };
  const std::array<Candidate,6> candidates = {{
      {std::abs(boundaryPosition.x() - hx), kPlusX},
      {std::abs(boundaryPosition.x() + hx), kMinusX},
      {std::abs(boundaryPosition.z() - hz), kPlusZ},
      {std::abs(boundaryPosition.z() + hz), kMinusZ},
      {std::abs(boundaryPosition.y() - hy), kPlusY},
      {std::abs(boundaryPosition.y() + hy), kMinusY}
  }};

  auto best = candidates.front();
  for (const auto& candidate : candidates) {
    if (candidate.residual < best.residual) best = candidate;
  }
  return best.face;
}

std::pair<G4double,G4double> DetectorConstruction::GetFaceLocalUV(
    G4int face, const G4ThreeVector& worldPosition) const {
  // These axes exactly match BuildFaceRotation(): local +z is the outward
  // face normal, while local x/y are the readout U/V axes.
  switch (face) {
    case kPlusX:  return {-worldPosition.z(),  worldPosition.y()};
    case kMinusX: return { worldPosition.z(),  worldPosition.y()};
    case kPlusZ:  return { worldPosition.x(),  worldPosition.y()};
    case kMinusZ: return {-worldPosition.x(),  worldPosition.y()};
    case kPlusY:  return { worldPosition.x(), -worldPosition.z()};
    case kMinusY: return { worldPosition.x(),  worldPosition.z()};
    default:      return {std::numeric_limits<G4double>::quiet_NaN(),
                          std::numeric_limits<G4double>::quiet_NaN()};
  }
}

const ChannelDescriptor* DetectorConstruction::GetChannelDescriptor(G4int id) const {
  if (id < 0 || id >= static_cast<G4int>(fChannels.size())) return nullptr;
  return &fChannels[id];
}

void DetectorConstruction::RegisterVolume(G4LogicalVolume* lv, G4int id,
                                          const G4String& name) {
  if (!lv) return;
  fVolumeIds[lv] = id;
  VolumeDescriptor d;
  d.volumeID = id;
  d.name = name;
  d.materialName = lv->GetMaterial() ? lv->GetMaterial()->GetName() : "Unknown";
  fVolumeDescriptors.push_back(d);
}

G4int DetectorConstruction::GetVolumeID(const G4LogicalVolume* lv) const {
  if (!lv) return -1;
  auto it = fVolumeIds.find(lv);
  return it == fVolumeIds.end() ? -1 : it->second;
}

G4ThreeVector DetectorConstruction::FaceNormal(G4int face) const {
  switch (face) {
    case kPlusX:  return G4ThreeVector( 1, 0, 0);
    case kMinusX: return G4ThreeVector(-1, 0, 0);
    case kPlusZ:  return G4ThreeVector( 0, 0, 1);
    case kMinusZ: return G4ThreeVector( 0, 0,-1);
    case kPlusY:  return G4ThreeVector( 0, 1, 0);
    case kMinusY: return G4ThreeVector( 0,-1, 0);
    default:      return G4ThreeVector();
  }
}

G4RotationMatrix* DetectorConstruction::BuildFaceRotation(G4int face) {
  auto* rot = new G4RotationMatrix();
  if (face == kPlusX)       rot->rotateY( 90*deg);
  else if (face == kMinusX) rot->rotateY(-90*deg);
  else if (face == kPlusZ)  {}
  else if (face == kMinusZ) rot->rotateY(180*deg);
  else if (face == kPlusY)  rot->rotateX(-90*deg);
  else if (face == kMinusY) rot->rotateX( 90*deg);
  fRotations.push_back(rot);
  return rot;
}

std::vector<std::pair<G4double,G4double>>
DetectorConstruction::BuildPlaneCoordinates(G4int face,
                                            G4double hx, G4double hy, G4double hz) const {
  const auto& p = fPlanes[face];
  std::vector<std::pair<G4double,G4double>> coords;
  if (!p.enabled) return coords;

  G4double halfU = 0.0;
  G4double halfV = 0.0;
  if (face == kPlusX || face == kMinusX) {
    halfU = hz;
    halfV = hy;
  } else if (face == kPlusZ || face == kMinusZ) {
    halfU = hx;
    halfV = hy;
  } else {
    halfU = hx;
    halfV = hz;
  }

  const G4double holeU = fTyvekThickness > 0.0 ? fTyvekHoleSize.x() : 0.0;
  const G4double holeV = fTyvekThickness > 0.0 ? fTyvekHoleSize.y() : 0.0;
  const G4double footprintU = std::max({fReadoutSize.x(), fCouplingSize.x(), holeU});
  const G4double footprintV = std::max({fReadoutSize.y(), fCouplingSize.y(), holeV});

  G4int nU = p.nU;
  G4int nV = p.nV;
  G4double pitchU = p.pitchU;
  G4double pitchV = p.pitchV;

  if (p.autoFill) {
    const G4double usableU = std::max(0.0, 2.0*halfU - 2.0*p.marginU);
    const G4double usableV = std::max(0.0, 2.0*halfV - 2.0*p.marginV);
    pitchU = footprintU + p.gapU;
    pitchV = footprintV + p.gapV;
    nU = pitchU > 0.0 ? static_cast<G4int>(std::floor((usableU + p.gapU) / pitchU)) : 0;
    nV = pitchV > 0.0 ? static_cast<G4int>(std::floor((usableV + p.gapV) / pitchV)) : 0;
  }

  if (nU <= 0 || nV <= 0) return coords;

  for (G4int iu = 0; iu < nU; ++iu) {
    const G4double u = (iu - 0.5*(nU-1))*pitchU + p.offsetU;
    for (G4int iv = 0; iv < nV; ++iv) {
      const G4double v = (iv - 0.5*(nV-1))*pitchV + p.offsetV;
      if (std::abs(u) + 0.5*footprintU > halfU + 1.0e-9*mm ||
          std::abs(v) + 0.5*footprintV > halfV + 1.0e-9*mm) {
        continue;
      }
      coords.emplace_back(u, v);
    }
  }
  return coords;
}


void DetectorConstruction::ValidateConfiguration() const {
  auto fatal = [](const char* code, const std::string& message) {
    G4Exception("DetectorConstruction::ValidateConfiguration", code,
                FatalException, message.c_str());
  };
  auto positive3 = [&](const G4ThreeVector& v, const char* name) {
    if (v.x() <= 0.0 || v.y() <= 0.0 || v.z() <= 0.0) {
      fatal("InvalidDimension", std::string(name) + " dimensions must all be > 0.");
    }
  };

  if (fCoreHalf.x() <= 0.0 || fCoreHalf.y() <= 0.0 || fCoreHalf.z() <= 0.0)
    fatal("InvalidCoreSize", "Core dimensions must all be > 0.");
  positive3(fScintSize, "Scintillator");
  positive3(fReadoutSize, "Readout element");
  positive3(fCouplingSize, "Coupling element");
  if (fOuterShellThickness < 0.0 || fTyvekThickness < 0.0 ||
      fScintGap < 0.0 || fSourceDist < 0.0)
    fatal("NegativeGeometryParameter", "Thicknesses, scintillator gap and source distance must be >= 0.");
  if (fUseOuterShell && fOuterShellThickness <= 0.0)
    fatal("ZeroOuterShell", "Outer shell is enabled but outer-shell thickness is not > 0.");
  if (fTyvekThickness > 0.0) {
    if (fTyvekHoleSize.x() <= 0.0 || fTyvekHoleSize.y() <= 0.0)
      fatal("InvalidTyvekHole", "Tyvek hole X/Y dimensions must be > 0 when Tyvek is enabled.");
    if (fTyvekHoleSize.x() + 1e-12*mm < std::max(fReadoutSize.x(), fCouplingSize.x()) ||
        fTyvekHoleSize.y() + 1e-12*mm < std::max(fReadoutSize.y(), fCouplingSize.y()))
      fatal("TyvekHoleTooSmall", "Tyvek hole must be at least as large as the readout/coupling footprint in X and Y.");
    if (fTyvekConfig == "Inside" &&
        (fTyvekThickness >= fCoreHalf.x() || fTyvekThickness >= fCoreHalf.y() || fTyvekThickness >= fCoreHalf.z()))
      fatal("TyvekTooThick", "Inside Tyvek thickness must be smaller than every core half-dimension.");
  }

  const G4double hx = fCoreHalf.x(), hy = fCoreHalf.y(), hz = fCoreHalf.z();
  const G4double holeU = fTyvekThickness > 0.0 ? fTyvekHoleSize.x() : 0.0;
  const G4double holeV = fTyvekThickness > 0.0 ? fTyvekHoleSize.y() : 0.0;
  const G4double footprintU = std::max({fReadoutSize.x(), fCouplingSize.x(), holeU});
  const G4double footprintV = std::max({fReadoutSize.y(), fCouplingSize.y(), holeV});

  for (G4int face=0; face<6; ++face) {
    const auto& cfg = fPlanes[face];
    if (!cfg.enabled) continue;
    G4double halfU=0.0, halfV=0.0;
    if (face==kPlusX || face==kMinusX) { halfU=hz; halfV=hy; }
    else if (face==kPlusZ || face==kMinusZ) { halfU=hx; halfV=hy; }
    else { halfU=hx; halfV=hz; }

    if (cfg.autoFill) {
      const G4double usableU = 2.0*halfU - 2.0*cfg.marginU;
      const G4double usableV = 2.0*halfV - 2.0*cfg.marginV;
      if (usableU + 1e-12*mm < footprintU || usableV + 1e-12*mm < footprintV) {
        fatal("AutoFillNoElement", std::string("AutoFill on face ") + FaceName(face) +
              " cannot fit even one readout element with the requested margins.");
      }
      continue;
    }

    if (cfg.nU <= 0 || cfg.nV <= 0)
      fatal("InvalidGridCount", std::string("Enabled fixed grid on face ") + FaceName(face) + " must have nU,nV > 0.");
    if ((cfg.nU > 1 && cfg.pitchU <= 0.0) || (cfg.nV > 1 && cfg.pitchV <= 0.0))
      fatal("InvalidGridPitch", std::string("Enabled fixed grid on face ") + FaceName(face) + " requires positive pitch when count > 1.");
    if ((cfg.nU > 1 && cfg.pitchU + 1e-12*mm < footprintU) ||
        (cfg.nV > 1 && cfg.pitchV + 1e-12*mm < footprintV))
      fatal("OverlappingGrid", std::string("Fixed-grid pitch on face ") + FaceName(face) +
            " is smaller than the readout/coupling/hole footprint and would overlap elements.");

    const G4double extentU = 0.5*(cfg.nU-1)*cfg.pitchU + 0.5*footprintU + std::abs(cfg.offsetU);
    const G4double extentV = 0.5*(cfg.nV-1)*cfg.pitchV + 0.5*footprintV + std::abs(cfg.offsetV);
    if (extentU > halfU + 1e-9*mm || extentV > halfV + 1e-9*mm)
      fatal("GridOutsideFace", std::string("Fixed grid on face ") + FaceName(face) +
            " extends beyond the radiator face. Reduce count/pitch/offset or use AutoFill.");
  }
}

G4double DetectorConstruction::ComputeScintTopY() const {
  const G4double hy = fCoreHalf.y();
  const G4double tShell = fUseOuterShell ? fOuterShellThickness : 0.0;
  const G4double sensorTop = hy + tShell + fCouplingSize.z() + fReadoutSize.z();
  const G4double tyvekTop = (fTyvekConfig == "Inside")
      ? hy : hy + tShell + fTyvekThickness;
  return std::max(sensorTop, tyvekTop) + fScintGap + 0.5*fScintSize.y();
}

G4double DetectorConstruction::GetGunPosY() const {
  return ComputeScintTopY() + 0.5*fScintSize.y() + fSourceDist;
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  ValidateConfiguration();

  for (auto* rot : fRotations) delete rot;
  fRotations.clear();
  fChannels.clear();
  fVolumeIds.clear();
  fVolumeDescriptors.clear();

  fLogicWorld = nullptr;
  fLogicCore = nullptr;
  fLogicOuterShell = nullptr;
  fLogicScintTop = nullptr;
  fLogicScintBottom = nullptr;
  fLogicTyvek = nullptr;
  fLogicCoupling = nullptr;
  fLogicSensor = nullptr;

  const G4double hx = fCoreHalf.x();
  const G4double hy = fCoreHalf.y();
  const G4double hz = fCoreHalf.z();
  const G4double tShell = fUseOuterShell ? fOuterShellThickness : 0.0;
  const G4double tTyv = std::max(0.0, fTyvekThickness);
  const G4double tGel = fCouplingSize.z();
  const G4double tSensor = fReadoutSize.z();

  const G4double maxReach = std::max({hx, hy, hz}) + tShell + tTyv + tGel + tSensor;
  const G4double worldSize = 3.0*maxReach + fScintGap + fScintSize.y() + fSourceDist + 100.0*mm;

  auto* solidWorld = new G4Box("World", worldSize, worldSize, worldSize);
  fLogicWorld = new G4LogicalVolume(solidWorld, fAir, "World");
  RegisterVolume(fLogicWorld, 2, "World");
  fLogicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), fLogicWorld,
                                      "World", nullptr, false, 0, fCheckOverlaps);

  const G4double scintY = ComputeScintTopY();
  auto* solidScint = new G4Box("ScintBox", 0.5*fScintSize.x(),
                               0.5*fScintSize.y(), 0.5*fScintSize.z());
  fLogicScintTop = new G4LogicalVolume(solidScint, fScintMat, "ScintTopLogic");
  fLogicScintBottom = new G4LogicalVolume(solidScint, fScintMat, "ScintBottomLogic");
  RegisterVolume(fLogicScintTop, 5, "ScintTop");
  RegisterVolume(fLogicScintBottom, 6, "ScintBottom");
  new G4PVPlacement(nullptr, G4ThreeVector(0, scintY, 0), fLogicScintTop,
                    "ScintTop", fLogicWorld, false, 0, fCheckOverlaps);
  new G4PVPlacement(nullptr, G4ThreeVector(0,-scintY, 0), fLogicScintBottom,
                    "ScintBottom", fLogicWorld, false, 1, fCheckOverlaps);

  auto* solidCore = new G4Box("CoreBox", hx, hy, hz);
  fLogicCore = new G4LogicalVolume(solidCore, fCoreMat, "CoreBox");
  RegisterVolume(fLogicCore, 0, "Core");
  auto* coreVis = new G4VisAttributes(G4Colour(0.0,0.4,1.0,0.2));
  coreVis->SetForceSolid(true);
  fLogicCore->SetVisAttributes(coreVis);
  auto* physCore = new G4PVPlacement(nullptr, G4ThreeVector(), fLogicCore,
                                     "CorePhysical", fLogicWorld, false, 0, fCheckOverlaps);

  G4VPhysicalVolume* physOuterShell = nullptr;
  if (fUseOuterShell) {
    auto* solidOuterShell = new G4SubtractionSolid(
        // R1 deliberately keeps Geant4 object names unchanged so any legacy
        // name-based debug/analysis path sees the exact pre-refactor geometry.
        "AcrylicShell",
        new G4Box("AcrOut", hx+tShell, hy+tShell, hz+tShell),
        new G4Box("AcrIn", hx, hy, hz));
    fLogicOuterShell = new G4LogicalVolume(solidOuterShell, fOuterShellMat, "AcrylicBox");
    // Legacy VolumeMap label is also retained until R2 migrates ROOT metadata.
    RegisterVolume(fLogicOuterShell, 1, "Acrylic");
    physOuterShell = new G4PVPlacement(nullptr, G4ThreeVector(), fLogicOuterShell,
                                       "AcrylicPhysical", fLogicWorld, false, 0, fCheckOverlaps);
  }

  G4VSolid* tyvekBase = nullptr;
  G4Box* punchHole = nullptr;
  G4OpticalSurface* tyvekSurface = nullptr;
  if (tTyv > 0.0) {
    G4ThreeVector tyvIn, tyvOut;
    if (fTyvekConfig == "Inside") {
      tyvIn  = G4ThreeVector(hx-tTyv, hy-tTyv, hz-tTyv);
      tyvOut = G4ThreeVector(hx-1.0*um, hy-1.0*um, hz-1.0*um);
    } else {
      tyvIn  = G4ThreeVector(hx+tShell, hy+tShell, hz+tShell);
      tyvOut = G4ThreeVector(hx+tShell+tTyv, hy+tShell+tTyv, hz+tShell+tTyv);
    }

    tyvekBase = new G4SubtractionSolid(
        "TyvekShellBase",
        new G4Box("TyvekOuter", tyvOut.x(), tyvOut.y(), tyvOut.z()),
        new G4Box("TyvekInner", tyvIn.x(), tyvIn.y(), tyvIn.z()));
    punchHole = new G4Box("ReadoutHole", 0.5*fTyvekHoleSize.x(),
                          0.5*fTyvekHoleSize.y(), tTyv+5.0*mm);

    tyvekSurface = new G4OpticalSurface("TyvekSurface");
    tyvekSurface->SetType(dielectric_metal);
    tyvekSurface->SetModel(unified);
    tyvekSurface->SetFinish(ground);
    std::vector<G4double> e, r;
    ReadOpticalData(DataPath(fTyvekDataset), e, r, true);
    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("REFLECTIVITY", e, r);
    tyvekSurface->SetMaterialPropertiesTable(mpt);
  }

  auto* sensorSurface = new G4OpticalSurface("SensorSurface");
  sensorSurface->SetType(dielectric_metal);
  sensorSurface->SetFinish(polished);
  sensorSurface->SetModel(glisur);
  std::vector<G4double> ePde, pde;
  ReadOpticalData(DataPath(fSensorPdeDataset), ePde, pde, true);
  std::vector<G4double> efficiency = pde;
  if (fSensorMode == "ConcentratorInput") {
    std::fill(efficiency.begin(), efficiency.end(), 0.0);
  }
  std::vector<G4double> reflectivity(ePde.size(), 0.0);
  auto* sensorMpt = new G4MaterialPropertiesTable();
  sensorMpt->AddProperty("EFFICIENCY", ePde, efficiency);
  sensorMpt->AddProperty("REFLECTIVITY", ePde, reflectivity);
  sensorSurface->SetMaterialPropertiesTable(sensorMpt);

  auto* gelInterface = new G4OpticalSurface("RadiatorGelSurface");
  gelInterface->SetType(dielectric_dielectric);
  gelInterface->SetModel(unified);
  gelInterface->SetFinish(polished);

  auto* solidCoupling = new G4Box("CouplingBox", 0.5*fCouplingSize.x(),
                                  0.5*fCouplingSize.y(), 0.5*tGel);
  fLogicCoupling = new G4LogicalVolume(solidCoupling, fCouplingMat, "CouplingBox");
  RegisterVolume(fLogicCoupling, 3, "Coupling");

  auto* solidSensor = new G4Box("SensorBox", 0.5*fReadoutSize.x(),
                                0.5*fReadoutSize.y(), 0.5*tSensor);
  auto* readoutBulk = (fSensorMode == "ConcentratorInput") ? fAir : fSensorMat;
  fLogicSensor = new G4LogicalVolume(solidSensor, readoutBulk, "SensorBox");
  RegisterVolume(fLogicSensor, 4, "Readout");

  auto* holes = (tTyv > 0.0) ? new G4MultiUnion("TyvekReadoutHoles") : nullptr;
  G4bool hasHole = false;
  G4int channelID = 0;

  for (G4int face = 0; face < 6; ++face) {
    const auto coords = BuildPlaneCoordinates(face, hx, hy, hz);

    G4double hAxis = 0.0;
    if (face == kPlusX || face == kMinusX) hAxis = hx;
    else if (face == kPlusZ || face == kMinusZ) hAxis = hz;
    else hAxis = hy;

    const G4double surfaceDistance = hAxis + tShell;
    const G4double tyvekCenter = (fTyvekConfig == "Inside")
        ? hAxis - 0.5*tTyv : hAxis + tShell + 0.5*tTyv;

    const auto& cfg = fPlanes[face];
    std::vector<G4double> uniqueU;
    std::vector<G4double> uniqueV;
    for (const auto& uv : coords) {
      if (std::find(uniqueU.begin(), uniqueU.end(), uv.first) == uniqueU.end()) uniqueU.push_back(uv.first);
      if (std::find(uniqueV.begin(), uniqueV.end(), uv.second) == uniqueV.end()) uniqueV.push_back(uv.second);
    }
    std::sort(uniqueU.begin(), uniqueU.end());
    std::sort(uniqueV.begin(), uniqueV.end());

    for (const auto& uv : coords) {
      const G4double u = uv.first;
      const G4double v = uv.second;

      auto* rot = BuildFaceRotation(face);
      const G4ThreeVector posHole = (*rot) * G4ThreeVector(u, v, tyvekCenter);
      const G4ThreeVector posGel  = (*rot) * G4ThreeVector(u, v, surfaceDistance + 0.5*tGel);
      const G4ThreeVector posSen  = (*rot) * G4ThreeVector(u, v, surfaceDistance + tGel + 0.5*tSensor);

      if (holes) {
        holes->AddNode(*punchHole, G4Transform3D(*rot, posHole));
        hasHole = true;
      }

      // IMPORTANT: BuildFaceRotation() is a direct rotation of the solid
      // (local -> world).  Use the G4Transform3D placement overload, whose
      // rotation has that convention.  The pointer+translation overload uses
      // the rotation of the daughter reference frame and would invert the
      // side-face local axes (+/-X,+/-Y), making local +z point inward.
      auto* physC = new G4PVPlacement(G4Transform3D(*rot, posGel), fLogicCoupling,
                                      "CouplingPhys", fLogicWorld, false,
                                      channelID, fCheckOverlaps);
      auto* physS = new G4PVPlacement(G4Transform3D(*rot, posSen), fLogicSensor,
                                      "SensorPhys", fLogicWorld, false,
                                      channelID, fCheckOverlaps);

      if (fUseOuterShell) {
        // Surface name retained for strict pre/post-R1 regression compatibility.
        new G4LogicalBorderSurface("AcrylicToGel_" + std::to_string(channelID),
                                   physOuterShell, physC, gelInterface);
      } else {
        new G4LogicalBorderSurface("CoreToGel_" + std::to_string(channelID),
                                   physCore, physC, gelInterface);
      }
      new G4LogicalBorderSurface("GelToSensor_" + std::to_string(channelID),
                                 physC, physS, sensorSurface);

      ChannelDescriptor desc;
      desc.channelID = channelID;
      desc.faceID = face;
      desc.localU = u;
      desc.localV = v;
      desc.center = posSen;
      desc.normal = FaceNormal(face);
      desc.entranceCenter = posSen - 0.5*tSensor*desc.normal;
      desc.sizeU = fReadoutSize.x();
      desc.sizeV = fReadoutSize.y();

      desc.gridU = static_cast<G4int>(std::distance(uniqueU.begin(),
          std::lower_bound(uniqueU.begin(), uniqueU.end(), u)));
      desc.gridV = static_cast<G4int>(std::distance(uniqueV.begin(),
          std::lower_bound(uniqueV.begin(), uniqueV.end(), v)));

      fChannels.push_back(desc);
      ++channelID;
    }
  }

  if (tyvekBase) {
    G4VSolid* tyvekSolid = tyvekBase;
    if (hasHole) {
      holes->Voxelize();
      tyvekSolid = new G4SubtractionSolid("TyvekShell", tyvekBase, holes);
    }
    fLogicTyvek = new G4LogicalVolume(tyvekSolid, fDummyTyvekMat, "TyvekLogic");
    RegisterVolume(fLogicTyvek, 7, "Tyvek");
    G4LogicalVolume* mother = (fTyvekConfig == "Inside") ? fLogicCore : fLogicWorld;
    new G4PVPlacement(nullptr, G4ThreeVector(), fLogicTyvek,
                      "TyvekPhysical", mother, false, 0, fCheckOverlaps);
    new G4LogicalSkinSurface("TyvekSkin", fLogicTyvek, tyvekSurface);
  }

  ++fGeometryVersion;
  return physWorld;
}

void DetectorConstruction::ConstructSDandField() {}

void DetectorConstruction::SetWaterRayleighFile(G4String filename) {
  // Validation/expert override for the Water_Default Rayleigh mean-free-path
  // table. Keep the profile identity, RINDEX and ABSLENGTH unchanged so
  // RunConfig provenance remains authoritative for the active test.
  auto* profile = FindMaterialProfile("Water_Default");
  if (!profile) {
    G4Exception("DetectorConstruction::SetWaterRayleighFile",
                "MissingWaterProfile", FatalException,
                "Water_Default profile is not registered.");
    return;
  }

  // Restrict this route to a data-set basename. This keeps RunConfig dataset
  // naming/fingerprinting consistent with the normal OPTICAL_DATA_DIR path.
  if (filename.find('/') != std::string::npos ||
      filename.find('\\') != std::string::npos) {
    G4Exception("DetectorConstruction::SetWaterRayleighFile",
                "WaterRayleighFileMustBeBasename", JustWarning,
                "Use a filename from the project data/ directory, not a path.");
    return;
  }

  const G4String resolved = DataPath(filename);
  std::vector<G4double> re, rl;
  ReadOpticalData(resolved, re, rl);
  for (auto& x : rl) x *= profile->rayleighLengthUnit;

  std::vector<G4double> ne, nv;
  ReadOpticalData(DataPath(profile->rindexDataset), ne, nv);
  std::vector<G4double> ae, al;
  ReadOpticalData(DataPath(profile->absDataset), ae, al);
  for (auto& x : al) x *= profile->absLengthUnit;

  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", ne, nv);
  mpt->AddProperty("ABSLENGTH", ae, al);
  mpt->AddProperty("RAYLEIGH", re, rl);

  profile->rayleighDataset = filename;
  profile->rayleighStatus = "RuntimeValidationOverride_via_waterRayleighFile";
  profile->material->SetMaterialPropertiesTable(mpt);

  // If Water_Default is the active Core, its G4Material pointer is unchanged
  // and therefore already sees the replaced MPT.
  if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

void DetectorConstruction::SetPmmaAbsFile(G4String filename) {
  // Validation/expert override for PMMA_SNO_A ABSLENGTH.  This is deliberately
  // analogous to /detector/waterRayleighFile: profile identity, RINDEX and
  // geometry are preserved; the active absorption dataset/status/policy change.
  auto* profile = FindMaterialProfile("PMMA_SNO_A");
  if (!profile) {
    G4Exception("DetectorConstruction::SetPmmaAbsFile",
                "MissingPmmaProfile", FatalException,
                "PMMA_SNO_A profile is not registered.");
    return;
  }

  // Keep this route inside the versioned project data/ directory so metadata
  // dataset names and the optical-data fingerprint remain reproducible.
  if (filename.find('/') != std::string::npos ||
      filename.find('\\') != std::string::npos) {
    G4Exception("DetectorConstruction::SetPmmaAbsFile",
                "PmmaAbsFileMustBeBasename", JustWarning,
                "Use a filename from the project data/ directory, not a path.");
    return;
  }

  // Validate/read the requested table before mutating the active profile.
  std::vector<G4double> testE, testL;
  ReadOpticalData(DataPath(filename), testE, testL);

  profile->absDataset = filename;
  profile->absStatus = "RuntimeValidationOverride_via_pmmaAbsFile";

  // The override table may contain explicitly modelled transport support
  // outside the profile's measured ABS domain.  Measured validity remains
  // authoritative in MaterialProfile::absMeasuredMinNm/MaxNm.
  profile->outOfRangePolicy =
      "ExplicitAbsTableOverride;"
      "MeasuredValidityFromProfile;"
      "NumericalSupportFromDataset;"
      "OutsideSupport_Geant4EndpointClamp";

  ApplyMaterialProfileOptics(*profile);

  // PMMA_SNO_A may be the Core, OuterShell, or both. All such roles share the
  // same curated G4Material pointer, so replacing its MPT updates the active
  // material without changing geometry.
  if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

void DetectorConstruction::SetCouplingFile(G4String filepath) {
  // Legacy expert override retained for backward compatibility.  The material
  // identity remains PDMS_Generic; only its run-local RINDEX dataset changes.
  // This preserves pre-R2D behavior while keeping provenance/fingerprint
  // getters synchronized with the active profile.
  auto* profile = FindMaterialProfile(fCouplingProfileName);
  if (!profile) {
    G4Exception("DetectorConstruction::SetCouplingFile",
                "MissingCouplingProfile", FatalException,
                "Active coupling profile is not registered.");
    return;
  }

  std::vector<G4double> e, n;
  ReadOpticalData(filepath, e, n);
  std::vector<G4double> ae, al;
  ReadOpticalData(DataPath(profile->absDataset), ae, al);
  for (auto& x : al) x *= profile->absLengthUnit;

  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", e, n);
  mpt->AddProperty("ABSLENGTH", ae, al);
  profile->rindexDataset = filepath;
  profile->rindexStatus = "RuntimeOverride_via_legacy_couplingFile";
  profile->material->SetMaterialPropertiesTable(mpt);
  fCouplingMat = profile->material;

  if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}
