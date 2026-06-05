#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"

#include <fstream>
#include <sstream>
#include <algorithm>

//--------------------------------------------------------------------
DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fCheckOverlaps(true),
    fWater(nullptr), fAcrylic(nullptr), fAir(nullptr), 
    fCouplingMat(nullptr), fLogicCoupling(nullptr)
{
  fMessenger = new DetectorMessenger(this);
  DefineMaterials();
}

DetectorConstruction::~DetectorConstruction() {
  delete fMessenger;
}

//--------------------------------------------------------------------
void DetectorConstruction::ReadOpticalData(const G4String& filename, 
                                           std::vector<G4double>& photonEnergies, 
                                           std::vector<G4double>& properties) 
{
  std::ifstream file(filename);
  if (!file.is_open()) {
    G4String msg = "Khong the mo file du lieu quang hoc: " + filename;
    G4Exception("DetectorConstruction::ReadOpticalData", "FileOpenError", FatalException, msg);
  }

  struct DataPoint { G4double energy; G4double value; };
  std::vector<DataPoint> data;

  G4String line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;

    std::stringstream ss(line);
    G4double wavelength_nm, value;
    
    std::replace(line.begin(), line.end(), ',', ' ');
    std::stringstream ss_clean(line);

    if (ss_clean >> wavelength_nm >> value) {
      G4double energy_eV = (1239.84193 / wavelength_nm) * eV;
      data.push_back({energy_eV, value});
    }
  }
  file.close();

  std::sort(data.begin(), data.end(), [](const DataPoint& a, const DataPoint& b) {
    return a.energy < b.energy;
  });

  photonEnergies.clear();
  properties.clear();
  for (const auto& p : data) {
    photonEnergies.push_back(p.energy);
    properties.push_back(p.value);
  }
}

//--------------------------------------------------------------------
void DetectorConstruction::DefineMaterials()
{
  G4NistManager* nist = G4NistManager::Instance();

  // ── 1. KHÔNG KHÍ (AIR) ───────────────────────────────────────────
  fAir = nist->FindOrBuildMaterial("G4_AIR");
  {
    std::vector<G4double> airEnergies = { 1.0*eV, 9.0*eV }; 
    std::vector<G4double> airRIndex = { 1.0000, 1.0000 }; 
    std::vector<G4double> airAbsLen = { 1.0e10*m, 1.0e10*m }; 

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX",    airEnergies, airRIndex);
    mpt->AddProperty("ABSLENGTH", airEnergies, airAbsLen);
    fAir->SetMaterialPropertiesTable(mpt);
  }

  // ── 2. NƯỚC ──────────────────────────────────
  fWater = nist->FindOrBuildMaterial("G4_WATER");
  {
    std::vector<G4double> energy_Rindex, rindex;
    std::vector<G4double> energy_Abs, absLength;

    ReadOpticalData("../data/Water_RINDEX.txt", energy_Rindex, rindex);
    ReadOpticalData("../data/Water_ABS.txt", energy_Abs, absLength);

    for(auto& val : absLength) { val *= m; }

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", energy_Rindex, rindex);
    mpt->AddProperty("ABSLENGTH", energy_Abs, absLength);
    mpt->AddConstProperty("RESOLUTIONSCALE", 1.0); 
    fWater->SetMaterialPropertiesTable(mpt);
  }

  // ── 3. ACRYLIC ───────────────────────────────
  fAcrylic = new G4Material("Acrylic", 1.19*g/cm3, 3);
  fAcrylic->AddElement(nist->FindOrBuildElement("C"), 5);
  fAcrylic->AddElement(nist->FindOrBuildElement("H"), 8);
  fAcrylic->AddElement(nist->FindOrBuildElement("O"), 2);
  {
    std::vector<G4double> energy_Rindex, rindex;
    std::vector<G4double> energy_Abs, absLength;

    ReadOpticalData("../data/Acrylic_RINDEX.txt", energy_Rindex, rindex);
    ReadOpticalData("../data/Acrylic_ABS.txt", energy_Abs, absLength);

    for(auto& val : absLength) { val *= mm; }

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", energy_Rindex, rindex);
    mpt->AddProperty("ABSLENGTH", energy_Abs, absLength);
    fAcrylic->SetMaterialPropertiesTable(mpt);
  }

  // ── 4. LỚP VỎ COUPLING LINH HOẠT (MẶC ĐỊNH LÀ KHÔNG KHÍ) ────────
  fCouplingMat = new G4Material("CouplingMat", 1.2*g/cm3, 1);
  fCouplingMat->AddElement(nist->FindOrBuildElement("C"), 1); // Giả lập chất hữu cơ
  {
    // Mặc định ban đầu là Không khí (n=1.0)
    std::vector<G4double> defE = { 1.0*eV, 9.0*eV }; 
    std::vector<G4double> defR = { 1.0, 1.0 }; 
    std::vector<G4double> defA = { 1.0e10*m, 1.0e10*m }; 
    auto* mptCoup = new G4MaterialPropertiesTable();
    mptCoup->AddProperty("RINDEX", defE, defR);
    mptCoup->AddProperty("ABSLENGTH", defE, defA);
    fCouplingMat->SetMaterialPropertiesTable(mptCoup);
  }
}

//--------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // ── KÍCH THƯỚC HÌNH HỌC ──────────────────────────────────────────
  G4double waterHalf    = 5.0*cm;              
  G4double acrylicThick = 2.0*mm;              
  G4double acrylicHalf  = waterHalf + acrylicThick; 
  
  G4double couplingThick= 2.0*mm;              // Lớp vỏ quang học mới
  G4double couplingHalf = acrylicHalf + couplingThick;
  
  G4double worldHalf    = 50.*cm;              

  // 1. WORLD (KHÔNG KHÍ)
  auto solidWorld = new G4Box("World", worldHalf, worldHalf, worldHalf);
  auto logicWorld = new G4LogicalVolume(solidWorld, fAir, "World");
  auto physWorld  = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, fCheckOverlaps);
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

  // 2. LỚP COUPLING BOX (Bọc trong World)
  auto solidCoupling = new G4Box("CouplingBox", couplingHalf, couplingHalf, couplingHalf);
  fLogicCoupling = new G4LogicalVolume(solidCoupling, fCouplingMat, "CouplingBox");
  auto* coupVis = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.3)); // Trực quan màu vàng
  coupVis->SetForceWireframe(true);
  fLogicCoupling->SetVisAttributes(coupVis);
  auto physCoupling = new G4PVPlacement(nullptr, G4ThreeVector(), fLogicCoupling, "CouplingBox", logicWorld, false, 0, fCheckOverlaps);

  // 3. BỂ ACRYLIC (Bọc trong CouplingBox)
  auto solidAcrylic = new G4Box("AcrylicBox", acrylicHalf, acrylicHalf, acrylicHalf);
  auto logicAcrylic = new G4LogicalVolume(solidAcrylic, fAcrylic, "AcrylicBox");
  auto* acrVis = new G4VisAttributes(G4Colour(0.5, 0.8, 1.0, 0.3));
  acrVis->SetForceWireframe(true);
  logicAcrylic->SetVisAttributes(acrVis);
  auto physAcrylic = new G4PVPlacement(nullptr, G4ThreeVector(), logicAcrylic, "AcrylicBox", fLogicCoupling, false, 0, fCheckOverlaps);

  // 4. LÕI NƯỚC (Bọc trong Acrylic)
  auto solidWater = new G4Box("WaterBox", waterHalf, waterHalf, waterHalf);
  auto logicWater = new G4LogicalVolume(solidWater, fWater, "WaterBox");
  auto* watVis = new G4VisAttributes(G4Colour(0.0, 0.4, 1.0, 0.5));
  watVis->SetForceSolid(true);
  logicWater->SetVisAttributes(watVis);
  auto physWater = new G4PVPlacement(nullptr, G4ThreeVector(), logicWater, "WaterBox", logicAcrylic, false, 0, fCheckOverlaps);

  // ── ÉP BUỘC GEANT4 NHẬN DIỆN BỀ MẶT ĐIỆN MÔI ─────────────
  auto* polishedSurface = new G4OpticalSurface("PolishedSurface");
  polishedSurface->SetType(dielectric_dielectric);
  polishedSurface->SetModel(unified);
  polishedSurface->SetFinish(polished); 

  new G4LogicalBorderSurface("W2A", physWater, physAcrylic, polishedSurface);
  new G4LogicalBorderSurface("A2W", physAcrylic, physWater, polishedSurface);
  // Sửa lại mặt biên ngoài thành Acrylic <-> Coupling
  new G4LogicalBorderSurface("A2Coup", physAcrylic, physCoupling, polishedSurface);
  new G4LogicalBorderSurface("Coup2A", physCoupling, physAcrylic, polishedSurface);
  
  return physWorld;
}

//--------------------------------------------------------------------
void DetectorConstruction::ConstructSDandField() {}

//--------------------------------------------------------------------
void DetectorConstruction::SetCouplingFile(G4String filepath) {
    std::vector<G4double> energies, rindex;
    
    // Đọc file txt để lấy dải RINDEX
    ReadOpticalData(filepath, energies, rindex); 

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", energies, rindex);

    std::vector<G4double> absEnergies = { energies.front(), energies.back() };
    std::vector<G4double> absLengths  = { 1.0e10*m, 1.0e10*m };
    mpt->AddProperty("ABSLENGTH", absEnergies, absLengths);

    if (fCouplingMat) {
        fCouplingMat->SetMaterialPropertiesTable(mpt);
        

        if (G4RunManager::GetRunManager()->GetUserDetectorConstruction() != nullptr) {
            G4RunManager::GetRunManager()->PhysicsHasBeenModified(); 
        }
        
        G4cout << "\n=================================================================" << G4endl;
        G4cout << "====> [FLEXIBLE OPTICS] DA NAP THANH CONG RINDEX TU: " << filepath << G4endl;
        G4cout << "=================================================================\n" << G4endl;
    } else {
        G4cout << "\n[ERROR] fCouplingMat is NULL! Khong the nap quang hoc!" << G4endl;
    }
}