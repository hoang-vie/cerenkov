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
    fCouplingMat(nullptr), fSensorMat(nullptr),
    fLogicCoupling(nullptr), fLogicSensor(nullptr)
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
    std::istringstream iss(line);
    G4double wavelength, val;
    if (iss >> wavelength >> val) {
      G4double energy = (1239.84193 * eV * nm) / (wavelength * nm);
      data.push_back({energy, val});
    }
  }

  std::sort(data.begin(), data.end(), [](const DataPoint& a, const DataPoint& b) {
    return a.energy < b.energy;
  });

  photonEnergies.clear();
  properties.clear();
  for (const auto& pt : data) {
    photonEnergies.push_back(pt.energy);
    properties.push_back(pt.value);
  }
}

//--------------------------------------------------------------------
void DetectorConstruction::DefineMaterials()
{
  auto nist = G4NistManager::Instance();
  fAir = nist->FindOrBuildMaterial("G4_AIR");
  fWater = nist->FindOrBuildMaterial("G4_WATER");
  fAcrylic = nist->FindOrBuildMaterial("G4_PLEXIGLASS"); 
  fCouplingMat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE"); // Dummy
  fSensorMat = nist->FindOrBuildMaterial("G4_Si");

  auto* mptCouplingDef = new G4MaterialPropertiesTable();
  std::vector<G4double> cWvl, cRindex;
  ReadOpticalData("../data/Gel_RINDEX.txt", cWvl, cRindex); // Mặc định dùng Gel
  mptCouplingDef->AddProperty("RINDEX", cWvl, cRindex);
  std::vector<G4double> cAbsE = { cWvl.front(), cWvl.back() };
  std::vector<G4double> cAbsL = { 1.0e10*m, 1.0e10*m };
  mptCouplingDef->AddProperty("ABSLENGTH", cAbsE, cAbsL);
  fCouplingMat->SetMaterialPropertiesTable(mptCouplingDef);

  fSensorMat = nist->FindOrBuildMaterial("G4_Si"); // Silicon (Đại diện cho khối SiPM)

  // ── 1. NƯỚC ──
  auto* mptWater = new G4MaterialPropertiesTable();
  std::vector<G4double> wWvl, wRindex, wAbsWvl, wAbs, wRayWvl, wRay;
  
  ReadOpticalData("../data/Water_RINDEX.txt", wWvl, wRindex);
  ReadOpticalData("../data/Water_ABS.txt", wAbsWvl, wAbs);
  ReadOpticalData("../data/Water_RAYLEIGH.txt", wRayWvl, wRay);
  
  for (auto& val : wAbs) val *= m; 
  for (auto& val : wRay) val *= m; 
  
  mptWater->AddProperty("RINDEX", wWvl, wRindex);
  mptWater->AddProperty("ABSLENGTH", wAbsWvl, wAbs);
  mptWater->AddProperty("RAYLEIGH", wRayWvl, wRay); 
  
  fWater->SetMaterialPropertiesTable(mptWater);

  // ── 2. NHỰA ACRYLIC ──
  auto* mptAcrylic = new G4MaterialPropertiesTable();
  std::vector<G4double> aWvl, aRindex, aAbsWvl, aAbs;
  ReadOpticalData("../data/Acrylic_RINDEX.txt", aWvl, aRindex);
  ReadOpticalData("../data/Acrylic_ABS.txt", aAbsWvl, aAbs);
  for (auto& val : aAbs) val *= mm; 
  mptAcrylic->AddProperty("RINDEX", aWvl, aRindex);
  mptAcrylic->AddProperty("ABSLENGTH", aAbsWvl, aAbs);
  fAcrylic->SetMaterialPropertiesTable(mptAcrylic);


  // ── DUMMY MPT CHO SENSOR 
  auto* mptSensorVol = new G4MaterialPropertiesTable();
  // Khai báo mảng năng lượng (từ 800nm đến 200nm)
  std::vector<G4double> sEne = { 1.55*eV, 6.20*eV }; 
  std::vector<G4double> sRindex = { 1.6, 1.6 }; // Chiết suất giả của Expoy
  std::vector<G4double> sAbs = { 1.0*nm, 1.0*nm }; // Hấp thụ ngay lập tức sau 1 nanomet
  mptSensorVol->AddProperty("RINDEX", sEne, sRindex);
  mptSensorVol->AddProperty("ABSLENGTH", sEne, sAbs);
  // Gán bảng thuộc tính này cho vật liệu khối Cảm biến
  fSensorMat->SetMaterialPropertiesTable(mptSensorVol);

  //tat cerenkov ring ngoai ko khi
  /*
  // ── 3. KHÔNG KHÍ (WORLD) ──
  auto* mptAir = new G4MaterialPropertiesTable();
  std::vector<G4double> airWvl, airRindex;
  ReadOpticalData("../data/Air_RINDEX.txt", airWvl, airRindex);
  mptAir->AddProperty("RINDEX", airWvl, airRindex);
  fAir->SetMaterialPropertiesTable(mptAir);*/
}

//--------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // ── KÍCH THƯỚC CÁC LỚP BỌC (Xếp lồng nhau) ──
  G4double waterHalf    = 50.0 * mm;
  G4double acrylicHalf  = 53.0 * mm;  // Bề dày nhựa = 3.0 mm
  G4double couplingHalf = 53.5 * mm;  // Bề dày lớp air/gel = 0.5 mm
  G4double sensorHalf   = 53.6 * mm; 

  // ── 1. WORLD ──
  auto solidWorld = new G4Box("World", 1.0*m, 1.0*m, 1.0*m);
  auto logicWorld = new G4LogicalVolume(solidWorld, fAir, "World");
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, fCheckOverlaps);

  // ── 2. CẢM BIẾN (PHYSICAL DETECTOR) ──
  auto solidSensor = new G4Box("SensorBox", sensorHalf, sensorHalf, sensorHalf);
  fLogicSensor = new G4LogicalVolume(solidSensor, fSensorMat, "SensorBox");
  fLogicSensor->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.5))); // Đỏ bán trong suốt
  auto physSensor = new G4PVPlacement(nullptr, G4ThreeVector(), fLogicSensor, "SensorBox", logicWorld, false, 0, fCheckOverlaps);

  // ── 3. LỚP COUPLING (MỠ / GEL / KHÔNG KHÍ) ──
  auto solidCoupling = new G4Box("CouplingBox", couplingHalf, couplingHalf, couplingHalf);
  fLogicCoupling = new G4LogicalVolume(solidCoupling, fCouplingMat, "CouplingBox");
  fLogicCoupling->SetVisAttributes(new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.3))); 
  auto physCoupling = new G4PVPlacement(nullptr, G4ThreeVector(), fLogicCoupling, "CouplingBox", fLogicSensor, false, 0, fCheckOverlaps);

  // ── 4. NHỰA ACRYLIC ──
  auto solidAcrylic = new G4Box("AcrylicBox", acrylicHalf, acrylicHalf, acrylicHalf);
  auto logicAcrylic = new G4LogicalVolume(solidAcrylic, fAcrylic, "AcrylicBox");
  logicAcrylic->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.2))); 
  auto physAcrylic = new G4PVPlacement(nullptr, G4ThreeVector(), logicAcrylic, "AcrylicBox", fLogicCoupling, false, 0, fCheckOverlaps);

  // ── 5. LÕI NƯỚC ──
  auto solidWater = new G4Box("WaterBox", waterHalf, waterHalf, waterHalf);
  auto logicWater = new G4LogicalVolume(solidWater, fWater, "WaterBox");
  auto* watVis = new G4VisAttributes(G4Colour(0.0, 0.4, 1.0, 0.5));
  watVis->SetForceSolid(true);
  logicWater->SetVisAttributes(watVis);
  auto physWater = new G4PVPlacement(nullptr, G4ThreeVector(), logicWater, "WaterBox", logicAcrylic, false, 0, fCheckOverlaps);


  // =========================================================================
  // ĐỊNH NGHĨA QUANG HỌC BỀ MẶT BẰNG G4LogicalBorderSurface
  // =========================================================================
  
  // A. Khúc xạ/Phản xạ điện môi (Nước <-> Nhựa <-> Mỡ)
  auto* polishedSurface = new G4OpticalSurface("PolishedSurface");
  polishedSurface->SetType(dielectric_dielectric);
  polishedSurface->SetModel(unified);
  polishedSurface->SetFinish(polished); 
  
  new G4LogicalBorderSurface("W2A", physWater, physAcrylic, polishedSurface);
  new G4LogicalBorderSurface("A2W", physAcrylic, physWater, polishedSurface);
  new G4LogicalBorderSurface("A2C", physAcrylic, physCoupling, polishedSurface);
  new G4LogicalBorderSurface("C2A", physCoupling, physAcrylic, polishedSurface);

  // B. Bề mặt Cảm biến (Coupling <-> Sensor) - dielectric_metal
  auto* sensorSurface = new G4OpticalSurface("SensorSurface");
  sensorSurface->SetType(dielectric_metal);
  sensorSurface->SetFinish(polished);
  sensorSurface->SetModel(glisur);

  std::vector<G4double> e_pde, pde;
  ReadOpticalData("../data/SiPM_PDE.txt", e_pde, pde);

  std::vector<G4double> reflectivity;
  for (double val : pde) { reflectivity.push_back(1.0 - val); }

  auto* mptSensor = new G4MaterialPropertiesTable();
  mptSensor->AddProperty("EFFICIENCY", e_pde, pde);
  mptSensor->AddProperty("REFLECTIVITY", e_pde, reflectivity);
  sensorSurface->SetMaterialPropertiesTable(mptSensor);

  // Gán bề mặt cho hạt đi từ Coupling chạm vào Sensor
  new G4LogicalBorderSurface("Coup2Sensor", physCoupling, physSensor, sensorSurface);

  return physWorld;
}

//--------------------------------------------------------------------
void DetectorConstruction::ConstructSDandField() {}

//--------------------------------------------------------------------
void DetectorConstruction::SetCouplingFile(G4String filepath) {
    std::vector<G4double> energies, rindex;
    ReadOpticalData(filepath, energies, rindex); 

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", energies, rindex);

    std::vector<G4double> absEnergies = { energies.front(), energies.back() };
    std::vector<G4double> absLengths  = { 1.0e10*m, 1.0e10*m }; // Coi mỡ trong suốt tuyệt đối
    mpt->AddProperty("ABSLENGTH", absEnergies, absLengths);

    if (fCouplingMat) {
        fCouplingMat->SetMaterialPropertiesTable(mpt);
        
        if (G4RunManager::GetRunManager()->GetUserDetectorConstruction() != nullptr) {
            G4RunManager::GetRunManager()->PhysicsHasBeenModified(); 
        }
        
        G4cout << "\n=================================================================" << G4endl;
        G4cout << "====> Set RINDEX Completed: " << filepath << G4endl;
        G4cout << "=================================================================\n" << G4endl;
    }
}