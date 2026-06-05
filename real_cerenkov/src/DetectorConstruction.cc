#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"

#include <fstream>
#include <sstream>
#include <algorithm>

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction(), fCheckOverlaps(true), fWater(nullptr), fAcrylic(nullptr), fAir(nullptr), fCouplingMat(nullptr), fSensorMat(nullptr), fLogicCoupling(nullptr), fLogicSensor(nullptr) {
  // === GIÁ TRỊ MẶC ĐỊNH ===
  fTyvekConfig = "Inside";
  fScintSize = G4ThreeVector(100.0*mm, 10.0*mm, 100.0*mm); 
  fSipmSize = G4ThreeVector(25.0*mm, 25.0*mm, 1.0*mm);
  
  fMessenger = new DetectorMessenger(this);
  DefineMaterials();
}
DetectorConstruction::~DetectorConstruction() { delete fMessenger; }

// --- CÁC HÀM SETTER TỪ MACRO ---
void DetectorConstruction::SetTyvekPosition(G4String pos) { fTyvekConfig = pos; }
void DetectorConstruction::SetScintSize(G4ThreeVector size) { fScintSize = size; }
void DetectorConstruction::SetSipmSize(G4ThreeVector size) { fSipmSize = size; }

void DetectorConstruction::ReadOpticalData(const G4String& filename, std::vector<G4double>& photonEnergies, std::vector<G4double>& properties) {
  std::ifstream file(filename);
  if (!file.is_open()) { G4Exception("DetectorConstruction::ReadOpticalData", "FileOpenError", FatalException, "Khong the mo file du lieu quang hoc."); }
  struct DataPoint { G4double energy; G4double value; };
  std::vector<DataPoint> data; G4String line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line); G4double wavelength, val;
    if (iss >> wavelength >> val) data.push_back({(1239.84193 * eV * nm) / (wavelength * nm), val});
  }
  std::sort(data.begin(), data.end(), [](const DataPoint& a, const DataPoint& b) { return a.energy < b.energy; });
  photonEnergies.clear(); properties.clear();
  for (const auto& pt : data) { photonEnergies.push_back(pt.energy); properties.push_back(pt.value); }
}

void DetectorConstruction::DefineMaterials() {
  auto nist = G4NistManager::Instance();
  fAir = nist->FindOrBuildMaterial("G4_AIR");
  fWater = nist->FindOrBuildMaterial("G4_WATER");
  fAcrylic = nist->FindOrBuildMaterial("G4_PLEXIGLASS"); 
  fCouplingMat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE"); 
  fSensorMat = nist->FindOrBuildMaterial("G4_Si");
  fScintMat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  auto* mptCoupling = new G4MaterialPropertiesTable();
  std::vector<G4double> cWvl, cRindex; ReadOpticalData("../data/Gel_RINDEX.txt", cWvl, cRindex); 
  mptCoupling->AddProperty("RINDEX", cWvl, cRindex);
  std::vector<G4double> cAbsE = { cWvl.front(), cWvl.back() }; std::vector<G4double> cAbsL = { 1.0e10*m, 1.0e10*m };
  mptCoupling->AddProperty("ABSLENGTH", cAbsE, cAbsL);
  fCouplingMat->SetMaterialPropertiesTable(mptCoupling);

  auto* mptWater = new G4MaterialPropertiesTable();
  std::vector<G4double> wWvl, wRindex, wAbsWvl, wAbs, wRayWvl, wRay;
  ReadOpticalData("../data/Water_RINDEX.txt", wWvl, wRindex); ReadOpticalData("../data/Water_ABS.txt", wAbsWvl, wAbs); ReadOpticalData("../data/Water_RAYLEIGH.txt", wRayWvl, wRay);
  for (auto& val : wAbs) val *= m; for (auto& val : wRay) val *= m; 
  mptWater->AddProperty("RINDEX", wWvl, wRindex); mptWater->AddProperty("ABSLENGTH", wAbsWvl, wAbs); mptWater->AddProperty("RAYLEIGH", wRayWvl, wRay); 
  fWater->SetMaterialPropertiesTable(mptWater);

  auto* mptAcrylic = new G4MaterialPropertiesTable();
  std::vector<G4double> aWvl, aRindex, aAbsWvl, aAbs;
  ReadOpticalData("../data/Acrylic_RINDEX.txt", aWvl, aRindex); ReadOpticalData("../data/Acrylic_ABS.txt", aAbsWvl, aAbs);
  for (auto& val : aAbs) val *= mm; 
  mptAcrylic->AddProperty("RINDEX", aWvl, aRindex); mptAcrylic->AddProperty("ABSLENGTH", aAbsWvl, aAbs);
  fAcrylic->SetMaterialPropertiesTable(mptAcrylic);

  auto* mptSensorVol = new G4MaterialPropertiesTable();
  std::vector<G4double> sEne = { 1.55*eV, 6.20*eV }; std::vector<G4double> sRindex = { 1.6, 1.6 }; std::vector<G4double> sAbs = { 1.0*nm, 1.0*nm }; 
  mptSensorVol->AddProperty("RINDEX", sEne, sRindex); mptSensorVol->AddProperty("ABSLENGTH", sEne, sAbs);
  fSensorMat->SetMaterialPropertiesTable(mptSensorVol);
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  G4double waterHalf    = 50.0 * mm;
  G4double acrylicHalf  = 53.0 * mm;  
  
  // Lấy dữ liệu kích thước từ Macro
  G4double scintHalfX = fScintSize.x() / 2.0;
  G4double scintHalfY = fScintSize.y() / 2.0;
  G4double scintHalfZ = fScintSize.z() / 2.0;
  G4double scintPosY  = 72.0 * mm; 

  G4double sipmSizeX = fSipmSize.x(); 
  G4double sipmSizeY = fSipmSize.y();
  G4double sensorHalfZ  = fSipmSize.z() / 2.0;
  G4double couplingHalfZ= 0.25 * mm; 
  
  std::vector<G4double> sipmPositions = {-25.0 * mm, 25.0 * mm}; 

  // 1. WORLD (DARK BOX)
  auto solidWorld = new G4Box("World", 1.0*m, 1.0*m, 1.0*m);
  auto logicWorld = new G4LogicalVolume(solidWorld, fAir, "World");
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, fCheckOverlaps);

  // 2. HODOSCOPE (TRIGGER)
  auto solidScint = new G4Box("ScintBox", scintHalfX, scintHalfY, scintHalfZ);
  auto logicScint = new G4LogicalVolume(solidScint, fScintMat, "ScintLogic");
  logicScint->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.5))); 
  new G4PVPlacement(nullptr, G4ThreeVector(0, scintPosY, 0), logicScint, "ScintTop", logicWorld, false, 0, fCheckOverlaps);
  new G4PVPlacement(nullptr, G4ThreeVector(0, -scintPosY, 0), logicScint, "ScintBottom", logicWorld, false, 1, fCheckOverlaps);

  // 3. TANK NHỰA ACRYLIC & NƯỚC
  auto solidAcrylic = new G4Box("AcrylicBox", acrylicHalf, acrylicHalf, acrylicHalf);
  auto logicAcrylic = new G4LogicalVolume(solidAcrylic, fAcrylic, "AcrylicBox");
  logicAcrylic->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.1))); 
  auto physAcrylic = new G4PVPlacement(nullptr, G4ThreeVector(), logicAcrylic, "AcrylicBox", logicWorld, false, 0, fCheckOverlaps);

  auto solidWater = new G4Box("WaterBox", waterHalf, waterHalf, waterHalf);
  auto logicWater = new G4LogicalVolume(solidWater, fWater, "WaterBox");
  auto* watVis = new G4VisAttributes(G4Colour(0.0, 0.4, 1.0, 0.2)); watVis->SetForceSolid(true); logicWater->SetVisAttributes(watVis);
  auto physWater = new G4PVPlacement(nullptr, G4ThreeVector(), logicWater, "WaterPhysical", logicAcrylic, false, 0, fCheckOverlaps);

  // =========================================================================
  // 4. KIẾN TRÚC RẼ NHÁNH: TYVEK INSIDE VÀ TYVEK OUTSIDE
  // =========================================================================
  G4VSolid* solidTyvek = nullptr;
  G4LogicalVolume* logicTyvek = nullptr;
  G4double tyvekThickness = 0.1 * mm;
  auto punchHole = new G4Box("Hole", sipmSizeX/2, sipmSizeY/2, 2.0*mm); // Lưỡi dao cắt (dài 2mm để đảm bảo thủng)

  if (fTyvekConfig == "Inside") {
      G4double tyvekHalfIn = waterHalf + tyvekThickness; // 50.1
      auto outerTyvek = new G4Box("OuterTyvek", tyvekHalfIn, tyvekHalfIn, tyvekHalfIn);
      auto innerTyvek = new G4Box("InnerTyvek", waterHalf, waterHalf, waterHalf);
      solidTyvek = new G4SubtractionSolid("TyvekShell_Init", outerTyvek, innerTyvek);
      
      int holeCount = 0;
      for (int face = 0; face < 4; face++) {
          for (G4double yPos : sipmPositions) {
              for (G4double hPos : sipmPositions) { 
                  G4ThreeVector pos; G4RotationMatrix* rot = new G4RotationMatrix();
                  if (face == 0)      { pos = G4ThreeVector(tyvekHalfIn, yPos, hPos); rot->rotateY(90*deg); }
                  else if (face == 1) { pos = G4ThreeVector(-tyvekHalfIn, yPos, hPos); rot->rotateY(-90*deg); }
                  else if (face == 2) { pos = G4ThreeVector(hPos, yPos, tyvekHalfIn); } 
                  else if (face == 3) { pos = G4ThreeVector(hPos, yPos, -tyvekHalfIn); }
                  solidTyvek = new G4SubtractionSolid("TyvekShell_" + std::to_string(++holeCount), solidTyvek, punchHole, rot, pos);
              }
          }
      }
      logicTyvek = new G4LogicalVolume(solidTyvek, fAcrylic, "TyvekLogic");
      logicTyvek->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.9))); 
      new G4PVPlacement(nullptr, G4ThreeVector(), logicTyvek, "TyvekPhysical", logicAcrylic, false, 0, fCheckOverlaps); // Nằm trong Acrylic
      
      G4cout << ">>> DETECTOR CONFIG: MANG TYVEK DUOC DAT TRONG BE NHUA (WAVEGUIDE ENABLED)" << G4endl;
  } 
  else if (fTyvekConfig == "Outside") {
      G4double tyvekHalfOut = acrylicHalf + tyvekThickness; // 53.1
      auto outerTyvek = new G4Box("OuterTyvek", tyvekHalfOut, tyvekHalfOut, tyvekHalfOut);
      auto innerTyvek = new G4Box("InnerTyvek", acrylicHalf, acrylicHalf, acrylicHalf);
      solidTyvek = new G4SubtractionSolid("TyvekShell_Init", outerTyvek, innerTyvek);
      
      int holeCount = 0;
      for (int face = 0; face < 4; face++) {
          for (G4double yPos : sipmPositions) {
              for (G4double hPos : sipmPositions) { 
                  G4ThreeVector pos; G4RotationMatrix* rot = new G4RotationMatrix();
                  if (face == 0)      { pos = G4ThreeVector(tyvekHalfOut, yPos, hPos); rot->rotateY(90*deg); }
                  else if (face == 1) { pos = G4ThreeVector(-tyvekHalfOut, yPos, hPos); rot->rotateY(-90*deg); }
                  else if (face == 2) { pos = G4ThreeVector(hPos, yPos, tyvekHalfOut); } 
                  else if (face == 3) { pos = G4ThreeVector(hPos, yPos, -tyvekHalfOut); }
                  solidTyvek = new G4SubtractionSolid("TyvekShell_" + std::to_string(++holeCount), solidTyvek, punchHole, rot, pos);
              }
          }
      }
      logicTyvek = new G4LogicalVolume(solidTyvek, fAcrylic, "TyvekLogic");
      logicTyvek->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.9))); 
      new G4PVPlacement(nullptr, G4ThreeVector(), logicTyvek, "TyvekPhysical", logicWorld, false, 0, fCheckOverlaps); // Bọc ngoài Acrylic
      
      G4cout << ">>> DETECTOR CONFIG: MANG TYVEK BOC NGOAI BE NHUA (WAVEGUIDE KILLED)" << G4endl;
  }

  // =========================================================================
  // QUANG HỌC BỀ MẶT
  // =========================================================================
  auto* tyvekSurface = new G4OpticalSurface("TyvekSurface");
  tyvekSurface->SetType(dielectric_metal); 
  tyvekSurface->SetModel(unified);
  tyvekSurface->SetFinish(ground); 
  std::vector<G4double> e_ref = {1.55*eV, 6.20*eV}, r_ref = {0.95, 0.95}; 
  auto* mptTyvek = new G4MaterialPropertiesTable(); mptTyvek->AddProperty("REFLECTIVITY", e_ref, r_ref);
  tyvekSurface->SetMaterialPropertiesTable(mptTyvek);
  
  // Phủ mặt Tyvek: Chỉ 1 dòng lệnh này sẽ áp dụng hiệu ứng phản xạ cho MỌI bề mặt của Tyvek (dù nó ở trong hay ngoài)
  if(logicTyvek) new G4LogicalSkinSurface("TyvekSkin", logicTyvek, tyvekSurface);

  auto* transparentSurface = new G4OpticalSurface("TransSurf");
  transparentSurface->SetType(dielectric_dielectric);
  transparentSurface->SetModel(unified);
  transparentSurface->SetFinish(polished); 
  new G4LogicalBorderSurface("W2A", physWater, physAcrylic, transparentSurface);
  new G4LogicalBorderSurface("A2W", physAcrylic, physWater, transparentSurface);

  auto* sensorSurface = new G4OpticalSurface("SensorSurface");
  sensorSurface->SetType(dielectric_metal); sensorSurface->SetFinish(polished); sensorSurface->SetModel(glisur);
  std::vector<G4double> e_pde, pde; ReadOpticalData("../data/SiPM_PDE.txt", e_pde, pde);
  std::vector<G4double> sipmReflectivity(pde.size(), 0.20); 
  auto* mptSensor = new G4MaterialPropertiesTable();
  mptSensor->AddProperty("EFFICIENCY", e_pde, pde); mptSensor->AddProperty("REFLECTIVITY", e_pde, sipmReflectivity);
  sensorSurface->SetMaterialPropertiesTable(mptSensor);

  // =========================================================================
  // 5. GẮN SIPM VÀ GEL (Luôn luôn nằm sát mặt ngoài vách Acrylic 53.0mm)
  // =========================================================================
  auto solidCoupling = new G4Box("CouplingBox", sipmSizeX/2, sipmSizeY/2, couplingHalfZ);
  fLogicCoupling = new G4LogicalVolume(solidCoupling, fCouplingMat, "CouplingBox");
  fLogicCoupling->SetVisAttributes(new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.5))); 
  auto solidSensor = new G4Box("SensorBox", sipmSizeX/2, sipmSizeY/2, sensorHalfZ);
  fLogicSensor = new G4LogicalVolume(solidSensor, fSensorMat, "SensorBox");
  fLogicSensor->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.8))); 

  G4int copyNo = 0; 
  for (int face = 0; face < 4; face++) {
      for (G4double yPos : sipmPositions) {
          for (G4double hPos : sipmPositions) { 
              G4ThreeVector posC, posS;
              G4RotationMatrix* rot = new G4RotationMatrix();
              if (face == 0) {
                  posC = G4ThreeVector(acrylicHalf + couplingHalfZ, yPos, hPos);      
                  posS = G4ThreeVector(acrylicHalf + 2*couplingHalfZ + sensorHalfZ, yPos, hPos); rot->rotateY(90*deg);
              } else if (face == 1) { 
                  posC = G4ThreeVector(-(acrylicHalf + couplingHalfZ), yPos, hPos);
                  posS = G4ThreeVector(-(acrylicHalf + 2*couplingHalfZ + sensorHalfZ), yPos, hPos); rot->rotateY(-90*deg);
              } else if (face == 2) { 
                  posC = G4ThreeVector(hPos, yPos, acrylicHalf + couplingHalfZ);
                  posS = G4ThreeVector(hPos, yPos, acrylicHalf + 2*couplingHalfZ + sensorHalfZ);
              } else if (face == 3) { 
                  posC = G4ThreeVector(hPos, yPos, -(acrylicHalf + couplingHalfZ));
                  posS = G4ThreeVector(hPos, yPos, -(acrylicHalf + 2*couplingHalfZ + sensorHalfZ)); rot->rotateY(180*deg);
              }

              auto physCoupling = new G4PVPlacement(rot, posC, fLogicCoupling, "CouplingBox", logicWorld, false, copyNo, fCheckOverlaps);
              auto physSensor = new G4PVPlacement(rot, posS, fLogicSensor, "SensorBox", logicWorld, false, copyNo, fCheckOverlaps);

              new G4LogicalBorderSurface("A2C_Window", physAcrylic, physCoupling, transparentSurface);
              new G4LogicalBorderSurface("Coup2Sensor", physCoupling, physSensor, sensorSurface);
              copyNo++; 
          }
      }
  }
  return physWorld;
}

void DetectorConstruction::ConstructSDandField() {}

void DetectorConstruction::SetCouplingFile(G4String filepath) {
    std::vector<G4double> energies, rindex; ReadOpticalData(filepath, energies, rindex); 
    auto* mpt = new G4MaterialPropertiesTable(); mpt->AddProperty("RINDEX", energies, rindex);
    std::vector<G4double> absEnergies = { energies.front(), energies.back() };
    std::vector<G4double> absLengths  = { 1.0e10*m, 1.0e10*m }; 
    mpt->AddProperty("ABSLENGTH", absEnergies, absLengths);
    if (fCouplingMat) {
        fCouplingMat->SetMaterialPropertiesTable(mpt);
        if (G4RunManager::GetRunManager()->GetUserDetectorConstruction() != nullptr) { G4RunManager::GetRunManager()->PhysicsHasBeenModified(); }
    }
}