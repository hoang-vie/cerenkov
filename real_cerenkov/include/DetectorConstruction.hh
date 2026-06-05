#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
  // Các hàm API để Messenger gọi vào
  void SetCouplingFile(G4String filepath);
  void SetTyvekPosition(G4String pos);
  void SetScintSize(G4ThreeVector size);
  void SetSipmSize(G4ThreeVector size);

protected:
  void DefineMaterials();
  G4bool fCheckOverlaps;

  G4Material* fWater;
  G4Material* fAcrylic;
  G4Material* fAir;
  G4Material* fCouplingMat; 
  G4Material* fSensorMat;   
  G4Material* fScintMat;

  G4LogicalVolume* fLogicCoupling; 
  G4LogicalVolume* fLogicSensor;

private:
  DetectorMessenger* fMessenger;
  void ReadOpticalData(const G4String& filename, 
                       std::vector<G4double>& photonEnergies, 
                       std::vector<G4double>& properties);

  // Biến lưu trữ cấu hình Macro
  G4String fTyvekConfig;
  G4ThreeVector fScintSize;
  G4ThreeVector fSipmSize;
};

#endif