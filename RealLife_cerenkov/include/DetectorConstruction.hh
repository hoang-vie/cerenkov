#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
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

  void SetCouplingFile(G4String filepath);

protected:
  void DefineMaterials();

  G4bool   fCheckOverlaps;

  // Materials
  G4Material* fWater;
  G4Material* fAcrylic;
  G4Material* fAir;
  G4Material* fCouplingMat; 
  G4Material* fSensorMat;   // Vật liệu cho vỏ Cảm biến (Si hoặc Glass)

  G4LogicalVolume* fLogicCoupling; 
  G4LogicalVolume* fLogicSensor;

private:
  DetectorMessenger* fMessenger;

  void ReadOpticalData(const G4String& filename, 
                       std::vector<G4double>& photonEnergies, 
                       std::vector<G4double>& properties);
};

#endif