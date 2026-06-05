#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class DetectorMessenger; // Forward declaration

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

  // Hàm nhận lệnh từ Macro để đổi file RINDEX
  void SetCouplingFile(G4String filepath);

protected:
  void DefineMaterials();

  G4bool   fCheckOverlaps;

  // Materials
  G4Material* fWater;
  G4Material* fAcrylic;
  G4Material* fAir;
  G4Material* fCouplingMat; // Vật liệu linh hoạt cho lớp vỏ SD2

  G4LogicalVolume* fLogicSD;
  G4LogicalVolume* fLogicCoupling; // Giữ con trỏ để update Material

private:
  DetectorMessenger* fMessenger;

  // Hàm tiện ích để đọc dữ liệu quang học từ file txt/csv
  void ReadOpticalData(const G4String& filename, 
                       std::vector<G4double>& photonEnergies, 
                       std::vector<G4double>& properties);
};

#endif