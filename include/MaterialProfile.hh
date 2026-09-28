#ifndef MaterialProfile_h
#define MaterialProfile_h 1

#include "globals.hh"

class G4Material;

// Curated bulk-material definition used by detector geometry roles.
// This is intentionally a data/provenance container rather than an arbitrary
// runtime material loader: every profile must be registered in code with
// explicit optical datasets and unit conventions.
struct MaterialProfile {
  G4String profileName;        // stable framework profile ID
  G4String materialLabel;      // human-readable physical material label
  G4Material* material = nullptr;

  G4String rindexDataset;
  G4String absDataset;
  G4String rayleighDataset;    // empty when no Rayleigh table is configured

  // Numeric length values in the text datasets are multiplied by these Geant4
  // units before being attached to the material-properties table.
  G4double absLengthUnit = 1.0;
  G4double rayleighLengthUnit = 1.0;

  // Empirical/measured validity of the absorption model.
  // These are intentionally distinct from the numerical wavelength support
  // of absDataset.  Negative values mean "fall back to dataset support".
  G4double absMeasuredMinNm = -1.0;
  G4double absMeasuredMaxNm = -1.0;

  G4String rindexStatus;
  G4String absStatus;
  G4String rayleighStatus;
  G4String outOfRangePolicy;
};

#endif
