# Framework configuration audit — general-v1.10

This document is the human-readable companion to `tools/check_config_metadata_contract.py`.
Its purpose is to keep the public configuration surface explicit and to ensure that ordinary detector studies remain reproducible from ROOT metadata rather than from a macro filename or an operator's memory.

## Audit result

- **51 executable project commands** are currently registered.
- **6 UI directories** group those commands (`/detector/`, `/readout/`, `/trigger/`, `/output/`, `/generator/`, `/CRY/`).
- Every registered executable command has an explicit metadata/reproducibility contract.
- Five commands are retained only as legacy aliases.
- Three commands are expert/validation optical-table overrides.
- `/output/fileName` is intentionally not duplicated in `RunConfig` because it changes storage naming, not simulated detector/physics state.

Run the mechanical check with:

```bash
python3 tools/check_config_metadata_contract.py
```

Adding a new project command without adding an explicit contract entry makes this check fail.

## Detector and geometry commands

| Command | Authoritative resolved metadata | Notes |
|---|---|---|
| `/detector/tyvekPosition` | `RunConfig.TyvekConfig` | `Inside` / `Outside`. |
| `/detector/coreMaterial` | `RunConfig.CoreMaterial`, `CoreProfile` | Legacy convenience selector; prefer `coreProfile`. |
| `/detector/coreProfile` | `RunConfig.CoreProfile` plus profile dataset/status fields | Canonical curated Core material selector. |
| `/detector/coreSize` | `RunConfig.CoreX/Y/Z_mm` | Full Core dimensions. |
| `/detector/scintSize` | `RunConfig.ScintX/Y/Z_mm` | Full scintillator dimensions. |
| `/detector/readoutSize` | `RunConfig.ReadoutX/Y/Z_mm`, `ChannelMap.SizeU/V_mm` | Canonical readout-element dimensions. |
| `/detector/sipmSize` | same as `readoutSize` | Legacy alias. |
| `/detector/sipmPitch` | `ReadoutPlaneConfig.PitchU/V_mm` | Legacy 2x2-uniform-grid helper. |
| `/detector/tyvekHoleSize` | `RunConfig.HoleX_mm`, `HoleY_mm` | Legacy 3-vector interface; Z is ignored by geometry. |
| `/detector/couplingSize` | `RunConfig.CouplingX/Y/Z_mm` | Coupling element dimensions. |
| `/detector/outerShellThickness` | `RunConfig.OuterShellThickness_mm` | Canonical shell thickness. |
| `/detector/useOuterShell` | `RunConfig.UseOuterShell` | Canonical shell enable flag. |
| `/detector/outerShellProfile` | `RunConfig.OuterShellProfile` plus profile fields | Currently curated `PMMA_SNO_A`. |
| `/detector/setAcrylicThickness` | `RunConfig.OuterShellThickness_mm` | Legacy alias. |
| `/detector/useAcrylicShell` | `RunConfig.UseOuterShell` | Legacy alias. |
| `/detector/setTyvekThickness` | `RunConfig.TyvekThickness_mm` | Tyvek disabled when thickness is zero. |
| `/detector/setSiPMOffsetY` | `ReadoutPlaneConfig.OffsetV_mm` | Legacy all-face offset helper. |
| `/detector/scintGap` | `RunConfig.ScintGap_mm` | Trigger-paddle gap from detector envelope. |
| `/detector/sourceDistance` | `RunConfig.SourceDistance_mm` | Used by automatic source-Y placement. |
| `/detector/couplingProfile` | `RunConfig.CouplingProfile` plus profile fields | Currently curated `PDMS_Generic`. |
| `/detector/sensorMode` | `RunConfig.SensorMode` | `DirectSiPM` / `ConcentratorInput` (`SiPM` retained as accepted alias). |

## Expert optical-table overrides

These commands exist for controlled validation/expert work and are not the preferred production material API.

| Command | Resolved metadata | Reproducibility note |
|---|---|---|
| `/detector/waterRayleighFile` | `WaterRayleighDataset`, profile Rayleigh dataset/status, optical-data fingerprint | Restricted to a versioned project `data/` basename. |
| `/detector/pmmaAbsFile` | `AcrylicAbsDataset`, active PMMA profile fields, optical-data fingerprint | Restricted to a versioned project `data/` basename. |
| `/detector/couplingFile` | coupling profile RINDEX dataset/status, optical-data fingerprint | Legacy expert override; unlike the two commands above it can still name an external path. The exact bytes are fingerprinted, but portable replay requires preserving that external file. |

The last behavior is intentionally documented rather than silently changed in v1.10, because it is a legacy expert API. New ordinary studies should select the curated coupling profile instead.

## Readout layout commands

| Command | Authoritative resolved metadata |
|---|---|
| `/readout/uniformGrid` | six `ReadoutPlaneConfig` rows plus actual `ChannelMap` |
| `/readout/planeGrid` | per-face `Enabled`, `NU/NV`, `PitchU/V`, `OffsetU/V`; actual `ChannelMap` |
| `/readout/planeAuto` | per-face `Enabled`, `AutoFill`, `GapU/V`, `MarginU/V`; resolved elements in `ChannelMap` |

`ReadoutPlaneConfig` records the requested policy. `ChannelMap` records what geometry was actually instantiated. Both are required for a self-describing flexible-layout run.

## Trigger commands

| Command | RunConfig field |
|---|---|
| `/trigger/scintThreshold` | `ScintThreshold_MeV` |
| `/trigger/coincidenceWindow` | `CoincidenceWindow_ns` |
| `/trigger/integrationWindow` | `IntegrationWindow_ns` |

## Output commands

| Command | RunConfig field / contract |
|---|---|
| `/output/eventMode` | `OutputEventMode` |
| `/output/incidentDetail` | `IncidentDetail` |
| `/output/photonBirthDetail` | `PhotonBirthDetail` |
| `/output/coreBoundaryMode` | `CoreBoundaryMode` |
| `/output/coreBoundaryDetail` | `CoreBoundaryDetail` |
| `/output/boundaryDebug` | `BoundaryDebug` |
| `/output/maxOpticalBounces` | `MaxOpticalBounces` |
| `/output/compressionLevel` | `CompressionLevel` |
| `/output/fileName` | storage destination/name; intentionally not physics metadata |

## CRY and generic generator commands

| Command | Authoritative resolved metadata |
|---|---|
| `/CRY/file` | full resolved `RunConfig.CRYConfig` file contents |
| `/CRY/input` | full resolved `RunConfig.CRYConfig` after `/CRY/update` |
| `/CRY/update` | action command; the resolved configuration is stored in `CRYConfig` |
| `/CRY/timeMode` | `CRYTimeMode` |
| `/generator/mode` | `GeneratorMode` |
| `/generator/particle` | `GunParticle` |
| `/generator/energy` | `GunEnergy_MeV` |
| `/generator/position` | `GunX/Y/Z_mm` |
| `/generator/direction` | `GunDirX/Y/Z` |
| `/generator/polarization` | `GunPolX/Y/Z` |
| `/generator/time` | `GunTime_ns` |
| `/generator/useAutoY` | `GunUseAutoY`; `GunY_mm` stores the resolved automatic Y when enabled |

## Simulation-affecting settings that are intentionally not project commands

Not every source-level constant should be exposed as a scan parameter. The following are intentionally fixed by the validated framework baseline or derived automatically:

- physics list: `FTFP_BERT + G4EmStandardPhysics_option4 + G4OpticalPhysics`;
- default production cut (runtime value is stored in `RunConfig`);
- optical process activation (runtime states are stored in `RunConfig`; Geant4's own UI may change them for validation);
- air, scintillator bulk, silicon readout bulk and dummy Tyvek bulk material identities;
- Tyvek optical-surface implementation (`dielectric_metal`, UNIFIED, ground) and sensor/interface surface models;
- effective SiPM PDE dataset and operating-point label;
- effective Tyvek total-reflectivity dataset;
- rectangular/cuboid detector topology and six axis-aligned face convention;
- world size, which is derived from the configured detector/source envelope with an internal safety margin;
- geometry overlap checking, which remains enabled.

These assumptions are part of the source revision identified by `GitCommit`/`ProjectRevision`. They should become user-facing commands only if a real study requirement appears; if exposed later, the framework completion rule requires metadata, validation and documentation at the same time.

## Runtime evidence for the v1.10 contract

The target Geant4 machine has exercised the following after applying the v1.10 patch:

1. **MT/output/provenance smoke** — 100 events, 2 threads, zero readout channels: `OUTPUT SANITY: PASS`; `FRAMEWORK METADATA VALIDATION: PASS`.
2. **Flexible per-face layout** — 10 events with a deliberately asymmetric 12-channel layout: `ChannelMap=12`, `ReadoutPlaneConfig=6`, `OUTPUT SANITY: PASS`; metadata PASS.
3. **Integrated DirectSiPM readout** — 100 events, 24 channels, 4 threads: `Photon=1,002,772`, `SensorIncident=21,384`, `DetectedSignal=5,728`, `ChannelEvent=2,400`; output and metadata PASS.
4. **Core-boundary mode/detail/TIR regression** — OFF/basic/phaseSpace/TIR cases all PASS.

These tests establish framework wiring/configuration behavior. They do not replace the dedicated physics validation suite.
