# ROOT tree / branch reference — general-v1.10
> Baseline: **general-v1.10 framework metadata for the current simulation model**. This file follows the ntuple definitions in `src/RunAction.cc` and complements `SCHEMA.md`.
## 1. How to think about the ROOT file
All ntuples belong to one simulation run. Use IDs for joins; do not rely on ROOT row order. The most common join keys are:

- event-level: `(RunID, EventID)`;
- non-optical track: `(RunID, EventID, TrackID)`;
- optical photon: `(RunID, EventID, PhotonID)`;
- readout channel: `(RunID, ChannelID)` via `ChannelMap`.

`Photon.ParentTrackID` is the immediate Geant4 parent. `Photon.SourceTrackID` is the propagated non-optical optical-source track. For primary-emitter selections, resolve `ParentTrackID` to `Track` and require `Track.ParentID==0` plus the requested PDG. For cascade/ancestor selections, use `SourceTrackID` and the source track's `RootPrimaryTrackID`.

## 2. Output tiers

- `thinAll`: thin trigger/scintillator/primary information for all events, detailed track/photon/readout rows only for triggered events.
- `triggered`: detailed information only for triggered events, while `RunStatistics` retains generated-event denominators.
- `fullAll`: full detail for every event; intended for validation/debug.

## 3. Type legend

- **Int** — integer branch (`CreateNtupleIColumn`).
- **Float** — single-precision floating branch (`CreateNtupleFColumn`).
- **Double** — double-precision branch (`CreateNtupleDColumn`).
- **String** — string column; inspect backend representation before binding manually in ROOT C++.

## `TriggerEvent` — 13 branches

Thin event and trigger summary.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `EventID` | Int | Event identifier. |
| 2 | `TriggerPassed` | Int | 1 if the configured top/bottom trigger condition passed, otherwise 0. |
| 3 | `EdepTop_MeV` | Float | Total event energy deposited in the top paddle; trigger thresholding uses individual grouped pulses. |
| 4 | `EdepBottom_MeV` | Float | Bottom trigger-paddle counterpart. |
| 5 | `TimeTop_ns` | Float | Representative top trigger-pulse time. |
| 6 | `TimeBottom_ns` | Float | Representative bottom trigger-pulse time. |
| 7 | `NPrimaryNonOptical` | Int | Number of non-optical primary particles in the event. |
| 8 | `NNonOpticalTracks` | Int | Number of non-optical tracks. |
| 9 | `NOpticalTracks` | Int | Number of optical-photon tracks. |
| 10 | `NIncident` | Int | Readout-incident count for the event. |
| 11 | `NDetected` | Int | Detected-signal count for the event. |
| 12 | `EdepCoreTotal_MeV` | Float | Total deposited energy in the Core for the event. |

## `ScintDeposit` — 6 branches

Sparse raw scintillator step deposits for offline trigger reconstruction.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `EventID` | Int | Event identifier. |
| 2 | `PaddleID` | Int | Trigger scintillator identifier (top/bottom mapping defined by detector code). |
| 3 | `TrackID` | Int | Track responsible for this non-zero scintillator step deposit. |
| 4 | `Time_ns` | Float | Pre-step global time of the deposit. |
| 5 | `Edep_MeV` | Float | Deposited energy in this scintillator step. |

## `Track` — 37 branches

Generic non-optical track truth.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `EventID` | Int | Event identifier. |
| 2 | `TrackID` | Int | Non-optical Geant4 TrackID. |
| 3 | `ParentID` | Int | Immediate Geant4 parent; 0 denotes primary. |
| 4 | `RootPrimaryTrackID` | Int | Resolved root non-optical primary TrackID for ancestry. |
| 5 | `PDG` | Int | PDG code. |
| 6 | `Charge_e` | Float | Charge in units of elementary charge. |
| 7 | `CreatorProcessType` | Int | Creator process type. |
| 8 | `CreatorProcessSubType` | Int | Creator process subtype. |
| 9 | `CreatorProcessName` | String | Creator process name. |
| 10 | `BirthVolumeID` | Int | Volume category at track birth. |
| 11 | `VertexKE_MeV` | Float | Kinetic energy at track vertex. |
| 12 | `VertexTime_ns` | Float | Track vertex global time. |
| 13 | `VertexX_mm` | Float | Vertex x. |
| 14 | `VertexY_mm` | Float | Vertex y. |
| 15 | `VertexZ_mm` | Float | Vertex z. |
| 16 | `DirX` | Float | Vertex direction x. |
| 17 | `DirY` | Float | Vertex direction y. |
| 18 | `DirZ` | Float | Vertex direction z. |
| 19 | `CoreEntrySet` | Int | 1 if a Core logical-volume entry was recorded. |
| 20 | `CoreEntryX_mm` | Float | Global x of first Core entry. |
| 21 | `CoreEntryY_mm` | Float | Global y of first Core entry. |
| 22 | `CoreEntryZ_mm` | Float | Global z of first Core entry. |
| 23 | `CoreExitSet` | Int | 1 if a Core logical-volume exit was recorded. |
| 24 | `CoreExitX_mm` | Float | Global x of last Core exit. |
| 25 | `CoreExitY_mm` | Float | Global y of last Core exit. |
| 26 | `CoreExitZ_mm` | Float | Global z of last Core exit. |
| 27 | `EndTime_ns` | Float | Track end global time. |
| 28 | `EndProperTime_ns` | Float | Track proper time at end. |
| 29 | `EndX_mm` | Float | Track end x. |
| 30 | `EndY_mm` | Float | Track end y. |
| 31 | `EndZ_mm` | Float | Track end z. |
| 32 | `EndVolumeID` | Int | Volume category at end. |
| 33 | `EndProcessType` | Int | Ending process type. |
| 34 | `EndProcessSubType` | Int | Ending process subtype. |
| 35 | `EndProcessName` | String | Ending process name. |
| 36 | `EndTrackStatus` | Int | Integer Geant4 track status at finalization. |

## `TrackVolume` — 8 branches

Per-track per-volume deposited energy and path length.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `EventID` | Int | Event identifier. |
| 2 | `TrackID` | Int | Non-optical track ID. |
| 3 | `VolumeID` | Int | Logical-volume category; resolve via VolumeMap. |
| 4 | `Edep_MeV` | Float | Sum of local deposited energy by this track in this volume. |
| 5 | `DeltaKE_MeV` | Float | Sum of KE_pre-KE_post for the track in this volume. |
| 6 | `Length_mm` | Float | Accumulated path length in the volume. |
| 7 | `NSteps` | Int | Number of contributing steps. |

## `Photon` — 36 branches

One compact summary row per optical photon. `PhotonID` is the optical Geant4 TrackID.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `EventID` | Int | Event identifier. |
| 2 | `PhotonID` | Int | Optical photon Geant4 TrackID. |
| 3 | `ParentTrackID` | Int | Immediate Geant4 parent TrackID. |
| 4 | `SourceTrackID` | Int | Propagated non-optical source track used across possible optical re-emission chains. |
| 5 | `CreatorProcessType` | Int | Geant4 creator-process type. |
| 6 | `CreatorProcessSubType` | Int | Geant4 creator-process subtype. |
| 7 | `OriginVolumeID` | Int | Volume category where the optical photon was created; resolve with VolumeMap. |
| 8 | `Wavelength_nm` | Float | Photon wavelength. |
| 9 | `BirthTime_ns` | Float | Photon creation global time. |
| 10 | `EmissionAngle_deg` | Float | Cherenkov emission angle when defined; NaN otherwise. |
| 11 | `BirthX_mm` | Float | Optional birth global x; may be NaN/omitted by detail policy. |
| 12 | `BirthY_mm` | Float | Optional birth global y. |
| 13 | `BirthZ_mm` | Float | Optional birth global z. |
| 14 | `BirthDirX` | Float | Optional birth direction x. |
| 15 | `BirthDirY` | Float | Optional birth direction y. |
| 16 | `BirthDirZ` | Float | Optional birth direction z. |
| 17 | `ReachedCoreWall` | Int | 1 after the first Core logical-volume wall encounter is recorded. |
| 18 | `FirstWallTime_ns` | Float | Time of first Core-wall encounter. |
| 19 | `FirstWallX_mm` | Float | Global x at first Core-wall encounter. |
| 20 | `FirstWallY_mm` | Float | Global y at first Core-wall encounter. |
| 21 | `FirstWallZ_mm` | Float | Global z at first Core-wall encounter. |
| 22 | `ExtractedToReadout` | Int | 1 once the photon is classified as extracted toward readout. |
| 23 | `ExtractionTime_ns` | Float | Time of first extraction-to-readout classification. |
| 24 | `NBoundaryReflections` | Int | Total counted optical boundary reflections. |
| 25 | `NTIR` | Int | Total counted total-internal-reflection occurrences. |
| 26 | `NTyvekReflections` | Int | Total counted Tyvek reflections. |
| 27 | `NRayleigh` | Int | Total Rayleigh scatters. |
| 28 | `EndReason` | Int | Photon end classification; see EndReason table below. |
| 29 | `EndVolumeID` | Int | Volume category where track ended. |
| 30 | `EndTime_ns` | Float | Photon end global time. |
| 31 | `KilledBySafetyCut` | Int | 1 if killed by configured max-optical-bounce safety cut. |
| 32 | `FirstRayleighTrackLength_mm` | Float | Track length at first Rayleigh scatter; used for MFP validation. |
| 33 | `EndTrackLength_mm` | Float | Total photon track length when it ended. |
| 34 | `AcrylicAbsOutOfRange` | Int | 1 if the photon traversed acrylic where absorption wavelength was outside the declared empirical absorption domain; an explicit extension or endpoint policy may apply. |
| 35 | `EmitterBeta` | Float | For Cherenkov photons, exact step-average emitter beta used by current G4Cerenkov sampling; NaN for non-Cherenkov photons. |

**`Photon.EndReason` codes**

| Code | Meaning |
|---:|---|
| 0 | Unknown / not classified |
| 1 | Detected |
| 2 | Bulk `OpAbsorption` |
| 3 | World escape |
| 4 | User safety bounce cut |
| 5 | Boundary absorption |
| 6 | Other stop |

## `SensorIncident` — 24 branches

Photon phase space at readout entrance.
One row per selected entrance crossing from coupling into a readout element. In `ConcentratorInput`, phase-space detail is the intended ROBAST/concentrator source.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `EventID` | Int | Event identifier. |
| 2 | `PhotonID` | Int | Optical photon ID; join to Photon on (RunID,EventID,PhotonID). |
| 3 | `ChannelID` | Int | Readout channel; resolve with ChannelMap. |
| 4 | `Wavelength_nm` | Float | Incident photon wavelength. |
| 5 | `Time_ns` | Float | Readout entrance time. |
| 6 | `LocalX_mm` | Float | Legacy local x relative to readout-volume centre. |
| 7 | `LocalY_mm` | Float | Legacy local y relative to readout-volume centre. |
| 8 | `LocalZ_mm` | Float | Legacy local z relative to readout-volume centre. |
| 9 | `LocalDirX` | Float | Incident direction x in channel local frame; full only in phaseSpace detail. |
| 10 | `LocalDirY` | Float | Incident direction y in local frame. |
| 11 | `LocalDirZ` | Float | Incident direction z in local frame; local +z is outward. |
| 12 | `LocalPolX` | Float | Incident polarization x in local frame. |
| 13 | `LocalPolY` | Float | Incident polarization y in local frame. |
| 14 | `LocalPolZ` | Float | Incident polarization z in local frame. |
| 15 | `IncidentAngle_deg` | Float | Signed/local incidence angle derived from the pre-boundary incident direction. |
| 16 | `NBoundaryReflections` | Int | Photon reflection counter at readout arrival. |
| 17 | `NTIR` | Int | TIR counter at readout arrival. |
| 18 | `NTyvekReflections` | Int | Tyvek-reflection counter at readout arrival. |
| 19 | `NRayleigh` | Int | Rayleigh counter at readout arrival. |
| 20 | `EntranceX_mm` | Float | Entrance-plane local x. |
| 21 | `EntranceY_mm` | Float | Entrance-plane local y. |
| 22 | `EntranceZ_mm` | Float | Entrance-plane local z; nominal entrance is z=0. |
| 23 | `Forward` | Int | 1 when LocalDirZ>0, i.e. ray points outward into the readout element. |

## `CoreBoundaryCrossing` — 22 branches

Optional first-successful-exit optical phase space at the Core logical-volume boundary. The tree is present but empty when `/output/coreBoundaryMode off`. In `firstExit`, one photon contributes at most one row, and reflections/TIR do not count as exits. Position, direction and polarization are the post-boundary state. Polarization is NaN in `basic` detail and stored in `phaseSpace`.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `EventID` | Int | Event identifier. |
| 2 | `PhotonID` | Int | Optical photon TrackID. |
| 3 | `ParentTrackID` | Int | Immediate Geant4 parent. |
| 4 | `SourceTrackID` | Int | Propagated non-optical source track. |
| 5 | `RootPrimaryTrackID` | Int | Resolved root primary for the optical photon. |
| 6 | `OriginVolumeID` | Int | Photon birth volume category. |
| 7 | `FaceID` | Int | `0=+X,1=-X,2=+Z,3=-Z,4=+Y,5=-Y`. |
| 8 | `X_mm` | Float | Global boundary x. |
| 9 | `Y_mm` | Float | Global boundary y. |
| 10 | `Z_mm` | Float | Global boundary z. |
| 11 | `LocalU_mm` | Float | Face-local U coordinate, consistent with readout-face axes. |
| 12 | `LocalV_mm` | Float | Face-local V coordinate. |
| 13 | `DirX` | Float | Post-boundary global direction x. |
| 14 | `DirY` | Float | Post-boundary global direction y. |
| 15 | `DirZ` | Float | Post-boundary global direction z. |
| 16 | `CosThetaNormal` | Float | Dot(post direction, outward face normal). |
| 17 | `Wavelength_nm` | Float | Photon wavelength at exit. |
| 18 | `Time_ns` | Float | Global first-exit time. |
| 19 | `PolX` | Float | Post-boundary global polarization x; NaN in basic detail. |
| 20 | `PolY` | Float | Post-boundary global polarization y; NaN in basic detail. |
| 21 | `PolZ` | Float | Post-boundary global polarization z; NaN in basic detail. |

## `DetectedSignal` — 5 branches

Detected optical photons in DirectSiPM mode.
Only meaningful in `DirectSiPM`; no fake detections are generated in `ConcentratorInput`.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `EventID` | Int | Event identifier. |
| 2 | `PhotonID` | Int | Detected photon ID; join to Photon. |
| 3 | `ChannelID` | Int | Detected readout channel. |
| 4 | `Time_ns` | Float | Detection time. |

## `ChannelEvent` — 7 branches

Sparse per-event readout-channel aggregate.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `EventID` | Int | Event identifier. |
| 2 | `ChannelID` | Int | Readout channel. |
| 3 | `NIncident` | Int | Number of incident photons in this event/channel aggregate. |
| 4 | `NDetected` | Int | Number of detections. |
| 5 | `FirstIncidentTime_ns` | Float | Earliest incident time for this event/channel. |
| 6 | `FirstDetectedTime_ns` | Float | Earliest detection time. |

## `RunConfig` — 160 branches

Run configuration and provenance. `general-v1.10` preserves the complete `general-v1.9` prefix (columns 0--156) unchanged and appends build/source provenance at columns 157--159: `GitCommit`, `GitBranch`, and `GitDirty`. `general-v1.9` had appended the Core-boundary recorder metadata at 154--156; the complete `general-v1.8` prefix (0--153) and older `general-v1.6` prefix (0--98) therefore remain backward-compatible.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `GeometryVersion` | Int | Geometry configuration version. |
| 2 | `CoreX_mm` | Float | Core x size. |
| 3 | `CoreY_mm` | Float | Core y size. |
| 4 | `CoreZ_mm` | Float | Core z size. |
| 5 | `ScintX_mm` | Float | Scintillator x size. |
| 6 | `ScintY_mm` | Float | Scintillator y size. |
| 7 | `ScintZ_mm` | Float | Scintillator z size. |
| 8 | `ReadoutX_mm` | Float | Readout-element x/footprint setting. |
| 9 | `ReadoutY_mm` | Float | Readout-element y/footprint setting. |
| 10 | `ReadoutZ_mm` | Float | Readout-element thickness/depth setting. |
| 11 | `CouplingX_mm` | Float | Coupling x size. |
| 12 | `CouplingY_mm` | Float | Coupling y size. |
| 13 | `CouplingZ_mm` | Float | Coupling thickness/depth. |
| 14 | `HoleX_mm` | Float | Configured optical/readout opening x. |
| 15 | `HoleY_mm` | Float | Configured optical/readout opening y. |
| 16 | `AcrylicThickness_mm` | Float | Legacy compatibility field; value equals generic `OuterShellThickness_mm` in this schema. |
| 17 | `TyvekThickness_mm` | Float | Configured Tyvek thickness. |
| 18 | `ScintGap_mm` | Float | Configured scintillator-to-detector gap. |
| 19 | `SourceDistance_mm` | Float | Configured source-plane/source distance. |
| 20 | `TyvekConfig` | String | Tyvek placement/configuration string. |
| 21 | `CoreMaterial` | String | Core material label. |
| 22 | `UseAcrylicShell` | Int | Legacy compatibility field; value equals generic `UseOuterShell`. |
| 23 | `SensorMode` | String | DirectSiPM or ConcentratorInput. |
| 24 | `ScintThreshold_MeV` | Float | Trigger scintillator energy threshold. |
| 25 | `CoincidenceWindow_ns` | Float | Top-bottom trigger coincidence window. |
| 26 | `IntegrationWindow_ns` | Float | Scintillator deposit integration window. |
| 27 | `OutputEventMode` | String | thinAll, triggered or fullAll. |
| 28 | `IncidentDetail` | String | auto, minimal or phaseSpace. |
| 29 | `PhotonBirthDetail` | String | basic or full. |
| 30 | `BoundaryDebug` | Int | 1 when debug optical-interaction logging is enabled. |
| 31 | `MaxOpticalBounces` | Int | Safety bounce limit; 0 means unlimited. |
| 32 | `ChannelCount` | Int | Actual generated readout-channel count. |
| 33 | `CompressionLevel` | Int | ROOT output compression setting. |
| 34 | `AcrylicRindexDataset` | String | Legacy PMMA-profile dataset field retained for v1.6 compatibility. |
| 35 | `AcrylicAbsDataset` | String | Legacy PMMA-profile dataset field retained for v1.6 compatibility. |
| 36 | `TyvekDataset` | String | Dataset/file identifier recorded for provenance. |
| 37 | `GeneratorMode` | String | Gun/CRY/etc generator mode. |
| 38 | `CRYTimeMode` | String | CRY time handling mode. |
| 39 | `CRYConfig` | String | Resolved CRY configuration string. |
| 40 | `Geant4Version` | String | Geant4 version string. |
| 41 | `PhysicsList` | String | Physics-list label. |
| 42 | `RNGSeed` | Double | Base random seed. |
| 43 | `AirRindexDataset` | String | Dataset/file identifier recorded for provenance. |
| 44 | `WaterRindexDataset` | String | Dataset/file identifier recorded for provenance. |
| 45 | `WaterAbsDataset` | String | Dataset/file identifier recorded for provenance. |
| 46 | `WaterRayleighDataset` | String | Dataset/file identifier recorded for provenance. |
| 47 | `CouplingRindexDataset` | String | Dataset/file identifier recorded for provenance. |
| 48 | `SensorPDEDataset` | String | Dataset/file identifier recorded for provenance. |
| 49 | `GunParticle` | String | Generic gun particle name. |
| 50 | `GunEnergy_MeV` | Float | Resolved gun kinetic energy. |
| 51 | `GunX_mm` | Float | Resolved gun x. |
| 52 | `GunY_mm` | Float | Resolved gun y. |
| 53 | `GunZ_mm` | Float | Resolved gun z. |
| 54 | `GunDirX` | Float | Gun direction x. |
| 55 | `GunDirY` | Float | Gun direction y. |
| 56 | `GunDirZ` | Float | Gun direction z. |
| 57 | `GunTime_ns` | Float | Gun time. |
| 58 | `GunUseAutoY` | Int | Whether automatic gun y placement was enabled. |
| 59 | `GunPolX` | Float | Gun optical polarization x. |
| 60 | `GunPolY` | Float | Gun optical polarization y. |
| 61 | `GunPolZ` | Float | Gun optical polarization z. |
| 62 | `SchemaVersion` | String | ROOT schema semantic version. |
| 63 | `NThreads` | Int | Configured worker-thread count. |
| 64 | `RNGEngine` | String | Random-number engine name. |
| 65 | `CRYVersion` | String | CRY version string. |
| 66 | `DefaultCut_mm` | Float | Default Geant4 production cut. |
| 67 | `CerenkovEnabled` | Int | Runtime activation state of Cherenkov process. |
| 68 | `ScintillationEnabled` | Int | Runtime activation state of scintillation. |
| 69 | `AbsorptionEnabled` | Int | Runtime activation state of optical absorption. |
| 70 | `RayleighEnabled` | Int | Runtime activation state of Rayleigh scattering. |
| 71 | `MieEnabled` | Int | Runtime activation state of Mie scattering. |
| 72 | `BoundaryEnabled` | Int | Runtime activation state of optical boundary process. |
| 73 | `WLSEnabled` | Int | Runtime activation state of WLS. |
| 74 | `WLS2Enabled` | Int | Runtime activation state of WLS2. |
| 75 | `AcrylicAbsValidMin_nm` | Float | Minimum documented/tabulated validity wavelength for this dataset. |
| 76 | `AcrylicAbsValidMax_nm` | Float | Maximum documented/tabulated validity wavelength for this dataset. |
| 77 | `AcrylicAbsOutOfRangePolicy` | String | Policy outside measured acrylic absorption range. |
| 78 | `CouplingAbsDataset` | String | Dataset/file identifier recorded for provenance. |
| 79 | `CouplingMaterialModel` | String | Coupling bulk-material model label. |
| 80 | `SensorModel` | String | Sensor model label. |
| 81 | `ProjectRevision` | String | Compact build revision (`<git-short-hash>` with `-dirty` suffix when applicable; `unknown` outside a Git checkout). |
| 82 | `AcrylicRindexValidMin_nm` | Float | Minimum documented/tabulated validity wavelength for this dataset. |
| 83 | `AcrylicRindexValidMax_nm` | Float | Maximum documented/tabulated validity wavelength for this dataset. |
| 84 | `CouplingRindexValidMin_nm` | Float | Minimum documented/tabulated validity wavelength for this dataset. |
| 85 | `CouplingRindexValidMax_nm` | Float | Maximum documented/tabulated validity wavelength for this dataset. |
| 86 | `CouplingAbsValidMin_nm` | Float | Minimum documented/tabulated validity wavelength for this dataset. |
| 87 | `CouplingAbsValidMax_nm` | Float | Maximum documented/tabulated validity wavelength for this dataset. |
| 88 | `WaterRindexValidMin_nm` | Float | Minimum documented/tabulated validity wavelength for this dataset. |
| 89 | `WaterRindexValidMax_nm` | Float | Maximum documented/tabulated validity wavelength for this dataset. |
| 90 | `WaterAbsValidMin_nm` | Float | Minimum documented/tabulated validity wavelength for this dataset. |
| 91 | `WaterAbsValidMax_nm` | Float | Maximum documented/tabulated validity wavelength for this dataset. |
| 92 | `WaterRayleighValidMin_nm` | Float | Minimum documented/tabulated validity wavelength for this dataset. |
| 93 | `WaterRayleighValidMax_nm` | Float | Maximum documented/tabulated validity wavelength for this dataset. |
| 94 | `WaterRayleighModelStatus` | String | Status/provenance description for water Rayleigh model. |
| 95 | `OpticalDataOutOfRangePolicy` | String | General optical-table out-of-range behavior. |
| 96 | `CoreEntryExitSemantic` | String | Definition of Core entry/exit fields. |
| 97 | `SensorPDEOperatingPoint` | String | Operating point associated with PDE table. |
| 98 | `OpticalDataFingerprint_FNV1a64` | String | Combined FNV-1a-64 fingerprint of exact optical data files. |
| 99 | `MaterialProfileRegistryVersion` | String | Curated bulk-profile registry contract; current value `curated-bulk-v2_absMeasuredRange`. |
| 100 | `LegacyRunConfigCompatibility` | String | Explicit compatibility statement for the preserved v1.6 prefix. |
| 101 | `CoreProfile` | String | Active Core `MaterialProfile` ID. |
| 102 | `CoreProfileG4Material` | String | Actual Geant4 material name used by the Core profile. |
| 103 | `CoreProfileRindexDataset` | String | Active Core RINDEX dataset. |
| 104 | `CoreProfileAbsDataset` | String | Active Core ABSLENGTH dataset. |
| 105 | `CoreProfileRayleighDataset` | String | Active Core RAYLEIGH dataset; empty when not configured. |
| 106 | `CoreProfileRindexValidMin_nm` | Float | Minimum tabulated wavelength for active Core RINDEX dataset. |
| 107 | `CoreProfileRindexValidMax_nm` | Float | Maximum tabulated wavelength for active Core RINDEX dataset. |
| 108 | `CoreProfileAbsValidMin_nm` | Float | Minimum tabulated wavelength for active Core absorption dataset. |
| 109 | `CoreProfileAbsValidMax_nm` | Float | Maximum tabulated wavelength for active Core absorption dataset. |
| 110 | `CoreProfileRayleighValidMin_nm` | Float | Minimum tabulated wavelength for active Core Rayleigh dataset, or -1 if absent. |
| 111 | `CoreProfileRayleighValidMax_nm` | Float | Maximum tabulated wavelength for active Core Rayleigh dataset, or -1 if absent. |
| 112 | `CoreProfileRindexStatus` | String | Provenance/model status of Core RINDEX. |
| 113 | `CoreProfileAbsStatus` | String | Provenance/model status of Core ABSLENGTH. |
| 114 | `CoreProfileRayleighStatus` | String | Provenance/model status of Core Rayleigh model. |
| 115 | `CoreProfileOutOfRangePolicy` | String | Geant4/profile policy outside Core optical-table ranges. |
| 116 | `UseOuterShell` | Int | 1 when the generic OuterShell geometry role is enabled. |
| 117 | `OuterShellThickness_mm` | Float | Generic OuterShell thickness. |
| 118 | `OuterShellProfile` | String | Selected OuterShell `MaterialProfile` ID. Recorded even when the shell is disabled so the configured selection is reproducible. |
| 119 | `OuterShellProfileG4Material` | String | Actual Geant4 material name selected for OuterShell. |
| 120 | `OuterShellProfileRindexDataset` | String | OuterShell RINDEX dataset. |
| 121 | `OuterShellProfileAbsDataset` | String | OuterShell ABSLENGTH dataset. |
| 122 | `OuterShellProfileRayleighDataset` | String | OuterShell RAYLEIGH dataset; empty when not configured. |
| 123 | `OuterShellProfileRindexValidMin_nm` | Float | Minimum tabulated wavelength for OuterShell RINDEX. |
| 124 | `OuterShellProfileRindexValidMax_nm` | Float | Maximum tabulated wavelength for OuterShell RINDEX. |
| 125 | `OuterShellProfileAbsValidMin_nm` | Float | Minimum tabulated wavelength for OuterShell absorption. |
| 126 | `OuterShellProfileAbsValidMax_nm` | Float | Maximum tabulated wavelength for OuterShell absorption. |
| 127 | `OuterShellProfileRayleighValidMin_nm` | Float | Minimum tabulated wavelength for OuterShell Rayleigh, or -1 if absent. |
| 128 | `OuterShellProfileRayleighValidMax_nm` | Float | Maximum tabulated wavelength for OuterShell Rayleigh, or -1 if absent. |
| 129 | `OuterShellProfileRindexStatus` | String | Provenance/model status of OuterShell RINDEX. |
| 130 | `OuterShellProfileAbsStatus` | String | Provenance/model status of OuterShell ABSLENGTH. |
| 131 | `OuterShellProfileRayleighStatus` | String | Provenance/model status of OuterShell Rayleigh. |
| 132 | `OuterShellProfileOutOfRangePolicy` | String | OuterShell optical-table out-of-range policy. |
| 133 | `CouplingProfile` | String | Active Coupling `MaterialProfile` ID. |
| 134 | `CouplingProfileG4Material` | String | Actual Geant4 material name used for Coupling. |
| 135 | `CouplingProfileRindexDataset` | String | Active Coupling RINDEX dataset. |
| 136 | `CouplingProfileAbsDataset` | String | Active Coupling ABSLENGTH dataset. |
| 137 | `CouplingProfileRayleighDataset` | String | Active Coupling RAYLEIGH dataset; empty when not configured. |
| 138 | `CouplingProfileRindexValidMin_nm` | Float | Minimum tabulated wavelength for Coupling RINDEX. |
| 139 | `CouplingProfileRindexValidMax_nm` | Float | Maximum tabulated wavelength for Coupling RINDEX. |
| 140 | `CouplingProfileAbsValidMin_nm` | Float | Minimum tabulated wavelength for Coupling absorption. |
| 141 | `CouplingProfileAbsValidMax_nm` | Float | Maximum tabulated wavelength for Coupling absorption. |
| 142 | `CouplingProfileRayleighValidMin_nm` | Float | Minimum tabulated wavelength for Coupling Rayleigh, or -1 if absent. |
| 143 | `CouplingProfileRayleighValidMax_nm` | Float | Maximum tabulated wavelength for Coupling Rayleigh, or -1 if absent. |
| 144 | `CouplingProfileRindexStatus` | String | Provenance/model status of Coupling RINDEX. |
| 145 | `CouplingProfileAbsStatus` | String | Provenance/model status of Coupling ABSLENGTH. |
| 146 | `CouplingProfileRayleighStatus` | String | Provenance/model status of Coupling Rayleigh. |
| 147 | `CouplingProfileOutOfRangePolicy` | String | Coupling optical-table out-of-range policy. |
| 148 | `CoreProfileAbsMeasuredMin_nm` | Float | Minimum wavelength of explicitly declared empirical ABS validity for the active Core profile; `-1` if undeclared. |
| 149 | `CoreProfileAbsMeasuredMax_nm` | Float | Maximum wavelength of explicitly declared empirical ABS validity for the active Core profile; `-1` if undeclared. |
| 150 | `OuterShellProfileAbsMeasuredMin_nm` | Float | Minimum wavelength of explicitly declared empirical ABS validity for the active OuterShell profile; `-1` if undeclared. |
| 151 | `OuterShellProfileAbsMeasuredMax_nm` | Float | Maximum wavelength of explicitly declared empirical ABS validity for the active OuterShell profile; `-1` if undeclared. |
| 152 | `CouplingProfileAbsMeasuredMin_nm` | Float | Minimum wavelength of explicitly declared empirical ABS validity for the active Coupling profile; `-1` if undeclared. |
| 153 | `CouplingProfileAbsMeasuredMax_nm` | Float | Maximum wavelength of explicitly declared empirical ABS validity for the active Coupling profile; `-1` if undeclared. |
| 154 | `CoreBoundaryMode` | String | `off` (default) or `firstExit`. |
| 155 | `CoreBoundaryDetail` | String | `basic` or `phaseSpace`. |
| 156 | `CoreBoundarySemantic` | String | Explicit first-successful-exit/post-boundary semantic and face convention. |
| 157 | `GitCommit` | String | Full Git commit hash embedded at CMake configure time, or `unknown` when unavailable. |
| 158 | `GitBranch` | String | Git branch embedded at CMake configure time, or `unknown` when unavailable. |
| 159 | `GitDirty` | Int | 1 when the source working tree was dirty at CMake configure time; 0 otherwise. |

## `ChannelMap` — 19 branches

Readout channel geometry.
Canonical actual readout geometry. Do not infer face/channel layout from arithmetic on ChannelID.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `ChannelID` | Int | Unique readout channel ID within the run. |
| 2 | `FaceID` | Int | Face code: 0=+X, 1=-X, 2=+Z, 3=-Z, 4=+Y, 5=-Y. |
| 3 | `FaceName` | String | Human-readable face label. |
| 4 | `GridU` | Int | Grid index U. |
| 5 | `GridV` | Int | Grid index V. |
| 6 | `LocalU_mm` | Float | Element centre coordinate on face-local U axis. |
| 7 | `LocalV_mm` | Float | Element centre coordinate on face-local V axis. |
| 8 | `CenterX_mm` | Float | Global readout-volume centre x. |
| 9 | `CenterY_mm` | Float | Global centre y. |
| 10 | `CenterZ_mm` | Float | Global centre z. |
| 11 | `NormalX` | Float | Global outward unit-normal x. |
| 12 | `NormalY` | Float | Global outward unit-normal y. |
| 13 | `NormalZ` | Float | Global outward unit-normal z. |
| 14 | `SizeU_mm` | Float | Element footprint along local U. |
| 15 | `SizeV_mm` | Float | Element footprint along local V. |
| 16 | `EntranceCenterX_mm` | Float | Global entrance-plane centre x. |
| 17 | `EntranceCenterY_mm` | Float | Global entrance-plane centre y. |
| 18 | `EntranceCenterZ_mm` | Float | Global entrance-plane centre z. |

**Face IDs**: `0=+X`, `1=-X`, `2=+Z`, `3=-Z`, `4=+Y`, `5=-Y`.

## `VolumeMap` — 4 branches

Stable volume category map for this schema.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `VolumeID` | Int | Stable/category volume ID used by output rows. |
| 2 | `VolumeName` | String | Logical category name. |
| 3 | `MaterialName` | String | Material assigned in this run. |

## `RunStatistics` — 6 branches

Per-worker run counters; sum rows for run totals.
In MT there is one row per worker. **Sum worker rows** for run totals.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `WorkerThreadID` | Int | Event-processing worker ID. Sum worker rows for run totals. |
| 2 | `GeneratedEvents` | Int | Events generated by this worker. |
| 3 | `TriggeredEvents` | Int | Events passing trigger on this worker. |
| 4 | `DetailedEvents` | Int | Events for which detailed tier was written. |
| 5 | `CRYSimulatedTime_s` | Double | CRY simulated-time contribution of this worker. |

## `BoundaryDebug` — 24 branches

Debug-only optical interaction log. Ordinary rows are optical boundary crossings.
For Stage-3D Rayleigh validation, the first Rayleigh scatter can also appear with
`BoundaryStatus=-1001` and `BoundaryStatusName=FirstRayleighScatter`, carrying
the same PRE/POST direction and polarization fields. Normally OFF in production.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `EventID` | Int | Event identifier. |
| 2 | `PhotonID` | Int | Optical photon ID. |
| 3 | `StepNumber` | Int | Geant4 optical track step number. |
| 4 | `PreVolumeID` | Int | Volume before the logged interaction. |
| 5 | `PostVolumeID` | Int | Volume after the logged interaction. |
| 6 | `BoundaryStatus` | Int | Integer G4OpBoundaryProcess status. |
| 7 | `Time_ns` | Float | Boundary-crossing time. |
| 8 | `X_mm` | Float | Global boundary x. |
| 9 | `Y_mm` | Float | Global boundary y. |
| 10 | `Z_mm` | Float | Global boundary z. |
| 11 | `PreDirX` | Float | Pre-boundary direction x. |
| 12 | `PreDirY` | Float | Pre-boundary direction y. |
| 13 | `PreDirZ` | Float | Pre-boundary direction z. |
| 14 | `PostDirX` | Float | Post-boundary direction x. |
| 15 | `PostDirY` | Float | Post-boundary direction y. |
| 16 | `PostDirZ` | Float | Post-boundary direction z. |
| 17 | `PrePolX` | Float | Pre-boundary polarization x. |
| 18 | `PrePolY` | Float | Pre-boundary polarization y. |
| 19 | `PrePolZ` | Float | Pre-boundary polarization z. |
| 20 | `PostPolX` | Float | Post-boundary polarization x. |
| 21 | `PostPolY` | Float | Post-boundary polarization y. |
| 22 | `PostPolZ` | Float | Post-boundary polarization z. |
| 23 | `BoundaryStatusName` | String | Readable G4OpBoundaryProcess status string. |

## `ReadoutPlaneConfig` — 15 branches

Configured readout layout per detector face.
Requested per-face readout policy; combine with ChannelMap to know what was requested and what geometry was actually generated.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `FaceID` | Int | Face code. |
| 2 | `FaceName` | String | Readable face name. |
| 3 | `Enabled` | Int | Whether this face readout is enabled. |
| 4 | `AutoFill` | Int | 1 for auto-fill policy, 0 for requested fixed grid. |
| 5 | `NU` | Int | Requested fixed-grid U count. |
| 6 | `NV` | Int | Requested fixed-grid V count. |
| 7 | `PitchU_mm` | Float | Requested U pitch. |
| 8 | `PitchV_mm` | Float | Requested V pitch. |
| 9 | `OffsetU_mm` | Float | Requested U offset. |
| 10 | `OffsetV_mm` | Float | Requested V offset. |
| 11 | `GapU_mm` | Float | Auto-fill U gap. |
| 12 | `GapV_mm` | Float | Auto-fill V gap. |
| 13 | `MarginU_mm` | Float | Auto-fill U margin. |
| 14 | `MarginV_mm` | Float | Auto-fill V margin. |

## `PrimaryTruth` — 12 branches

Lightweight non-optical primary truth aligned with the thin event tier.
Lightweight non-optical primaries, particularly important for unbiased `thinAll` production samples.

| # | Branch | Type | Meaning / unit |
|---:|---|---|---|
| 0 | `RunID` | Int | Run identifier. |
| 1 | `EventID` | Int | Event identifier. |
| 2 | `TrackID` | Int | Non-optical primary TrackID. |
| 3 | `PDG` | Int | PDG code. |
| 4 | `VertexKE_MeV` | Float | Primary kinetic energy at generation. |
| 5 | `VertexTime_ns` | Float | Primary generation time. |
| 6 | `VertexX_mm` | Float | Primary vertex x. |
| 7 | `VertexY_mm` | Float | Primary vertex y. |
| 8 | `VertexZ_mm` | Float | Primary vertex z. |
| 9 | `DirX` | Float | Primary direction x. |
| 10 | `DirY` | Float | Primary direction y. |
| 11 | `DirZ` | Float | Primary direction z. |

## 4. Essential joins

### Photon -> non-optical emitter

```text
Photon.(RunID,EventID,ParentTrackID)
    -> Track.(RunID,EventID,TrackID)
```

Use this for immediate-emitter and primary-emitter selections.

### Photon -> propagated source/ancestor

```text
Photon.SourceTrackID -> Track.TrackID -> Track.RootPrimaryTrackID
```

Use for primary-initiated cascade selections and future optical re-emission chains.

### Incident/detection -> photon

```text
SensorIncident.(RunID,EventID,PhotonID) -> Photon
DetectedSignal.(RunID,EventID,PhotonID) -> Photon
```

### Readout record -> geometry

```text
SensorIncident.ChannelID / DetectedSignal.ChannelID / ChannelEvent.ChannelID
    -> ChannelMap.ChannelID
```

### Volume ID -> volume/material

```text
Track.BirthVolumeID, TrackVolume.VolumeID, Photon.OriginVolumeID, ...
    -> VolumeMap.VolumeID
```

## 5. DirectSiPM versus ConcentratorInput

### DirectSiPM

- `SensorIncident` records photons at the effective sensor entrance before the PDE decision.
- `DetectedSignal` records successful boundary detections.
- `ChannelEvent` aggregates incidents/detections per event/channel.

### ConcentratorInput

- `SensorIncident` is the important output.
- automatic `incidentDetail` is `phaseSpace`; local direction/polarization/wavelength/time are retained.
- `EntranceZ_mm` is defined as 0 on the nominal entrance plane.
- local `+z` points outward on all six faces and valid forward entrance rays should have `Forward==1`.
- `DetectedSignal` must remain empty by construction.

## 6. Important semantic cautions

1. `CoreEntry`/`CoreExit` are crossings of the **Core logical volume**, not necessarily crossings of the entire detector envelope.
2. `Edep_MeV` and `DeltaKE_MeV` are different observables; do not use them interchangeably.
3. `Photon.EmissionAngle_deg` is only physically defined for processes that set it; non-Cherenkov photons can contain NaN.
4. `Photon.EmitterBeta` is only defined for Cherenkov photons.
5. Some high-detail photon birth or incident phase-space fields can be NaN depending on output policy.
6. `thinAll` production files intentionally omit detailed rows for rejected events; use `PrimaryTruth`, `TriggerEvent` and `RunStatistics` for unbiased denominators.
7. Do not assume a fixed 24-channel detector. `ChannelMap` is the source of truth.
8. Do not interpret PMMA/Tyvek/grease dataset names as proof that the final laboratory material has identical optical properties; they document the model used in the simulation.

## 7. Which trees are normally used for common tasks?

| Task | Main ntuples |
|---|---|
| Primary cosmic spectrum | `PrimaryTruth`, `TriggerEvent`, `RunStatistics` |
| Track/path/dE/dx | `Track`, `TrackVolume`, `VolumeMap` |
| Cherenkov generation | `Photon`, `Track`, `TrackVolume`, `VolumeMap` |
| Optical transport history | `Photon`; `BoundaryDebug` only for controlled debugging |
| Readout phase space | `SensorIncident`, `ChannelMap`, `Photon` |
| DirectSiPM detections | `DetectedSignal`, `SensorIncident`, `Photon`, `ChannelMap` |
| Per-channel occupancy/timing | `ChannelEvent`, `ChannelMap` |
| Trigger reconstruction/QA | `TriggerEvent`, `ScintDeposit` |
| Run provenance | `RunConfig`, `RunStatistics`, `ReadoutPlaneConfig`, `ChannelMap`, `VolumeMap` |
| ROBAST concentrator source | `SensorIncident` in `ConcentratorInput` phase-space mode |
