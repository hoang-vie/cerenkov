# General-v1 redesign summary

## Core goals

- Particle-general truth: no muon/electron/pion-specific simulation trees.
- Flexible readout geometry: independent per-face grids, auto-fill, arbitrary channel count.
- Two readout-response modes now: `DirectSiPM` and `ConcentratorInput`.
- Objective joins and process/volume metadata instead of analysis-specific truth labels.
- Small default production output while retaining unbiased primary/trigger denominators.

## Main schema

`TriggerEvent`, `ScintDeposit`, `PrimaryTruth`, `Track`, `TrackVolume`, `Photon`,
`SensorIncident`, `DetectedSignal`, `ChannelEvent`, `RunConfig`, `ChannelMap`,
`VolumeMap`, `RunStatistics`, `ReadoutPlaneConfig`; `BoundaryDebug` is debug-only.

`Photon.ParentTrackID` identifies the immediate Geant4 parent. `SourceTrackID`
propagates the non-optical optical-source ancestry through re-emission chains, so
analysis can distinguish direct primary emission from a primary-initiated cascade.

## Output policy

- `thinAll` (default): thin trigger/pulse/primary rows for every generated event;
  full track/photon/readout detail only for triggered events.
- `triggered`: only triggered-event data, with run-level generated denominator.
- `fullAll`: validation/debug mode for small samples.
- Full optical boundary crossing logging is OFF by default.
- Optical safety bounce limit is OFF by default (`0 = unlimited`).

## Physics/geometry changes

- Hardware-like trigger uses a scintillator integration window, Edep threshold and
  top-bottom coincidence window; no same-TrackID condition.
- CRY default timing is relative within a CRY event; absolute and zero modes remain.
- CRY latitude examples use magnetic latitude. For the Quy Nhon 2026 example the bundled IGRF-14 tilted-dipole helper gives about 4.63 deg; geographic 13.78 deg is not inserted directly into CRY.
- Tyvek baseline is an effective diffuse-dominated UNIFIED `dielectric_metal/ground`
  surface with wavelength-dependent reflectivity.
- PMMA attenuation baseline uses the published SNO UV-transmitting Sample-A Table-II
  data. The broad project RINDEX remains the production RINDEX because the exact
  published SNO RTT' table included here ends at 405 nm.
- Scintillation and WLS/WLS2 are schema-ready but OFF by default until calibrated
  material/fibre data and geometry are supplied.

## Flexible readout

Each of the six detector faces may be independently enabled and configured as a
fixed `NU x NV` grid or auto-filled from element footprint, gap and margins.
`ChannelMap` is the source of truth for face/position/normal; analysis never infers
a face from `ChannelID/4` or assumes 24 sensors.

## Analysis

`analysis/analyze_general.C` accepts a selected PDG and source mode:

- `primaryEmitter`: photon emitted directly by a selected primary particle;
- `emitter`: photon emitted directly by any selected-species track;
- `ancestor`: photon belongs to a cascade initiated by a selected primary species.

The same report can therefore be produced for muons, electrons, pions, protons,
kaons available through the Geant4 gun, or other PDG selections without changing
the simulation schema.

## Validation/status

Static checks performed in the development container:

- ntuple creation-order / fill-column index and type check: PASS;
- independent Cherenkov/TIR/attenuation script: PASS;
- source/header delimiter balance: PASS;
- hard-coded 24-channel / old TIR label / 17.38 latitude / 100-bounce audit: PASS
  outside the retained legacy analysis reference.

The development container does not contain a Geant4 package configuration or ROOT,
so this revision has NOT been compiled, linked, or run against Geant4/ROOT here.

## Hotfix 1 — build compatibility

- CRY 1.7 `CRYParticle` getters such as `t()` are non-const. Changed the CRY event minimum-time loop from `const auto*` to `auto*` so `GeneratePrimaries()` compiles against the bundled CRY 1.7 header.
- Replaced deprecated `G4AccumulableManager::RegisterAccumulable(...)` calls with `Register(...)` for current Geant4 releases.
- Re-ran `tools/check_ntuple_schema.py`: all 15 ntuples pass creation-order and fill-column type/index checks.

## Hotfix 2 — Geant4 UI range parameter registration
- Fixed `/output/maxOpticalBounces`: register parameter name before `SetRange`.
- Fixed `/output/compressionLevel`: register parameter name before `SetRange`.
- These commands now accept values used by the default macro (`0` and `6`, respectively).

## Hotfix 3 — logic/physics audit

- Added `TrackVolume.DeltaKE_MeV` so local Edep and track kinetic-energy loss are not conflated in stopping-power validation.
- `ScintDeposit` now stores raw sparse per-step deposits, allowing offline re-clustering of integration/threshold/coincidence.
- Tightened extraction to true transmission boundary statuses and detection to Coupling->Readout only.
- `SensorIncident` phase-space uses the incident (pre-boundary) direction/polarization.
- `ConcentratorInput` readout bulk is air rather than silicon; the optical border remains the terminal virtual aperture.
- AutoFill fixed-grid count/pitch metadata is neutralized; realized geometry is represented by `ChannelMap`.
- Added non-negative generic gun energy UI range.
- Added CRY simulated-time run normalization and magnetic-latitude helper.
- Added Frank-Tamm/PDE validation macros and ROOT output inspectors.

## Hotfix3 final audit

- Explicitly clean `G4LogicalBorderSurface` and `G4LogicalSkinSurface` tables before sequential Idle-state `ReinitializeGeometry()` for compatibility with Geant4 releases before the 11.4 kernel fix.
- Fixed ROOT `inspect_output.C` branch binding for RunStatistics integer columns (`int`, not `Long64_t`).
- `uniformGrid` now resets per-plane offsets to zero so a uniform command cannot inherit stale custom offsets.
- Added CRY ConcentratorInput and CRY DirectSiPM validation macros. The earlier geometry-rebuild macro was later replaced in hotfix4 because multi-event dynamic geometry is not a valid MT validation path.
- Re-audited photon extraction/detection semantics, ancestry, trigger pulse selection, geometry-version refresh and ntuple fill indices.


## Hotfix 4 — MT output and validation audit

- Fixed the end-of-run ROOT warning `Analysis_W001: Ntuple 11 does not exist`. `RunStatistics` is no longer filled from the MT master analysis manager. Each event-processing worker writes one local statistics row; analysis sums those rows for global totals.
- Expanded `RunStatistics` with `WorkerThreadID` and updated `inspect_output.C` to sum worker rows.
- Strengthened `inspect_output.C`: required-tree checks, exactly one `RunConfig` row, generated/channel expectations, output-mode row-count checks, PrimaryTruth multiplicity check, readout/detection consistency.
- Replaced the invalid multi-event MT geometry-rebuild validation with `validation_flexible_layout.mac`, which validates an asymmetric 12-channel layout in one geometry construction. Dynamic geometry remains restricted at Idle in an MT application.
- Added `validation_output_mt.mac` as the first runtime test for ROOT ntuple merging.
- Added `validate_dedx.C` and `validate_muon_decay.C`; tightened `validate_frank_tamm.C` so future scintillation/WLS photons are not counted as Cherenkov yield.
- Renamed ambiguous TriggerEvent multiplicities to `NPrimaryNonOptical`, `NNonOpticalTracks`, and `NOpticalTracks`.
- Added static `check_analysis_schema.py` and `check_macros.py` so analysis branch drift and unregistered project UI commands are caught before runtime.
- Schema version bumped to `general-v1.4`.
- `RunConfig.NThreads` is now read from the MT master run manager when metadata is written on worker 0; a worker-local `G4RunManager` must not be used as the source of the configured worker count.
- `Write()`/`CloseFile()` are executed unconditionally at end of run so metadata/statistics from a worker that processed zero events are still eligible for ROOT ntuple merging.
- `analyze_general.C` now distinguishes incident attempts from unique photons reaching the readout; transport/sensor efficiencies use unique-photon denominators while occupancy/incident spectra retain attempt counts.

## Hotfix 6 — physics audit 1

- Replaced SiO2 coupling alias with explicit generic PDMS optical coupling.
- Updated default water RINDEX/absorption provenance; PMMA broad RINDEX and
  measured-valid absorption ranges are now explicit.
- Hardened optical table parsing and added exact-data fingerprint metadata.
- Added `TrackingAction` so generic non-optical truth is not step-1 dependent.
- Added entrance-plane coordinates and signed-forward validation to
  `SensorIncident`; appended entrance centres to `ChannelMap`.
- Added photon track-length diagnostics and detailed pre/post boundary vectors.
- Removed duplicate Tyvek subtraction-solid construction.
- Added deterministic Geant4 validation macros/ROOT validators for six-face
  transforms, absorption, Rayleigh, Fresnel/Snell/TIR and Tyvek reflectivity.
- Schema version bumped to `general-v1.5`.
- This revision is an audit baseline, not yet `physics_validated_v1`; runtime
  validation on Geant4+ROOT is required.


## Hotfix6 audit 4 — exact Cherenkov emitter beta

- Schema bumped to `general-v1.6`.
- Appended `Photon.EmitterBeta` without renumbering existing Photon columns.
- For Cerenkov secondaries, stores `0.5*(beta_pre+beta_post)`, matching the beta convention used by `G4Cerenkov::PostStepDoIt` for photon-angle sampling.
- `validate_frank_tamm.C` now evaluates the primary-muon angle residual from the stored per-photon beta rather than reconstructing beta from vertex kinetic energy.
- Added a 1000-event full-chain monoenergetic-muon macro for exercising `analyze_general.C`.

## Framework metadata/API normalization — general-v1.10

- Added canonical `/detector/coreProfile Water_Default|PMMA_SNO_A`; retained `/detector/coreMaterial` as a convenience/legacy alias.
- Added CMake-configured source/build provenance and appended `RunConfig.GitCommit`, `GitBranch`, and `GitDirty`; `ProjectRevision` is now derived automatically from Git rather than a manually edited hotfix string.
- Preserved the complete general-v1.9 RunConfig prefix (columns 0--156) and extended the static schema checker with a v1.9-prefix hash guard.
- Clarified that `/detector/tyvekHoleSize` is a legacy 3-vector interface whose X/Y values define the opening while Z is ignored by geometry.
- Updated user-facing example macros to prefer canonical profile/generic geometry commands while retaining legacy commands in regression macros.
- Synchronized the framework/schema documentation with the current PMMA transport policy and added `FRAMEWORK_STATUS.md` to define scope, configuration/metadata contracts, and completion criteria.

### Framework contract audit (post-v1.10 runtime gates)
- Corrected the static macro checker summary to distinguish 51 executable UI commands from 6 UI directories.
- Added `tools/check_config_metadata_contract.py`: registering a new project command without an explicit metadata/reproducibility contract now fails the static suite.
- Added `FRAMEWORK_CONFIG_AUDIT.md` with the complete command-to-metadata map, legacy/expert exceptions, intentionally fixed framework assumptions, and target-machine runtime evidence.
- No Geant4 physics, optical tables, geometry semantics, ROOT schema, or event data model changed in this audit step.

### User-facing macro organization / F5
- Added a top-level `README.md` quick start and `macros/README.md` intent catalog so a new user can distinguish ordinary examples from validation, study, legacy, and visualization macros.
- Added `example_gun_direct_sipm.mac` as a deterministic generic-gun starting point alongside the existing CRY `DirectSiPM` and `ConcentratorInput` examples.
- Added `tools/check_macro_catalog.py`; every flat `macros/*.mac` file must have a known intent category, user-facing examples must be documented, use canonical rather than legacy detector commands, and contain the minimum run/source/readout/output structure.
- Updated framework/runbook documentation to direct ordinary studies to `example_*.mac` rather than adapting controlled validation macros.
- No Geant4 physics, optical data, geometry semantics, ROOT schema, or event data model changed in this step.

## Framework V1 freeze harness
- Added `tools/run_framework_static_checks.sh` as the one-command static contract gate.
- Added `tools/run_framework_freeze_regression.sh` for clean-build release regression and exact Git provenance checking.
- Added `FRAMEWORK_FREEZE_CHECKLIST.md` defining the framework freeze/tag procedure and post-freeze change policy.
- Updated `validation_final_integrated_freeze_1000.mac` to canonical `coreProfile` / `readoutSize` commands and explicit CoreBoundary OFF/basic metadata. No physics model or ROOT schema change.
