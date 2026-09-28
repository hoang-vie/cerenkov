# Cherenkov simulation framework status — general-v1.10

## Purpose

This repository is a **simulation framework for detector studies**, not a frozen detector-optimization result. The framework goal is to provide a validated Geant4 physics baseline plus a configurable detector/readout/generator/output surface so another user can change a study configuration without editing C++ for ordinary scans. Optimization of a particular sensor placement, face subset, or concentrator is a downstream use case, not a framework-completion criterion.

## Supported topology

The current geometry family is intentionally finite: a rectangular Core; optional uniform outer shell; optional Tyvek reflector; top/bottom scintillator trigger paddles; rectangular coupling/readout elements on any of six axis-aligned faces. Sizes, materials/profiles and per-face readout layouts are configurable, but arbitrary CAD/cylindrical/spherical geometries are outside this framework version.

## Canonical configuration surface

Prefer these commands in new macros:

- Core: `/detector/coreProfile`, `/detector/coreSize`.
- Outer shell: `/detector/useOuterShell`, `/detector/outerShellProfile`, `/detector/outerShellThickness`.
- Tyvek: `/detector/tyvekPosition`, `/detector/setTyvekThickness`, `/detector/tyvekHoleSize`. The legacy hole command is a 3-vector; only X/Y define the opening and Z is ignored.
- Coupling/readout: `/detector/couplingProfile`, `/detector/couplingSize`, `/detector/readoutSize`, `/detector/sensorMode`.
- Per-face layout: `/readout/uniformGrid`, `/readout/planeGrid`, `/readout/planeAuto`.
- Trigger: `/trigger/scintThreshold`, `/trigger/integrationWindow`, `/trigger/coincidenceWindow`.
- Generator: `/generator/*` for the generic gun and `/CRY/*` for CRY.
- Output: `/output/eventMode`, `/output/incidentDetail`, `/output/photonBirthDetail`, `/output/coreBoundaryMode`, `/output/coreBoundaryDetail`, `/output/boundaryDebug`, `/output/maxOpticalBounces`, `/output/fileName`, `/output/compressionLevel`.

Legacy aliases (`coreMaterial`, `sipmSize`, `sipmPitch`, `setAcrylicThickness`, `useAcrylicShell`, `setSiPMOffsetY`) remain supported for old macros/regression but should not be used as the primary public interface in new examples.

## Configuration -> metadata contract

The framework records the resolved configuration rather than relying on macro filenames:

| Configuration area | Authoritative output metadata | Status |
|---|---|---|
| Core/scintillator/readout/coupling dimensions | `RunConfig` | GREEN |
| Core/OuterShell/Coupling curated material profiles and optical datasets | `RunConfig` profile-aware fields | GREEN |
| Tyvek placement/thickness/opening and dataset | `RunConfig` | GREEN (hole Z is legacy/ignored) |
| Per-face enabled/grid/auto-fill/pitch/offset/gap/margin policy | `ReadoutPlaneConfig` | GREEN |
| Actual generated readout elements/positions/normals | `ChannelMap` | GREEN |
| Trigger thresholds/windows | `RunConfig` | GREEN |
| Gun or CRY mode/configuration | `RunConfig` | GREEN |
| RNG seed, engine, actual thread count | `RunConfig` | GREEN |
| Geant4 version, physics-list label, production cut, optical-process activation | `RunConfig` | GREEN |
| Output detail/safety/compression policy | `RunConfig` | GREEN |
| Core first-exit diagnostic policy | `RunConfig` + `CoreBoundaryCrossing` | GREEN |
| Source/build revision | `ProjectRevision`, `GitCommit`, `GitBranch`, `GitDirty` | GREEN from general-v1.10 |

`/output/fileName` is intentionally not duplicated in `RunConfig`: it is a storage destination/name, not a physics or detector configuration parameter.

## Curated-material policy

Ordinary users select validated curated profiles rather than arbitrary material files. Current role eligibility is deliberately conservative:

- Core: `Water_Default` or `PMMA_SNO_A`.
- OuterShell: `PMMA_SNO_A`.
- Coupling: `PDMS_Generic`.

Expert dataset overrides exist for controlled validation. Adding a new production profile should include provenance and validation rather than silently loading an arbitrary table.

## Physics baseline

The validated baseline remains `FTFP_BERT + G4EmStandardPhysics_option4 + G4OpticalPhysics`, with the documented default production cut and optical-process activation. Framework work after the `physics_validated_v1` tag is additive configuration/metadata/diagnostic work unless a change is explicitly revalidated as a physics change.

## Schema/provenance policy

RunConfig is append-only across the recent schema series:

- v1.6 prefix: columns 0--98.
- v1.8 prefix: columns 0--153.
- v1.9: Core first-exit metadata at 154--156.
- v1.10: Git source/build provenance at 157--159.

CMake captures Git state at **configure time**. If source files change afterwards, rerun `cmake ..` before rebuilding so `GitDirty`/revision metadata reflect the current source tree. A build from a source archive without `.git` remains supported and records `unknown` Git fields.


## User-facing macro policy

The flat `macros/` directory is now explicitly catalogued in `macros/README.md`. Ordinary users start from `example_*.mac`; `validation_*.mac` are controlled regression samples and `study_*.mac` are downstream analysis use cases rather than framework defaults. `tools/check_macro_catalog.py` statically enforces that every macro belongs to a known category, every public example is documented, public examples use canonical rather than legacy detector commands, and examples contain the minimum detector/readout/output/generator/run structure.

Current public starting points are a deterministic Gun+DirectSiPM example, CRY+DirectSiPM, and CRY+ConcentratorInput. These are examples of the configuration surface, not claims of an optimized detector design.

## Framework completion rule

A new user-facing parameter/feature is complete only when it has:

1. a stable macro/config command (or a documented reason to remain fixed);
2. resolved metadata in the ROOT output;
3. a static/runtime validation that checks its semantics;
4. documentation and at least one representative example when appropriate.

Detector-performance studies (best face, best sensor placement, optimal concentrator, etc.) are examples built **on** this framework and should not drive the core framework roadmap unless they expose a missing simulation capability or an implementation error.

## Verification status of the general-v1.10 framework patch

Static checks run on the patched source:

- `python3 tools/check_ntuple_schema.py` — PASS (16 ntuples; `RunConfig` 160 columns).
- `python3 tools/check_analysis_schema.py` — PASS.
- `python3 tools/check_macros.py` — PASS (51 executable project commands; 6 UI directories).
- `python3 tools/check_config_metadata_contract.py` — PASS (all 51 commands have an explicit metadata/reproducibility contract).
- `python3 validate_physics.py` — PASS; the script now reads the active production `PMMA_SNO_A` transport table rather than the historical measured-only table.
- `python3 tools/cry_magnetic_latitude.py` — PASS.
- CMake Git-provenance generation was exercised independently and produced commit/branch/dirty/revision fields correctly.

Target-machine runtime gates have now been exercised after the v1.10 patch:

- MT/output/provenance smoke: PASS (`general-v1.10`, Git provenance present and internally consistent).
- Deliberately asymmetric 12-channel flexible layout: PASS (`ChannelMap=12`, six `ReadoutPlaneConfig` rows).
- Integrated 24-channel `DirectSiPM` optical run: PASS with non-zero `Photon`, `SensorIncident`, `DetectedSignal`, and `ChannelEvent`.
- Core-boundary OFF/basic/phaseSpace/TIR regression: PASS.
- User-facing macro catalog/static policy: PASS (`tools/check_macro_catalog.py`).

See `FRAMEWORK_CONFIG_AUDIT.md` for the complete command-to-metadata contract and the exact runtime evidence. The dedicated physics validation suite remains a separate freeze gate.
