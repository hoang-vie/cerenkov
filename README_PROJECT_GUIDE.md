# Project guide

## Code structure

| File or directory | Responsibility |
|---|---|
| `cherenkov.cc` | Run manager, random seed, physics initialization and macro execution. |
| `DetectorConstruction` | Geometry, optical material profiles, surfaces and channel mapping. |
| `DetectorMessenger` | Detector/readout/trigger/output commands. |
| `PrimaryGeneratorAction`, `PrimaryGeneratorMessenger` | Generic Gun and CRY sources. |
| `TrackingAction` | Start/end truth for non-optical tracks. |
| `SteppingAction` | Step deposits, optical transport bookkeeping and readout boundary handling. |
| `EventAction` | Event buffers, pulse grouping, coincidence trigger and event output. |
| `RunAction` | Ntuple definitions, run metadata, per-worker statistics and file lifecycle. |
| `analysis/` | ROOT inspection, dedicated validators and general reporting. |
| `data/` | Versioned optical tables, model alternatives and provenance comments. |
| `tools/` | Static consistency checks and supporting utilities. |
| `macros/` | Configuration examples, validation samples and visualization. |

Signal collection is implemented in stepping/event actions. `ConstructSDandField()` is empty; the current model does not use separate SensitiveDetector classes.

## Event lifecycle

1. The configured source creates primary vertices.
2. Geant4 transports particles and creates secondaries, including Cherenkov photons.
3. Tracking and stepping actions collect truth, energy deposits and optical history.
4. Readout arrival is recorded separately from the PDE detection decision.
5. At event end, deposits in each trigger paddle are sorted and grouped using the configured integration window. Qualifying top/bottom pulses are paired within the coincidence window.
6. The selected output policy determines which event records are written.

The integration window is anchored to the start of each deposit group. Trigger timing is an energy-deposition approximation, not a digitized electronics waveform.

## Configure a study

Copy an `example_*.mac` and change one group of parameters at a time:

- `/detector/*`: dimensions, material profiles, shell, Tyvek and coupling.
- `/readout/*`: independent fixed grids or auto-fill layouts on six faces.
- `/generator/*`, `/CRY/*`: source settings.
- `/trigger/*`: deposit threshold and integration/coincidence windows.
- `/output/*`: filename, event policy and optional detail.

Set geometry before `/run/initialize`. In a multithreaded application, geometry and optical-table changes at Idle are rejected by the project messenger; use a separate process for each such configuration. Generator commands are normally placed after initialization, as in the examples.

`/generator/position` disables automatic source-Y placement. Put `/generator/useAutoY true` after it if automatic placement is intended. For CRY, call `/CRY/update` after accumulating `/CRY/input` lines.

Inspect `ReadoutPlaneConfig` for requested settings and `ChannelMap` for actual channels. Face IDs are `+X,-X,+Z,-Z,+Y,-Y` in order 0–5. Core volume ID is 0; resolve other volume categories through `VolumeMap`.

## Reproducible work

Pass a positive seed as the second executable argument. Record the macro, source revision, dependency versions, active optical tables and analysis selection. Git metadata is captured at CMake configure time: rerun CMake after relevant source/Git changes before building an archived result. A dirty flag identifies local modifications but does not store their content; retain the actual diff/source snapshot.

Use IDs to join ROOT trees. ROOT row order can change with multithreaded execution. Separate processes may reuse both RunID and EventID, so retain a file/sample identifier when combining datasets.

## Reading order

Start with [README](README.md), then [General model](README_GENERAL.md), [Validation](README_VALIDATION.md), [Schema](SCHEMA.md) and [branch reference](README_TREE_INFO.md). Use the [macro catalog](macros/README.md) to choose a starting configuration.
