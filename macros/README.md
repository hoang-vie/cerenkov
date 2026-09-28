# Macro catalog and user workflow

All macros configure the same `./cherenkov` executable. There is no separate validation simulator or production simulator. The categories below describe **intent**, not different physics engines.

## 1. User-facing examples — copy these for a new study

### `example_gun_direct_sipm.mac`
Deterministic 3 GeV downward `mu-` gun with six-face DirectSiPM readout. Use this when learning the configuration surface, debugging a geometry/readout change, or building a reproducible mono-particle study.

Run from `build/`:

```bash
./cherenkov example_gun_direct_sipm.mac 12345
```

### `example_direct_sipm.mac`
CRY sea-level muon example using DirectSiPM response. Use this as the starting point for detector-response/cosmic studies where the configured effective SiPM PDE is part of the simulation.

```bash
./cherenkov example_direct_sipm.mac 12345
```

### `example_concentrator_input.mac`
CRY example using `ConcentratorInput`. It records readout-entrance optical phase space instead of applying SiPM PDE. Use this when the downstream optical concentrator is to be studied separately (for example with ROBAST).

```bash
./cherenkov example_concentrator_input.mac 12345
```

These examples are **starting configurations, not optimized detector designs**. Copy one, change the parameters needed by the study, and rely on ROOT metadata rather than the macro filename as the authoritative run record.

## 2. Validation macros — do not use as production templates

Files named `validation_*.mac` create controlled samples for implementation/physics/regression checks. Their unusual geometry, disabled processes, monoenergetic optical photons, forced directions, or special output settings may be intentional.

Use them only with the corresponding validator/runbook step. See `VALIDATION_RUNBOOK.md`.

## 3. Study macros — downstream framework use cases

Files named `study_*.mac` are retained study configurations, currently including CRY/CoreBoundary examples. They demonstrate how the framework can be used for analysis questions, but their settings are not framework defaults and they are not required for framework correctness.

## 4. Legacy/historical macros

`r1_shell_*.mac` are retained to exercise compatibility with earlier shell APIs. New user macros should prefer the canonical commands documented in `FRAMEWORK_STATUS.md` and `FRAMEWORK_CONFIG_AUDIT.md`.

## 5. Visualization macros

- `init_vis.mac` — interactive visualization startup used when `./cherenkov` is launched without a macro argument.
- `vis.mac` — visualization helper commands.

They are not batch-production templates.

## 6. What to edit in a copied example

Typical study parameters are grouped as follows:

- execution: `/run/numberOfThreads`;
- Core/shell/Tyvek: `/detector/coreProfile`, `/detector/coreSize`, `/detector/useOuterShell`, `/detector/outerShellProfile`, `/detector/outerShellThickness`, `/detector/tyvekPosition`, `/detector/setTyvekThickness`;
- coupling/readout element: `/detector/couplingProfile`, `/detector/couplingSize`, `/detector/readoutSize`, `/detector/sensorMode`;
- readout layout: `/readout/uniformGrid`, `/readout/planeGrid`, `/readout/planeAuto`;
- trigger: `/trigger/*`;
- output: `/output/*`;
- source: `/generator/*` or `/CRY/*`.

For normal geometry scans in an MT build, run a **new process/configuration for each geometry** rather than repeatedly rebuilding geometry between multi-event MT runs.

## 7. Initialization order

A safe ordinary pattern is:

```text
# detector/readout/output configuration
/run/initialize
# generator/CRY configuration
/run/beamOn N
```

CRY additionally requires `/CRY/update` after its `/CRY/input ...` lines (or use `/CRY/file`).

## 8. Reproducibility

For a run intended to be compared, archived, or published:

- pass an explicit positive RNG seed as the second executable argument;
- rerun `cmake ..` before rebuilding after source/Git-state changes so build provenance is current;
- keep the ROOT file, because it stores the resolved configuration and source/build provenance;
- record external study context that is not part of the simulator itself (analysis selection, dataset grouping, scan labels, etc.).

## 9. Optional diagnostics

`/output/coreBoundaryMode firstExit` enables the non-invasive `CoreBoundaryCrossing` diagnostic. Use `basic` for position/direction/wavelength/time studies and `phaseSpace` when polarization is required. Leave it `off` for ordinary production when that checkpoint is not needed, because it can substantially increase file size.
