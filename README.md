# Cherenkov detector simulation framework

This repository is a Geant4 framework for **Cherenkov detector studies**. Its goal is not to prescribe one optimized detector configuration, but to provide a validated physics baseline and a configurable detector/readout/generator/output interface that other users can run without editing C++ for ordinary scans.

Current schema: `general-v1.10`.

## Quick start

Build:

```bash
mkdir -p build
cd build
cmake ..
cmake --build . -j
```

Run a deterministic gun example:

```bash
./cherenkov example_gun_direct_sipm.mac 12345
```

Run a CRY + DirectSiPM example:

```bash
./cherenkov example_direct_sipm.mac 12345
```

Run a CRY + ConcentratorInput example:

```bash
./cherenkov example_concentrator_input.mac 12345
```

The optional second argument is the RNG seed. The resolved run configuration, actual seed, thread count, material/profile provenance, channel geometry, and Git source/build provenance are written into the ROOT output.

## Where to start reading

- `macros/README.md` — which macro to copy for a new study and how macros are categorized.
- `README_PROJECT_GUIDE.md` — repository architecture and validation philosophy.
- `FRAMEWORK_STATUS.md` — supported topology, public configuration surface, fixed assumptions, and completion rules.
- `SCHEMA.md` / `README_TREE_INFO.md` — ROOT data model and branch semantics.
- `VALIDATION_RUNBOOK.md` — ordered validation/freeze procedure.
- `FRAMEWORK_CONFIG_AUDIT.md` — command-to-metadata contract.

## Public study modes

The ordinary user workflow is:

1. copy one `example_*.mac` file;
2. change detector/readout/generator/output configuration through macro commands;
3. choose an explicit seed for reproducible studies;
4. run `./cherenkov <macro> <seed>`;
5. use `RunConfig`, `ReadoutPlaneConfig`, `ChannelMap`, and the event/photon/readout trees as the authoritative record of what was simulated.

The two readout response modes are:

- `DirectSiPM` — readout entrance applies the configured effective SiPM PDE and writes `DetectedSignal`;
- `ConcentratorInput` — readout entrance records optical phase space for downstream concentrator/ray-tracing studies and terminates the photon without applying SiPM PDE.

`CoreBoundaryCrossing` is an optional diagnostic checkpoint for first successful optical escape from the Core; it is not a physical detector element.

## Macro categories

Do not treat every `.mac` file as a production template:

- `example_*.mac` — user-facing starting points for ordinary studies;
- `validation_*.mac` — controlled correctness/regression tests;
- `study_*.mac` — analysis-specific studies retained as examples of framework use, not framework defaults;
- `r1_shell_*.mac` — historical/legacy compatibility examples;
- `init_vis.mac`, `vis.mac` — interactive visualization support.

See `macros/README.md` for the full catalog and expected use.

## Validation before release/freeze

At minimum run:

```bash
python3 tools/check_ntuple_schema.py
python3 tools/check_analysis_schema.py
python3 tools/check_macros.py
python3 tools/check_config_metadata_contract.py
python3 tools/check_macro_catalog.py
python3 validate_physics.py
python3 tools/cry_magnetic_latitude.py
```

Then follow `VALIDATION_RUNBOOK.md` for Geant4/ROOT runtime gates.

## Scope

The supported detector family is intentionally finite: rectangular Core, optional uniform outer shell, optional Tyvek reflector, top/bottom scintillator trigger paddles, and rectangular coupling/readout elements on any of six axis-aligned faces. Arbitrary CAD/cylindrical/spherical geometry is outside the current framework version.

Material and sensor models are curated profiles with documented provenance and limitations. A validation PASS demonstrates that the implementation reproduces the configured model; it does not by itself prove that every input table is representative of a future physical detector.

## Freeze / release validation

For a framework release candidate, use `FRAMEWORK_FREEZE_CHECKLIST.md`. The one-command gates are:

```bash
tools/run_framework_static_checks.sh
tools/run_framework_freeze_regression.sh
```

The runtime freeze gate requires a clean Git tree so the produced ROOT metadata can prove `GitDirty=0` and match the exact release-candidate commit.
