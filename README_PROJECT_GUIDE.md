# Project guide — Geant4 Cherenkov detector

> Framework baseline documented here: **general-v1.10 on top of the `physics_validated_v1` physics freeze**. The schema remains append-only; source/build provenance and the Core first-exit diagnostic are additive.
>
> This document explains **what each project document/tool is for**, what is considered simulation, validation, or analysis, and in what order a new user should read the repository.

## 1. What this project is

This repository is a Geant4 framework for studying Cherenkov-light production and transport in a configurable detector. The current detector model can contain a water Core, optional acrylic/PMMA shell, Tyvek reflector, optical coupling, scintillator trigger paddles, and dynamically configured readout elements on any of the six detector faces.

Two readout modes are supported:

- **`DirectSiPM`** — a readout element represents an effective SiPM entrance. A wavelength-dependent PDE is applied at the readout boundary and successful detections are written to `DetectedSignal`.
- **`ConcentratorInput`** — a readout element represents the entrance aperture of an optical concentrator. No fake SiPM detection is applied; full local phase space is written to `SensorIncident` for later ROBAST/concentrator studies.

The project is deliberately split into three logical layers:

```text
                     configuration / macro
                            |
                            v
                  +--------------------+
                  | Geant4 simulation  |
                  |    ./cherenkov     |
                  +---------+----------+
                            |
                         ROOT file
                            |
               +------------+-------------+
               |                          |
               v                          v
        dedicated validation       production analysis
          validate_*.C             analyze_general.C
               |                          |
       compare with theory/          physics plots,
       external expectation          efficiencies, QA
```

A validation macro and a production macro therefore use the **same simulation executable and the same physics implementation**, but usually generate different event samples because their configurations are intentionally different.

## 2. Validation data versus production data

This distinction is central to understanding the repository.

### Validation data

Validation samples are deliberately simple and controlled. Examples:

- one monoenergetic 3 GeV muon through 10 cm water;
- stopped `mu+` for lifetime validation;
- monochromatic optical photons at one wavelength;
- known incidence angles and polarizations at one optical interface;
- one known ray on each of the six readout faces.

The purpose is to isolate one physical process and compare the Geant4 result to an independent expectation.

Examples:

- Cherenkov angle: `cos(theta_C) = 1 / (beta*n(lambda))`;
- Frank–Tamm yield;
- optical absorption: `P(x) = exp(-x/L_abs)`;
- Rayleigh mean free path;
- Fresnel coefficients and Snell's law;
- SiPM `Ndet/Ninc = PDE(lambda)` in a controlled mono-wavelength test.

### Production / ordinary simulation data

Production-like samples may use CRY, random cosmic rays, realistic geometry, all enabled optical processes, trigger selection, and many particle species. They are not expected to reproduce one fixed validation number such as `817 photons/cm` or `PDE=0.400` globally, because their energies, wavelengths, angles and track lengths vary event by event.

The **physical rules** validated in the controlled samples must still hold. For example, a Cherenkov photon in a production sample must still satisfy the Cherenkov-angle relation for its own emitter beta and wavelength.

### Analyzer versus validator

- `analysis/validate_*.C` files are **correctness tests**. They read ROOT data and check a narrowly defined process against theory/reference values.
- `analysis/analyze_general.C` is a **general production-analysis/reporting tool**. It selects species/ancestry, constructs spectra, efficiencies and detector-level plots.

Neither one generates Geant4 physics. Both consume ROOT files produced by `./cherenkov`.

The dedicated validation suite should be considered separately from `analyze_general.C`. A bug in plotting or selection in `analyze_general.C` does not automatically mean the Geant4 simulation is wrong.

## 3. Recommended reading order

For a new user, read the repository in this order:

1. **`README_PROJECT_GUIDE.md`** — this file. Understand the architecture and document map.
2. **`FRAMEWORK_STATUS.md`** — current framework scope, canonical configuration surface, configuration-to-metadata contract, and completion rule.
3. **`FRAMEWORK_CONFIG_AUDIT.md`** — complete project-command to ROOT-metadata contract and intentionally fixed assumptions.
4. **`README_VALIDATION.md`** — practical commands for building and validating the project.
5. **`README_TREE_INFO.md`** — ROOT ntuple/tree and branch reference.
6. **`README_GENERAL.md`** — high-level design of the general-v1 simulation architecture.
7. **`SCHEMA.md`** — formal output semantics, identities, joins, modes and ancestry definitions.
8. **`PHYSICS_VALIDATION.md`** — physics reasoning behind the validation strategy and what each validation actually proves.
9. **`VALIDATION_RUNBOOK.md`** — current validation command sequence maintained during development.
10. **`CHANGELOG_GENERAL_V1.md`** — chronological implementation changes.
11. **`AUDIT_NOTES.md`** and the `WATER_DEFAULT_*`, `PMMA_SNO_A_*`, `R2E_*`, and stage-specific audit files — historical development evidence; they describe the schema/state that existed at that audit stage and are not the current user-facing specification.

Use the current user-facing documents above as the source of truth for the framework interface. Stage/audit notes remain valuable provenance for why a physical model or implementation choice was made, but their historical schema labels should not be read as the current schema.

## 4. Existing documentation and what each file means

### User-facing architecture/reference

#### `FRAMEWORK_STATUS.md`
Defines the current supported detector topology, canonical macro commands, configuration-to-metadata contract, curated-material policy, and the rule for considering a new framework feature complete.

Use it when asking: **“What is the framework supposed to support, and what is still framework work versus a downstream detector study?”**

#### `FRAMEWORK_CONFIG_AUDIT.md`
Lists all executable project UI commands, their authoritative ROOT metadata, legacy/expert exceptions, and intentionally fixed simulation assumptions. It is backed by `tools/check_config_metadata_contract.py`, which fails if a newly registered command has no declared reproducibility contract.

Use it when asking: **“If I change this macro command, where is the resolved value recorded, and can I reproduce the run later?”**

#### `README_GENERAL.md`
Overview of the general-v1 architecture: particle-general truth, flexible readout geometry, output policies, trigger model, DirectSiPM/ConcentratorInput modes and species-general analysis.

Use it when asking: **“What can the simulation do?”**

#### `SCHEMA.md`
Canonical description of ROOT output semantics: ntuples, identities, joins, ancestry, output modes, channel/volume maps and readout modes.

Use it when asking: **“What does this field/tree actually mean?”**

#### `README_TREE_INFO.md`
User-friendly branch-by-branch companion to `SCHEMA.md`, including branch types, units and common joins.

Use it when writing ROOT/Python analysis or debugging a ROOT file.

### Physics validation

#### `PHYSICS_VALIDATION.md`
Explains the validation philosophy and physics expectations. The key idea is that validation is intentionally separated from CRY production: simple deterministic tests establish correctness before complex cosmic-ray data are interpreted.

Use it when asking: **“Why is this test valid and what does PASS actually prove?”**

#### `README_VALIDATION.md`
Operational guide: build commands, macro/validator pairings, expected PASS messages, generated plots, and interpretation of failures.

Use it when asking: **“What commands do I run?”**

#### `VALIDATION_RUNBOOK.md`
Development-era runbook containing the validation sequence and historical additions. It overlaps with `README_VALIDATION.md`; keep it as provenance and a lower-level checklist.

### Audit/history

#### `HOTFIX6_PHYSICS_AUDIT.md`
Documents the major physics/material audit: optical-material provenance, generic PDMS coupling, PMMA limitations, metadata improvements and runtime validation additions.

#### `HOTFIX6_AUDIT2_NOTES.md`
Documents fixes found by runtime testing, including placement/rotation and trigger-macro corrections.

#### `HOTFIX6_AUDIT3_NOTES.md`
Documents the six-face validator correction: rare opposite-face hits caused by real Fresnel reflection must not be classified as coordinate-transform failures.

#### `HOTFIX6_AUDIT4_NOTES.md`
Documents `Photon.EmitterBeta` and the exact per-photon Cherenkov-angle validation/report additions.

#### `CHANGELOG_GENERAL_V1.md`
Long-form chronological change history.

#### `AUDIT_NOTES.md`
Earlier audit notes and intentionally unresolved items. Do not treat an old open item as necessarily still open without checking later audit notes/changelog.

## 5. Important code directories

### `src/` and `include/`
The Geant4 simulation implementation.

Important responsibilities:

- `DetectorConstruction.*` — geometry, materials, optical surfaces, readout placement.
- `DetectorMessenger.*` — detector/readout configuration commands.
- `PrimaryGeneratorAction.*` / `PrimaryGeneratorMessenger.*` — gun/CRY/optical primary configuration.
- `TrackingAction.*` — non-optical track registration/final truth.
- `SteppingAction.*` — step-level transport bookkeeping, optical boundary/readout logic, Cherenkov photon information.
- `EventAction.*` — event aggregation, trigger logic and per-event output preparation.
- `RunAction.*` — ROOT ntuple schema, run metadata and output writing.
- `ActionInitialization.*` — installs user actions.
- `SimulationMetadata.hh` — shared project/version metadata.

### `macros/`
Geant4 configuration files. They do not contain a second simulation engine; they configure the same `./cherenkov` executable. The directory intentionally remains flat because CMake exposes `macros/*.mac` in the build directory. `macros/README.md` is the user-facing catalog.

Intent categories:

- `example_*.mac` — user-facing starting points for ordinary studies;
- `validation_*.mac` — controlled correctness/regression samples, not production templates;
- `study_*.mac` — downstream analysis-specific use cases, not framework defaults;
- `r1_shell_*.mac` — historical/legacy compatibility macros;
- `init_vis.mac` / `vis.mac` — interactive visualization support.

New users should copy an `example_*.mac` file rather than adapting a validation macro. `tools/check_macro_catalog.py` enforces this catalog and rejects legacy public commands in user-facing examples.

### `analysis/`
ROOT post-processing.

- `validate_*.C` — dedicated correctness validators;
- `inspect_output.C` — structural/schema sanity checks;
- `analyze_general.C` — general production/report analysis;
- `analyze_legacy_reference.C` — retained historical reference, not the canonical general-v1 analysis.

### `data/`
Optical/material tables and provenance notes. These are part of the physical model and should be versioned/fingerprinted with production results.

### `tools/`
Static consistency checks and utility calculations.

## 6. Current validation status and freeze policy

At the hotfix6/audit4 stage, dedicated runtime tests have demonstrated consistency of the implemented model for:

- ROOT schema and MT merging;
- deterministic geometry/path length;
- six-face readout rotations and world-to-local phase space;
- flexible channel placement and unit normals;
- trigger logic in the controlled through-going muon test;
- muon lifetime and Michel spectrum;
- Cherenkov Frank–Tamm yield;
- per-photon Cherenkov angle using `EmitterBeta`;
- PMMA absorption at the tested measured point;
- configured water Rayleigh mean free path;
- Fresnel s/p reflection, Snell refraction and TIR;
- effective Tyvek total reflectivity at the tested wavelength;
- effective S13360-3050CS PDE at 450 nm;
- DirectSiPM and ConcentratorInput structural behavior.

This does **not** mean every physical input is a measured model of the final detector. Current model limitations/assumptions include generic PDMS coupling, effective Tyvek rather than a measured detector-specific BRDF, effective-PDE SiPM response, and limited measured PMMA attenuation coverage.

A sensible freeze policy is:

- freeze simulation core + ROOT truth semantics + dedicated validators only after the full dedicated suite is reproducibly PASS;
- keep production QA and `analyze_general.C` on a separate development track until their selections/denominators/plots are audited;
- after freeze, any change to material tables, physics configuration, schema or detector transport code must create a new revision and rerun the affected validation tests.

## 7. Key principle

A PASS should always be interpreted narrowly.

For example:

- PMMA absorption PASS means Geant4 reproduces the supplied attenuation length in the dedicated test. It does not prove an arbitrary future PMMA sheet has that attenuation.
- Tyvek PASS means the effective model reproduces the configured total reflectivity. It does not validate a real detector-specific BRDF.
- PDE PASS means the effective boundary implementation reproduces the configured PDE under the controlled test condition. It does not model dark counts, afterpulsing or full SiPM electronics.

The repository should therefore distinguish three questions:

1. **Is the implementation internally correct?** — dedicated validation.
2. **Are the chosen physical input data representative of the real detector?** — material/sensor measurements and systematic studies.
3. **What does a production/cosmic sample imply for detector performance?** — production QA and `analyze_general.C`.
