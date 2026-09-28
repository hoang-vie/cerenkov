> **hotfix6_physics_audit7_final_integrated_freeze status:** physics-validation
> baseline complete. Deterministic Geant4+ROOT validation, schema/static checks,
> material/surface/sensor validation, charged-particle transport, muon decay,
> ancestry/bookkeeping, and the final integrated production-like regression pass.
> Remaining Tyvek, SiPM-PDE, and PDMS limitations are documented input-model
> caveats rather than unresolved implementation failures.

# Cherenkov general simulation redesign

This version generalizes particle truth, optical-photon ancestry, readout geometry, trigger logic, and output policy while keeping large optical output compact.

## Build

```bash
mkdir build && cd build
cmake ..
cmake --build . -j
```

CRY 1.7 is fetched by the existing ExternalProject rule. Optical-data and CRY-data paths are compiled as absolute paths, avoiding dependence on the current working directory.
 Git provenance is captured at CMake configure time; rerun `cmake ..` after committing or changing source files before the build you intend to archive/publish.

## Ordinary study examples

Start from the catalogued user examples rather than a validation macro:

```bash
./cherenkov example_gun_direct_sipm.mac 1234567
./cherenkov example_direct_sipm.mac 1234567
./cherenkov example_concentrator_input.mac 1234567
```

See `macros/README.md` for the intent of each category and which file to copy for a new study.

The optional second argument is an explicit RNG seed. Without it a high-resolution/random-device seed is generated; the actual seed is written to `RunConfig`. `RunConfig` also stores the schema version, actual thread count, RNG engine, CRY version, Geant4 version, default production cut and the activation state of all registered optical processes. From `general-v1.10` it additionally stores the Git commit, branch, dirty-tree flag, and a compact `ProjectRevision`, so a ROOT file identifies the source revision used to build the executable rather than relying on a manually edited hotfix label.

## Generator

CRY:
```text
/generator/mode CRY
/CRY/timeMode relative
... CRY input lines ...
/CRY/update
```

Generic validation/source gun:
```text
/generator/mode Gun
/generator/particle pi+
/generator/energy 1 GeV
/generator/position 0 60 0 mm
/generator/direction 0 -1 0
# Needed for primary optical-photon boundary tests:
/generator/polarization 1 0 0
```

The gun can therefore run the same detector for muons, electrons, pions, protons, kaons, etc. without changing the output schema. The bundled CRY configuration should be limited to species supported by CRY 1.7 (muons, neutrons, protons, electrons, photons and pions); other Geant4 species are available through the generic gun.

## Optical processes

The executable registers Geant4 optical physics with Cerenkov ON, absorption/Rayleigh/boundary ON, Scintillation OFF and WLS OFF. Geant4's `/process/optical/processActivation` command can activate registered processes later. The photon schema already distinguishes optical creator process type/subtype and propagates non-optical source ancestry through optical chains. The actual activation state (including WLS2) is written into `RunConfig` rather than inferred from a macro filename.

The schema is therefore ready for future scintillation/WLS output, but the present detector materials do **not** yet define a calibrated scintillation emission/yield model or a WLS-fibre material/geometry. Turning those processes on is not, by itself, a validated experimental model; their material property tables must be added from the actual detector/fibre data.

## Readout

`/detector/sensorMode DirectSiPM` applies PDE. See `macros/example_direct_sipm.mac`.

`/detector/sensorMode ConcentratorInput` records readout-entrance phase space and absorbs the virtual-aperture photon without counting a detection. The interactive `init_vis.mac` uses this mode for the current 10 x 10 mm concentrator-entrance study.

Readout element count is not fixed. Every face can be disabled, use an independent fixed `NU x NV` grid, or auto-fill from element size/gap/margins. `ChannelMap` stores the actual generated elements and `ReadoutPlaneConfig` stores the requested layout. See `SCHEMA.md` and `macros/example_concentrator_input.mac`.

For new macros, use the curated profile commands (`/detector/coreProfile`, `/detector/outerShellProfile`, `/detector/couplingProfile`) and the generic geometry commands (`/detector/readoutSize`, `/detector/outerShellThickness`, `/detector/useOuterShell`). Legacy aliases such as `/detector/coreMaterial`, `/detector/sipmSize`, and `/detector/setAcrylicThickness` remain supported for old validation macros.

## Analysis

ROOT example:
```text
root -l
.L analysis/analyze_general.C+
analyze_general("output_cerenkov_run0.root",13,true,"primaryEmitter");
analyze_general("output_cerenkov_run0.root",11,true,"primaryEmitter");
analyze_general("output_cerenkov_run0.root",211,false,"primaryEmitter");
analyze_general("output_cerenkov_run0.root",-211,false,"primaryEmitter");
# selectedPDG=0 selects all optical photons, useful for optical-primary validation:
analyze_general("output_cerenkov_run0.root",0,true,"emitter");
```

The report uses the lightweight `PrimaryTruth` tier when present so selected-primary energy/zenith distributions are not trigger-biased in `thinAll`. It also includes dE/dx, generated spectrum, Cherenkov angle and angle-vs-wavelength, transport funnel, incident spectrum/angle/aperture map, detected spectrum/timing, transport/sensor/total efficiencies, **species-selected** channel occupancies, reflection/Rayleigh summaries, trigger QA, and muon decay/Michel plots when applicable. Detected timing is plotted relative to the selected top-scintillator trigger pulse. Free-lifetime auto-fitting is intentionally restricted to an explicit `mu+` selection (`selectedPDG=-13, combineCharge=false`).

## SNO PMMA reference

Exact published SNO UV-transmitting acrylic Sample-A values are retained as the empirical anchor for PMMA absorption in the 300--500 nm measured domain. The active production transport table is `PMMA_SNO_UVT_SampleA_ABS_Measured300to500_Model200to900.txt`: 300--500 nm is explicitly tagged as measured Sample-A data, while 200--300 nm and 500--900 nm are explicitly labelled transport-policy extensions rather than direct SNO measurements. `RunConfig` stores both numerical support (`*AbsValidMin/Max_nm`) and empirical validity (`*AbsMeasuredMin/Max_nm`) together with the exact dataset/profile provenance. The PMMA RINDEX uses the documented broad Cauchy approximation profile rather than pretending the limited published SNO refractive-index measurements cover the full optical transport band.

The SNO material is a detector-grade UV-transmitting reference, not evidence that an arbitrary laboratory PMMA sheet has identical properties.

## Validation macros

- `validation_muon_decay.mac`: stopped mu+ proper-time / Michel-electron test.
- `validation_dedx_muon.mac`: 3 GeV mu- dE/dx with optical production disabled.
- `validate_physics.py`: independent table-based Cherenkov/TIR/attenuation checks.

## Build-status note

The simulator needs a Geant4 development installation; the supplied ROOT analysis macro additionally needs ROOT. The development container used to prepare this revision has neither a Geant4 package configuration nor `root-config`, so the code can be statically audited here but cannot be linked/run against Geant4 or ROOT here.


## Runtime validation

Follow `VALIDATION_RUNBOOK.md`. Run `validation_output_mt.mac` first; it specifically checks ROOT/Geant4 MT ntuple merging before any physics plots are interpreted.
