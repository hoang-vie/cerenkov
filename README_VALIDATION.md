# Validation guide — Geant4 Cherenkov project

> Framework baseline: **general-v1.10 on top of the `physics_validated_v1` physics freeze**. Physics validation remains the same baseline; v1.9/v1.10 changes are observer/metadata/framework additions.
>
> This guide is operational: it explains what to run, what data each test generates, which ROOT validator reads it, what a PASS means, and when to stop.

## 1. The most important concept

Validation and production use the **same Geant4 executable and simulation code**, but usually use different macros and therefore different event samples.

```text
validation_xxx.mac ------------------+
                                     |
production/example/cosmic macro ----+--> ./cherenkov --> ROOT data

validation ROOT --> validate_xxx.C       (correctness test)
production ROOT --> analyze_general.C    (physics/report analysis)
```

A validator is not a second physics engine. It is a ROOT program that reads simulation output and compares a narrow observable with an analytical/reference expectation.

Do not demand that a random cosmic sample reproduce a fixed number from a controlled validation sample. Example: the 450-nm PDE validation expects about 0.400, but a broadband cosmic sample has a wavelength-weighted detection probability that can be very different.

For `general-v1.10`, the first MT/output smoke test also runs `analysis/validate_framework_metadata.C` to verify that source/build provenance (`ProjectRevision`, `GitCommit`, `GitBranch`, `GitDirty`) is present and internally consistent.

## 2. Before running anything

From the project root, confirm Geant4, ROOT and compiler/CMake are available:

```bash
geant4-config --version
root-config --version
cmake --version
g++ --version
```

If Geant4 or ROOT is installed but not in the current shell environment, source the corresponding setup scripts before configuring the build.

## 3. Static checks

Run from the source root:

```bash
python3 tools/check_ntuple_schema.py
python3 tools/check_analysis_schema.py
python3 tools/check_macros.py
python3 tools/check_config_metadata_contract.py
python3 validate_physics.py
python3 tools/cry_magnetic_latitude.py
```

Meaning:

- `check_ntuple_schema.py` — checks ntuple creation/fill indices and types;
- `check_analysis_schema.py` — checks analysis branch bindings against the simulation schema;
- `check_macros.py` — checks project-specific executable UI commands used by macros;
- `check_config_metadata_contract.py` — requires every registered project command to have an explicit resolved-metadata/reproducibility contract and verifies referenced ROOT branches exist;
- `validate_physics.py` — independent numerical/table sanity calculations; it does **not** run Geant4;
- `cry_magnetic_latitude.py` — converts geographic location to the magnetic-latitude input convention used for the CRY setup.

If one of the first four framework/static checks fails, do not interpret ROOT physics results yet.

## 4. Clean build

Recommended when changing C++ source/schema or switching between important revisions:

```bash
rm -rf build
mkdir build
cd build
cmake ..
cmake --build . -j8
```

CMake copies `macros/*.mac` into `build/`. If only a source macro changed, rerun at least:

```bash
cmake ..
```

before executing the copied macro from `build/`.

## 5. Validation order

Run tests in this order. The order isolates cheap structural failures before expensive optical transport.

### 5.1 ROOT/MT output sanity

```bash
./cherenkov validation_output_mt.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_output_mt_run0.root",100,0)'
```

Expected final line:

```text
OUTPUT SANITY: PASS
```

This verifies the ROOT file is structurally usable, run statistics merge correctly in MT, expected ntuples exist and event/channel counts are consistent.

**If this fails: stop.** Do not trust later analysis.

### 5.2 Six-face geometry and coordinate transform

```bash
./cherenkov validation_phase_space_6faces.mac 12345
root -l -b -q '../analysis/validate_phase_space_6faces.C()'
```

Expected:

```text
PHASE-SPACE 6-FACE VALIDATION: PASS
```

Checks all six faces (`+X`, `-X`, `+Z`, `-Z`, `+Y`, `-Y`) using known rays.

For accepted entrance hits:

- `EntranceZ_mm` must be approximately 0;
- local `+z` must point outward;
- `Forward == 1` means `LocalDirZ > 0`;
- known normal-incidence rays should have local direction approximately `(0,0,1)` and incident angle near 0 deg.

Rare hits on the exact opposite face can be physical Fresnel-reflected photons and are accounted for by the audit3 validator.

Plots:

```text
build/validation_plots/phase_space_6faces_validation.{png,pdf}
```

### 5.3 Trigger logic

```bash
./cherenkov validation_trigger_muon.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_trigger_muon_run0.root",1000,0)'
```

This is the dedicated trigger test. A straight 3 GeV muon starts above the top paddle and crosses top scintillator -> radiator -> bottom scintillator.

With the configured low threshold and wide coincidence window, trigger efficiency should be close to unity (the current deterministic test has produced 1000/1000).

Do **not** use `validation_dedx_muon.mac` trigger fraction as a trigger validation observable; that macro positions the source specifically to isolate the 10-cm Core path.

### 5.4 Charged-particle path length and dE/dx bookkeeping

```bash
./cherenkov validation_dedx_muon.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_dedx_muon_run0.root",10000,0)'
root -l -b -q '../analysis/validate_dedx.C("validation_dedx_muon_run0.root",13,true)'
```

Checks:

- straight central Core path approximately 10 cm for the current validation geometry;
- `Edep/L` — local deposited energy per path;
- `DeltaKE/L` — primary kinetic-energy loss per path.

These two quantities are not identical because secondary particles can carry energy out of the volume.

Plot:

```text
build/validation_plots/dedx_validation.{png,pdf}
```

Interpret stopping power only against a reference with the same particle, material, energy and production-cut convention.

### 5.5 Free `mu+` lifetime and Michel spectrum

```bash
./cherenkov validation_muon_decay.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_muon_decay_run0.root",100000,0)'
root -l -b -q '../analysis/validate_muon_decay.C("validation_muon_decay_run0.root")'
```

Expected free-muon lifetime near 2.197 microseconds and a Michel positron spectrum with endpoint near 52.3 MeV kinetic energy.

Use stopped `mu+`; do not interpret stopped `mu-` in matter as a free-lifetime measurement because nuclear capture competes with decay.

Plot:

```text
build/validation_plots/muon_decay_validation.{png,pdf}
```

### 5.6 Cherenkov generation — Frank–Tamm and exact angle

First calculate the independent expectation from the active refractive-index tables:

```bash
cd ..
python3 validate_physics.py
cd build
```

Then run:

```bash
./cherenkov validation_frank_tamm_water.mac 12345
root -l -b -q '../analysis/validate_frank_tamm.C("validation_frank_tamm_water_run0.root")'
```

The validation macro isolates Cherenkov production; optical boundary/absorption/Rayleigh transport is not needed for this test.

Checks:

1. Frank–Tamm generated yield per Core path;
2. generated wavelength spectrum;
3. exact per-photon Cherenkov-angle relation using `Photon.EmitterBeta`:

```text
cos(theta_C) = 1 / (EmitterBeta * n(lambda))
```

`EmitterBeta` stores the step-average beta used by the current Geant4 Cherenkov process for that photon.

A successful current run has shown approximately 816.5 photons/cm versus an independent expectation near 817 photons/cm, and an angle residual centered extremely close to zero. These numbers are configuration-specific reference examples, not universal production expectations.

Plot:

```text
build/validation_plots/frank_tamm_validation.{png,pdf}
```

### 5.7 Bulk optical absorption and Rayleigh

Run both samples:

```bash
./cherenkov validation_optical_absorption_pmma.mac 12345
./cherenkov validation_optical_rayleigh_water.mac 12345
root -l -b -q '../analysis/validate_optical_bulk.C()'
```

Expected final line:

```text
OPTICAL BULK PROCESS VALIDATION: PASS
```

Interpretation:

- PMMA test checks that Geant4 obeys the configured attenuation length at the chosen tested wavelength;
- Rayleigh test checks that the first-scatter path follows the configured exponential mean-free-path model.

It does **not** prove that an arbitrary real PMMA/water sample has exactly the configured material constants.

Plot:

```text
build/validation_plots/optical_bulk_validation.{png,pdf}
```

### 5.8 Fresnel, Snell and TIR

```bash
./cherenkov validation_fresnel_tir_pmma_air.mac 12345
root -l -b -q '../analysis/validate_fresnel_tir.C()'
```

Expected final line:

```text
FRESNEL/SNELL/TIR VALIDATION: PASS
```

Checks controlled s/p polarization reflection probabilities, transmitted angle from Snell's law, and total internal reflection above the critical angle.

Plot:

```text
build/validation_plots/fresnel_tir_validation.{png,pdf}
```


#### PDMS_Generic coupling-material validation

The generic coupling profile is intentionally split into two levels of
evidence.

Its RINDEX table,
`OpticalGrease_PDMS_RINDEX_Roychowdhury2018_Cauchy3Anchor.txt`, is a
three-anchor Cauchy approximation. The real PMMA<->PDMS detector interface
has been validated in both directions against analytical Fresnel and Snell
predictions, including polarization dependence, Brewster behavior and the
PMMA->PDMS TIR threshold.

Its ABS table,
`OpticalGrease_PDMS_ABS_Generic10m.txt`, is different: the constant
10-m attenuation length is an explicit low-loss placeholder, not a vendor
measurement. Consequently the profile reports numerical support over
200--900 nm but empirical ABS validity as `(-1,-1)`.

The nominal 10-m implementation was tested directly in the real 1-mm
coupling geometry. A 0.9999-mm 400-nm path produced 101 losses in
1,000,000 photons versus 99.985 expected (`z = 0.102`), confirming the
runtime `m` unit conversion and absorption process.

Analytical single-pass sensitivity shows that for a 1-mm coupling a
1-m attenuation length still gives only about 0.10% loss at normal
incidence and 0.29% at 70 deg, while a much shorter 0.1-m attenuation
length could give approximately 1--3% loss. These are transport
sensitivities, not detector-level efficiency uncertainties.

`PDMS_Generic` is therefore frozen only as a generic detector-design
coupling model. Its ABS assumption must be replaced/revalidated when the
actual optical grease or coupling batch is selected.

### 5.9 Effective Tyvek reflectivity

```bash
./cherenkov validation_tyvek_reflectivity_450nm.mac 12345
root -l -b -q '../analysis/validate_tyvek_reflectivity.C()'
```

Expected final line:

```text
TYVEK REFLECTIVITY VALIDATION: PASS
```

This validates the **effective configured total reflectivity** at the test wavelength. It does not validate a detector-specific measured Tyvek BRDF.

Plot:

```text
build/validation_plots/tyvek_reflectivity_validation.{png,pdf}
```

### 5.10 Effective DirectSiPM PDE

```bash
./cherenkov validation_pde_450nm.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_pde_450nm_run0.root",100000,1)'
root -l -b -q '../analysis/validate_pde.C("validation_pde_450nm_run0.root",0.400)'
```

Checks a controlled monochromatic 450-nm sample against the configured S13360-3050CS effective PDE point.

Expected final line:

```text
DIRECT-SIPM PDE VALIDATION: PASS
```

Plot:

```text
build/validation_plots/pde_450nm_validation.{png,pdf}
```

Do not expect a broadband production sample's integrated `NDetected/NIncident` to equal 0.400. Production efficiency must be compared wavelength-by-wavelength or with a wavelength-weighted expectation.

### 5.11 Flexible readout geometry

```bash
./cherenkov validation_flexible_layout.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_flexible_layout_run0.root",10,12)'
root -l -b -q '../analysis/validate_flexible_layout.C("validation_flexible_layout_run0.root")'
```

Current asymmetric test expects per-face counts:

```text
+X  1
-X  4
+Z  3
-Z  0
+Y  2
-Y  2
Total = 12
```

Checks channel counts, local placements and unit normals.

Plot:

```text
build/validation_plots/flexible_layout_validation.{png,pdf}
```

### 5.12 CRY -> ConcentratorInput structural smoke test

```bash
./cherenkov validation_cry_concentrator.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_cry_concentrator_run0.root",100,24)'
```

Expected:

- 24 configured aperture channels in the current smoke-test geometry;
- `SensorIncident` can be non-zero;
- `DetectedSignal = 0` by construction;
- automatic incident detail is phase-space mode.

This is a smoke test, not yet a final ROBAST/concentrator validation.

### 5.13 CRY -> DirectSiPM structural smoke test

```bash
./cherenkov validation_cry_direct_sipm.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_cry_direct_sipm_run0.root",100,24)'
```

`DetectedSignal` can be non-zero. A 100-event cosmic sample is too small to estimate final detector efficiency reliably.

## 6. General integrated analysis test

`validation_full_muon_analysis_1000.mac` is useful for exercising many already-validated processes together and testing `analyze_general.C`:

```bash
./cherenkov validation_full_muon_analysis_1000.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_full_muon_analysis_1000_run0.root",1000,24)'
root -l -b -q '../analysis/analyze_general.C("validation_full_muon_analysis_1000_run0.root",13,false,"primaryEmitter")'
```

This test asks a different question from the dedicated validators:

- dedicated validator: **is one process correct against theory/reference?**
- integrated analysis: **do the processes and bookkeeping remain coherent when used together?**

At the current development stage, `analyze_general.C` should still be treated as an analysis layer under audit. A misleading plot/selection there does not override a dedicated validator result.

## 7. What validation data can and cannot prove

A controlled validation ROOT file and a production ROOT file share the same schema and are generated by the same simulation code, but their event distributions differ.

Some relations should hold on any suitable sample, including production:

- valid unit directions;
- valid channel IDs;
- local readout phase-space conventions;
- per-Cherenkov-photon angle relation when all required truth fields are present;
- non-negative path length/deposited energy;
- `NDetected <= NIncident`.

Other validation expectations are sample-specific:

- `817 photons/cm` is specific to the tested particle/material/energy/table range;
- `PDE=0.400` is specific to the controlled 450-nm point;
- one fixed Fresnel reflection probability is specific to the tested indices, angle and polarization.

Therefore do not run a dedicated validator blindly on arbitrary production data unless its assumptions match that data.

## 8. Validation output policy

Use `fullAll` for small validation/debug samples so required truth is retained.

Recommended normal production defaults remain:

- `thinAll`;
- `photonBirthDetail basic`;
- `incidentDetail auto`;
- `boundaryDebug false`;
- `maxOpticalBounces 0` unless a documented safety study justifies a positive limit.

`BoundaryDebug` can generate an enormous number of rows and should only be enabled in small controlled boundary tests.

## 9. Reproducibility

For every important validation result, record:

- project revision;
- `RunConfig.SchemaVersion`;
- seed;
- Geant4 version;
- active physics processes;
- material dataset names and validity ranges;
- `OpticalDataFingerprint_FNV1a64`;
- thread count;
- macro used.

Sequential fixed-seed validation should be reproducible. Do not rely on MT output row order as a reproducibility identifier; join using IDs.

## 10. Freeze rule

A validation/core revision is ready to freeze only when:

1. static checks pass;
2. all applicable dedicated runtime validators pass;
3. failures are understood rather than hidden by looser tolerances;
4. material/model limitations are documented;
5. schema and material provenance are recorded;
6. a representative integrated run passes structural QA.

After freeze, changes to simulation physics, geometry semantics, optical datasets, ROOT truth schema or dedicated validators require a new revision and rerunning the affected validation subset.
