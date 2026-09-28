# Validation runbook — general-v1.10 framework baseline

This order is intentional: validate ROOT/MT output first, then non-optical charged-particle physics, then optical production, sensor response, flexible geometry, and finally CRY.

Run from the project source root for static checks and from `build/` for Geant4 runs.

## 0. Rebuild and static checks

```bash
python3 tools/check_ntuple_schema.py
python3 tools/check_analysis_schema.py
python3 tools/check_macros.py
python3 tools/check_config_metadata_contract.py
python3 tools/check_macro_catalog.py
python3 validate_physics.py
python3 tools/cry_magnetic_latitude.py

cd build
cmake ..
make -j8
```

`cmake ..` is important after adding/removing macros because CMake copies `macros/*.mac` into `build/`.

## 1. MT output/merging smoke test — run this first

This specifically tests the bug that previously produced `Analysis_W001: Ntuple 11 does not exist` at end-of-run.

```bash
./cherenkov validation_output_mt.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_output_mt_run0.root",100,0)'
root -l -b -q '../analysis/validate_framework_metadata.C("validation_output_mt_run0.root")'
```

Expected:
- no `Analysis_W001`;
- `GeneratedEvents` summed across `RunStatistics` rows = 100;
- `TriggerEvent` has 100 rows because the macro uses `fullAll`;
- `ChannelCount = 0`;
- final line from `inspect_output.C`: `OUTPUT SANITY: PASS`;
- final line from `validate_framework_metadata.C`: `FRAMEWORK METADATA VALIDATION: PASS`, with `SchemaVersion=general-v1.10` and non-empty build/source provenance fields.

If this test fails, do not interpret any later ROOT physics plots yet.

## 2. Muon dE/dx / stopping-power bookkeeping

```bash
./cherenkov validation_dedx_muon.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_dedx_muon_run0.root",10000,0)'
root -l -b -q '../analysis/validate_dedx.C("validation_dedx_muon_run0.root",13,true)'
root -l -b -q '../analysis/analyze_general.C("validation_dedx_muon_run0.root",13,true,"primaryEmitter")'
```

The dedicated validator prints both:
- `Edep/L`: local deposited-energy observable;
- `DeltaKE/L`: kinetic-energy loss of the primary per Core path, the preferred quantity for a clean stopping-power comparison.

The 10-cm Core path should be obvious in the reported path-length distribution for the straight central gun geometry. Do not treat `Edep/L` and `DeltaKE/L` as identical definitions.

## 3. Muon free lifetime / Michel decay validation

```bash
./cherenkov validation_muon_decay.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_muon_decay_run0.root",100000,0)'
root -l -b -q '../analysis/validate_muon_decay.C("validation_muon_decay_run0.root")'
root -l -b -q '../analysis/analyze_general.C("validation_muon_decay_run0.root",-13,false,"primaryEmitter")'
```

This macro uses stopped `mu+` in a large water volume. Use `mu+` for the free-lifetime validation; negative muons in matter have capture physics in addition to decay.

## 4. Frank–Tamm Cherenkov yield in water

First print the independent expectation from the actual `Water_RINDEX.txt` table:

```bash
cd ..
python3 validate_physics.py
cd build
```

Then run Geant4:

```bash
./cherenkov validation_frank_tamm_water.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_frank_tamm_water_run0.root",100,0)'
root -l -b -q '../analysis/validate_frank_tamm.C("validation_frank_tamm_water_run0.root")'
```

`validate_frank_tamm.C` counts only Core-origin photons with a finite Cherenkov emission angle, so future scintillation/WLS photons are not silently included. Compare the printed mean photons/cm with `validate_physics.py` for the 3-GeV kinetic-muon beta and the same 200–900 nm optical table range.

## 5. DirectSiPM PDE at 450 nm

```bash
./cherenkov validation_pde_450nm.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_pde_450nm_run0.root",100000,1)'
root -l -b -q '../analysis/validate_pde.C("validation_pde_450nm_run0.root",0.400)'
```

Expected structure:
- one +Z readout channel;
- approximately 100000 `SensorIncident` rows;
- detection fraction near the 450-nm PDE table value, 0.400;
- no need for CRY or scintillator trigger to validate this boundary model because the macro uses `fullAll`.

## 6. Flexible readout layout

The MT application intentionally does not change geometry at Idle. Instead validate an arbitrary asymmetric layout in one geometry build:

```bash
./cherenkov validation_flexible_layout.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_flexible_layout_run0.root",10,12)'
```

Expected `ChannelMap` / `RunConfig.ChannelCount = 12` from the configured per-face counts:
`1 + 4 + 3 + 0 + 2 + 2`.

For geometry scans in normal MT production, start a new process/configuration for each geometry rather than changing the geometry between multi-event MT runs.

## 7. CRY -> ConcentratorInput smoke test

```bash
./cherenkov validation_cry_concentrator.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_cry_concentrator_run0.root",100,24)'
root -l -b -q '../analysis/analyze_general.C("validation_cry_concentrator_run0.root",13,true,"primaryEmitter")'
```

Expected:
- 24 aperture channels (`2 x 2 x 6`);
- `DetectedSignal = 0` by construction in `ConcentratorInput`;
- `SensorIncident` may be non-zero and contains phase space in automatic incident-detail mode;
- `PrimaryTruth`/`TriggerEvent` remain available for all events because the macro uses `thinAll`.

## 8. CRY -> DirectSiPM smoke test

```bash
./cherenkov validation_cry_direct_sipm.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_cry_direct_sipm_run0.root",100,24)'
root -l -b -q '../analysis/analyze_general.C("validation_cry_direct_sipm_run0.root",13,true,"primaryEmitter")'
```

`DetectedSignal` can be non-zero. A 100-event CRY smoke sample is too small to interpret zero detections as a physics failure by itself; increase statistics only after all structural checks pass.

## 9. Species-general reports

The same analysis can be rerun without regenerating the file when the species exists in the stored detailed sample:

```bash
# directly emitted by primary electrons
root -l -b -q '../analysis/analyze_general.C("FILE.root",11,true,"primaryEmitter")'

# directly emitted by pi+
root -l -b -q '../analysis/analyze_general.C("FILE.root",211,false,"primaryEmitter")'

# all photons whose propagated non-optical ancestry belongs to a primary muon cascade
root -l -b -q '../analysis/analyze_general.C("FILE.root",13,true,"ancestor")'
```

## 10. Only then start production

For a new study, copy one of the catalogued `example_*.mac` files described in `macros/README.md`. Do not start from a `validation_*.mac`: validation macros may intentionally use artificial geometry, process activation, optical primaries, or output settings. `study_*.mac` files are retained downstream use cases and likewise are not framework defaults.


Recommended production defaults:
- `thinAll`;
- `photonBirthDetail basic`;
- `incidentDetail auto`;
- `boundaryDebug false`;
- `maxOpticalBounces 0` unless a documented safety study justifies a positive limit.

Keep `fullAll` for controlled validation/debug samples because optical rows dominate file size.

## Audit-1 corrections (hotfix5_audit1)

Before optical validation, run the dedicated trigger test:

```bash
./cherenkov validation_trigger_muon.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_trigger_muon_run0.root",1000,0)'
```

The gun uses `useAutoY true`, so the straight 3 GeV muon starts above the top paddle and crosses top scintillator, radiator and bottom scintillator. With the 0.1 MeV threshold and 100 ns coincidence window, the trigger fraction should be close to unity; a large inefficiency is then a real trigger/geometry issue rather than the source-position artifact in `validation_dedx_muon.mac`.

`validation_dedx_muon.mac` intentionally starts at y=60 mm to isolate the 10 cm radiator path; its trigger fraction is therefore not a trigger validation observable.

## Hotfix6 physics-audit runtime suite

After the existing static checks and a clean rebuild, run the deterministic
suite **before CRY production**.  All commands below assume the build directory
and that CMake has copied the current macros there.

```bash
./cherenkov validation_phase_space_6faces.mac 12345
root -l -b -q '../analysis/validate_phase_space_6faces.C()'

./cherenkov validation_optical_absorption_pmma.mac 12345
./cherenkov validation_optical_rayleigh_water.mac 12345
root -l -b -q '../analysis/validate_optical_bulk.C()'

./cherenkov validation_fresnel_tir_pmma_air.mac 12345
root -l -b -q '../analysis/validate_fresnel_tir.C()'

./cherenkov validation_tyvek_reflectivity_450nm.mac 12345
root -l -b -q '../analysis/validate_tyvek_reflectivity.C()'
```

Expected final lines:

- `PHASE-SPACE 6-FACE VALIDATION: PASS`
- `OPTICAL BULK PROCESS VALIDATION: PASS`
- `FRESNEL/SNELL/TIR VALIDATION: PASS`
- `TYVEK REFLECTIVITY VALIDATION: PASS`

Interpretation is deliberately narrow:

- PMMA absorption test proves Geant4 obeys the supplied measured SNO 300-nm
  attenuation point; it does not validate extrapolation below 300 nm.
- Water Rayleigh test now exercises the production Zhang-Hu PMH pure-water
  mean-free-path table. The mean-free-path normalization/spectrum are validated;
  the residual `delta=0.039` vs Geant4 ideal-dipole angular-law difference is a
  documented small systematic for long-path/polarization studies.
- Tyvek test validates first-hit total reflectivity for the effective surface;
  it does not validate a real Tyvek BRDF.

If any test fails, stop before CRY/production and inspect that subsystem first.

## Validation plots add-on

The dedicated ROOT validators now save visual summaries under `build/validation_plots/` as both PNG and PDF.  The flexible-layout geometry check additionally uses:

```bash
./cherenkov validation_flexible_layout.mac 12345
root -l -b -q '../analysis/validate_flexible_layout.C()'
```

The general analysis macro already writes a multipage PDF report and remains the preferred high-level visualization for CRY/production samples.

## CB1 Core-boundary first-exit recorder

`CoreBoundaryCrossing` is an observer-only, optional checkpoint at the Core
logical-volume boundary.  It is OFF by default and records at most one row per
optical photon in `firstExit` mode.  A row is created only for a successful
Core -> non-Core optical transmission (`FresnelRefraction` or `Transmission`),
not for a reflection or TIR wall encounter.

After the static checks and a clean rebuild, run both deterministic validators:

```bash
./cherenkov validation_core_boundary_first_exit_6faces.mac 91001
root -l -b -q '../analysis/validate_core_boundary_first_exit_6faces.C()'

./cherenkov validation_core_boundary_modes.mac 91002
root -l -b -q '../analysis/validate_core_boundary_modes.C()'
```

Expected final lines:

- `CORE-BOUNDARY FIRST-EXIT 6-FACE VALIDATION: PASS`
- `CORE-BOUNDARY MODE/DETAIL/TIR VALIDATION: PASS`

The first validator checks all six face IDs, the shared readout-face U/V
convention, post-boundary direction, wavelength/time, polarization, metadata,
and the one-row-per-photon rule.  The second checks three semantic gates:
`mode=off` produces zero rows; `basic` keeps the polarization branches but
stores NaN; `phaseSpace` stores finite polarization; and a TIR-only Water->Air
trajectory may set `Photon.ReachedCoreWall` while producing zero
`CoreBoundaryCrossing` rows.

For an unbiased escape-efficiency measurement, use `fullAll`; this tree follows
the same detailed-event retention policy as the other detailed optical trees.
