# Validation and known limitations

## What a check establishes

Static checks inspect source/schema/configuration consistency. Runtime validators read controlled Geant4 samples and test selected observables. Neither establishes that all optical inputs describe a particular real detector. General cosmic production samples do not necessarily satisfy the assumptions of monoenergetic validation samples.

Keep validation evidence with its source revision, Geant4/ROOT/compiler versions, seed, macro, thread count and log. This document lists the tools included in the repository; it does not certify that they have all passed on your machine or current checkout.

## Static checks

From the source root:

```bash
bash tools/run_framework_static_checks.sh
```

This checks ntuple indices/types, recognized analysis branch names, registered macro commands, configuration-to-metadata references and macro organization, and runs independent optical-table calculations. These checks do not compile C++, verify every analysis binding type or validate numerical physics against Geant4. `validate_physics.py` prints calculations; its successful exit alone is not a quantitative regression acceptance test.

## Build and first runtime check

Build as described in [README](README.md). From `build/`:

```bash
./cherenkov validation_output_mt.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_output_mt_run0.root",100,0)'
root -l -b -q '../analysis/validate_framework_metadata.C("validation_output_mt_run0.root")'
```

Resolve simulation warnings, missing trees, count discrepancies and metadata failures before interpreting physics plots.

## Available runtime samples

Run each configuration with `./cherenkov <macro> 12345` from `build/`, then invoke its ROOT validator from `../analysis/`.

| Configuration | Validator | Scope |
|---|---|---|
| `validation_output_mt.mac` | `inspect_output.C`, `validate_framework_metadata.C` | ROOT structure, worker counts and metadata consistency. |
| `validation_trigger_muon.mac` | `inspect_output.C` | Through-going trigger sample; additionally inspect passed/generated counts. Structural PASS does not enforce a trigger-efficiency target. |
| `validation_dedx_muon.mac` | `validate_dedx.C`, `validate_mcs_muon.C` | Path/deposition/loss reporting and approximate multiple-scattering scale. |
| `validation_muon_decay.mac` | `validate_muon_decay.C` | Stopped positive-muon lifetime and positron endpoint. |
| `validation_frank_tamm_water.mac` | `validate_frank_tamm.C` | Generated yield; angle residual is plotted/reported. |
| `validation_optical_absorption_pmma.mac` and `validation_optical_rayleigh_water.mac` | `validate_optical_bulk.C` | Supplied absorption/scattering lengths in controlled samples. |
| `validation_fresnel_tir_pmma_air.mac` | `validate_fresnel_tir.C` | Selected reflection/refraction/polarization/TIR cases. |
| `validation_tyvek_reflectivity_450nm.mac` | `validate_tyvek_reflectivity.C` | Configured first-hit total reflectivity at 450 nm. |
| `validation_pde_450nm.mac` | `validate_pde.C` | Configured effective PDE at 450 nm. |
| `validation_flexible_layout.mac` | `validate_flexible_layout.C` | Expected asymmetric channel counts and unit normal magnitudes. |
| `validation_core_boundary_modes.mac` | `validate_core_boundary_modes.C` | Recorder modes, polarization detail and TIR exclusion. |
| `validation_core_boundary_first_exit_6faces.mac` | `validate_core_boundary_first_exit_6faces.C` | Controlled first-exit face/local-coordinate cases. |
| `validation_electron_brems.mac` | `validate_electron_brems.C` | Process-activity check; see known issue below. |
| `validation_final_integrated_freeze_1000.mac` | `inspect_output.C`, `validate_framework_metadata.C`, `analyze_general.C` | Combined transport/output sample, 1000 events and 24 channels. Filename is a historical label, not independent certification. |

For example:

```bash
./cherenkov validation_frank_tamm_water.mac 12345
root -l -b -q '../analysis/validate_frank_tamm.C("validation_frank_tamm_water_run0.root")'
./cherenkov validation_final_integrated_freeze_1000.mac 12345
root -l -b -q '../analysis/inspect_output.C("validation_final_integrated_freeze_1000_run0.root",1000,24)'
```

Two input files are needed for `validate_optical_bulk.C()`. The Fresnel and CoreBoundary validators consume several run files produced by their respective macros. Preserve those names or pass the supported filename/prefix arguments.

## Known implementation and analysis issues

These observations apply to the public source reviewed at `1be2c97b871cdf8da10568113a7baebb135f123d`. Recheck them if your local source differs. This documentation update does not modify simulation or analysis code.

- `validate_electron_brems.C` selects VolumeID 1 for its primary-Core statistics, whereas Core is 0 and the paired macro disables the shell. Its PASS does not require valid Core statistics. Correct this before relying on that part of the report.
- The ROOT angle analyzers interpolate refractive index linearly in wavelength; the simulation tables and independent Python calculation use photon-energy interpolation. Align these before interpreting an exact angle residual. The Frank–Tamm PASS currently tests yield, not the angle residual.
- General-analysis trigger percentage uses retained TriggerEvent rows. In triggered mode this describes retained rows, not generated-event trigger efficiency. Use summed RunStatistics counters for the latter.
- General-analysis efficiency histograms do not provide calculated statistical uncertainties. Add suitable uncertainties for quantitative comparisons.
- The general analyzer loads tables through `../data/` and assumes one run/file identity. Run from `build/`; combining outputs from separate processes requires a sample/file key.
- Many validators print FAIL or return early without setting a failing process exit status. The executable ignores the macro ApplyCommand return code. Inspect logs as well as shell/LSF job status.
- Many times are stored as Float nanoseconds; absolute timestamps can lose fine timing precision.

The dedicated dE/dx macro reports distributions without an automatic reference acceptance criterion. The multiple-scattering comparison is an approximate scale test. A PDE or Tyvek test at one wavelength does not validate its full spectrum or angular response.

## Physical model limitations

PMMA outside its measured absorption domain uses explicit extensions. PDMS attenuation is a generic assumption. Tyvek uses effective surface properties and SiPM uses effective PDE. Detector-specific optical measurements, systematic variations and experimental comparisons are needed for absolute detector-performance claims.

Use fullAll for small unbiased optical validation samples. Leave boundary debugging off for ordinary production; it can create very large files. A finite maxOpticalBounces is an artificial transport cut whose impact must be assessed.
