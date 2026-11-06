# Cherenkov detector simulation

A configurable Geant4 application for studying Cherenkov-light production, optical transport, and effective readout response in cuboid detectors.

The detector model supports a water or PMMA Core, an optional PMMA shell, a Tyvek-like reflector, PDMS optical coupling, readout elements on six faces, and two energy-deposition trigger paddles. Gun and CRY sources use the same simulation implementation.

## Scope and status

This repository provides simulation code, configurable examples, ROOT analysis macros, and dedicated checks of selected processes and output bookkeeping. It is intended for detector-design studies. Detector-specific measurements and comparison with experimental data are needed before interpreting its output as an absolute prediction for a particular instrument.

The current output schema is `general-v1.10`. The presence of validation tools does not establish that every configuration or dependency version has been tested. Record the source revision, software versions, seed, macro and test results for each reference run. See [Validation](README_VALIDATION.md).

## Readout modes

| Mode | Behaviour |
|---|---|
| `DirectSiPM` | Records incident photons and applies an effective wavelength-dependent PDE at the readout surface. |
| `ConcentratorInput` | Records entrance position, direction, polarization, wavelength and time for a downstream optical study; no SiPM detections are generated. |

A concentrator geometry and a complete ROBAST workflow are not included. The SiPM model covers effective PDE; it does not include a complete electronics/noise/saturation model.

## Dependencies and build

- CMake 3.22 or newer and a C++17 compiler.
- A Geant4 development installation; the default build requests UI and visualization components.
- ROOT for the supplied analysis and validation macros. The executable writes ROOT through Geant4 analysis and does not directly link a separate ROOT package in CMake.
- Python 3 for static checks.
- Network access for the current CMake rule to download and build CRY 1.7, including when using only Gun mode.

Set up the Geant4/ROOT environment for your installation, then run:

```bash
cmake -S . -B build
cmake --build build -j8
```

For a batch build without requesting Geant4 UI/visualization drivers:

```bash
cmake -S . -B build -DWITH_GEANT4_UIVIS=OFF
cmake --build build -j8
```

Keep the source and build directories available: optical and CRY data paths are compiled into the executable. The current install rule installs only the executable and is not a standalone relocatable distribution.

## First run

Run from `build/`; CMake exposes the source macros there as symbolic links:

```bash
cd build
./cherenkov example_gun_direct_sipm.mac 12345
root -l -b -q '../analysis/inspect_output.C("example_gun_direct_sipm_run0.root",100,24)'
```

The example generates 100 downward 3 GeV muons with 24 readout channels. Check the simulation log as well as the ROOT structural report. See [macro examples](macros/README.md) for CRY and concentrator-input configurations.

## Output and analysis

The output filename is `<configured-name>_run<RunID>.root`. Important trees include primary and track truth, photon histories, readout incidents/detections, trigger deposits, channel maps and run metadata.

| Event mode | Retained events |
|---|---|
| `fullAll` | Full detail for all events; useful for small validation samples. |
| `thinAll` | Thin event/primary/deposit information for all events; detailed data only for triggered events. |
| `triggered` | Triggered events only; generated-event totals remain in `RunStatistics`. |

Detailed photon distributions in `thinAll` and `triggered` are trigger-selected. Use appropriate denominators and state the selection when reporting efficiencies. `ChannelEvent` is sparse: absent event/channel rows do not imply a missing detector channel.

A general report can be produced from `build/`:

```bash
root -l -b -q '../analysis/analyze_general.C("example_gun_direct_sipm_run0.root",13,true,"primaryEmitter")'
```

This analyzer is a starting point for analysis and has known limitations listed in [Validation](README_VALIDATION.md). Inspect `RunConfig` before interpreting any result.

## Physical input assumptions

- Water optical properties use the datasets identified in `RunConfig` and the comments in `data/`.
- The PMMA absorption profile distinguishes the measured 300–500 nm domain from model extensions outside it. Its broad-band refractive index is an approximation.
- PDMS absorption length is a generic 10 m assumption.
- Tyvek reflectivity and SiPM PDE curves are effective approximations, not complete measured models of a particular assembled detector.
- Scintillation and WLS are disabled by default. Trigger paddles are evaluated through deposited energy; calibrated scintillation/light-readout models are not supplied.

## Documentation

- [Project guide](README_PROJECT_GUIDE.md): architecture and configuration workflow.
- [General model](README_GENERAL.md): physics, sources, readout and analysis semantics.
- [Validation guide](README_VALIDATION.md): available tests and limitations.
- [ROOT branch reference](README_TREE_INFO.md): branches, units and joins.
- [Schema](SCHEMA.md): output contract and selection semantics.
- [Macro catalog](macros/README.md): examples and controlled test configurations.
