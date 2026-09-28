# Detector-level stopped-muon validation — stage 1

This patch adds the **first detector-level muon-lifetime prerequisite**:
find a controlled `mu+` kinetic-energy window that stops and decays inside
the current Water Core while keeping the current detector geometry.

It intentionally does **not** replace the existing simple `validation_muon_decay.mac`.
That simple test remains the process-level check of Geant4 decay/proper-time bookkeeping.

## Added files

- `macros/validation_detector_muon_stopping_scan.mac`
- `analysis/validate_detector_muon_stopping_scan.C`

## Run

From source root, reconfigure so the macro is copied:

```bash
cd build
cmake ..
make -j8
./cherenkov validation_detector_muon_stopping_scan.mac 12345
root -l -b -q '../analysis/validate_detector_muon_stopping_scan.C()'
```

Expected outputs:

```text
validation_detector_muon_stopping_scan_run0.root
...
validation_detector_muon_stopping_scan_run9.root

validation_plots/detector_muon_stopping_scan.png
validation_plots/detector_muon_stopping_scan.pdf
```

Energy map:

```text
run0  = 10 MeV
run1  = 20 MeV
run2  = 30 MeV
run3  = 40 MeV
run4  = 50 MeV
run5  = 60 MeV
run6  = 80 MeV
run7  = 100 MeV
run8  = 150 MeV
run9  = 200 MeV
```

The validator selects `ParentID==0 && PDG==-13`, resolves the Core ID from
`VolumeMap`, and counts events where the primary `mu+` ends by `Decay`
inside the Core.

**Do not call the best coarse point the final lifetime configuration.**
Use it as the centre of a finer scan. The next stage is then a high-statistics
detector-level lifetime macro that compares:

1. `Track.EndProperTime_ns` (truth lifetime),
2. delayed Michel-positron activity relative to the prompt muon time,
3. later, reconstructed detector timing if the signal model is enabled.
