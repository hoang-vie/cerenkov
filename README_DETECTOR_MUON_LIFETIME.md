# Detector-geometry mu+ lifetime validation

This is **Stage 2** after the detector stopping-energy scan.

The coarse scan produced:

```text
30 MeV -> 2000/2000 primary mu+ decayed in the Water Core
mean EndY ~= +16.6 mm
```

Therefore the baseline lifetime run uses **30 MeV mu+** in the same detector geometry.

## What this test measures

It produces two lifetime observables from the same ROOT file.

### 1. Truth proper-time lifetime

Select:

```text
ParentID == 0
PDG == -13
EndVolumeID == Core
EndProcessName contains "Decay"
```

and histogram:

```text
Track.EndProperTime_ns / 1000
```

Then fit:

```text
A * exp(-t / tau)
```

over `0.1--10 us`.

### 2. Detector-like delayed timing

For each selected primary mu+, find a direct daughter positron:

```text
PDG == -11
ParentID == primary mu+ TrackID
```

and form:

```text
Delta t = positron VertexTime_ns - muon VertexTime_ns
```

This is not yet an electronics/readout reconstruction. It is a **geometry-level delayed
secondary observable**. It avoids reading `EndProperTime_ns` directly and should recover
the same decay constant because the positron is born at the muon decay time.

The validator also plots:

- Michel positron kinetic-energy spectrum;
- primary mu+ decay-position `EndY_mm` in the Water Core.

## Run

From the project root:

```bash
patch -p1 < add_detector_muon_lifetime_validation.patch
cd build
cmake ..
make -j8
./cherenkov validation_detector_muon_lifetime.mac 12345
root -l -b -q '../analysis/validate_detector_muon_lifetime.C()'
```

Expected ROOT output:

```text
validation_detector_muon_lifetime_run0.root
```

Expected plot:

```text
validation_plots/detector_muon_lifetime.png
validation_plots/detector_muon_lifetime.pdf
```

Expected physics:

```text
tau(mu+) ~ 2.197 us
Michel positron kinetic endpoint ~52.3 MeV
Core-decay fraction >98% for this selected 30 MeV baseline
```

## Important limitation

`Delta t = positron birth time - muon birth time` is **not yet a full detector
measurement**. It is the next intermediate level:

```text
process truth
    -> real detector geometry
    -> stopped mu+
    -> delayed Michel daughter timing
```

A later stage can replace daughter-birth time with actual scintillator/SiPM timing
after those signal models are declared ready for this validation.
