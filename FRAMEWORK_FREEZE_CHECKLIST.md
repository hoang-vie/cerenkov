# Framework V1 freeze checklist

This checklist freezes the **simulation framework**, not an optimized detector configuration. A release candidate is acceptable when the configuration surface, metadata/provenance, physics baseline, output schema, user examples, and regression gates are mutually consistent.

## 1. Prepare one clean release-candidate commit

Before the freeze regression, decide which source-level study utilities are intentionally part of the repository. Commit intended C++/Python/ROOT macros, documentation, and `.mac` files; move or discard scratch files. ROOT outputs and the build directory remain generated artifacts.

The freeze regression intentionally refuses to run from a dirty Git tree because `RunConfig.GitDirty=0` is part of the release provenance gate.

## 2. Run the static contract

From the repository root:

```bash
tools/run_framework_static_checks.sh
```

Required final line:

```text
FRAMEWORK STATIC CHECKS: PASS
```

This covers ntuple/schema consistency, analysis/schema compatibility, the 51-command UI inventory, command-to-metadata contract, macro catalog, physics-table sanity, CRY magnetic-latitude helper, and `git diff --check`.

## 3. Run the clean-build freeze regression

From the repository root, with Geant4/ROOT/CRY/ROBAST environment available:

```bash
tools/run_framework_freeze_regression.sh
```

Optional environment variables:

```bash
JOBS=8 SEED_BASE=26000 tools/run_framework_freeze_regression.sh
```

An alternate build directory may be given as the first argument.

The script re-runs CMake on the exact clean commit, builds, and validates:

1. the public deterministic Gun + DirectSiPM example;
2. asymmetric per-face readout layout;
3. CoreBoundary OFF/basic/phaseSpace/TIR semantics;
4. deterministic six-face first-exit geometry/local coordinates;
5. the 1000-event, 24-channel, MT integrated detector chain;
6. ROOT build provenance against the exact Git `HEAD`, including `GitDirty=0`.

Required final line:

```text
FRAMEWORK FREEZE REGRESSION: PASS
```

Logs are retained under `build/framework_freeze_logs/`.

## 4. Inspect the release-candidate metadata

The final integrated file is:

```text
build/validation_final_integrated_freeze_1000_run0.root
```

Its `RunConfig` must report:

- `SchemaVersion = general-v1.10`;
- `ChannelCount = 24`;
- `SensorMode = DirectSiPM`;
- `GitCommit` equal to `git rev-parse HEAD`;
- `GitDirty = 0`.

The exact event/optical counts are sample-dependent and are not themselves the API contract. Physics-specific deterministic checks remain in the validation suite/runbook.

## 5. Freeze/tag only after both gates pass

Recommended sequence:

```bash
git status --short
# must be empty

git tag -a framework_v1 -m "Freeze validated Cherenkov simulation framework v1 (general-v1.10)"
git show --stat framework_v1
```

`physics_validated_v1` remains the historical frozen physics baseline; `framework_v1` marks the later reusable framework/API/provenance/documentation freeze.

## 6. Change policy after freeze

For a new simulation-affecting feature or command, require all of:

```text
configuration/API
    -> metadata/provenance
    -> validation/regression
    -> documentation/example when user-facing
```

Append schema fields rather than reordering existing columns when possible. A physics/material-model change must be treated separately from an interface/documentation-only change and revalidated accordingly.
