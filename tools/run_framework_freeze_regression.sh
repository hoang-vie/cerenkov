#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${1:-$ROOT_DIR/build}"
JOBS="${JOBS:-4}"
SEED_BASE="${SEED_BASE:-26000}"
LOG_DIR="$BUILD_DIR/framework_freeze_logs"

is_dirty=0
if ! git -C "$ROOT_DIR" diff --quiet; then is_dirty=1; fi
if ! git -C "$ROOT_DIR" diff --cached --quiet; then is_dirty=1; fi
if [[ -n "$(git -C "$ROOT_DIR" ls-files --others --exclude-standard)" ]]; then is_dirty=1; fi
if [[ "$is_dirty" -ne 0 ]]; then
  echo "ERROR: freeze regression requires a clean Git working tree." >&2
  echo "Commit/move intentional source files first, then rerun CMake so GitDirty can be 0." >&2
  git -C "$ROOT_DIR" status --short >&2 || true
  exit 2
fi

"$ROOT_DIR/tools/run_framework_static_checks.sh"

# Reconfigure on the exact clean commit so build-time Git provenance is current.
cmake -S "$ROOT_DIR" -B "$BUILD_DIR"
cmake --build "$BUILD_DIR" -j"$JOBS"
mkdir -p "$LOG_DIR"
cd "$BUILD_DIR"

run_root_check() {
  local label="$1"
  local expression="$2"
  local required="$3"
  local log="$LOG_DIR/${label}.log"
  echo
  echo "=== $label ==="
  root -l -b -q "$expression" 2>&1 | tee "$log"
  if ! grep -Fq "$required" "$log"; then
    echo "ERROR: $label did not report: $required" >&2
    exit 3
  fi
}

run_sim() {
  local label="$1"
  local macro="$2"
  local seed="$3"
  local log="$LOG_DIR/${label}.log"
  echo
  echo "=== $label ==="
  ./cherenkov "$macro" "$seed" 2>&1 | tee "$log"
}

# 1) Public deterministic example: ordinary user-facing DirectSiPM chain.
run_sim "example_gun_direct_sipm" "example_gun_direct_sipm.mac" "$((SEED_BASE+1))"
run_root_check "example_gun_inspect" \
  '../analysis/inspect_output.C("example_gun_direct_sipm_run0.root",100,24)' \
  'OUTPUT SANITY: PASS'
run_root_check "example_gun_metadata" \
  '../analysis/validate_framework_metadata.C("example_gun_direct_sipm_run0.root")' \
  'FRAMEWORK METADATA VALIDATION: PASS'

# 2) Per-face flexible layout contract.
run_sim "flexible_layout" "validation_flexible_layout.mac" "$((SEED_BASE+2))"
run_root_check "flexible_layout_validator" \
  '../analysis/validate_flexible_layout.C("validation_flexible_layout_run0.root")' \
  'FLEXIBLE LAYOUT VALIDATION: PASS'
run_root_check "flexible_layout_inspect" \
  '../analysis/inspect_output.C("validation_flexible_layout_run0.root",10,12)' \
  'OUTPUT SANITY: PASS'

# 3) CoreBoundary OFF/basic/phaseSpace/TIR semantics.
run_sim "core_boundary_modes" "validation_core_boundary_modes.mac" "$((SEED_BASE+3))"
run_root_check "core_boundary_modes_validator" \
  '../analysis/validate_core_boundary_modes.C()' \
  'CORE-BOUNDARY MODE/DETAIL/TIR VALIDATION: PASS'

# 4) Deterministic six-face first-exit geometry/local-coordinate contract.
run_sim "core_boundary_6faces" "validation_core_boundary_first_exit_6faces.mac" "$((SEED_BASE+4))"
run_root_check "core_boundary_6faces_validator" \
  '../analysis/validate_core_boundary_first_exit_6faces.C()' \
  'CORE-BOUNDARY FIRST-EXIT 6-FACE VALIDATION: PASS'

# 5) Full normal detector chain, 24 channels, MT, 1000 events.
run_sim "final_integrated" "validation_final_integrated_freeze_1000.mac" "$((SEED_BASE+5))"
run_root_check "final_integrated_inspect" \
  '../analysis/inspect_output.C("validation_final_integrated_freeze_1000_run0.root",1000,24)' \
  'OUTPUT SANITY: PASS'
run_root_check "final_integrated_metadata" \
  '../analysis/validate_framework_metadata.C("validation_final_integrated_freeze_1000_run0.root")' \
  'FRAMEWORK METADATA VALIDATION: PASS'

# Freeze provenance must describe the clean commit used to build.
head_commit="$(git -C "$ROOT_DIR" rev-parse HEAD)"
meta_log="$LOG_DIR/final_integrated_metadata.log"
if ! grep -Fq "GitCommit     : $head_commit" "$meta_log"; then
  echo "ERROR: final ROOT GitCommit does not match repository HEAD $head_commit" >&2
  exit 4
fi
if ! grep -Fq 'GitDirty      : 0' "$meta_log"; then
  echo "ERROR: final ROOT was not produced from a clean build provenance state" >&2
  exit 4
fi

echo
echo "FRAMEWORK FREEZE REGRESSION: PASS"
echo "Commit: $head_commit"
echo "Logs  : $LOG_DIR"
