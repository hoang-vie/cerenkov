#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "=== Framework static checks ==="
python3 tools/check_ntuple_schema.py
python3 tools/check_analysis_schema.py
python3 tools/check_macros.py
python3 tools/check_config_metadata_contract.py
python3 tools/check_macro_catalog.py
python3 validate_physics.py
python3 tools/cry_magnetic_latitude.py

echo "=== Whitespace / patch sanity ==="
git diff --check

echo "FRAMEWORK STATIC CHECKS: PASS"
