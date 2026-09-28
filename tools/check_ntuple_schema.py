#!/usr/bin/env python3
"""Static consistency check for G4AnalysisManager ntuple IDs and fill-column types.
No Geant4 or ROOT installation is required.
"""
from pathlib import Path
import collections
import hashlib
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
run = (ROOT / "src" / "RunAction.cc").read_text()
event = (ROOT / "src" / "EventAction.cc").read_text()
header = (ROOT / "include" / "RunAction.hh").read_text()

enums = {name: int(value) for name, value in re.findall(r"(k\w+)\s*=\s*(\d+)", header)}
blocks = []
for actual_id, m in enumerate(re.finditer(
        r'am->CreateNtuple\("([^"]+)".*?am->FinishNtuple\((k\w+)\);', run, re.S)):
    name, enum_name = m.group(1), m.group(2)
    cols = [(c.group(1), c.group(2)) for c in re.finditer(
        r'am->CreateNtuple([IFDS])Column\("([^"]+)"\)', m.group(0))]
    expected = enums.get(enum_name)
    if expected != actual_id:
        raise SystemExit(
            f"Ntuple ID mismatch: {name} is created as {actual_id}, "
            f"but {enum_name}={expected}")
    blocks.append((enum_name, name, cols))

# general-v1.10 is intentionally additive for RunConfig: the complete v1.6,
# v1.8 and v1.9 prefixes must remain identical in name, type and order.
runconfig_cols = next(cols for enum_name, name, cols in blocks if name == "RunConfig")
legacy_prefix = "\n".join(f"{t}:{n}" for t, n in runconfig_cols[:99]).encode()
legacy_prefix_sha256 = hashlib.sha256(legacy_prefix).hexdigest()
expected_v16_prefix_sha256 = "8b2db2c5579b3d711e2f5ae3f97192d0f181a58bd2b990ce3f37ed6a6f16e204"
if legacy_prefix_sha256 != expected_v16_prefix_sha256:
    raise SystemExit(
        "RunConfig legacy v1.6 prefix (columns 0..98) changed unexpectedly: "
        f"sha256={legacy_prefix_sha256}")

v18_prefix = "\n".join(f"{t}:{n}" for t, n in runconfig_cols[:154]).encode()
v18_prefix_sha256 = hashlib.sha256(v18_prefix).hexdigest()
expected_v18_prefix_sha256 = "54a63ed8d1155c8b5ddfbdbed07d2114928fa822bb8be8b2477ef95fb7f7401b"
if v18_prefix_sha256 != expected_v18_prefix_sha256:
    raise SystemExit(
        "RunConfig legacy v1.8 prefix (columns 0..153) changed unexpectedly: "
        f"sha256={v18_prefix_sha256}")
v19_prefix = "\n".join(f"{t}:{n}" for t, n in runconfig_cols[:157]).encode()
v19_prefix_sha256 = hashlib.sha256(v19_prefix).hexdigest()
expected_v19_prefix_sha256 = "53fc0c5a66f3aef0e26f26488ab0803dc50beb36678cbd5e4a8e80dbd389f743"
if v19_prefix_sha256 != expected_v19_prefix_sha256:
    raise SystemExit(
        "RunConfig legacy v1.9 prefix (columns 0..156) changed unexpectedly: "
        f"sha256={v19_prefix_sha256}")
if len(runconfig_cols) != 160:
    raise SystemExit(
        f"general-v1.10 RunConfig must contain 160 columns, found {len(runconfig_cols)}")

for source_name, text in (("EventAction.cc", event), ("RunAction.cc", run)):
    fills = collections.defaultdict(list)
    for m in re.finditer(
        r'am->FillNtuple([IFDS])Column\(\s*(k\w+)\s*,\s*(\d+)\s*,',
        text):
        fills[m.group(2)].append((int(m.group(3)), m.group(1)))

    for enum_name, ntuple_name, cols in blocks:
        if enum_name not in fills:
            continue
        seen = collections.defaultdict(set)
        for idx, type_code in fills[enum_name]:
            seen[idx].add(type_code)
        missing = [idx for idx in range(len(cols)) if idx not in seen]
        errors = []
        for idx, used_types in seen.items():
            if idx >= len(cols):
                errors.append(f"column {idx} is out of range")
            elif cols[idx][0] not in used_types:
                errors.append(
                    f"column {idx} ({cols[idx][1]}) created as {cols[idx][0]} "
                    f"but filled as {sorted(used_types)}")
        if missing or errors:
            raise SystemExit(
                f"{source_name}: {ntuple_name}: missing={missing}, errors={errors}")


# Catch accidental immediately repeated output/allocation statements that can be
# semantically harmless (duplicate Fill) or expensive (duplicate allocation).
for source in (ROOT / "src").glob("*.cc"):
    lines = [ln.strip() for ln in source.read_text().splitlines() if ln.strip()]
    for a, b in zip(lines, lines[1:]):
        if a == b and ("FillNtuple" in a or "new PhotonInfo" in a):
            raise SystemExit(f"{source.name}: immediate duplicate statement: {a}")

print("Ntuple creation order and fill-column type/index checks: OK")
for enum_name, name, cols in blocks:
    print(f"  {enums[enum_name]:2d}  {name:22s}  {len(cols):2d} columns")
