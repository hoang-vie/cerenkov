#!/usr/bin/env python3

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

NOMINAL = DATA / \
    "PMMA_SNO_UVT_SampleA_ABS_Measured300to500_Model200to900_Candidate.txt"

UV_OPAQUE = DATA / \
    "PMMA_SNO_UVT_SampleA_ABS_UVOpaqueBelow300_Model200to900_Candidate.txt"

VISIBLE_TRANSPARENT = DATA / \
    "PMMA_SNO_UVT_SampleA_ABS_VisibleTransparentAbove500_Model200to900_Candidate.txt"

# Effectively opaque UV bound.
UV_OPAQUE_L_MM = 1.0e-3

# alpha = 1e-4 cm^-1 -> L = 10 / alpha = 100000 mm = 100 m.
VISIBLE_ALPHA_CM = 1.0e-4
VISIBLE_L_MM = 10.0 / VISIBLE_ALPHA_CM


def read_table(path):
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        cols = line.split()
        try:
            wl = float(cols[0])
            L = float(cols[1])
        except (ValueError, IndexError):
            continue

        rows.append((wl, L))

    rows.sort()
    return rows


nominal = read_table(NOMINAL)

if not nominal:
    raise SystemExit("ERROR: nominal Stage 3D candidate is empty")

if abs(nominal[0][0] - 200.0) > 1e-9:
    raise SystemExit("ERROR: nominal support must start at 200 nm")

if abs(nominal[-1][0] - 900.0) > 1e-9:
    raise SystemExit("ERROR: nominal support must end at 900 nm")


# ------------------------------------------------------------
# Variant B: UV opaque below 300 nm.
#
# Dense 1-nm sampling prevents a long interpolation interval from
# smearing the deliberately sharp conservative UV bound.
# >=300 nm follows nominal exactly.
# ------------------------------------------------------------

nominal_map = {wl: L for wl, L in nominal}

uv_rows = []

for wl in range(200, 300):
    uv_rows.append((float(wl), UV_OPAQUE_L_MM))

# Keep the conservative UV bound opaque essentially all the way to
# the measured 300-nm endpoint.  The near-endpoint node prevents
# Geant4 energy-linear interpolation from smearing the transition
# over the entire 299-300 nm interval.
uv_rows.append((299.999, UV_OPAQUE_L_MM))

for wl, L in nominal:
    if wl >= 300.0:
        uv_rows.append((wl, L))


# ------------------------------------------------------------
# Variant C: visible-transparent above 500 nm.
#
# <=500 nm follows nominal exactly.
# 501..900 nm use the high-transparency proxy.
# ------------------------------------------------------------

vis_rows = []

for wl, L in nominal:
    if wl <= 500.0:
        vis_rows.append((wl, L))

# Sharp systematic transition immediately above the last measured node.
vis_rows.append((500.001, VISIBLE_L_MM))

for wl in range(501, 901):
    vis_rows.append((float(wl), VISIBLE_L_MM))


def write_table(path, rows, description):
    rows = sorted(dict(rows).items())

    with path.open("w") as f:
        f.write("# PMMA Stage 3D systematic sensitivity candidate\n")
        f.write(f"# {description}\n")
        f.write("#\n")
        f.write("# IMPORTANT:\n")
        f.write("#   300-500 nm retains the SNO Sample-A measured-domain candidate.\n")
        f.write("#   Values outside 300-500 nm are explicit systematic models.\n")
        f.write("#\n")
        f.write("# wavelength_nm   absorption_length_mm\n")

        for wl, L in rows:
            f.write(f"{wl:10.4f}  {L:16.9g}\n")


write_table(
    UV_OPAQUE,
    uv_rows,
    "UV conservative bound: lambda < 300 nm effectively opaque; nominal elsewhere."
)

write_table(
    VISIBLE_TRANSPARENT,
    vis_rows,
    "Visible transparent bound: lambda > 500 nm alpha=1e-4 cm^-1; nominal elsewhere."
)

print("=== PMMA STAGE 3D SENSITIVITY VARIANTS ===")
print(f"NOMINAL             : {NOMINAL.name}")
print(f"UV_OPAQUE           : {UV_OPAQUE.name}")
print(f"VISIBLE_TRANSPARENT : {VISIBLE_TRANSPARENT.name}")
print()
print(f"UV opaque L_abs     = {UV_OPAQUE_L_MM:g} mm")
print(f"visible proxy alpha = {VISIBLE_ALPHA_CM:g} cm^-1")
print(f"visible proxy L_abs = {VISIBLE_L_MM:g} mm")
print()
print("PMMA STAGE 3D SENSITIVITY VARIANTS: CREATED")
