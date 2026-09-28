#!/usr/bin/env python3

from pathlib import Path
import math

ROOT = Path(__file__).resolve().parents[1]

SRC = (
    ROOT / "data" /
    "PMMA_SNO_UVT_SampleA_ABS_TableI_RT_plus_TableII_300to500_Candidate.txt"
)

OUT = (
    ROOT / "data" /
    "PMMA_SNO_UVT_SampleA_ABS_Measured300to500_Model200to900_Candidate.txt"
)

UV_MODEL_MIN_NM = 250.0
TRANSPORT_MIN_NM = 200.0
TRANSPORT_MAX_NM = 900.0

# Local log(alpha) slope obtained in Stage 3D-A from the
# same Sample-A measured nodes over 300-320 nm.
UV_K_PER_NM = 0.09526416


def read_table(path):
    rows = []

    for line in path.read_text().splitlines():
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        parts = line.replace(",", " ").split()

        vals = []
        for x in parts:
            try:
                vals.append(float(x))
            except ValueError:
                pass

        if len(vals) >= 2:
            rows.append((vals[0], vals[1]))

    rows.sort()
    return rows


measured = read_table(SRC)

if not measured:
    raise SystemExit("No measured PMMA data found.")

if abs(measured[0][0] - 300.0) > 1e-9:
    raise SystemExit("Expected first measured wavelength = 300 nm.")

if abs(measured[-1][0] - 500.0) > 1e-9:
    raise SystemExit("Expected last measured wavelength = 500 nm.")


def nearest(wl):
    return min(measured, key=lambda p: abs(p[0] - wl))


wl300, L300 = nearest(300.0)
wl500, L500 = nearest(500.0)

alpha300 = 10.0 / L300
alpha500 = 10.0 / L500


def uv_alpha(wl_nm):
    # Endpoint-anchored local UV-edge continuation:
    #
    # alpha(lambda) =
    # alpha(300) * exp[k * (300-lambda)]
    #
    # This is explicitly a policy model, NOT measured SNO data.
    return alpha300 * math.exp(UV_K_PER_NM * (300.0 - wl_nm))


alpha250 = uv_alpha(UV_MODEL_MIN_NM)
L250 = 10.0 / alpha250


out = []

# Deep-UV policy:
# below 250 nm keep the already effectively-opaque 250-nm value.
out.append((TRANSPORT_MIN_NM, L250, "MODEL_UV_OPAQUE_HOLD250"))
out.append((UV_MODEL_MIN_NM, L250, "MODEL_UV_EDGE"))

# 1-nm grid between 251 and 299 nm.
# Dense sampling suppresses interpolation artifacts because Geant4
# interpolates material properties in photon energy.
for wl in range(251, 300):
    alpha = uv_alpha(float(wl))
    Lmm = 10.0 / alpha
    out.append((float(wl), Lmm, "MODEL_UV_EDGE"))

# Directly validated measured domain.
for wl, Lmm in measured:
    out.append((wl, Lmm, "MEASURED_SNO_SAMPLE_A"))

# Conservative visible continuation.
# alpha held at the measured 500-nm endpoint.
out.append((TRANSPORT_MAX_NM, L500, "MODEL_VISIBLE_HOLD500"))

# Remove duplicate wavelengths, measured data win.
priority = {
    "MODEL_UV_OPAQUE_HOLD250": 0,
    "MODEL_UV_EDGE": 1,
    "MODEL_VISIBLE_HOLD500": 1,
    "MEASURED_SNO_SAMPLE_A": 2,
}

merged = {}

for wl, Lmm, status in out:
    old = merged.get(wl)

    if old is None or priority[status] > priority[old[1]]:
        merged[wl] = (Lmm, status)


with OUT.open("w") as f:
    f.write(
        "# PMMA_SNO_A ABS transport policy candidate\n"
        "#\n"
        "# Empirical measured/validated domain:\n"
        "#   300-500 nm\n"
        "#   SNO Sample A\n"
        "#\n"
        "# Model extensions:\n"
        "#   250-300 nm : endpoint-anchored local UV-edge exponential\n"
        "#                slope from same-sample 300-320 nm data\n"
        "#   200-250 nm : hold L_abs(250); material is already effectively\n"
        "#                opaque across a 3-mm shell at this point\n"
        "#   500-900 nm : conservative constant-alpha continuation using\n"
        "#                alpha(500 nm)\n"
        "#\n"
        "# IMPORTANT:\n"
        "#   Values outside 300-500 nm are POLICY MODEL values, not direct\n"
        "#   measured SNO Sample-A absorption data.\n"
        "#\n"
        "# wavelength_nm   absorption_length_mm\n"
    )

    for wl in sorted(merged):
        Lmm, status = merged[wl]
        f.write(
            f"{wl:10.4f}  {Lmm:16.9g}"
            f"  # {status}\n"
        )


print("=== PMMA STAGE 3D-B POLICY CANDIDATE ===")
print(f"source   = {SRC.name}")
print(f"output   = {OUT.name}")
print()
print(f"alpha300 = {alpha300:.9g} cm^-1")
print(f"L300     = {L300:.9g} mm")
print(f"UV k     = {UV_K_PER_NM:.9g} nm^-1")
print()
print(f"alpha250 = {alpha250:.9g} cm^-1")
print(f"L250     = {L250:.9g} mm")

T250 = math.exp(-3.0 / L250)

print(f"T250 through 3 mm = {100*T250:.9g}%")
print()
print(f"alpha500 = {alpha500:.9g} cm^-1")
print(f"L500     = {L500:.9g} mm")
print()
print("Measured validity = [300,500] nm")
print("Transport support = [200,900] nm")
print()
print("PMMA STAGE 3D-B POLICY CANDIDATE: CREATED")
