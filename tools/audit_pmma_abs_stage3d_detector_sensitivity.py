#!/usr/bin/env python3

from pathlib import Path
import math

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

HC = 1239.841984332  # eV nm

FILES = {
    "NOMINAL":
        "PMMA_SNO_UVT_SampleA_ABS_Measured300to500_Model200to900_Candidate.txt",

    "UV_OPAQUE":
        "PMMA_SNO_UVT_SampleA_ABS_UVOpaqueBelow300_Model200to900_Candidate.txt",

    "VISIBLE_TRANSPARENT":
        "PMMA_SNO_UVT_SampleA_ABS_VisibleTransparentAbove500_Model200to900_Candidate.txt",
}

PDE_FILE = DATA / "SiPM_PDE.txt"

PATHS_MM = [3.0, 6.0, 10.0]

STEP_NM = 0.1


def read_table(path):
    rows = []

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        cols = line.split()

        try:
            x = float(cols[0])
            y = float(cols[1])
        except (IndexError, ValueError):
            continue

        rows.append((x, y))

    rows.sort()
    return rows


def interp_linear(rows, x):
    if x <= rows[0][0]:
        return rows[0][1]

    if x >= rows[-1][0]:
        return rows[-1][1]

    lo = 0
    hi = len(rows) - 1

    while hi - lo > 1:
        mid = (lo + hi) // 2
        if rows[mid][0] <= x:
            lo = mid
        else:
            hi = mid

    x0, y0 = rows[lo]
    x1, y1 = rows[hi]

    t = (x - x0) / (x1 - x0)

    return y0 + t*(y1-y0)


def abs_energy_table(path):
    wl_rows = read_table(path)

    energy_rows = [
        (HC / wl, L)
        for wl, L in wl_rows
    ]

    energy_rows.sort()

    return energy_rows


def abs_length_geant4_like(energy_rows, wl_nm):
    E = HC / wl_nm
    return interp_linear(energy_rows, E)


pde = read_table(PDE_FILE)

tables = {
    name: abs_energy_table(DATA / filename)
    for name, filename in FILES.items()
}


def spectral_bin(wl):
    if wl < 300.0:
        return "<300"
    if wl <= 500.0:
        return "300-500"
    return ">500"


def integrate(policy, path_mm):
    table = tables[policy]

    totals = {
        "<300": 0.0,
        "300-500": 0.0,
        ">500": 0.0,
    }

    transmitted = {
        "<300": 0.0,
        "300-500": 0.0,
        ">500": 0.0,
    }

    wl = 200.0

    while wl < 900.0:
        wm = wl + 0.5*STEP_NM

        PDE = max(0.0, interp_linear(pde, wm))

        # Frank-Tamm spectral dependence, with the PDE factor.
        w = PDE / (wm*wm) * STEP_NM

        L = abs_length_geant4_like(table, wm)

        if L <= 0.0:
            T = 0.0
        else:
            T = math.exp(-path_mm/L)

        b = spectral_bin(wm)

        totals[b] += w
        transmitted[b] += w*T

        wl += STEP_NM

    total_in = sum(totals.values())
    total_out = sum(transmitted.values())

    return totals, transmitted, total_in, total_out


print("=== PMMA STAGE 3D DETECTOR-LEVEL SENSITIVITY AUDIT ===")

for path_mm in PATHS_MM:

    print()
    print(f"--- PMMA path length = {path_mm:.1f} mm ---")

    results = {}

    for policy in FILES:

        total, trans, total_in, total_out = integrate(policy, path_mm)

        results[policy] = total_out

        survival = total_out / total_in

        print(
            f"{policy:20s} "
            f"PDE-weighted survival={100*survival:9.5f}% "
            f"loss={100*(1-survival):9.5f}%"
        )

        for b in ("<300", "300-500", ">500"):
            if total[b] > 0:
                bin_survival = trans[b] / total[b]
            else:
                bin_survival = 0.0

            print(
                f"    {b:>7}: "
                f"input={100*total[b]/total_in:7.3f}% "
                f"T={100*bin_survival:9.5f}%"
            )

    nominal = results["NOMINAL"]

    print()
    print("  Yield change relative to NOMINAL:")

    for policy in ("UV_OPAQUE", "VISIBLE_TRANSPARENT"):
        delta = 100.0 * (results[policy]/nominal - 1.0)

        print(
            f"    {policy:20s}: "
            f"{delta:+.6f}%"
        )

print()
print("Interpretation:")
print("  3 mm = normal shell crossing.")
print("  6 and 10 mm are oblique/multiple-path stress scales.")
print("  This is PDE/lambda^2 weighted spectral transport, not yet a full event MC.")
print()
print("PMMA STAGE 3D DETECTOR SENSITIVITY AUDIT: COMPLETE")
