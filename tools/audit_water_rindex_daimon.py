#!/usr/bin/env python3
"""Independent audit of data/Water_RINDEX_Daimon2007_21p5C.txt.

Reference model:
  M. Daimon and A. Masumura, Applied Optics 46, 3811-3820 (2007),
  DOI 10.1364/AO.46.003811.  Four-term Sellmeier fit for HPLC distilled
  water at 21.5 C, published validity 0.182-1.129 um.

This script intentionally does NOT derive its expected n(lambda) values from
project data.  The published Sellmeier coefficients are hard-coded below.
It also emulates the interpolation relevant to Geant4: material-property values
are linearly interpolated versus photon energy, not wavelength.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

HC_EV_NM = 1239.84193
ALPHA = 1.0 / 137.035999084
MUON_MASS_MEV = 105.6583755

# Daimon & Masumura 2007, 21.5 C. lambda is in micrometres.
B = (5.689093832e-1, 1.719708856e-1, 2.062501582e-2, 1.123965424e-1)
C = (5.110301794e-3, 1.825180155e-2, 2.624158904e-2, 1.067505178e1)
PUBLISHED_MIN_NM = 182.0
PUBLISHED_MAX_NM = 1129.0


def sellmeier_n(wavelength_nm: float) -> float:
    lam_um = wavelength_nm * 1.0e-3
    lam2 = lam_um * lam_um
    n2 = 1.0
    for b, c in zip(B, C):
        n2 += b * lam2 / (lam2 - c)
    return math.sqrt(n2)


def read_table(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    for lineno, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) != 2:
            raise ValueError(f'{path}:{lineno}: expected 2 columns')
        wl, n = map(float, parts)
        rows.append((wl, n))
    if len(rows) < 2:
        raise ValueError(f'{path}: need at least two rows')
    if any(rows[i+1][0] <= rows[i][0] for i in range(len(rows)-1)):
        raise ValueError(f'{path}: wavelength column must be strictly increasing')
    return rows


def interp_linear(x: float, xs: list[float], ys: list[float]) -> float:
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    f = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] + f * (ys[hi] - ys[lo])


def g4_interpolated_n(wavelength_nm: float, rows: list[tuple[float, float]]) -> float:
    # DetectorConstruction converts wavelength to photon energy and sorts by energy.
    points = sorted((HC_EV_NM / wl, n) for wl, n in rows)
    energies = [p[0] for p in points]
    values = [p[1] for p in points]
    return interp_linear(HC_EV_NM / wavelength_nm, energies, values)


def cherenkov_angle_deg(beta: float, n: float) -> float:
    if beta * n <= 1.0:
        return float('nan')
    return math.degrees(math.acos(1.0 / (beta * n)))


def integrate_frank_tamm(rows: list[tuple[float, float]], beta: float, use_g4_interp: bool) -> float:
    lo, hi = rows[0][0], rows[-1][0]
    step_nm = 0.01
    count = int(round((hi - lo) / step_nm))
    total = 0.0
    prev = None
    prev_wl = None
    for i in range(count + 1):
        wl = lo + i * step_nm
        n = g4_interpolated_n(wl, rows) if use_g4_interp else sellmeier_n(wl)
        term = max(0.0, 1.0 - 1.0 / (beta * beta * n * n)) / (wl * wl)
        if prev is not None:
            total += 0.5 * (prev + term) * (wl - prev_wl)
        prev, prev_wl = term, wl
    # lambda in nm: integral d(lambda_nm)/lambda_nm^2 -> 1e7 per cm.
    return 2.0 * math.pi * ALPHA * 1.0e7 * total


def max_interp_error(rows: list[tuple[float, float]], lo: float, hi: float, step_nm: float = 0.001):
    n_steps = int(round((hi - lo) / step_nm))
    max_dn = (0.0, lo)
    max_dtheta = (0.0, lo)
    for i in range(n_steps + 1):
        wl = lo + i * step_nm
        nt = sellmeier_n(wl)
        ng = g4_interpolated_n(wl, rows)
        dn = abs(ng - nt)
        da = abs(cherenkov_angle_deg(1.0, ng) - cherenkov_angle_deg(1.0, nt))
        if dn > max_dn[0]:
            max_dn = (dn, wl)
        if da > max_dtheta[0]:
            max_dtheta = (da, wl)
    return max_dn, max_dtheta


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    path = root / 'data' / 'Water_RINDEX_Daimon2007_21p5C.txt'
    rows = read_table(path)

    wl_min, wl_max = rows[0][0], rows[-1][0]
    domain_ok = wl_min >= PUBLISHED_MIN_NM and wl_max <= PUBLISHED_MAX_NM

    residuals = [(abs(n - sellmeier_n(wl)), wl, n, sellmeier_n(wl)) for wl, n in rows]
    max_node = max(residuals)

    (max_dn_all, max_da_all) = max_interp_error(rows, wl_min, wl_max)
    vis_lo, vis_hi = max(300.0, wl_min), min(600.0, wl_max)
    (max_dn_vis, max_da_vis) = max_interp_error(rows, vis_lo, vis_hi)

    gamma = 1.0 + 3000.0 / MUON_MASS_MEV
    beta_3gev = math.sqrt(1.0 - 1.0 / (gamma * gamma))
    y_cont = integrate_frank_tamm(rows, beta_3gev, use_g4_interp=False)
    y_g4 = integrate_frank_tamm(rows, beta_3gev, use_g4_interp=True)
    rel_ppm = (y_g4 / y_cont - 1.0) * 1.0e6

    print('=== WATER RINDEX: DAIMON 2007 INDEPENDENT AUDIT ===')
    print(f'file: {path}')
    print(f'rows={len(rows)} range={wl_min:g}-{wl_max:g} nm')
    print(f'published Sellmeier validity={PUBLISHED_MIN_NM:g}-{PUBLISHED_MAX_NM:g} nm: '
          f'{"PASS" if domain_ok else "FAIL"}')
    print(f'max node |n_table-n_Sellmeier|={max_node[0]:.3e} at {max_node[1]:g} nm')
    print(f'max Geant4 energy-linear interpolation |dn| over table range='
          f'{max_dn_all[0]:.3e} at {max_dn_all[1]:.3f} nm')
    print(f'max beta=1 Cherenkov-angle interpolation error over table range='
          f'{max_da_all[0]:.6f} deg at {max_da_all[1]:.3f} nm')
    print(f'max |dn| over 300-600 nm={max_dn_vis[0]:.3e} at {max_dn_vis[1]:.3f} nm')
    print(f'max beta=1 angle error over 300-600 nm={max_da_vis[0]:.6f} deg')
    print(f'3 GeV mu beta={beta_3gev:.9f}')
    print(f'Frank-Tamm 200-900 nm continuous Sellmeier={y_cont:.6f} photons/cm')
    print(f'Frank-Tamm 200-900 nm sampled+G4 interpolation={y_g4:.6f} photons/cm')
    print(f'sampling/interpolation yield bias={y_g4-y_cont:+.6f} photons/cm ({rel_ppm:+.1f} ppm)')

    # The node tolerance checks whether the project table was actually sampled
    # from the cited formula.  The interpolation limits are QA limits, not claims
    # about the experimental uncertainty of Daimon & Masumura.
    passed = (
        domain_ok
        and max_node[0] < 5e-9
        and max_dn_all[0] < 4e-4
        and max_dn_vis[0] < 2e-5
        and abs(rel_ppm) < 100.0
    )
    print('WATER RINDEX DAIMON AUDIT:', 'PASS' if passed else 'FAIL')
    return 0 if passed else 1


if __name__ == '__main__':
    sys.exit(main())
