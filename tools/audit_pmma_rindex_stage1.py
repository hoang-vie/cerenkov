#!/usr/bin/env python3
"""
PMMA_SNO_A Stage-1 audit: refractive-index model.

This is deliberately an audit of the CURRENT project model, not a replacement.

Current project RINDEX:
  data/PMMA_RINDEX_Roychowdhury2020_Cauchy3Anchor.txt

The project table is a 3-term Cauchy curve constrained to three refractive-index
anchors quoted by Roychowdhury et al. (2020):
  n(300 nm)=1.528, n(500 nm)=1.492, n(1000 nm)=1.479.

Important scope:
Roychowdhury et al. did NOT publish this 3-term Cauchy curve as their material
model. Their optical functions were obtained from spectroscopic ellipsometry /
transmission modelling using a Cody-Lorentz oscillator plus Gaussian
oscillators. Therefore this audit can prove:
  - the project table really implements its declared 3-anchor approximation;
  - Geant4 interpolation error relative to that approximation is small;
  - the approximation is in the range of other published PMMA measurements.
It cannot prove the full 200-900 nm project curve equals Roychowdhury's raw
spectral n(lambda) point-by-point.

Independent PMMA cross-checks:
  Beadie et al., Applied Optics 54, F139-F143 (2015), Table 2:
    472.9 nm -> 1.498795
    780.4 nm -> 1.486417
  These are direct independent PMMA measurements adjusted to 20.1 C.

  Sultanova et al., Acta Physica Polonica A 116, 585-587 (2009):
    PMMA Sellmeier fit at 20 C over 0.4368-1.052 um,
      n^2 - 1 = 1.1819 lambda^2 / (lambda^2 - 0.011313)
    with lambda in um.
The independent samples are NOT the same acrylic batch as Roychowdhury or SNO.
Their difference is treated as a material/sample systematic, not a node error.
"""

from __future__ import annotations
from pathlib import Path
import bisect
import math
import sys

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "PMMA_RINDEX_Roychowdhury2020_Cauchy3Anchor.txt"

# Project-declared 3-anchor Cauchy model, lambda in um.
A = 1.47495512821
B = 0.00397275641026
C = 7.21153846154e-05

ROY_ANCHORS = {
    300.0: 1.528,
    500.0: 1.492,
    1000.0: 1.479,
}

# Beadie et al. 2015, Table 2: independent PMMA index measurements at 20.1 C.
BEADIE = {
    472.9: 1.498795,
    780.4: 1.486417,
}

HC_EV_NM = 1239.8419843320026
ALPHA = 1.0 / 137.035999084
MUON_MASS_MEV = 105.6583755


def cauchy_n(lam_nm: float) -> float:
    u = lam_nm * 1e-3
    return A + B / (u*u) + C / (u**4)


def sultanova_n(lam_nm: float) -> float:
    u = lam_nm * 1e-3
    return math.sqrt(1.0 + 1.1819*u*u/(u*u - 0.011313))


def read_table(path: Path):
    rows = []
    for raw in path.read_text().splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        a, b = s.split()[:2]
        rows.append((float(a), float(b)))
    rows.sort()
    if not rows:
        raise RuntimeError(f"no numeric rows in {path}")
    return rows


def energy_interp(rows, lam_nm: float) -> float:
    # Emulate G4MaterialPropertyVector interpolation after wavelength->energy
    # conversion and ascending-energy sorting.
    erows = sorted((HC_EV_NM/lam, n) for lam, n in rows)
    xs = [x for x, _ in erows]
    e = HC_EV_NM/lam_nm
    if e <= xs[0]:
        return erows[0][1]
    if e >= xs[-1]:
        return erows[-1][1]
    i = bisect.bisect_right(xs, e)
    e0, n0 = erows[i-1]
    e1, n1 = erows[i]
    f = (e-e0)/(e1-e0)
    return n0 + f*(n1-n0)


def theta_deg(n: float, beta: float = 1.0) -> float:
    x = 1.0/(beta*n)
    if x >= 1.0:
        return float("nan")
    return math.degrees(math.acos(x))


def integrate_frank_tamm(nfunc, beta: float, lo_nm: float, hi_nm: float,
                         step_nm: float = 0.05) -> float:
    # photons / cm
    total = 0.0
    prev = None
    nsteps = int(round((hi_nm-lo_nm)/step_nm))
    for i in range(nsteps+1):
        lam = lo_nm + i*step_nm
        n = nfunc(lam)
        term = max(0.0, 1.0 - 1.0/(beta*beta*n*n))
        lam_cm = lam*1e-7
        y = 2.0*math.pi*ALPHA*term/(lam_cm*lam_cm)*1e-7
        if prev is not None:
            total += 0.5*(prev+y)*step_nm
        prev = y
    return total


def max_scan(func, lo, hi, step=0.01):
    best = (-1.0, None)
    x = lo
    while x <= hi + 0.5*step:
        y = func(x)
        if y > best[0]:
            best = (y, x)
        x += step
    return best


def main() -> int:
    rows = read_table(DATA)
    wl_min, wl_max = rows[0][0], rows[-1][0]

    node_res = max((abs(n-cauchy_n(lam)), lam) for lam, n in rows)

    anchor_res = []
    for lam, ref in ROY_ANCHORS.items():
        calc = cauchy_n(lam)
        anchor_res.append((abs(calc-ref), lam, calc, ref))

    dn_all = max_scan(lambda l: abs(energy_interp(rows,l)-cauchy_n(l)),
                      wl_min, wl_max)
    da_all = max_scan(lambda l: abs(theta_deg(energy_interp(rows,l)) -
                                    theta_deg(cauchy_n(l))),
                      wl_min, wl_max)
    dn_vis = max_scan(lambda l: abs(energy_interp(rows,l)-cauchy_n(l)),
                      max(300.0,wl_min), min(600.0,wl_max))
    da_vis = max_scan(lambda l: abs(theta_deg(energy_interp(rows,l)) -
                                    theta_deg(cauchy_n(l))),
                      max(300.0,wl_min), min(600.0,wl_max))

    gamma = 1.0 + 3000.0/MUON_MASS_MEV
    beta = math.sqrt(1.0 - 1.0/(gamma*gamma))
    y_cont = integrate_frank_tamm(cauchy_n, beta, wl_min, wl_max)
    y_g4 = integrate_frank_tamm(lambda l: energy_interp(rows,l), beta,
                                wl_min, wl_max)
    ppm = (y_g4/y_cont-1.0)*1e6

    # Independent direct Beadie points.
    beadie_rows = []
    for lam, ref in sorted(BEADIE.items()):
        cur = cauchy_n(lam)
        beadie_rows.append((lam,cur,ref,cur-ref,100*(cur/ref-1.0)))

    # Independent Sultanova PMMA curve over the overlap with our 200-900 nm table.
    sul_lo, sul_hi = max(436.8,wl_min), min(900.0,wl_max)
    sul_max = max_scan(lambda l: abs(cauchy_n(l)-sultanova_n(l)),
                       sul_lo, sul_hi)
    sul_at = [(lam,cauchy_n(lam),sultanova_n(lam),
               cauchy_n(lam)-sultanova_n(lam),
               theta_deg(cauchy_n(lam))-theta_deg(sultanova_n(lam)))
              for lam in (450.0,500.0,600.0,800.0)
              if sul_lo <= lam <= sul_hi]

    y_cur_common = integrate_frank_tamm(cauchy_n,beta,sul_lo,sul_hi)
    y_sul_common = integrate_frank_tamm(sultanova_n,beta,sul_lo,sul_hi)

    print("=== PMMA_SNO_A RINDEX STAGE-1 AUDIT ===")
    print(f"file: {DATA}")
    print(f"rows={len(rows)} range={wl_min:g}-{wl_max:g} nm")
    print()
    print("Declared Roychowdhury 3-anchor Cauchy approximation:")
    for _,lam,calc,ref in sorted(anchor_res, key=lambda x:x[1]):
        print(f"  {lam:6.1f} nm: Cauchy={calc:.9f} quoted anchor={ref:.9f} "
              f"delta={calc-ref:+.3e}")
    print(f"max table-node |n_table-n_Cauchy|={node_res[0]:.3e} at {node_res[1]:g} nm")
    print()
    print("Geant4 energy-linear interpolation relative to continuous Cauchy:")
    print(f"  max |dn| 200-900 nm = {dn_all[0]:.3e} at {dn_all[1]:.3f} nm")
    print(f"  max beta=1 angle error = {da_all[0]:.6f} deg")
    print(f"  max |dn| 300-600 nm = {dn_vis[0]:.3e} at {dn_vis[1]:.3f} nm")
    print(f"  max beta=1 angle error 300-600 nm = {da_vis[0]:.6f} deg")
    print(f"  3-GeV mu Frank-Tamm continuous Cauchy = {y_cont:.6f} photons/cm")
    print(f"  sampled+G4 interpolation = {y_g4:.6f} photons/cm")
    print(f"  interpolation yield bias = {y_g4-y_cont:+.6f} photons/cm ({ppm:+.1f} ppm)")
    print()
    print("Independent PMMA sample cross-check: Beadie et al. 2015 Table 2")
    for lam,cur,ref,dn,pct in beadie_rows:
        print(f"  {lam:6.1f} nm: current={cur:.9f} Beadie={ref:.9f} "
              f"delta={dn:+.6f} ({pct:+.3f}%)")
    print()
    print("Independent PMMA sample cross-check: Sultanova et al. 2009 Sellmeier")
    print(f"  overlap={sul_lo:g}-{sul_hi:g} nm")
    print(f"  max |current-Sultanova|={sul_max[0]:.6f} at {sul_max[1]:.1f} nm")
    for lam,cur,ref,dn,dth in sul_at:
        print(f"  {lam:6.1f} nm: current={cur:.9f} Sultanova={ref:.9f} "
              f"delta_n={dn:+.6f} delta_theta(beta=1)={dth:+.4f} deg")
    print(f"  3-GeV mu Frank-Tamm over common {sul_lo:g}-{sul_hi:g} nm:")
    print(f"    current={y_cur_common:.6f} /cm Sultanova={y_sul_common:.6f} /cm "
          f"delta={100*(y_cur_common/y_sul_common-1.0):+.3f}%")
    print()
    print("Interpretation:")
    print("  - The PROJECT TABLE accurately implements its declared 3-anchor Cauchy curve.")
    print("  - The Roychowdhury paper used a much richer ellipsometry oscillator model;")
    print("    matching only n(300), n(500), n(1000) does NOT establish point-by-point")
    print("    equivalence to the published Roychowdhury spectral curve.")
    print("  - Other measured PMMA samples are ~0.003-0.005 higher in n over much of")
    print("    the visible/NIR. This is a material/sample systematic, not a table bug.")
    print("  - PMMA_SNO_A is therefore a HYBRID reference profile: generic-commercial-PMMA")
    print("    RINDEX plus SNO Sample-A ABSLENGTH; its RINDEX is not SNO-Sample-A-specific.")

    passed = (
        wl_min >= 200.0 and wl_max >= 900.0
        and max(a[0] for a in anchor_res) < 1e-10
        and node_res[0] < 1e-8
        and dn_vis[0] < 3e-5
        and abs(ppm) < 50.0
    )
    print()
    print("PMMA RINDEX 3-ANCHOR IMPLEMENTATION AUDIT: " + ("PASS" if passed else "FAIL"))
    print("PMMA RINDEX FULL ROYCHOWDHURY SPECTRAL EQUIVALENCE: NOT_ESTABLISHED")
    print("PMMA_SNO_A RINDEX MATERIAL STATUS: GENERIC_REFERENCE_NOT_SNO_SPECIFIC")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
