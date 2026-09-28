#!/usr/bin/env python3
"""
Independent physical audit of the current Water_Default Rayleigh table.

This tool intentionally does NOT modify the simulation. It asks two questions:

1) What spectral law is actually encoded by Water_RAYLEIGH_Lambda4_Baseline.txt?
2) Is that law consistent with a pure-fresh-water molecular scattering reference?

References used for the numerical cross-check:
- Morel (1974) pure-water values as reproduced in Twardowski et al.,
  Biogeosciences 4, 1041-1058 (2007), Table 2.
- Buiteveld, Hakvoort & Donze (1994) pure-water values as reproduced in the
  same Table 2.
- The Morel pure-SEAWATER parameterization widely used in ocean optics:
    b_sw(lambda) = 0.00288 * (lambda/500 nm)^(-4.32) [m^-1].

The purpose is provenance/physics auditing, not to declare a replacement model.
"""

from __future__ import annotations
from pathlib import Path
import math
import sys

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "Water_RAYLEIGH_Lambda4_Baseline.txt"

# Fresh-water total molecular scattering coefficients b_w [m^-1].
# Values are from Twardowski et al. (2007), Table 2, which tabulates
# Morel (1974) and Buiteveld et al. (1994).
MOREL_FRESH = {
    350.0: 103.9e-4,
    375.0: 77.1e-4,
    400.0: 58.3e-4,
    425.0: 45.0e-4,
    450.0: 35.0e-4,
    475.0: 27.8e-4,
    500.0: 22.2e-4,
    525.0: 18.0e-4,
    550.0: 14.9e-4,
    575.0: 12.5e-4,
    600.0: 10.9e-4,
}

BUITEVELD_FRESH = {
    350.0: 92.0e-4,
    375.0: 69.2e-4,
    400.0: 53.0e-4,
    425.0: 41.2e-4,
    450.0: 32.6e-4,
    475.0: 26.0e-4,
    500.0: 21.1e-4,
    525.0: 17.2e-4,
    550.0: 14.2e-4,
    575.0: 11.8e-4,
    600.0: 9.9e-4,
}

def read_table(path: Path):
    rows = []
    for raw in path.read_text().splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        a, b = s.split()[:2]
        rows.append((float(a), float(b)))
    if not rows:
        raise RuntimeError(f"no numeric rows in {path}")
    return rows

def interp_lambda(rows, lam):
    rows = sorted(rows)
    if lam <= rows[0][0]:
        return rows[0][1]
    if lam >= rows[-1][0]:
        return rows[-1][1]
    for (x0, y0), (x1, y1) in zip(rows, rows[1:]):
        if x0 <= lam <= x1:
            f = (lam - x0) / (x1 - x0)
            return y0 + f * (y1 - y0)
    raise AssertionError

def fit_power_law(rows, lo=200.0, hi=800.0):
    # Fit b(lambda)=A*(lambda/500)^(-p) in log space.
    pts = [(lam, 1.0/L) for lam, L in rows if lo <= lam <= hi]
    xs = [math.log(lam/500.0) for lam, _ in pts]
    ys = [math.log(b) for _, b in pts]
    xm = sum(xs)/len(xs)
    ym = sum(ys)/len(ys)
    cov = sum((x-xm)*(y-ym) for x, y in zip(xs, ys))
    var = sum((x-xm)**2 for x in xs)
    slope = cov/var
    intercept = ym - slope*xm
    return math.exp(intercept), -slope

def morel_seawater_b(lam):
    return 0.00288 * (lam/500.0)**(-4.32)

rows = read_table(DATA)
A_fit, p_fit = fit_power_law(rows)

# How closely does the current table follow the canonical Morel pure-seawater
# power law over the range where the table tracks it almost exactly?
resid = []
for lam, L in rows:
    if 200.0 <= lam <= 800.0:
        b = 1.0/L
        ref = morel_seawater_b(lam)
        resid.append((abs(b/ref - 1.0), lam, b, ref))
max_sw_resid, max_sw_lam, _, _ = max(resid)

print("=== WATER RAYLEIGH STAGE-3 PHYSICAL AUDIT ===")
print(f"file: {DATA}")
print(f"rows={len(rows)} range={rows[0][0]:g}-{rows[-1][0]:g} nm")
print()
print("Encoded spectral law (fit over 200-800 nm):")
print(f"  b(lambda) = {A_fit:.7f} * (lambda/500 nm)^(-{p_fit:.5f}) m^-1")
print("Canonical Morel pure-SEAWATER law:")
print("  b_sw(lambda) = 0.0028800 * (lambda/500 nm)^(-4.32000) m^-1")
print(f"  max relative residual current-vs-seawater law over 200-800 nm = "
      f"{100*max_sw_resid:.3f}% at {max_sw_lam:g} nm")
print()
print("Fresh-water reference comparison:")
print(" lambda   current b      Morel fresh   ratio    Buiteveld fresh  ratio")
rat_morel = []
rat_buit = []
for lam in sorted(MOREL_FRESH):
    L = interp_lambda(rows, lam)
    b = 1.0/L
    rm = b/MOREL_FRESH[lam]
    rb = b/BUITEVELD_FRESH[lam]
    rat_morel.append((rm, lam))
    rat_buit.append((rb, lam))
    print(f" {lam:5.0f}  {b:11.6g}   {MOREL_FRESH[lam]:11.6g}  "
          f"{rm:6.3f}    {BUITEVELD_FRESH[lam]:11.6g}  {rb:6.3f}")

# Direct-path impact at two representative wavelengths.
print()
print("Direct-path scattering probability P=1-exp(-d/L):")
for lam in (400.0, 500.0):
    Lcur = interp_lambda(rows, lam)
    Lmf = 1.0/MOREL_FRESH[lam]
    for d in (0.10, 1.0, 5.0):
        pcur = 1.0 - math.exp(-d/Lcur)
        pref = 1.0 - math.exp(-d/Lmf)
        print(f"  {lam:3.0f} nm, d={d:>4.1f} m: current={100*pcur:.4f}%  "
              f"Morel-fresh={100*pref:.4f}%  delta={100*(pcur-pref):+.4f} pp")

print()
comment_claim = "lambda^-4"
print("Audit interpretation:")
print(f"  - File comment claims approximately {comment_claim}, but the numerical table")
print(f"    fits exponent p={p_fit:.4f} over 200-800 nm.")
print("  - Its normalization and slope are numerically consistent with the classic")
print("    Morel pure-SEAWATER parameterization, not with the pure-fresh-water")
print("    reference values tabulated by Morel/Buiteveld.")
print("  - Therefore the current table remains useful as a historical regression")
print("    baseline, but it is NOT accepted as an absolute Water_Default pure-water")
print("    Rayleigh model.")

# The audit itself passes if it successfully identifies the current historical
# model. Physical acceptance is deliberately a separate status.
audit_ok = (
    abs(A_fit/0.00288 - 1.0) < 0.01
    and abs(p_fit - 4.32) < 0.02
    and max_sw_resid < 0.01
    and (1.0/interp_lambda(rows, 500.0))/MOREL_FRESH[500.0] > 1.20
)

print()
print("WATER RAYLEIGH MODEL IDENTIFICATION AUDIT: " + ("PASS" if audit_ok else "FAIL"))
print("LEGACY RAYLEIGH PHYSICAL STATUS: REJECTED_AS_PURE_WATER")

sys.exit(0 if audit_ok else 1)
