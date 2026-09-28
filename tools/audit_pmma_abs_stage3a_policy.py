#!/usr/bin/env python3
"""
PMMA_SNO_A Stage-3A: out-of-range absorption policy sensitivity.

This is a sensitivity study, NOT a production-physics replacement.

Measured production ABS domain:
    300-500 nm, SNO Sample-A blank specimen, Zwinkels et al. 1990 Table II.

The current Geant4 behavior endpoint-holds:
    lambda < 300 nm -> L_abs(300 nm)
    lambda > 500 nm -> L_abs(500 nm)

We quantify how much a 3-mm shell's PDE-weighted optical upper scale changes
under several deliberately labelled stress variants.

Stress variants are NOT all physical models:
- UV-transparent: no absorption below 300 nm (upper-transmission bound).
- UV-local-edge: continue the local log(alpha) slope defined by same-paper
  Table-I RT points at 300 and 310 nm. This is an extrapolation sensitivity,
  not measured below 300 nm.
- UV-opaque: zero transmission below 300 nm (lower-transmission bound).
- VIS-transparent: no absorption above 500 nm.
- VIS-Zhang proxies: above 500 nm use constant alpha equal to the single
  587.6-nm absorption coefficient of two independent commercial PMMA samples
  from Zhang et al. 2020. These are sample-spread proxies, NOT SNO extensions.

The spectral weight is the Water_Default beta=1 Frank-Tamm spectrum multiplied
by the configured SiPM PDE. It is an optical upper-scale diagnostic, not a
full detector hit prediction.
"""

from __future__ import annotations
from pathlib import Path
import bisect
import math
import sys

ROOT = Path(__file__).resolve().parents[1]
ABSFILE = ROOT / "data" / "PMMA_SNO_UVT_SampleA_ABS_TableII.txt"
WATER_N = ROOT / "data" / "Water_RINDEX_Daimon2007_21p5C.txt"
PDEFILE = ROOT / "data" / "SiPM_PDE.txt"

HC_EV_NM = 1239.8419843320026
THICK_MM = 3.0

# Same-paper Table-I RT absorption coefficients [cm^-1].
ALPHA_300 = 0.410
ALPHA_310 = 0.091

# Independent generic-commercial-PMMA single-point benchmarks at 587.6 nm.
# Source: Zhang et al., Applied Optics 59, 2337-2344 (2020), as distributed
# by refractiveindex.info. These are NOT SNO Sample-A values.
ZHANG_MITSUBISHI_5876_ALPHA_CM = 0.038786
ZHANG_TOMSON_5876_ALPHA_CM = 0.049384


def read2(path: Path):
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


def interp_energy(nodes, lam_nm):
    en = sorted((HC_EV_NM/w, y) for w, y in nodes)
    xs = [e for e, _ in en]
    e = HC_EV_NM/lam_nm
    if e <= xs[0]:
        return en[0][1]
    if e >= xs[-1]:
        return en[-1][1]
    i = bisect.bisect_right(xs, e)
    e0, y0 = en[i-1]
    e1, y1 = en[i]
    return y0 + (y1-y0)*(e-e0)/(e1-e0)


def interp_lambda(nodes, lam_nm):
    xs = [x for x, _ in nodes]
    if lam_nm <= xs[0]:
        return nodes[0][1]
    if lam_nm >= xs[-1]:
        return nodes[-1][1]
    i = bisect.bisect_right(xs, lam_nm)
    x0, y0 = nodes[i-1]
    x1, y1 = nodes[i]
    return y0 + (y1-y0)*(lam_nm-x0)/(x1-x0)


def integrate(f, lo=200.0, hi=900.0, step=0.05):
    x = lo
    prev = f(x)
    s = 0.0
    while x < hi - 1e-12:
        x2 = min(hi, x+step)
        cur = f(x2)
        s += 0.5*(prev+cur)*(x2-x)
        x, prev = x2, cur
    return s


def main():
    absnodes = read2(ABSFILE)
    water = read2(WATER_N)
    pde = read2(PDEFILE)

    if absnodes[0][0] != 300.0 or absnodes[-1][0] != 500.0:
        print("FAIL: expected SNO production ABS domain 300-500 nm")
        return 1

    def weight(lam):
        n = interp_energy(water, lam)
        ft = max(0.0, 1.0 - 1.0/(n*n))/(lam*lam)
        return ft*interp_lambda(pde, lam)

    def current_T(lam, thick=THICK_MM):
        L = interp_energy(absnodes, lam)
        return math.exp(-thick/L)

    # Continue the local log(alpha) slope from same-paper RT values at 300/310.
    log_slope = (math.log(ALPHA_310)-math.log(ALPHA_300))/10.0

    def uv_edge_alpha(lam):
        return math.exp(math.log(ALPHA_300) + log_slope*(lam-300.0))

    def uv_edge_T(lam, thick=THICK_MM):
        if lam < 300.0:
            alpha = uv_edge_alpha(lam)
            return math.exp(-alpha*(thick/10.0))
        return current_T(lam, thick)

    def uv_transparent_T(lam, thick=THICK_MM):
        return 1.0 if lam < 300.0 else current_T(lam, thick)

    def uv_opaque_T(lam, thick=THICK_MM):
        return 0.0 if lam < 300.0 else current_T(lam, thick)

    def vis_transparent_T(lam, thick=THICK_MM):
        return 1.0 if lam > 500.0 else current_T(lam, thick)

    def vis_proxy_T(lam, alpha_cm, thick=THICK_MM):
        if lam > 500.0:
            return math.exp(-alpha_cm*(thick/10.0))
        return current_T(lam, thick)

    base = integrate(weight)
    current = integrate(lambda l: weight(l)*current_T(l))

    policies = [
        ("current endpoint clamp", lambda l: current_T(l)),
        ("UV transparent bound", lambda l: uv_transparent_T(l)),
        ("UV local-edge continuation", lambda l: uv_edge_T(l)),
        ("UV opaque bound", lambda l: uv_opaque_T(l)),
        ("VIS transparent bound", lambda l: vis_transparent_T(l)),
        ("VIS Zhang-Mitsubishi proxy", lambda l: vis_proxy_T(l, ZHANG_MITSUBISHI_5876_ALPHA_CM)),
        ("VIS Zhang-Tomson proxy", lambda l: vis_proxy_T(l, ZHANG_TOMSON_5876_ALPHA_CM)),
    ]

    print("=== PMMA_SNO_A ABS OUT-OF-RANGE STAGE-3A SENSITIVITY ===")
    print(f"production ABS file: {ABSFILE}")
    print("measured/tabulated ABS domain: 300-500 nm")
    print(f"baseline shell thickness: {THICK_MM:g} mm")
    print()
    print("PDE-weighted optical upper-scale transmission:")
    for name, T in policies:
        val = integrate(lambda l: weight(l)*T(l))
        print(f"  {name:29s}: Trel(no-shell)={100*val/base:8.3f}%  "
              f"delta vs current={100*(val/current-1.0):+7.3f}%")

    print()
    print("Current endpoint-clamp post-shell spectral composition:")
    for lo, hi in ((200,300),(300,500),(500,700),(700,900)):
        part = integrate(lambda l: weight(l)*current_T(l), lo, hi)/current
        print(f"  {lo:3.0f}-{hi:3.0f} nm: {100*part:7.3f}%")
    outside = (
        integrate(lambda l: weight(l)*current_T(l), 200,300)
        + integrate(lambda l: weight(l)*current_T(l), 500,900)
    )/current
    print(f"  outside 300-500 nm total: {100*outside:.3f}%")

    print()
    print("UV local-edge stress continuation derived from SNO Table-I RT 300/310 nm:")
    for lam in (299,295,290,285,280,270):
        a = uv_edge_alpha(lam)
        T = math.exp(-a*(THICK_MM/10.0))
        print(f"  {lam:3.0f} nm: alpha_extrap={a:10.4g} cm^-1  T(3mm)={100*T:8.4f}%")
    print("  NOTE: values below 300 nm are extrapolated sensitivity points, NOT SNO measurements.")

    print()
    print("Independent generic-PMMA visible benchmark at 587.6 nm:")
    t_sno = math.exp(-0.007*(THICK_MM/10.0))  # current >500 endpoint alpha
    t_m = math.exp(-ZHANG_MITSUBISHI_5876_ALPHA_CM*(THICK_MM/10.0))
    t_t = math.exp(-ZHANG_TOMSON_5876_ALPHA_CM*(THICK_MM/10.0))
    print(f"  current SNO 500-nm endpoint alpha=0.007 cm^-1 -> T(3mm)={100*t_sno:.4f}%")
    print(f"  Zhang Mitsubishi alpha=0.038786 cm^-1       -> T(3mm)={100*t_m:.4f}%")
    print(f"  Zhang Tomson alpha=0.049384 cm^-1           -> T(3mm)={100*t_t:.4f}%")
    print("  These different PMMA samples are context for material spread, not SNO replacements.")

    # Geometry scaling: show how the same policy uncertainty grows for thick PMMA.
    print()
    print("Geometry scaling: UV policy sensitivity for a 25-mm PMMA path:")
    cur25 = integrate(lambda l: weight(l)*current_T(l,25.0))
    edge25 = integrate(lambda l: weight(l)*uv_edge_T(l,25.0))
    opaque25 = integrate(lambda l: weight(l)*uv_opaque_T(l,25.0))
    print(f"  current clamp             : {100*cur25/base:.3f}% of no-shell upper scale")
    print(f"  UV local-edge continuation: {100*edge25/base:.3f}% "
          f"({100*(edge25/cur25-1):+.3f}% vs current)")
    print(f"  UV opaque bound           : {100*opaque25/base:.3f}% "
          f"({100*(opaque25/cur25-1):+.3f}% vs current)")

    print()
    print("Decision logic:")
    print("  - Above 500 nm, the current endpoint is already almost transparent for 3 mm;")
    print("    visible sample-spread proxies change the integrated upper scale only at the")
    print("    few-0.1% level in this deliberately crude test.")
    print("  - Below 300 nm, continuing the measured UV-edge trend changes the integrated")
    print("    upper scale by several percent, so the low-UV policy dominates.")
    print("  - Therefore the current endpoint clamp below 300 nm must NOT be promoted as")
    print("    physically validated. A production UV extension/cut policy is required.")
    print("  - Do not adopt the local-edge extrapolation as truth; use it only as a sensitivity")
    print("    until a numeric broad-range detector-grade acrylic source is selected.")

    # Numerical regression gate on the sensitivity tool itself.
    current_ratio = current/base
    edge_ratio = integrate(lambda l: weight(l)*uv_edge_T(l))/base
    passed = (
        0.97 < current_ratio < 1.0
        and edge_ratio < current_ratio
        and outside > 0.20
    )
    print()
    print("PMMA ABS OUT-OF-RANGE SENSITIVITY AUDIT: " + ("PASS" if passed else "FAIL"))
    print("PMMA <300 NM PRODUCTION POLICY: NEEDS_PHYSICAL_EXTENSION_OR_EXPLICIT_CUT")
    print("PMMA >500 NM PRODUCTION POLICY: LOWER_PRIORITY_BUT_NOT_SNO_MEASURED")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
