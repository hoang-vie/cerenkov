#!/usr/bin/env python3
"""Independent audit of Water_ABS_Fewell2019_Pure23C.txt.

Reference values are the adopted pure-water absorption coefficients from
Fewell & von Trojan, Applied Optics 58, 2408-2421 (2019), Table 1, at 23 C.
The project stores ABSLENGTH = 1/alpha in metres.

This script deliberately distinguishes:
  (1) node-by-node provenance/conversion correctness;
  (2) published coefficient uncertainty;
  (3) interpolation-policy sensitivity between tabulated nodes.

It does not claim detector water equals the paper's pure-water reference.
"""
from __future__ import annotations

import math
from pathlib import Path
import sys

HC_EV_NM = 1239.8419843320026

# wavelength_nm, alpha_m^-1, 1-sigma uncertainty_m^-1
# Fewell & von Trojan 2019, Table 1. Parenthesized uncertainties in the paper
# are expanded here explicitly.
REF = [
    (200.0,0.602,0.030),(205.0,0.234,0.026),(210.0,0.160,0.024),(215.0,0.137,0.017),
    (220.0,0.101,0.014),(225.0,0.075,0.011),(230.0,0.070,0.010),(235.0,0.065,0.010),
    (240.0,0.057,0.009),(245.0,0.042,0.006),(250.0,0.059,0.018),(260.0,0.052,0.012),
    (270.0,0.044,0.008),(280.0,0.0223,0.0033),(290.0,0.0094,0.0020),
    (300.0,0.0047,0.0010),(310.0,0.0024,0.0006),(320.0,0.0015,0.0004),
    (330.0,0.0011,0.0004),(340.0,0.00085,0.00026),(350.0,0.00089,0.00033),
    (360.0,0.00115,0.00033),(370.0,0.00124,0.00033),(380.0,0.00143,0.00029),
    (390.0,0.00170,0.00030),(400.0,0.00222,0.00027),(410.0,0.00266,0.00032),
    (420.0,0.00312,0.00030),(430.0,0.00376,0.00031),(440.0,0.0052,0.0004),
    (450.0,0.0081,0.0005),(460.0,0.0091,0.0005),(470.0,0.0103,0.0004),
    (480.0,0.0121,0.0005),(490.0,0.0146,0.0005),(500.0,0.0207,0.0007),
    (510.0,0.0330,0.0012),(520.0,0.0392,0.0013),(530.0,0.0424,0.0012),
    (540.0,0.0475,0.0014),(550.0,0.0563,0.0018),(560.0,0.0619,0.0010),
    (570.0,0.0694,0.0011),(580.0,0.0896,0.0012),(590.0,0.1354,0.0012),
    (600.0,0.2232,0.0017),(610.0,0.2653,0.0019),(620.0,0.2760,0.0025),
    (630.0,0.2917,0.0027),(640.0,0.3106,0.0028),(650.0,0.340,0.003),
    (660.0,0.410,0.004),(670.0,0.439,0.004),(680.0,0.465,0.004),(690.0,0.516,0.004),
    (700.0,0.624,0.006),(710.0,0.828,0.007),(720.0,1.236,0.011),(725.0,1.496,0.006),
    (729.9,1.97,0.05),(739.6,2.79,0.04),(759.9,2.87,0.04),(770.4,2.82,0.04),
    (780.0,2.69,0.04),(789.9,2.46,0.04),(800.0,2.24,0.04),(825.1,2.62,0.05),
    (850.3,4.21,0.05),(875.7,5.04,0.05),(899.3,6.35,0.06),
]


def read2(path: Path):
    out=[]
    for raw in path.read_text().splitlines():
        s=raw.strip()
        if not s or s.startswith('#'):
            continue
        p=s.split()
        if len(p) < 2:
            continue
        out.append((float(p[0]), float(p[1])))
    return out


def lin(x, x0, y0, x1, y1):
    if x1 == x0:
        return y0
    return y0 + (y1-y0)*(x-x0)/(x1-x0)


def interp_energy(wl_nm: float, nodes):
    # Geant4 MPT interpolation variable is photon energy. Sort ascending energy.
    e=HC_EV_NM/wl_nm
    en=sorted((HC_EV_NM/w, y) for w,y in nodes)
    if e <= en[0][0]: return en[0][1]
    if e >= en[-1][0]: return en[-1][1]
    for (e0,y0),(e1,y1) in zip(en,en[1:]):
        if e0 <= e <= e1:
            return lin(e,e0,y0,e1,y1)
    raise RuntimeError('interpolation bracket failure')


def main():
    root=Path(__file__).resolve().parents[1]
    path=root/'data'/'Water_ABS_Fewell2019_Pure23C.txt'
    if not path.exists():
        print(f'ERROR: missing {path}')
        return 2
    tab=read2(path)
    print('=== WATER ABSORPTION: FEWELL 2019 INDEPENDENT AUDIT ===')
    print(f'file: {path}')
    print(f'rows={len(tab)} range={tab[0][0]:g}-{tab[-1][0]:g} nm')

    ok=True
    if len(tab) != len(REF):
        print(f'FAIL: expected {len(REF)} Table-1 nodes in project range, found {len(tab)}')
        ok=False

    max_rel=0.0; max_abs=0.0; max_w=None
    for i,(wr,a,sig) in enumerate(REF):
        if i >= len(tab): break
        wt,L=tab[i]
        if abs(wt-wr)>1e-9:
            print(f'FAIL: row {i} wavelength project={wt} ref={wr}')
            ok=False
            continue
        expected=1.0/a
        rel=abs(L/expected-1.0)
        ab=abs(L-expected)
        if rel>max_rel:
            max_rel,max_abs,max_w=rel,ab,wr
    print(f'max node reciprocal mismatch={max_rel:.3e} relative ({max_abs:.6g} m) at {max_w:g} nm')
    if max_rel > 1e-7:
        print('FAIL: project ABSLENGTH nodes are not reciprocal Table-1 coefficients within rounding tolerance')
        ok=False
    else:
        print('node-by-node Table-1 reciprocal conversion: PASS')

    vis=[r for r in REF if 300.0 <= r[0] <= 600.0]
    worst=max(vis,key=lambda r:r[2]/r[1])
    w,a,sig=worst
    print(f'largest published relative 1-sigma alpha uncertainty in 300-600 nm={100*sig/a:.1f}% at {w:g} nm')
    print(f'  L_abs={1/a:.6g} m; same first-order relative uncertainty applies to L_abs')

    # For a direct 10-cm path, propagate coefficient uncertainty to survival P=exp(-alpha*x).
    x=0.10
    vals=[]
    for w,a,sig in vis:
        p=math.exp(-a*x)
        sigp=p*x*sig
        vals.append((sigp,w,p,a,sig))
    sigp,w,p,a,sig=max(vals)
    print(f'max 1-sigma survival uncertainty for direct 10 cm path in 300-600 nm={sigp:.3e} '
          f'({100*sigp:.4f} percentage points) at {w:g} nm')

    # Interpolation-policy sensitivity: current Geant4 behavior interpolates stored L_abs
    # linearly in photon energy. Compare to a different, equally plausible table use:
    # interpolate the published alpha values linearly in photon energy, then invert.
    # This is a sensitivity estimate, not a claim that one interpolation is uniquely correct.
    lnodes=[(w,L) for w,L in tab]
    anodes=[(w,a) for w,a,_ in REF]
    def compare_range(lo,hi,path_m):
        n=50000
        best_rel=(-1,None,None,None)
        best_p=(-1,None,None,None)
        for i in range(n+1):
            q=lo+(hi-lo)*i/n
            Lg=interp_energy(q,lnodes)
            ai=interp_energy(q,anodes)
            La=1.0/ai
            rel=abs(Lg/La-1.0)
            if rel>best_rel[0]: best_rel=(rel,q,Lg,La)
            pg=math.exp(-path_m/Lg); pa=math.exp(-path_m/La)
            dp=abs(pg-pa)
            if dp>best_p[0]: best_p=(dp,q,pg,pa)
        return best_rel,best_p

    rel10,p10=compare_range(300.0,600.0,0.10)
    rel1m,p1m=compare_range(300.0,600.0,1.0)
    print('interpolation-policy sensitivity in 300-600 nm:')
    print(f'  max |Delta L|/L = {100*rel10[0]:.2f}% at {rel10[1]:.3f} nm')
    print(f'  max direct-path survival difference at 10 cm = {p10[0]:.6f} '
          f'({100*p10[0]:.3f} percentage points) at {p10[1]:.3f} nm')
    print(f'  max direct-path survival difference at 1 m  = {p1m[0]:.6f} '
          f'({100*p1m[0]:.3f} percentage points) at {p1m[1]:.3f} nm')
    print('  NOTE: this is an interpolation-model systematic, not a node-value error.')

    print('WATER ABS FEWELL AUDIT: ' + ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1

if __name__ == '__main__':
    sys.exit(main())
