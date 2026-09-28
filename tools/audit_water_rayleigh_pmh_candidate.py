#!/usr/bin/env python3
from __future__ import annotations
import math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
CAND=ROOT/"data"/"Water_RAYLEIGH_ZhangHu2009_PMH_21p5C.txt"
LEGACY=ROOT/"data"/"Water_RAYLEIGH_Lambda4_Baseline.txt"

KB=1.380649e-23
DELTA=0.039
TC=21.5

# Daimon & Masumura 2007, 21.5 C, published four-term Sellmeier coefficients.
B=(5.689093832e-1,1.719708856e-1,2.062501582e-2,1.123965424e-1)
C=(5.110301794e-3,1.825180155e-2,2.624158904e-2,1.067505178e1)

# Zhang & Hu 2009 Table 1 PMH beta(90) values at 20 C, delta=0.039,
# in units m^-1 sr^-1.
ZH_TABLE1_PMH={
    366.0:4.500e-4,
    405.0:2.903e-4,
    436.0:2.120e-4,
    546.0:0.824e-4,
    578.0:0.650e-4,
}
MOREL_MEAS={
    366.0:4.53e-4,
    405.0:2.90e-4,
    436.0:2.12e-4,
    546.0:0.835e-4,
    578.0:0.660e-4,
}

def n_daimon(lam_nm):
    u=lam_nm*1e-3
    l2=u*u
    return math.sqrt(1.0+sum(b*l2/(l2-c) for b,c in zip(B,C)))

def kappa_kell_kpa_inv(t):
    return (
        50.88496e-8 + 6.163813e-9*t + 1.459187e-11*t*t
        +20.08438e-14*t**3 -58.47727e-17*t**4 +410.4110e-20*t**5
    )/(1+19.67348e-3*t)

def pmh_d(n):
    n2=n*n
    return (n2-1.0)*(1.0+(2.0/3.0)*(n2+2.0)*((n2-1.0)/(3.0*n))**2)

def beta90(lam_nm,tc=TC,delta=DELTA):
    n=n_daimon(lam_nm)
    d=pmh_d(n)
    f=(6+6*delta)/(6-7*delta)
    kap=kappa_kell_kpa_inv(tc)*1e-3
    lam=lam_nm*1e-9
    return math.pi**2/(2*lam**4)*d*d*KB*(tc+273.15)*kap*f

def b_total(lam_nm,tc=TC,delta=DELTA):
    return (8*math.pi/3)*beta90(lam_nm,tc,delta)*(2+delta)/(1+delta)

def read(path):
    rows=[]
    for raw in path.read_text().splitlines():
        s=raw.strip()
        if not s or s.startswith("#"): continue
        a,b=s.split()[:2]
        rows.append((float(a),float(b)))
    return rows

def interp(rows,x):
    rows=sorted(rows)
    if x<=rows[0][0]: return rows[0][1]
    if x>=rows[-1][0]: return rows[-1][1]
    for (x0,y0),(x1,y1) in zip(rows,rows[1:]):
        if x0<=x<=x1:
            f=(x-x0)/(x1-x0)
            return y0+f*(y1-y0)
    raise RuntimeError

cand=read(CAND)
legacy=read(LEGACY)

print("=== WATER RAYLEIGH STAGE-3B CANDIDATE AUDIT ===")
print(f"candidate={CAND}")
print(f"rows={len(cand)} range={cand[0][0]:g}-{cand[-1][0]:g} nm")
print(f"T={TC:g} C, delta={DELTA:g}")
print()

# Candidate file consistency
max_rel=0.0; max_lam=None
for lam,L in cand:
    expected=1.0/b_total(lam)
    rel=abs(L/expected-1.0)
    if rel>max_rel: max_rel,max_lam=rel,lam
print(f"candidate file vs analytic PMH formula max relative mismatch={max_rel:.3e} at {max_lam:g} nm")

# Reproduce Zhang-Hu Table 1 at 20 C.
print()
print("Zhang-Hu 2009 Table-1 PMH reproduction at 20 C (using Daimon n as optical input):")
max_table=0.0
for lam,ref in ZH_TABLE1_PMH.items():
    calc=beta90(lam,20.0,DELTA)
    rel=(calc/ref-1.0)
    max_table=max(max_table,abs(rel))
    meas=MOREL_MEAS[lam]
    print(f" {lam:5.0f} nm: calc beta90={calc*1e4:7.4f}e-4  "
          f"ZH-PMH={ref*1e4:7.4f}e-4  diff={100*rel:+.3f}%  "
          f"Morel meas={meas*1e4:7.4f}e-4")
print(f"max |calc-ZH Table1|/ZH = {100*max_table:.3f}%")

# Legacy-vs-candidate
print()
print("Legacy baseline vs candidate total scattering:")
for lam in (300.0,350.0,400.0,450.0,500.0,550.0,600.0,700.0):
    Lold=interp(legacy,lam)
    Lnew=interp(cand,lam)
    bold=1/Lold; bnew=1/Lnew
    print(f" {lam:5.0f} nm: legacy L={Lold:9.3f} m  candidate L={Lnew:9.3f} m  "
          f"legacy/candidate b={bold/bnew:6.3f}")

print()
print("Detector-scale direct scattering probabilities:")
for lam in (400.0,500.0):
    Lold=interp(legacy,lam); Lnew=interp(cand,lam)
    for d in (0.10,1.0,5.0):
        po=1-math.exp(-d/Lold); pn=1-math.exp(-d/Lnew)
        print(f" {lam:3.0f} nm d={d:3.1f} m: legacy={100*po:.4f}%  candidate={100*pn:.4f}%  "
              f"delta={100*(po-pn):+.4f} pp")

print()
print("IMPORTANT angular-model caveat:")
print("  G4OpRayleigh uses its built-in polarized Rayleigh angular sampling and does")
print("  not consume a material depolarization-ratio property. Therefore this candidate")
print("  fixes the mean-free-path normalization/spectrum only; delta=0.039 is not yet")
print("  independently represented in the Geant4 scattering phase function.")

ok=(max_rel<1e-8 and max_table<0.002)
print()
print("WATER RAYLEIGH PMH CANDIDATE NUMERICAL AUDIT: "+("PASS" if ok else "FAIL"))
print("PRODUCTION DEFAULT SWITCH: PROMOTED_IN_STAGE3E")
sys.exit(0 if ok else 1)
