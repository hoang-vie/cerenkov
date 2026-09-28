#!/usr/bin/env python3
"""Quantify the angular systematic from G4OpRayleigh's delta=0 dipole law
versus a Cabannes pure-water phase function with depolarization delta=0.039.

This is an analytical systematic estimate; it does not alter Geant4.
"""
from __future__ import annotations
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TABLE=ROOT/"data"/"Water_RAYLEIGH_ZhangHu2009_PMH_21p5C.txt"
DELTA=0.039


def read(path):
    rows=[]
    for raw in path.read_text().splitlines():
        s=raw.strip()
        if not s or s.startswith('#'): continue
        a,b=s.split()[:2]; rows.append((float(a),float(b)))
    return rows


def interp(rows,x):
    if x<=rows[0][0]: return rows[0][1]
    if x>=rows[-1][0]: return rows[-1][1]
    for (x0,y0),(x1,y1) in zip(rows,rows[1:]):
        if x0<=x<=x1:
            t=(x-x0)/(x1-x0); return y0+t*(y1-y0)
    raise RuntimeError

# Relative-to-uniform-sphere density for linearly polarized incidence, where
# x = n_out . epsilon_in.  G4 is delta=0.  Cabannes delta>0 fills the dipole
# null and changes normalization.
def rel_density(x,delta):
    return 3.0/(2.0+delta)*(1.0-(1.0-delta)*x*x)

# Numerical integrations in x are sufficient because the phase law depends
# only on x=n.epsilon for a fixed linear input polarization.
def integrate(fn,a=-1.0,b=1.0,n=400000):
    h=(b-a)/n
    s=0.5*(fn(a)+fn(b))
    for i in range(1,n): s+=fn(a+i*h)
    return s*h

# Uniform-sphere expectation has dx/2, hence total variation = 1/4 integral.
tv=0.25*integrate(lambda x: abs(rel_density(x,DELTA)-rel_density(x,0.0)))

# Unpolarized polar law p(u) proportional 1+f u^2.
f=(1-DELTA)/(1+DELTA)
def cdf_u(u,ff):
    norm=2.0*(1.0+ff/3.0)
    return ((u+u**3*ff/3.0)-(-1.0-ff/3.0))/norm
max_cdf=0.0; at_u=0.0
for i in range(200001):
    u=-1.0+2.0*i/200000.0
    d=abs(cdf_u(u,f)-cdf_u(u,1.0))
    if d>max_cdf: max_cdf,at_u=d,u

# Probability in the two cones around +/- input polarization axis.
def double_cone_fraction(alpha_deg,delta):
    c=math.cos(math.radians(alpha_deg))
    # two cones and uniform sphere dx/2 -> integral c..1 of relative density
    return integrate(lambda x: rel_density(x,delta),c,1.0,n=100000)

rows=read(TABLE)
print("=== WATER RAYLEIGH ANGULAR DELTA SYSTEMATIC ===")
print(f"delta_water={DELTA}")
print("G4OpRayleigh angular sampler corresponds to delta=0 dipole scattering.")
print(f"full polarized phase-function total-variation distance = {100*tv:.3f}%")
print(f"unpolarized polar-CDF max difference = {100*max_cdf:.3f} percentage points at cos(theta)={at_u:.3f}")
for a in (5,10,20):
    p0=double_cone_fraction(a,0.0); pd=double_cone_fraction(a,DELTA)
    print(f"two cones within {a:2d} deg of +/- polarization axis: G4={100*p0:.5f}%  water-delta={100*pd:.5f}%  delta={100*(pd-p0):+.5f} pp")

print("\nDetector-path scale (candidate MFP):")
for lam in (400.0,500.0):
    L=interp(rows,lam)
    for d in (0.10,1.0,5.0):
        psc=1.0-math.exp(-d/L)
        # TV-weighted fraction is a useful scale for how many *all input photons*
        # can be redistributed solely by the phase-law difference after one scatter.
        scale=psc*tv
        print(f" {lam:3.0f} nm d={d:3.1f} m: Pscatter={100*psc:.4f}%  Pscatter*TV={100*scale:.6f}%")

print("\nInterpretation:")
print("  - delta=0.039 is a real water-polarization correction; Geant4 does not model it explicitly.")
print("  - The largest relative difference occurs near the ideal dipole null, but the integrated")
print("    linearly-polarized phase-function difference is only about 1.1% of scattered photons.")
print("  - Combined with the very small molecular-scattering probability over a 10-cm detector")
print("    path, this is a tiny baseline detector systematic; retain it explicitly for long-path")
print("    timing/polarization studies rather than silently claiming exact Cabannes angular physics.")
print("WATER RAYLEIGH ANGULAR DELTA SYSTEMATIC: DOCUMENTED")
