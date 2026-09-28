#!/usr/bin/env python3
"""Independent numerical sanity checks for the CherenkovWater input tables.
No Geant4 dependency. Run from the project root.
"""
from pathlib import Path
import math

DATA = Path(__file__).resolve().parent / "data"
HC_EV_NM = 1239.841984
ALPHA = 1/137.035999177

MUON_MASS_MEV = 105.6583755
MUON_KINETIC_MEV = 3000.0
MUON_GAMMA_3GEV = 1.0 + MUON_KINETIC_MEV / MUON_MASS_MEV
MUON_BETA_3GEV = math.sqrt(1.0 - 1.0/(MUON_GAMMA_3GEV*MUON_GAMMA_3GEV))


def read2(name):
    out=[]
    for line in (DATA/name).read_text().splitlines():
        line=line.strip()
        if not line or line.startswith('#'): continue
        a,b=line.split()[:2]
        out.append((float(a),float(b)))
    return out


def interp(tab, wavelength_nm):
    """Match Geant4 G4MaterialPropertyVector interpolation: linear in photon energy.
    Input tables are wavelength/value; Geant4 converts wavelength to energy first.
    Values outside the tabulated energy range are held at the nearest endpoint.
    """
    etab=sorted((HC_EV_NM/lam, val) for lam,val in tab)
    e=HC_EV_NM/wavelength_nm
    if e <= etab[0][0]: return etab[0][1]
    if e >= etab[-1][0]: return etab[-1][1]
    for (e0,y0),(e1,y1) in zip(etab,etab[1:]):
        if e0 <= e <= e1:
            f=(e-e0)/(e1-e0)
            return y0+f*(y1-y0)
    raise RuntimeError


def theta_c_deg(n,beta=1.0):
    q=1/(beta*n)
    return math.degrees(math.acos(q)) if q <= 1 else float('nan')


def crit_deg(n1,n2):
    if n1 <= n2: return None
    return math.degrees(math.asin(n2/n1))


def frank_tamm_per_cm(rindex, beta=1.0, lo_nm=200, hi_nm=900):
    # integrate dN/dx = 2*pi*alpha * integral (1-1/(beta^2 n^2))/lambda^2 dlambda
    # with lambda in cm. Trapezoidal integration on a 0.1 nm grid.
    step=0.1
    nsteps=int((hi_nm-lo_nm)/step)
    total=0.0
    prev=None
    for i in range(nsteps+1):
        lam=lo_nm+i*step
        n=interp(rindex,lam)
        term=max(0.0,1.0-1.0/(beta*beta*n*n))
        lam_cm=lam*1e-7
        y=2*math.pi*ALPHA*term/(lam_cm*lam_cm) * 1e-7 # per cm per nm
        if prev is not None: total += 0.5*(prev+y)*step
        prev=y
    return total


def weighted_pde_per_cm(rindex,pde,beta=1.0,lo_nm=200,hi_nm=900):
    step=0.1; total=0.0; prev=None
    nsteps=int((hi_nm-lo_nm)/step)
    for i in range(nsteps+1):
        lam=lo_nm+i*step
        n=interp(rindex,lam)
        term=max(0.0,1.0-1.0/(beta*beta*n*n))
        lam_cm=lam*1e-7
        y=2*math.pi*ALPHA*term/(lam_cm*lam_cm)*1e-7*interp(pde,lam)
        if prev is not None: total += 0.5*(prev+y)*step
        prev=y
    return total

# Keep these names synchronized with DetectorConstruction.hh defaults.
water_n=read2('Water_RINDEX_Daimon2007_21p5C.txt')
acr_n=read2('PMMA_RINDEX_Roychowdhury2020_Cauchy3Anchor.txt')
gel_n=read2('OpticalGrease_PDMS_RINDEX_Roychowdhury2018_Cauchy3Anchor.txt')
pde=read2('SiPM_PDE.txt')
water_abs=read2('Water_ABS_Fewell2019_Pure23C.txt')
acr_abs_legacy=read2('Acrylic_ABS.txt')
pmma_abs=read2('PMMA_SNO_UVT_SampleA_ABS_Measured300to500_Model200to900.txt')
ray=read2('Water_RAYLEIGH_ZhangHu2009_PMH_21p5C.txt')
ray_legacy=read2('Water_RAYLEIGH_Lambda4_Baseline.txt')
grease_abs=read2('OpticalGrease_PDMS_ABS_Generic10m.txt')

print('=== Cherenkov angles, beta=1 ===')
for lam in (300,350,400,420,450,500,600):
    print(f'{lam:4.0f} nm : water {theta_c_deg(interp(water_n,lam)):7.3f} deg | acrylic {theta_c_deg(interp(acr_n,lam)):7.3f} deg')

print('\n=== Frank-Tamm yield, beta=1, 200-900 nm ===')
print(f'Water   : {frank_tamm_per_cm(water_n):.1f} photons/cm')
print(f'Acrylic : {frank_tamm_per_cm(acr_n):.1f} photons/cm')
print(f'3 GeV kinetic muon beta = {MUON_BETA_3GEV:.8f}')
print(f'Water, 3 GeV kinetic muon: {frank_tamm_per_cm(water_n, MUON_BETA_3GEV):.1f} photons/cm')
print('\nPDE-weighted upper-scale (NOT a detector hit prediction):')
print(f'Water   : {weighted_pde_per_cm(water_n,pde):.1f} /cm')
print(f'Acrylic : {weighted_pde_per_cm(acr_n,pde):.1f} /cm')

print('\n=== Critical angles (from normal) ===')
for lam in (350,400,450,500,600):
    nw,na,ng=interp(water_n,lam),interp(acr_n,lam),interp(gel_n,lam)
    def fmt(v): return 'none' if v is None else f'{v:.2f} deg'
    print(f'{lam:4.0f} nm : Acrylic->Gel {fmt(crit_deg(na,ng))}; Gel->Water {fmt(crit_deg(ng,nw))}; Gel->Acrylic {fmt(crit_deg(ng,na))}')

print('\n=== Selected attenuation sanity checks ===')
print('NOTE: PMMA_SNO_A uses measured SNO Sample-A values at 300-500 nm;')
print('      200-300 and 500-900 nm are explicit transport-policy extensions, not direct measurements.')
for lam in (250,300,320,330,340,350,375,400,420,450,500,600,700,850):
    wa=interp(water_abs,lam) # m
    aa_legacy=interp(acr_abs_legacy,lam)   # mm
    aa_pmma=interp(pmma_abs,lam)           # mm; active production transport table
    rl=interp(ray,lam)                     # m, production Zhang-Hu PMH
    rl_legacy=interp(ray_legacy,lam)       # m, historical seawater-like baseline
    ga=interp(grease_abs,lam)               # m, explicit baseline assumption
    pwater=math.exp(-0.05/wa)
    plegacy=math.exp(-50.0/aa_legacy)
    ppm=math.exp(-50.0/aa_pmma)
    print(f'{lam:4.0f} nm : Water ABS={wa:10.4g} m, survival/5cm={pwater:.6f}; '
          f'legacy PMMA ABS={aa_legacy:10.4g} mm, surv={plegacy:.6f}; '
          f'PMMA_SNO_A ABS={aa_pmma:9.4g} mm, surv={ppm:.6f}; '
          f'RayleighPMH={rl:9.4g} m; legacyRay={rl_legacy:9.4g} m; grease ABS={ga:g} m')
