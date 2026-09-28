#!/usr/bin/env python3
"""
PMMA_SNO_A Stage-3B: densify only the measured 300-500 nm SNO absorption domain.

No production physics is changed.

This candidate intentionally solves one problem only: the current production
table is sparse, so Geant4 interpolates ABSLENGTH between widely separated
Table-II points. The same Zwinkels-Davidson-Dodd paper provides denser Sample-A
Table-I RT absorption coefficients over 300-405 nm. We merge those exact
published nodes with the compatible Table-II blank-specimen nodes.

We do NOT extend below 300 nm or above 500 nm in this stage.
"""

from pathlib import Path
import bisect, math, sys

ROOT=Path(__file__).resolve().parents[1]
CURRENT=ROOT/"data"/"PMMA_SNO_UVT_SampleA_ABS_TableII.txt"
CAND=ROOT/"data"/"PMMA_SNO_UVT_SampleA_ABS_TableI_RT_plus_TableII_300to500_Candidate.txt"
WATER=ROOT/"data"/"Water_RINDEX_Daimon2007_21p5C.txt"
PDE=ROOT/"data"/"SiPM_PDE.txt"
HC=1239.8419843320026

REF=[
    (300.0,0.410,"TableI+II"),
    (310.0,0.091,"TableI_RT"),
    (320.0,0.061,"TableI_RT"),
    (325.0,0.055,"TableII_blank"),
    (335.0,0.038,"TableI_RT"),
    (350.0,0.025,"TableI+II"),
    (365.0,0.017,"TableI_RT"),
    (385.0,0.010,"TableI_RT"),
    (400.0,0.009,"TableII_blank"),
    (405.0,0.009,"TableI_RT"),
    (450.0,0.009,"TableII_blank"),
    (500.0,0.007,"TableII_blank"),
]

def read2(p):
    out=[]
    for raw in p.read_text().splitlines():
        s=raw.strip()
        if not s or s.startswith("#"): continue
        a,b=s.split()[:2]
        out.append((float(a),float(b)))
    out.sort()
    return out

def e_interp(tab,wl):
    en=sorted((HC/x,y) for x,y in tab)
    xs=[x for x,_ in en]; e=HC/wl
    if e<=xs[0]: return en[0][1]
    if e>=xs[-1]: return en[-1][1]
    i=bisect.bisect_right(xs,e)
    e0,y0=en[i-1]; e1,y1=en[i]
    return y0+(y1-y0)*(e-e0)/(e1-e0)

def l_interp(tab,wl):
    xs=[x for x,_ in tab]
    if wl<=xs[0]: return tab[0][1]
    if wl>=xs[-1]: return tab[-1][1]
    i=bisect.bisect_right(xs,wl)
    x0,y0=tab[i-1]; x1,y1=tab[i]
    return y0+(y1-y0)*(wl-x0)/(x1-x0)

def integ(f,a,b,step=.05):
    x=a; prev=f(x); s=0.0
    while x<b-1e-12:
        x2=min(b,x+step); cur=f(x2)
        s += .5*(prev+cur)*(x2-x)
        x,prev=x2,cur
    return s

def main():
    cur=read2(CURRENT); cand=read2(CAND)
    water=read2(WATER); pde=read2(PDE)

    print("=== PMMA_SNO_A ABS STAGE-3B MEASURED-DOMAIN DENSIFICATION ===")
    print(f"current : {CURRENT.name} ({len(cur)} rows, {cur[0][0]:g}-{cur[-1][0]:g} nm)")
    print(f"candidate: {CAND.name} ({len(cand)} rows, {cand[0][0]:g}-{cand[-1][0]:g} nm)")
    print()

    ok=True
    if len(cand)!=len(REF):
        print(f"FAIL: candidate rows={len(cand)} expected={len(REF)}")
        ok=False

    max_rel=0.0; max_w=None
    print("Exact published-node conversion:")
    for i,(wl,a,src) in enumerate(REF):
        if i>=len(cand): break
        wc,L=cand[i]
        expected=10.0/a
        rel=abs(L/expected-1)
        max_rel=max(max_rel,rel)
        if rel==max_rel: max_w=wl
        if abs(wc-wl)>1e-9: ok=False
        print(f" {wl:5.0f} nm {src:13s}: alpha={a:.6f} cm^-1  "
              f"L_expected={expected:9.4f} mm  L_candidate={L:9.4f} mm")
    print(f"max node conversion mismatch={max_rel:.3e} at {max_w:g} nm")

    # Current sparse vs candidate, both interpreted exactly as Geant4 does:
    # ABSLENGTH interpolation linear in photon energy.
    max_alpha=(-1,None,None,None)
    max_T3=(-1,None,None,None)
    max_T25=(-1,None,None,None)
    wl=300.0
    while wl<=500.000001:
        Lc=e_interp(cur,wl); Ld=e_interp(cand,wl)
        ac=10.0/Lc; ad=10.0/Ld
        r=abs(ac/ad-1)
        t3c=math.exp(-3/Lc); t3d=math.exp(-3/Ld)
        t25c=math.exp(-25/Lc); t25d=math.exp(-25/Ld)
        if r>max_alpha[0]: max_alpha=(r,wl,ac,ad)
        if abs(t3c-t3d)>max_T3[0]: max_T3=(abs(t3c-t3d),wl,t3c,t3d)
        if abs(t25c-t25d)>max_T25[0]: max_T25=(abs(t25c-t25d),wl,t25c,t25d)
        wl += .01

    print()
    print("Current sparse Table-II vs densified same-paper candidate:")
    print(f" max effective-alpha difference={100*max_alpha[0]:.2f}% at {max_alpha[1]:.2f} nm")
    print(f" max direct T difference, 3 mm ={100*max_T3[0]:.3f} pp at {max_T3[1]:.2f} nm")
    print(f" max direct T difference, 25 mm={100*max_T25[0]:.3f} pp at {max_T25[1]:.2f} nm")

    def weight(wl):
        n=e_interp(water,wl)
        return max(0.0,1.0-1.0/(n*n))/(wl*wl)*l_interp(pde,wl)
    w=integ(weight,300,500)
    Tc=integ(lambda l:weight(l)*math.exp(-3/e_interp(cur,l)),300,500)
    Td=integ(lambda l:weight(l)*math.exp(-3/e_interp(cand,l)),300,500)
    print()
    print("300-500 nm PDE-weighted optical upper-scale, 3-mm shell:")
    print(f" current sparse / no-shell  = {100*Tc/w:.6f}%")
    print(f" densified / no-shell       = {100*Td/w:.6f}%")
    print(f" densified vs current        = {100*(Td/Tc-1):+.6f}%")

    print()
    print("Interpretation:")
    print("  - The candidate contains only published same-paper Sample-A nodes.")
    print("  - It removes the known 310/320/... sparse-interpolation mismatch without")
    print("    inventing any absorption values outside the measured 300-500 nm domain.")
    print("  - Detector-scale impact for a 3-mm shell is tiny in the integrated 300-500")
    print("    nm PDE-weighted upper scale, even though local UV transmission can shift")
    print("    by ~0.75 percentage point.")
    print("  - This candidate is suitable for runtime validation and eventual promotion")
    print("    inside the measured domain.")
    print("  - It deliberately does NOT solve the <300 nm or >500 nm production policy.")

    passed = ok and max_rel<1e-8 and max_T3[0]<0.01 and abs(Td/Tc-1)<0.002
    print()
    print("PMMA SNO MEASURED-DOMAIN DENSIFICATION AUDIT: "+("PASS" if passed else "FAIL"))
    print("PMMA 300-500 NM CANDIDATE STATUS: READY_FOR_RUNTIME_VALIDATION")
    print("PMMA OUT-OF-RANGE STATUS: STILL_OPEN")
    return 0 if passed else 1

if __name__=="__main__":
    sys.exit(main())
