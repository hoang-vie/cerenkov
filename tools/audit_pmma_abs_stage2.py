#!/usr/bin/env python3
"""
PMMA_SNO_A Stage-2 audit: SNO Sample-A absorption table, interpolation,
immersion/scattering systematics, and out-of-range exposure.

No production physics is changed by this tool.

Primary reference:
J. C. Zwinkels, W. F. Davidson, C. X. Dodd,
"Optical properties of UV transmitting acrylics for use in a heavy water
Cerenkov detector", Applied Optics 29 (1990) 3240-3248,
DOI 10.1364/AO.29.003240.

Published Table II, Sample-A blank specimen:
 wavelength_nm, alpha_abs [cm^-1], scattering_coeff [cm^-1]
 300  0.410  0.00119
 325  0.055  0.00084
 350  0.025  0.00054
 400  0.009  0.00023
 450  0.009  0.00023
 500  0.007  0.00020

Project conversion:
  ABSLENGTH_mm = 10 / alpha_cm^-1

Independent in-paper cross-check:
Table I gives denser Sample-A absorption coefficients for the RT and RTT'
methods over 300-405 nm. Table-II blank alpha follows the RT result at the
shared wavelengths (e.g. 300, 350, ~400 nm). We use the Table-I RT points to
quantify the distortion introduced by storing sparse ABSLENGTH nodes and then
letting Geant4 interpolate ABSLENGTH linearly in photon energy.

Table II also gives a D2O-immersed specimen, used here only as a material
condition systematic; production remains the blank-specimen reference.
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

# Table II: blank specimen, alpha absorption coefficient [cm^-1],
# and scattering coefficient [cm^-1].
TABLE2_BLANK = [
    (300.0, 0.410, 0.00119),
    (325.0, 0.055, 0.00084),
    (350.0, 0.025, 0.00054),
    (400.0, 0.009, 0.00023),
    (450.0, 0.009, 0.00023),
    (500.0, 0.007, 0.00020),
]

# Table II: Sample-A specimen immersed in D2O; absorption coefficients [cm^-1].
TABLE2_IMMERSED_ALPHA = {
    300.0: 0.488,
    325.0: 0.099,
    350.0: 0.043,
    400.0: 0.013,
    450.0: 0.009,
    500.0: 0.006,
}

# Table I: Sample-A absorption coefficient alpha [cm^-1].
# RT is the method underlying the blank-specimen absorption curve summarized
# in Table II; RTT' is included to expose method dependence.
TABLE1_RT = {
    300.0: 0.410,
    310.0: 0.091,
    320.0: 0.061,
    335.0: 0.038,
    350.0: 0.025,
    365.0: 0.017,
    385.0: 0.010,
    405.0: 0.009,
}
TABLE1_RTT = {
    300.0: 0.424,
    310.0: 0.093,
    320.0: 0.059,
    335.0: 0.037,
    350.0: 0.022,
    365.0: 0.014,
    385.0: 0.007,
    405.0: 0.006,
}


def read2(path: Path):
    rows=[]
    for raw in path.read_text().splitlines():
        s=raw.strip()
        if not s or s.startswith("#"):
            continue
        a,b=s.split()[:2]
        rows.append((float(a),float(b)))
    rows.sort()
    if not rows:
        raise RuntimeError(f"no numeric rows in {path}")
    return rows


def interp_energy(nodes, lam_nm):
    """Linear interpolation in photon energy, with endpoint clamp."""
    en=sorted((HC_EV_NM/w,y) for w,y in nodes)
    xs=[e for e,_ in en]
    e=HC_EV_NM/lam_nm
    if e <= xs[0]:
        return en[0][1]
    if e >= xs[-1]:
        return en[-1][1]
    i=bisect.bisect_right(xs,e)
    e0,y0=en[i-1]; e1,y1=en[i]
    f=(e-e0)/(e1-e0)
    return y0+f*(y1-y0)


def interp_lambda(nodes, lam_nm):
    xs=[x for x,_ in nodes]
    if lam_nm <= xs[0]: return nodes[0][1]
    if lam_nm >= xs[-1]: return nodes[-1][1]
    i=bisect.bisect_right(xs,lam_nm)
    x0,y0=nodes[i-1]; x1,y1=nodes[i]
    f=(lam_nm-x0)/(x1-x0)
    return y0+f*(y1-y0)


def survival_from_Lmm(Lmm, thickness_mm):
    return math.exp(-thickness_mm/Lmm)


def survival_from_alpha_cm(alpha, thickness_mm):
    return math.exp(-alpha*thickness_mm/10.0)


def integrate(func, lo, hi, step=0.05):
    x=lo
    prev=func(x)
    s=0.0
    while x < hi-1e-12:
        x2=min(hi,x+step)
        cur=func(x2)
        s += 0.5*(prev+cur)*(x2-x)
        x,prev=x2,cur
    return s


def main():
    tab=read2(ABSFILE)
    water=read2(WATER_N)
    pde=read2(PDEFILE)

    print("=== PMMA_SNO_A ABSORPTION STAGE-2 AUDIT ===")
    print(f"file: {ABSFILE}")
    print(f"rows={len(tab)} range={tab[0][0]:g}-{tab[-1][0]:g} nm")
    print()

    ok=True
    if len(tab) != len(TABLE2_BLANK):
        print(f"FAIL: expected {len(TABLE2_BLANK)} Table-II nodes, found {len(tab)}")
        ok=False

    max_rel=0.0; max_w=None
    for i,(wr,alpha,scat) in enumerate(TABLE2_BLANK):
        if i >= len(tab): break
        wt,L=tab[i]
        if abs(wt-wr)>1e-9:
            print(f"FAIL row {i}: wavelength project={wt:g} ref={wr:g}")
            ok=False
            continue
        expected=10.0/alpha
        rel=abs(L/expected-1.0)
        if rel>max_rel:
            max_rel,max_w=rel,wr
        print(f"{wr:5.0f} nm: alpha={alpha:.6g} cm^-1  "
              f"expected Labs={expected:.6g} mm  project={L:.6g} mm  "
              f"scat={scat:.6g} cm^-1")
    print(f"max Table-II reciprocal-conversion mismatch={max_rel:.3e} at {max_w:g} nm")
    if max_rel < 1e-7:
        print("Table-II blank alpha -> ABSLENGTH conversion: PASS")
    else:
        print("Table-II blank alpha -> ABSLENGTH conversion: FAIL")
        ok=False

    # Cross-check sparse L interpolation against denser same-paper RT values.
    print()
    print("Same-paper Table-I RT cross-check of current Geant4 ABSLENGTH interpolation:")
    max_alpha_rel=(-1,None,None,None)
    max_t3=(-1,None,None,None)
    max_t25=(-1,None,None,None)
    for lam,alpha_ref in TABLE1_RT.items():
        Lg=interp_energy(tab,lam)
        alpha_g4=10.0/Lg
        rel=abs(alpha_g4/alpha_ref-1.0)
        d3=abs(survival_from_Lmm(Lg,3.0)-survival_from_alpha_cm(alpha_ref,3.0))
        d25=abs(survival_from_Lmm(Lg,25.0)-survival_from_alpha_cm(alpha_ref,25.0))
        if rel>max_alpha_rel[0]: max_alpha_rel=(rel,lam,alpha_g4,alpha_ref)
        if d3>max_t3[0]: max_t3=(d3,lam,alpha_g4,alpha_ref)
        if d25>max_t25[0]: max_t25=(d25,lam,alpha_g4,alpha_ref)
        print(f" {lam:5.0f} nm: alpha_G4eff={alpha_g4:.6f}  "
              f"TableI_RT={alpha_ref:.6f}  delta={100*(alpha_g4/alpha_ref-1):+.2f}%")
    print(f"max |alpha_G4eff/TableI_RT-1|={100*max_alpha_rel[0]:.2f}% "
          f"at {max_alpha_rel[1]:g} nm")
    print(f"max direct transmission difference: 3 mm={100*max_t3[0]:.3f} pp "
          f"at {max_t3[1]:g} nm; 25 mm={100*max_t25[0]:.3f} pp "
          f"at {max_t25[1]:g} nm")
    print("NOTE: this is a sparse-table/interpolation systematic, not a Table-II node error.")

    # Method spread and immersion condition.
    print()
    print("Material/measurement-condition systematics from the same SNO paper:")
    max_method=0.0
    for lam in TABLE1_RT:
        r=TABLE1_RT[lam]; rp=TABLE1_RTT[lam]
        max_method=max(max_method,abs(rp/r-1.0))
    print(f"  max Table-I RT vs RTT' alpha difference over 300-405 nm = {100*max_method:.1f}%")

    max_imm=(-1,None,None,None)
    for lam,alpha_blank,_ in TABLE2_BLANK:
        alpha_imm=TABLE2_IMMERSED_ALPHA[lam]
        Tb=survival_from_alpha_cm(alpha_blank,3.0)
        Ti=survival_from_alpha_cm(alpha_imm,3.0)
        d=abs(Ti-Tb)
        if d>max_imm[0]: max_imm=(d,lam,Tb,Ti)
        print(f" {lam:5.0f} nm, 3-mm shell: blank T={Tb:.6f}  "
              f"D2O-immersed T={Ti:.6f}  delta={100*(Ti-Tb):+.3f} pp")
    print(f"  max blank-vs-immersed 3-mm transmission difference = "
          f"{100*max_imm[0]:.3f} pp at {max_imm[1]:g} nm")

    # Omitted bulk scattering from Table II.
    print()
    print("Published Sample-A scattering coefficient omitted from current PMMA bulk model:")
    max_ps=(-1,None)
    for lam,alpha,scat in TABLE2_BLANK:
        p=1.0-math.exp(-scat*3.0/10.0)
        if p>max_ps[0]: max_ps=(p,lam)
        print(f" {lam:5.0f} nm: s={scat:.6g} cm^-1  "
              f"P(scatter in 3 mm)={100*p:.4f}%  s/alpha={100*scat/alpha:.2f}%")
    print(f"  largest 3-mm scattering probability = {100*max_ps[0]:.4f}% "
          f"at {max_ps[1]:g} nm")
    print("  This is small for the current 3-mm baseline shell, but it is a documented omission.")

    # Endpoint-clamp exposure. This is NOT an error estimate; it only reports
    # how much of the optical spectrum lives outside the measured ABS domain.
    def water_n(lam):
        return interp_energy(water,lam)
    def pde_lam(lam):
        # Diagnostic weighting only; wavelength-linear interpolation matches the
        # existing validate_physics.py style and is not claimed as detector truth.
        return interp_lambda(pde,lam)
    def ft_weight(lam):
        n=water_n(lam)
        return max(0.0,1.0-1.0/(n*n))/(lam*lam)
    def ft_pde_weight(lam):
        return ft_weight(lam)*pde_lam(lam)
    def T_current(lam):
        return survival_from_Lmm(interp_energy(tab,lam),3.0)

    raw_total=integrate(ft_weight,200,900)
    pde_total=integrate(ft_pde_weight,200,900)
    raw_parts=[]
    pde_parts=[]
    post_parts=[]
    for lo,hi in ((200,300),(300,500),(500,700),(700,900)):
        raw_parts.append((lo,hi,integrate(ft_weight,lo,hi)/raw_total))
        pde_parts.append((lo,hi,integrate(ft_pde_weight,lo,hi)/pde_total))
    post_total=integrate(lambda l: ft_pde_weight(l)*T_current(l),200,900)
    for lo,hi in ((200,300),(300,500),(500,700),(700,900)):
        post_parts.append((lo,hi,integrate(lambda l: ft_pde_weight(l)*T_current(l),lo,hi)/post_total))

    print()
    print("Out-of-range exposure diagnostic for current 3-mm endpoint clamp:")
    print("  Fractions below are spectral-weight exposure, NOT a claim about true extrapolation error.")
    print("  interval       raw Cherenkov     PDE-weighted upper-scale   after current 3-mm clamp")
    for r,p,q in zip(raw_parts,pde_parts,post_parts):
        print(f"  {r[0]:3.0f}-{r[1]:3.0f} nm      {100*r[2]:7.3f}%              "
              f"{100*p[2]:7.3f}%                 {100*q[2]:7.3f}%")
    outside_post=post_parts[0][2]+post_parts[2][2]+post_parts[3][2]
    print(f"  current post-shell PDE-weighted fraction outside measured 300-500 nm = "
          f"{100*outside_post:.3f}%")
    print("  Therefore endpoint handling is NOT a negligible bookkeeping detail.")
    print("  An 'opaque outside range' variant would be only a stress-test bound, not a physical model.")

    print()
    print("Interpretation:")
    print("  - Published Table-II blank absorption nodes are correctly converted to ABSLENGTH.")
    print("  - Existing 300-nm Geant4 OpAbsorption runtime validation remains a valid implementation check.")
    print("  - Sparse ABSLENGTH interpolation is not identical to the denser Table-I absorption curve,")
    print("    but for a 3-mm shell the largest checked direct-transmission difference is sub-percent.")
    print("  - Blank-vs-immersed acrylic condition is a larger ~percent-level UV systematic.")
    print("  - The paper's bulk scattering coefficient is currently omitted; its 3-mm probability is tiny.")
    print("  - The major unresolved issue is the endpoint clamp outside 300-500 nm, because a nontrivial")
    print("    fraction of the PDE-weighted optical spectrum lies outside that tabulated ABS domain.")

    passed = (
        ok
        and max_rel < 1e-7
        and max_t3[0] < 0.01
        and max_ps[0] < 0.001
    )
    print()
    print("PMMA SNO TABLE-II ABS NODE/IMPLEMENTATION AUDIT: " + ("PASS" if passed else "FAIL"))
    print("PMMA SNO ABS INTERPOLATION SYSTEMATIC: DOCUMENTED")
    print("PMMA SNO ABS OUT-OF-RANGE POLICY: REQUIRES_DECISION")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
