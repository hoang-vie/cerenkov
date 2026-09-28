#!/usr/bin/env python3

from pathlib import Path
import bisect
import math

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

HC = 1239.8419843320026
WL = 400.0

ANGLES = [0.,30.,34.,40.,41.5,41.7,41.8,50.]

def read_table(path):
    out=[]
    for raw in path.read_text().splitlines():
        s=raw.strip()
        if not s or s.startswith("#"):
            continue
        a,b=s.split()[:2]
        out.append((float(a),float(b)))
    return out

def interp_energy(rows, wl):
    e=sorted((HC/x,y) for x,y in rows)
    xs=[x for x,_ in e]
    E=HC/wl

    if E<=xs[0]:
        return e[0][1]

    if E>=xs[-1]:
        return e[-1][1]

    i=bisect.bisect_right(xs,E)

    E0,y0=e[i-1]
    E1,y1=e[i]

    f=(E-E0)/(E1-E0)

    return y0+f*(y1-y0)

n1=interp_energy(
    read_table(DATA/"PMMA_RINDEX_Roychowdhury2020_Cauchy3Anchor.txt"),
    WL
)

n2=interp_energy(
    read_table(DATA/"Air_RINDEX.txt"),
    WL
)

critical=math.degrees(math.asin(n2/n1))
brewster=math.degrees(math.atan(n2/n1))

print("=== PMMA -> AIR @ 400 NM ===")
print(f"n_PMMA    = {n1:.9f}")
print(f"n_air     = {n2:.9f}")
print(f"Brewster  = {brewster:.9f} deg")
print(f"critical  = {critical:.9f} deg")
print()

for deg in ANGLES:

    ti=math.radians(deg)
    st=(n1/n2)*math.sin(ti)

    if st >= 1.0:
        print(
            f"{deg:6.1f} deg  "
            f"TIR  R=100%"
        )
        continue

    tt=math.asin(st)

    ci=math.cos(ti)
    ct=math.cos(tt)

    Rs=((n1*ci-n2*ct)/(n1*ci+n2*ct))**2
    Rp=((n2*ci-n1*ct)/(n2*ci+n1*ct))**2

    print(
        f"{deg:6.1f} deg  "
        f"theta_t={math.degrees(tt):10.6f}  "
        f"Rs={100*Rs:10.6f}%  "
        f"Rp={100*Rp:10.6f}%"
    )

print()
print("PMMA STAGE 3E PMMA->AIR ANALYTICAL AUDIT: COMPLETE")
