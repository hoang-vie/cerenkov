#!/usr/bin/env python3

import math

THICKNESS_MM = 1.0

angles = [0., 30., 45., 60., 70.]

cases = [
    ("0.1 m",   100.0),
    ("1 m",    1000.0),
    ("10 m",  10000.0),
    ("100 m",100000.0),
]

print("=== PDMS STAGE 4D ABS SENSITIVITY ===")
print()
print("Coupling physical thickness = 1.0 mm")
print("10 m is the current generic placeholder, not measured data.")
print()

for deg in angles:

    theta = math.radians(deg)

    path = THICKNESS_MM / math.cos(theta)

    print(
        f"--- incidence = {deg:5.1f} deg, "
        f"path = {path:.6f} mm ---"
    )

    for name,L in cases:

        T = math.exp(-path/L)

        print(
            f"{name:>7s}: "
            f"T={100*T:10.6f}%  "
            f"loss={100*(1-T):10.6f}%"
        )

    print(
        f"{'transparent':>7s}: "
        f"T={100.0:10.6f}%  "
        f"loss={0.0:10.6f}%"
    )

    print()

print("PDMS STAGE 4D ABS SENSITIVITY AUDIT: COMPLETE")
