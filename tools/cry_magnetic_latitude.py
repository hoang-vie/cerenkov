#!/usr/bin/env python3
"""Approximate CRY magnetic latitude with the IGRF-14 centered tilted dipole.

CRY 1.7 defines its latitude parameter relative to the magnetic equator/poles
and applies a centered-dipole geomagnetic cutoff.  This helper therefore uses
only the degree-1 IGRF-14 coefficients, matching that level of approximation.
The predictive coefficients below are valid for 2025.0--2030.0.
"""
import argparse
import datetime as dt
import math

# IGRF-14, epoch 2025.0 [nT] and 2025-2030 secular variation [nT/year].
G10_2025, G11_2025, H11_2025 = -29350.0, -1410.3, 4545.5
DG10, DG11, DH11 = 12.6, 10.0, -21.5


def decimal_year(d: dt.date) -> float:
    start = dt.date(d.year, 1, 1)
    end = dt.date(d.year + 1, 1, 1)
    return d.year + (d - start).days / (end - start).days


def tilted_dipole_latitude(lat_deg: float, lon_deg: float, date: dt.date) -> float:
    year = decimal_year(date)
    if not (2025.0 <= year <= 2030.0):
        raise ValueError("This embedded IGRF-14 predictive helper is limited to 2025-2030.")
    dy = year - 2025.0
    g10 = G10_2025 + DG10 * dy
    g11 = G11_2025 + DG11 * dy
    h11 = H11_2025 + DH11 * dy

    # Geographic unit vector and geomagnetic north dipole axis.  The overall
    # sign is irrelevant to CRY's cos^4 cutoff, but this convention returns
    # the conventional northern geomagnetic pole near 81 N, 73 W.
    axis = (-g11, -h11, -g10)
    an = math.sqrt(sum(v*v for v in axis))
    axis = tuple(v/an for v in axis)

    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    r = (math.cos(lat)*math.cos(lon),
         math.cos(lat)*math.sin(lon),
         math.sin(lat))
    dot = max(-1.0, min(1.0, sum(a*b for a,b in zip(axis, r))))
    return math.degrees(math.asin(dot))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lat", type=float, default=13.78, help="geographic latitude [deg]")
    ap.add_argument("--lon", type=float, default=109.22, help="geographic longitude east [deg]")
    ap.add_argument("--date", default="2026-06-18", help="YYYY-MM-DD")
    args = ap.parse_args()
    date = dt.date.fromisoformat(args.date)
    mlat = tilted_dipole_latitude(args.lat, args.lon, date)
    cutoff = 14.8 * math.cos(math.radians(mlat))**4
    print(f"geographic: lat={args.lat:.6f} deg, lon={args.lon:.6f} deg, date={date}")
    print(f"IGRF-14 tilted-dipole magnetic latitude: {mlat:.6f} deg")
    print(f"CRY latitude input (rounded): {mlat:.2f}")
    print(f"CRY dipole proton cutoff approximation: {cutoff:.6f} GeV")


if __name__ == "__main__":
    main()
