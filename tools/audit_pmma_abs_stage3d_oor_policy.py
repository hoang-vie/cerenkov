#!/usr/bin/env python3

from pathlib import Path
import math

ROOT = Path(__file__).resolve().parents[1]

ABS_FILE = ROOT / "data" / \
    "PMMA_SNO_UVT_SampleA_ABS_TableI_RT_plus_TableII_300to500_Candidate.txt"

PDE_FILE = ROOT / "data" / "SiPM_PDE.txt"

SHELL_MM = 3.0


def read_table(path):
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        cols = line.replace(",", " ").split()

        nums = []
        for c in cols:
            try:
                nums.append(float(c))
            except ValueError:
                pass

        if len(nums) >= 2:
            rows.append((nums[0], nums[1]))

    rows.sort()
    return rows


def linear_fit(xs, ys):
    n = len(xs)
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x*x for x in xs)
    sxy = sum(x*y for x, y in zip(xs, ys))

    den = n*sxx - sx*sx
    if den == 0:
        raise RuntimeError("Degenerate fit")

    m = (n*sxy - sx*sy) / den
    b = (sy - m*sx) / n
    return m, b


def fit_log_alpha(rows, lo, hi):
    pts = []

    for wl, L_mm in rows:
        if lo <= wl <= hi:
            alpha_cm = 10.0 / L_mm
            pts.append((wl, math.log(alpha_cm)))

    if len(pts) < 2:
        raise RuntimeError(f"Not enough points in {lo}-{hi} nm")

    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]

    m, b = linear_fit(xs, ys)
    return m, b, len(pts)


def alpha_from_fit(fit, wl):
    m, b, _ = fit
    return math.exp(m*wl + b)


def transmission(alpha_cm, path_mm):
    # path_mm -> cm
    return math.exp(-alpha_cm * path_mm / 10.0)


abs_rows = read_table(ABS_FILE)

if not abs_rows:
    raise SystemExit("No PMMA ABS data read")

wmin = min(w for w, _ in abs_rows)
wmax = max(w for w, _ in abs_rows)

print("=== PMMA ABS STAGE 3D OUT-OF-RANGE POLICY AUDIT ===")
print(f"dataset = {ABS_FILE.name}")
print(f"measured candidate range = [{wmin:.6f}, {wmax:.6f}] nm")

if abs(wmin - 300.0) > 1e-9 or abs(wmax - 500.0) > 1e-9:
    raise SystemExit("ERROR: expected measured range exactly 300-500 nm")

windows = [
    (300.0, 320.0),
    (300.0, 335.0),
    (310.0, 350.0),
]

fits = {}

print("\n--- UV Urbach-like log(alpha) fits ---")

for lo, hi in windows:
    fit = fit_log_alpha(abs_rows, lo, hi)
    fits[(lo, hi)] = fit
    m, b, n = fit

    print(
        f"fit {lo:.0f}-{hi:.0f} nm: "
        f"N={n}  ln(alpha_cm^-1) = {m:.8f}*lambda_nm + {b:.8f}"
    )

targets = [299.0, 290.0, 280.0, 270.0, 250.0]

print("\n--- UV extrapolation spread ---")

for wl in targets:
    vals = []

    print(f"{wl:.1f} nm:")

    for win, fit in fits.items():
        alpha = alpha_from_fit(fit, wl)
        L_mm = 10.0 / alpha
        T = transmission(alpha, SHELL_MM)

        vals.append(alpha)

        print(
            f"  fit {win[0]:.0f}-{win[1]:.0f}: "
            f"alpha={alpha:.6g} cm^-1  "
            f"L={L_mm:.6g} mm  "
            f"T_3mm={100*T:.3f}%"
        )

    spread = max(vals) / min(vals)

    print(f"  alpha spread max/min = {spread:.3f}x")


# exact endpoint
L500 = min(abs_rows, key=lambda x: abs(x[0] - 500.0))[1]
alpha500 = 10.0 / L500

# Literature proxy quoted/discussed by the SNO paper for purer optical-grade PMMA.
alpha_visible_transparent_proxy = 1.0e-4

print("\n--- Visible-side continuation sensitivity ---")
print(
    f"500-nm endpoint: alpha={alpha500:.6g} cm^-1, "
    f"L={L500:.6g} mm"
)

for path_mm in (3.0, 10.0):
    t_hold = transmission(alpha500, path_mm)
    t_trans = transmission(alpha_visible_transparent_proxy, path_mm)

    print(
        f"path={path_mm:.1f} mm: "
        f"500-endpoint constant-alpha T={100*t_hold:.5f}% ; "
        f"alpha=1e-4 proxy T={100*t_trans:.5f}% ; "
        f"delta={100*(t_trans-t_hold):.5f} percentage-point"
    )


# Simple Cherenkov x PDE relevance proxy.
if PDE_FILE.exists():
    pde = read_table(PDE_FILE)

    if len(pde) >= 2:
        print("\n--- SiPM spectral relevance proxy ---")
        print(
            f"PDE table range = "
            f"[{pde[0][0]:.3f}, {pde[-1][0]:.3f}] nm"
        )

        bins = {
            "<300": 0.0,
            "300-500": 0.0,
            ">500": 0.0,
        }

        total = 0.0

        for (x0, y0), (x1, y1) in zip(pde[:-1], pde[1:]):
            xm = 0.5*(x0+x1)
            ym = 0.5*(y0+y1)
            dx = x1-x0

            # Frank-Tamm spectral factor ~1/lambda^2;
            # this is only a relevance proxy, not detector efficiency.
            weight = max(0.0, ym) * dx / (xm*xm)

            total += weight

            if xm < 300.0:
                bins["<300"] += weight
            elif xm <= 500.0:
                bins["300-500"] += weight
            else:
                bins[">500"] += weight

        if total > 0:
            for key in ("<300", "300-500", ">500"):
                print(
                    f"{key:>7}: "
                    f"{100*bins[key]/total:.3f}% "
                    f"of PDE/lambda^2 proxy"
                )

print("\nInterpretation:")
print("  - 300-500 nm remains measured/runtime-validated.")
print("  - UV fits are sensitivity models only; compare their spread before")
print("    adopting any <300 nm production continuation.")
print("  - >500 nm endpoint holding may be tested as an explicitly named")
print("    conservative visible continuation, not labeled measured data.")

print("\nPMMA STAGE 3D OOR POLICY AUDIT: READY_FOR_DECISION")
