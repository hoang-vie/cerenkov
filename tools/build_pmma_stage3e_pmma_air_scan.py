#!/usr/bin/env python3

from pathlib import Path
import math

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "macros" / \
    "validation_pmma_interface_stage3e_pmma_to_air_scan.mac"

# Around Brewster (~33.65 deg) and critical (~41.73 deg).
angles = [0.0, 30.0, 34.0, 40.0, 41.5, 41.7, 41.8, 50.0]

lines = r'''# PMMA Stage 3E-D
# PMMA -> Air Fresnel/Snell/TIR scan at 400 nm.
#
# Water core half-size = 50 mm.
# PMMA shell thickness = 3 mm.
# Outer PMMA/Air boundary = x = +53 mm.
#
# Photon starts inside PMMA at x = 52.5 mm.
#
# Incidence plane = X-Y.
# Incident direction = (+cos(theta), +sin(theta), 0).
# s polarization = +Z.
# p polarization = (-sin(theta), +cos(theta), 0).

/control/useDoublePrecision true
/control/verbose 1
/run/verbose 1
/run/numberOfThreads 1

/detector/coreMaterial Water
/detector/coreSize 100 100 100 mm
/detector/useOuterShell true
/detector/outerShellProfile PMMA_SNO_A
/detector/outerShellThickness 3 mm
/detector/setTyvekThickness 0 mm
/readout/uniformGrid 0 0 0 0

/output/fileName validation_pmma_interface_stage3e_pmma_to_air_scan
/output/eventMode fullAll
/output/photonBirthDetail basic
/output/incidentDetail minimal
/output/boundaryDebug true
/output/maxOpticalBounces 0

/process/optical/processActivation Cerenkov false
/process/optical/processActivation Scintillation false
/process/optical/processActivation OpAbsorption false
/process/optical/processActivation OpRayleigh false
/process/optical/processActivation OpMieHG false
/process/optical/processActivation OpBoundary true
/process/optical/processActivation OpWLS false
/process/optical/processActivation OpWLS2 false

/run/initialize

/generator/mode Gun
/generator/particle opticalphoton
/generator/energy 3.09960496083 eV
/generator/useAutoY false
/generator/position 52.5 0 0 mm
/generator/time 0 ns
/run/printProgress 20000

'''

run = 0

for deg in angles:
    th = math.radians(deg)

    dx = math.cos(th)
    dy = math.sin(th)

    # s polarization
    lines += f'''
# run{run}: theta_i={deg:g} deg, s
/generator/direction {dx:.15f} {dy:.15f} 0
/generator/polarization 0 0 1
/run/beamOn 100000
'''
    run += 1

    # p polarization, in X-Y incidence plane.
    px = -math.sin(th)
    py =  math.cos(th)

    lines += f'''
# run{run}: theta_i={deg:g} deg, p
/generator/direction {dx:.15f} {dy:.15f} 0
/generator/polarization {px:.15f} {py:.15f} 0
/run/beamOn 100000
'''
    run += 1

OUT.write_text(lines)

print(f"Wrote {OUT}")
print(f"runs = {run}")
