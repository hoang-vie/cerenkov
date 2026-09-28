#!/usr/bin/env python3

from pathlib import Path
import math

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "macros" / \
    "validation_pmma_interface_stage3e_pmma_to_water_scan.mac"

angles = [0.0, 30.0, 42.0, 60.0, 63.0, 63.3, 63.5, 70.0]

lines = r'''# PMMA Stage 3E-C
# PMMA -> Water Fresnel/Snell/TIR scan at 400 nm.
#
# Photon starts in the +X PMMA shell at x=50.5 mm and travels inward.
# Interface Water/PMMA is x=50 mm.
#
# Incidence plane = X-Y.
# Incident direction = (-cos(theta), +sin(theta), 0).
# s polarization = +Z.
# p polarization = (+sin(theta), +cos(theta), 0).

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

/output/fileName validation_pmma_interface_stage3e_pmma_to_water_scan
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

# 0.5 mm inside the PMMA shell.
/generator/position 50.5 0 0 mm

/generator/time 0 ns
/run/printProgress 20000

'''

run = 0

for deg in angles:
    th = math.radians(deg)

    dx = -math.cos(th)
    dy =  math.sin(th)

    # s
    lines += f'''
# run{run}: theta_i={deg:g} deg, s
/generator/direction {dx:.15f} {dy:.15f} 0
/generator/polarization 0 0 1
/run/beamOn 100000
'''
    run += 1

    # p: perpendicular to (-cos,sin,0), within X-Y
    px = math.sin(th)
    py = math.cos(th)

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
