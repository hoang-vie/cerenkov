#!/usr/bin/env python3

from pathlib import Path
import math

ROOT = Path(__file__).resolve().parents[1]

OUT = ROOT / "macros" / \
    "validation_pdms_interface_stage4c_pmma_to_pdms_scan.mac"

angles = [
    0.0,
    30.0,
    45.0,
    60.0,
    70.0,
    70.5,
    70.67,
    70.68,
]

lines = r'''# PDMS Stage 4C-A
# PMMA -> PDMS interface Fresnel/Snell/TIR validation at 400 nm.
#
# Real detector geometry:
#   Water core half-size       = 50 mm
#   PMMA shell thickness       = 3 mm
#   PMMA/PDMS boundary on +X   = x = 53 mm
#   PDMS coupling              = x = 53..54 mm
#
# One centered +X coupling/readout channel.
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

/detector/couplingSize 25 25 1 mm
/detector/readoutSize 25 25 1 mm

/readout/uniformGrid 0 0 0 0
/readout/planeGrid +X true 1 1 0 0 0 0

/detector/couplingProfile PDMS_Generic
/detector/sensorMode ConcentratorInput

/output/fileName validation_pdms_interface_stage4c_pmma_to_pdms_scan
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

# 400 nm
/generator/energy 3.09960496083 eV

/generator/useAutoY false

# 0.5 mm inside PMMA before its outer +X face.
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
# run{run}: theta_i={deg:g} deg, s polarization
/generator/direction {dx:.15f} {dy:.15f} 0
/generator/polarization 0 0 1
/run/beamOn 100000
'''
    run += 1

    # p polarization, in incidence plane and perpendicular to direction.
    px = -math.sin(th)
    py =  math.cos(th)

    lines += f'''
# run{run}: theta_i={deg:g} deg, p polarization
/generator/direction {dx:.15f} {dy:.15f} 0
/generator/polarization {px:.15f} {py:.15f} 0
/run/beamOn 100000
'''
    run += 1

OUT.write_text(lines)

print(f"Wrote {OUT}")
print(f"runs = {run}")
