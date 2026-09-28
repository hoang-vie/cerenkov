#!/usr/bin/env python3

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAC = ROOT / "macros"

policies = {
    "nominal":
        "PMMA_SNO_UVT_SampleA_ABS_Measured300to500_Model200to900_Candidate.txt",

    "uvopaque":
        "PMMA_SNO_UVT_SampleA_ABS_UVOpaqueBelow300_Model200to900_Candidate.txt",

    "vistransparent":
        "PMMA_SNO_UVT_SampleA_ABS_VisibleTransparentAbove500_Model200to900_Candidate.txt",
}

energies = [
    (280.0, 4.42800708690000),
    (400.0, 3.09960496083000),
    (600.0, 2.06640330722000),
]

for tag, dataset in policies.items():

    prefix = f"validation_pmma_abs_stage3d_d2_{tag}"
    path = MAC / f"{prefix}.mac"

    lines = f"""# PMMA Stage 3D-D2: 3-mm absorption transport check
/control/useDoublePrecision true
/control/verbose 1
/run/verbose 1
/run/numberOfThreads 1

/detector/coreMaterial Acrylic
/detector/coreSize 6 100 100 mm
/detector/useOuterShell false
/detector/setTyvekThickness 0 mm
/readout/uniformGrid 0 0 0 0

/detector/pmmaAbsFile {dataset}

/output/fileName {prefix}
/output/eventMode fullAll
/output/photonBirthDetail basic
/output/incidentDetail minimal
/output/boundaryDebug false
/output/maxOpticalBounces 0

/process/optical/processActivation Cerenkov false
/process/optical/processActivation Scintillation false
/process/optical/processActivation OpAbsorption true
/process/optical/processActivation OpRayleigh false
/process/optical/processActivation OpMieHG false
/process/optical/processActivation OpBoundary false
/process/optical/processActivation OpWLS false
/process/optical/processActivation OpWLS2 false

/run/initialize

/generator/mode Gun
/generator/particle opticalphoton
/generator/useAutoY false
/generator/position 0 0 0 mm
/generator/direction 1 0 0
/generator/polarization 0 1 0
/generator/time 0 ns
/run/printProgress 10000

"""

    for i, (wl, E) in enumerate(energies):
        lines += f"""# run{i}: {wl:.0f} nm
/generator/energy {E:.14f} eV
/run/beamOn 50000

"""

    path.write_text(lines)
    print(path.name)

print("PMMA STAGE 3D-D2 MACROS: CREATED")
