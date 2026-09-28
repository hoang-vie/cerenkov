#!/usr/bin/env python3
"""Static audit of the project UI-command -> output-metadata contract.

This is intentionally a framework-maintenance check, not a physics validator.
It fails when a new project-owned command is registered without an explicit
reproducibility contract, or when the contract references a ROOT branch that
no longer exists in RunAction.
"""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
MESSENGERS = list((ROOT / "src").glob("*Messenger.cc"))
RUN_ACTION = (ROOT / "src" / "RunAction.cc").read_text()
SOURCES = "\n".join(p.read_text() for p in MESSENGERS)

registered = set(re.findall(r'new\s+G4UIcmd\w*\("(/[^"]+)"', SOURCES))

# Values are tuples of "Tree:Branch" references.  Special tokens beginning
# with "@" document intentionally non-branch semantics.
CONTRACT = {
    "/detector/tyvekPosition": ("RunConfig:TyvekConfig",),
    "/detector/coreMaterial": ("RunConfig:CoreMaterial", "RunConfig:CoreProfile"),
    "/detector/coreProfile": ("RunConfig:CoreProfile",),
    "/detector/coreSize": ("RunConfig:CoreX_mm", "RunConfig:CoreY_mm", "RunConfig:CoreZ_mm"),
    "/detector/scintSize": ("RunConfig:ScintX_mm", "RunConfig:ScintY_mm", "RunConfig:ScintZ_mm"),
    "/detector/readoutSize": ("RunConfig:ReadoutX_mm", "RunConfig:ReadoutY_mm", "RunConfig:ReadoutZ_mm", "ChannelMap:SizeU_mm", "ChannelMap:SizeV_mm"),
    "/detector/sipmSize": ("RunConfig:ReadoutX_mm", "RunConfig:ReadoutY_mm", "RunConfig:ReadoutZ_mm", "@legacy-alias"),
    "/detector/sipmPitch": ("ReadoutPlaneConfig:PitchU_mm", "ReadoutPlaneConfig:PitchV_mm", "@legacy-alias"),
    "/detector/tyvekHoleSize": ("RunConfig:HoleX_mm", "RunConfig:HoleY_mm", "@legacy-z-ignored"),
    "/detector/couplingSize": ("RunConfig:CouplingX_mm", "RunConfig:CouplingY_mm", "RunConfig:CouplingZ_mm"),
    "/detector/outerShellThickness": ("RunConfig:OuterShellThickness_mm",),
    "/detector/useOuterShell": ("RunConfig:UseOuterShell",),
    "/detector/outerShellProfile": ("RunConfig:OuterShellProfile",),
    "/detector/setAcrylicThickness": ("RunConfig:OuterShellThickness_mm", "@legacy-alias"),
    "/detector/useAcrylicShell": ("RunConfig:UseOuterShell", "@legacy-alias"),
    "/detector/setTyvekThickness": ("RunConfig:TyvekThickness_mm",),
    "/detector/setSiPMOffsetY": ("ReadoutPlaneConfig:OffsetV_mm", "@legacy-alias"),
    "/detector/scintGap": ("RunConfig:ScintGap_mm",),
    "/detector/sourceDistance": ("RunConfig:SourceDistance_mm",),
    "/detector/couplingProfile": ("RunConfig:CouplingProfile",),
    "/detector/couplingFile": ("RunConfig:CouplingProfileRindexDataset", "RunConfig:CouplingProfileRindexStatus", "RunConfig:OpticalDataFingerprint_FNV1a64", "@expert-override"),
    "/detector/waterRayleighFile": ("RunConfig:WaterRayleighDataset", "RunConfig:CoreProfileRayleighDataset", "RunConfig:OpticalDataFingerprint_FNV1a64", "@expert-override"),
    "/detector/pmmaAbsFile": ("RunConfig:AcrylicAbsDataset", "RunConfig:OuterShellProfileAbsDataset", "RunConfig:OpticalDataFingerprint_FNV1a64", "@expert-override"),
    "/detector/sensorMode": ("RunConfig:SensorMode",),
    "/readout/uniformGrid": ("ReadoutPlaneConfig:NU", "ReadoutPlaneConfig:NV", "ReadoutPlaneConfig:PitchU_mm", "ReadoutPlaneConfig:PitchV_mm", "ChannelMap:ChannelID"),
    "/readout/planeGrid": ("ReadoutPlaneConfig:Enabled", "ReadoutPlaneConfig:NU", "ReadoutPlaneConfig:NV", "ReadoutPlaneConfig:PitchU_mm", "ReadoutPlaneConfig:PitchV_mm", "ReadoutPlaneConfig:OffsetU_mm", "ReadoutPlaneConfig:OffsetV_mm", "ChannelMap:ChannelID"),
    "/readout/planeAuto": ("ReadoutPlaneConfig:Enabled", "ReadoutPlaneConfig:AutoFill", "ReadoutPlaneConfig:GapU_mm", "ReadoutPlaneConfig:GapV_mm", "ReadoutPlaneConfig:MarginU_mm", "ReadoutPlaneConfig:MarginV_mm", "ChannelMap:ChannelID"),
    "/trigger/scintThreshold": ("RunConfig:ScintThreshold_MeV",),
    "/trigger/coincidenceWindow": ("RunConfig:CoincidenceWindow_ns",),
    "/trigger/integrationWindow": ("RunConfig:IntegrationWindow_ns",),
    "/output/eventMode": ("RunConfig:OutputEventMode",),
    "/output/incidentDetail": ("RunConfig:IncidentDetail",),
    "/output/photonBirthDetail": ("RunConfig:PhotonBirthDetail",),
    "/output/coreBoundaryMode": ("RunConfig:CoreBoundaryMode",),
    "/output/coreBoundaryDetail": ("RunConfig:CoreBoundaryDetail",),
    "/output/boundaryDebug": ("RunConfig:BoundaryDebug",),
    "/output/maxOpticalBounces": ("RunConfig:MaxOpticalBounces",),
    "/output/fileName": ("@storage-destination-not-physics-metadata",),
    "/output/compressionLevel": ("RunConfig:CompressionLevel",),
    "/CRY/file": ("RunConfig:CRYConfig",),
    "/CRY/input": ("RunConfig:CRYConfig",),
    "/CRY/update": ("RunConfig:CRYConfig", "@action-command-resolved-config-recorded"),
    "/CRY/timeMode": ("RunConfig:CRYTimeMode",),
    "/generator/mode": ("RunConfig:GeneratorMode",),
    "/generator/particle": ("RunConfig:GunParticle",),
    "/generator/energy": ("RunConfig:GunEnergy_MeV",),
    "/generator/position": ("RunConfig:GunX_mm", "RunConfig:GunY_mm", "RunConfig:GunZ_mm"),
    "/generator/direction": ("RunConfig:GunDirX", "RunConfig:GunDirY", "RunConfig:GunDirZ"),
    "/generator/polarization": ("RunConfig:GunPolX", "RunConfig:GunPolY", "RunConfig:GunPolZ"),
    "/generator/time": ("RunConfig:GunTime_ns",),
    "/generator/useAutoY": ("RunConfig:GunUseAutoY", "RunConfig:GunY_mm"),
}

errors = []
missing = sorted(registered - CONTRACT.keys())
stale = sorted(CONTRACT.keys() - registered)
if missing:
    errors.append("Registered commands missing from metadata contract:\n  " + "\n  ".join(missing))
if stale:
    errors.append("Contract entries with no registered command:\n  " + "\n  ".join(stale))

created = {
    "RunConfig": set(),
    "ChannelMap": set(),
    "ReadoutPlaneConfig": set(),
}
current = None
for line in RUN_ACTION.splitlines():
    m = re.search(r'CreateNtuple\("(RunConfig|ChannelMap|ReadoutPlaneConfig)"', line)
    if m:
        current = m.group(1)
        continue
    if current and "FinishNtuple" in line:
        current = None
        continue
    if current:
        m = re.search(r'CreateNtuple(?:I|F|D|S)Column\("([^"]+)"\)', line)
        if m:
            created[current].add(m.group(1))

for cmd, refs in sorted(CONTRACT.items()):
    for ref in refs:
        if ref.startswith("@"):
            continue
        tree, branch = ref.split(":", 1)
        if branch not in created.get(tree, set()):
            errors.append(f"{cmd}: contract references missing branch {ref}")

if errors:
    print("\n".join(errors))
    sys.exit(1)

legacy = sum(any(x == "@legacy-alias" for x in refs) for refs in CONTRACT.values())
expert = sum(any(x == "@expert-override" for x in refs) for refs in CONTRACT.values())
print(
    "Configuration metadata contract: OK "
    f"({len(registered)} commands; {legacy} legacy aliases; {expert} expert overrides)"
)
