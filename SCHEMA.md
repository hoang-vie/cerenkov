# General production schema

The simulation stores particle/track facts rather than particle-specific analysis variables. No ntuple is tied to muons, electrons, pions, or a fixed number of readout channels.

## Output modes

- `thinAll` (default): small `TriggerEvent`/`ScintDeposit` rows plus lightweight `PrimaryTruth` rows for every generated event; detailed rows only for triggered events.
- `triggered`: only triggered events are written.
- `fullAll`: detailed rows for every event. Use for validation/debug samples, not large CRY production.

High-multiplicity optional detail:
- `/output/photonBirthDetail basic|full`: `basic` stores no birth XYZ/direction; `full` does.
- `/output/incidentDetail auto|minimal|phaseSpace`: `auto` is `minimal` for DirectSiPM and `phaseSpace` for ConcentratorInput.
- `/output/coreBoundaryMode off|firstExit`: OFF by default. `firstExit` stores at most one `CoreBoundaryCrossing` row per photon, only after the first successful optical transmission from Core to a non-Core volume.
- `/output/coreBoundaryDetail basic|phaseSpace`: both modes store first-exit position/direction/time/wavelength; `phaseSpace` additionally stores post-boundary polarization. Polarization branches remain present and are NaN in `basic` mode.
- `/output/boundaryDebug false|true`: full boundary-crossing rows are OFF by default.
- `/output/maxOpticalBounces 0`: 0 means unlimited. A positive value is a safety cut and is recorded per photon.

## Identity and joins

Within one event, `(EventID, TrackID)` identifies a Geant4 track. Optical photon `PhotonID` is its Geant4 track ID. Each output file is run-specific (`*_run<RunID>.root`), and `RunID` is also stored in every ntuple.

`Photon.ParentTrackID` is the immediate Geant4 parent. If it resolves directly to a row in `Track`, that row is the immediate non-optical emitter. This is the field used for **primary-emitter** selections: the immediate parent must have `ParentID==0` and the selected PDG.

`Photon.SourceTrackID` is the non-optical source track propagated through an optical chain. For a direct Cherenkov/scintillation photon it is the immediate non-optical parent; for a later WLS photon it is inherited from the absorbed optical parent. Following the source track's `RootPrimaryTrackID` gives the **primary-initiated cascade** selection.

Example: `primary muon -> delta electron -> Cherenkov photon -> WLS photon`. The Cherenkov photon's immediate parent is the electron; the WLS photon's immediate parent is optical, while `SourceTrackID` remains the electron and that source track's `RootPrimaryTrackID` points back to the primary muon.

These definitions are intentionally separate. `RootPrimaryTrackID==0` or `SourceTrackID==0` means ancestry could not be resolved; the simulation does not substitute a plausible-looking ancestor ID.

## Ntuples

### `RunConfig`
One metadata row per run. `general-v1.10` preserves the complete `general-v1.9` RunConfig prefix (columns 0--156) unchanged and appends `GitCommit`, `GitBranch`, and `GitDirty` at columns 157--159. `general-v1.9` had appended `CoreBoundaryMode`, `CoreBoundaryDetail`, and `CoreBoundarySemantic` at 154--156; the complete `general-v1.8` prefix (0--153) and older `general-v1.6` prefix (0--98) therefore remain unchanged for backward analysis compatibility. Columns 99--147 contain authoritative profile-aware metadata for the active `Core`, `OuterShell`, and `Coupling` bulk-material roles. Columns 148--153 contain explicitly declared empirical absorption-validity ranges (`*AbsMeasuredMin_nm`, `*AbsMeasuredMax_nm`) for those three roles. `*AbsValidMin_nm` / `*AbsValidMax_nm` describe numerical/tabulated transport support, whereas `*AbsMeasuredMin_nm` / `*AbsMeasuredMax_nm` describe direct empirical validity; an undeclared measured range is stored as `(-1,-1)` nm.

The legacy `Acrylic*`, `Water*`, and older coupling fields remain present in columns 0--98 so existing analysis continues to work. For new code, the appended role/profile fields are authoritative. R2E is a metadata/schema refactor only: it does not broaden which profiles are permitted for OuterShell or Coupling and does not alter any material table, geometry, or optical process. The pre-existing exact-file `OpticalDataFingerprint_FNV1a64` remains an optical-data-bundle fingerprint; profile selection is recorded explicitly by the new role fields.

RunConfig also contains geometry, trigger settings, output policy, sensor model/operating point, source/build revision (`ProjectRevision`, full Git commit, branch and dirty flag), generator mode/CRY configuration, resolved generic-gun configuration, Geant4 version, physics-list label, compression, channel count, RNG seed, actual thread count, RNG engine, CRY version, default production cut, and the run-time activation state of Cerenkov, Scintillation, absorption, Rayleigh, Mie, boundary, WLS, and WLS2 optical processes.

### `RunStatistics`
One row per event-processing worker: `RunID`, `WorkerThreadID`, generated, triggered, detailed event counts, and CRY simulated time. Sum all rows to obtain run totals. This layout is compatible with ROOT ntuple merging in Geant4 MT and preserves the trigger-efficiency denominator even in `triggered` output mode.

### `TriggerEvent`
Small event row for trigger QA, denominator accounting, and event-level summaries:
`EventID`, pass flag, top/bottom Edep, top/bottom pulse time, `NPrimaryNonOptical`, `NNonOpticalTracks`, `NOpticalTracks`, incident/detected counts, and total Core Edep. The explicit names avoid implying that optical-primary validation tracks are rows in the generic non-optical `Track` table.

### `ScintDeposit`
Sparse raw non-zero scintillator step deposits for events written to the thin tier. Each row stores `PaddleID`, `TrackID`, pre-step global time and deposited energy. Because clustering is not frozen into this table, integration-window, energy-threshold and top/bottom coincidence definitions can be rescanned offline without rerunning Geant4.

### `PrimaryTruth`
Lightweight **non-optical** primary rows aligned with the thin event tier: EventID, TrackID, PDG, vertex kinetic energy/time/position and direction. In `thinAll` this preserves unbiased CRY/generator primary spectra for rejected as well as triggered events without writing secondary or optical detail.

### `Track`
One row per non-optical Geant4 track.  Since hotfix6 registration occurs in
`TrackingAction::PreUserTrackingAction` and final state is updated in
`PostUserTrackingAction`, the truth row no longer depends on a track reaching
its first transport step:
- identity/ancestry: TrackID, ParentID, RootPrimaryTrackID, PDG, charge;
- creator process type/subtype/name;
- birth volume and vertex KE/time/position/direction;
- first Core entry and last Core exit.  These are crossings of the Core
  **logical volume**; they are not detector-envelope crossings when an inside
  daughter such as Tyvek is present;
- end time, **proper time**, end position/volume/process/status.

Muon lifetime validation is derived from this generic table; no dedicated muon tree is required.

### `TrackVolume`
Sparse `(TrackID, VolumeID)` aggregate for every non-optical track: local deposited energy `Edep`, summed kinetic-energy change `DeltaKE = sum(KE_pre-KE_post)`, path length and step count. `Edep/L` is a local/restricted-deposition observable; for clean through-going charged-particle stopping-power validation, `DeltaKE/L` is the more direct quantity.

### `Photon`
Exactly one compact summary row per optical photon:
- ParentTrackID and SourceTrackID;
- optical creator process type/subtype;
- origin volume, wavelength, birth time, emission angle (NaN when not defined);
- optional birth XYZ/direction;
- first Core-wall encounter;
- first extraction toward readout;
- objective transport counters: boundary reflections, TIR, Tyvek reflections, Rayleigh scatters;
- end reason/volume/time and safety-cut flag;
- first-Rayleigh track length and final track length for runtime mean-free-path
  validation.

EndReason codes:
0 unknown/not classified; 1 detected; 2 bulk `OpAbsorption`; 3 world escape; 4 user safety bounce cut; 5 boundary absorption; 6 other stop.


`EmitterBeta` is NaN for non-Cerenkov photons. For Cerenkov photons it stores the exact step-average beta, `0.5*(beta_pre+beta_post)`, used by the current `G4Cerenkov::PostStepDoIt` implementation when sampling `cos(theta)=1/(beta*n(lambda))`. It is appended after the pre-existing Photon columns to preserve their indices.

### `SensorIncident`
One row each time a selected photon enters a readout element from its coupling
volume.  Always stores wavelength/time, legacy local XYZ relative to the
readout-volume centre, entrance-plane XYZ (`EntranceZ_mm=0` at the nominal
entrance), signed local incident angle, `Forward`, and transport counters.
Local +z points outward on all six faces; `Forward==1` therefore requires
`LocalDirZ>0`.  Full local direction and polarization are stored in
`phaseSpace` mode and are NaN in `minimal` mode.

For `ConcentratorInput`, `phaseSpace` is the automatic default and this ntuple
is the intended concentrator input.

### `CoreBoundaryCrossing`
Optional first-exit optical phase space at the **Core logical-volume boundary**. The tree always exists in the schema but has zero rows when `/output/coreBoundaryMode off` (default). In `firstExit` mode, each optical photon contributes at most one row: the first boundary step satisfying Core -> non-Core plus an optical transmission status (`FresnelRefraction` or `Transmission`). A mere wall encounter, Fresnel reflection, or TIR does **not** create a row, so this is intentionally stricter than `Photon.ReachedCoreWall/FirstWall*`.

The recorded state is the **post-boundary** state, i.e. the phase space with which the photon begins downstream optical transport after leaving Core. Columns are RunID/EventID; PhotonID; ParentTrackID/SourceTrackID/RootPrimaryTrackID; OriginVolumeID; FaceID; global XYZ; face-local UV; global post-boundary direction; `CosThetaNormal` relative to the outward Core-face normal; wavelength; global time; and optional post-boundary polarization XYZ. Face IDs reuse the readout convention `0=+X`, `1=-X`, `2=+Z`, `3=-Z`, `4=+Y`, `5=-Y`. Local U/V axes are exactly the same axes generated by `BuildFaceRotation` for the corresponding readout face.

This tree provides the intermediate checkpoint needed to separate Core escape from downstream collection and sensor detection: `Photon` (generated) -> `CoreBoundaryCrossing` (escaped Core) -> `SensorIncident` (reached readout) -> `DetectedSignal` (detected).
Like `Photon`, `SensorIncident`, and other detailed optical trees, `CoreBoundaryCrossing` follows the existing detailed-event tier: `fullAll` writes all events, while `thinAll`/`triggered` retain detailed rows only for triggered events. Use `fullAll` for unbiased Core-escape efficiency studies.

### `DetectedSignal`
DirectSiPM only: compact `(EventID, PhotonID, ChannelID, Time)` detection record. Wavelength, emitter species, angle, and history are joined from `Photon`/`SensorIncident` instead of duplicated.

### `ChannelEvent`
Sparse event-channel aggregate for active channels: NIncident, NDetected, first incident/detected times. There are no fixed `Hit_SiPM_0...23` columns.

### `ChannelMap`
One row per dynamically generated readout element: face, grid indices, local
UV, global readout-volume centre, global entrance-plane centre, normal and
element size.

Face IDs: `0=+X`, `1=-X`, `2=+Z`, `3=-Z`, `4=+Y`, `5=-Y`.

### `ReadoutPlaneConfig`
Six tiny metadata rows per run, one per detector face. Stores enabled/auto-fill flags plus requested `NU/NV`, pitch, offset, gap and margin parameters. Together with `ChannelMap`, this records both the requested layout policy and the actual generated channels.

### `VolumeMap`
One row per registered logical-volume category that actually exists in the run. Core keeps stable ID 0 for current analysis; other existing detector categories use reserved IDs. The registry is generated from the geometry rather than a fixed output table, so future WLS/concentrator volumes can be added without changing the analysis schema.

### `BoundaryDebug`
Debug-only optical interaction log. Ordinary rows are optical boundary crossings
with integer/readable boundary status plus pre/post direction and polarization.
Stage 3D also permits a tagged first-Rayleigh row (`BoundaryStatus=-1001`,
`BoundaryStatusName=FirstRayleighScatter`) when boundary debug is enabled, so
the G4OpRayleigh angular/polarization law can be validated without changing the
ROOT schema. OFF by default; use only for controlled validation/debug runs.

## Flexible readout

Each of the six faces is independent.

Fixed grid:
`/readout/planeGrid <face> <enabled> nU nV pitchU_mm pitchV_mm offsetU_mm offsetV_mm`

Automatic fill:
`/readout/planeAuto <face> <enabled> gapU_mm gapV_mm marginU_mm marginV_mm`

Uniform fixed grid on all faces:
`/readout/uniformGrid nU nV pitchU_mm pitchV_mm`

The same readout layout is used by both response modes:
- `DirectSiPM`: the element is the effective SiPM entrance and PDE is applied at the boundary.
- `ConcentratorInput`: the element is the concentrator entrance aperture; no fake PDE detection is generated and full local phase space is recorded.

## Species-general analysis

The simulation does not create `Muon`, `Electron`, or `Pion` trees. Analysis selects by PDG and role. `analysis/analyze_general.C` supports:
- `primaryEmitter`: immediate `Photon.ParentTrackID` resolves to a primary non-optical track of the selected species;
- `emitter`: immediate `Photon.ParentTrackID` resolves to any non-optical track of the selected species;
- `ancestor`: use propagated `SourceTrackID`, then select by that source track's `RootPrimaryTrackID`. This includes optical re-emission chains such as future WLS.

Examples:
- combined mu-/mu+: `selectedPDG=13, combineCharge=true`;
- pi+ only: `selectedPDG=211, combineCharge=false`;
- pi- only: `selectedPDG=-211, combineCharge=false`.

For dedicated optical-primary validation, `analysis/analyze_general.C` treats `selectedPDG=0` as an all-optical-photon selection because a primary optical photon has no non-optical parent row in `Track`.
