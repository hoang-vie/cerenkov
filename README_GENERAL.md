# General simulation model

## Geometry and physics

The supported topology is a cuboid radiator with optional uniform PMMA shell and Tyvek-like reflector, rectangular coupling/readout elements on six axis-aligned faces, and top/bottom trigger paddles. It is configurable within this family; arbitrary CAD geometry is not provided through macros.

The executable uses `FTFP_BERT`, replaces its electromagnetic component with `G4EmStandardPhysics_option4`, and registers `G4OpticalPhysics`. The default production cut is 0.1 mm. Cherenkov, absorption, Rayleigh and optical-boundary processing are enabled; scintillation, Mie and WLS/WLS2 are disabled. Macro commands can change optical activation states, which are recorded in `RunConfig`.

Enabling a process alone does not supply calibrated emission/yield/material properties for a new scintillator or WLS system.

## Sources

Gun mode configures particle, kinetic energy, position, direction, time and optical-photon polarization. It is useful for controlled samples.

CRY mode maps CRY coordinates into the detector's Y-vertical convention and can create multiple primaries per event. Supported project time modes are `relative`, `absolute` and `zero`; default `relative` subtracts the earliest primary time in each CRY event. Many output times are single-precision nanoseconds: use relative timing for fine time differences, and assess precision before using absolute timestamps.

The supplied CRY examples select muons. Other supported CRY species require corresponding configuration and validation. The magnetic-latitude helper is an approximate centered-dipole calculation with its stated date range.

## Optical inputs

Material profiles associate Geant4 materials with optical tables and provenance. Core supports Water_Default and PMMA_SNO_A through the public interface; OuterShell currently supports PMMA_SNO_A and Coupling PDMS_Generic.

Tables are converted from wavelength to photon energy before attachment to Geant4. Data outside a table's numerical support use the endpoint behaviour described by the model. Numerical table support and empirical measurement validity are distinct concepts.

The bundled historical/candidate tables are retained for comparisons. A filename in `data/` does not make it active; inspect the selected profile and dataset fields in `RunConfig`.

## Readout

DirectSiPM uses an absorbing effective surface with zero configured reflectivity and a wavelength-dependent detection efficiency. This is a photon-counting response model, without full SiPM microcell/electronics behaviour.

ConcentratorInput uses a virtual absorbing entrance with zero detection efficiency. `SensorIncident` is the interface to a downstream optical study. Its directions and polarization are the incident pre-boundary state expressed in the channel's local frame. Local +z points outward; the nominal entrance plane has EntranceZ_mm=0.

## Photon origin and transport

The output distinguishes immediate parent, propagated non-optical source and root primary. For example, a muon can produce a secondary electron that emits a photon. Selecting primary-muon emission differs from selecting the entire muon-initiated cascade.

`analyze_general.C` provides `primaryEmitter`, `emitter` and `ancestor` selections. `selectedPDG=0` selects all optical photons, including optical-primary test samples.

`Photon.ReachedCoreWall` denotes an encounter. Optional `CoreBoundaryCrossing` denotes the first successful transmission out of Core; reflection and total internal reflection do not count. These are logical-volume boundaries. `SensorIncident` and `DetectedSignal` then distinguish readout arrival and detection.

## Analysis quantities

`TrackVolume.Edep_MeV` is local energy deposition; `DeltaKE_MeV` is the accumulated kinetic-energy change of a track. Secondary transport means these need not agree. State the observable and selection when comparing stopping-power results.

An efficiency must specify its population and denominator: generated photons, unique photons reaching a readout, detected photons, generated events or triggered events. In thinAll/triggered, detailed optical output is conditional on the trigger. Use RunStatistics for generated-event totals, especially in triggered mode.

See [Validation](README_VALIDATION.md) for analysis caveats and [Schema](SCHEMA.md) for exact output semantics.
