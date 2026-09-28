Cherenkov Simulation Framework v1
================================

Version
-------
Framework tag   : framework_v1
Schema          : general-v1.10
Validated commit: e3c429f8a682d2fe6307a8c0508994d983f8f69e

Overview
--------
This project is a reusable Geant4 framework for simulating a cuboid Cherenkov detector.

The framework is designed so that most detector studies can be performed by editing
Geant4 macro files rather than modifying the C++ source code.

Typical simulation chain:

    Gun / CRY cosmic-ray generator
                |
                v
        Charged-particle transport
                |
                v
          Cherenkov production
                |
                v
            Water Core
                |
                v
       Optical photon transport
                |
                v
       Acrylic / Tyvek / Coupling
                |
                v
      DirectSiPM or ConcentratorInput
                |
                v
             ROOT output


1. Detector Geometry
====================

The supported detector family is based on a rectangular/cuboid radiator.

Conceptual geometry:

                    +Y

                [ Scintillator ]

                +-------------+
                |   Readout   |
                +-------------+
                     |   |
            -----------------------
           |                       |
      -X   |      WATER CORE       |   +X
           |                       |
            -----------------------
                     |   |
                +-------------+
                |   Readout   |
                +-------------+

                    -Y

The +Z and -Z faces can also contain readout elements.

Main detector components:

- Core
  Main Cherenkov radiator.

- Outer shell
  Optional acrylic/PMMA shell surrounding the Core.

- Tyvek
  Optional reflective layer.

- Coupling
  Optical coupling layer between detector material and readout.

- Readout
  SiPM or concentrator-input elements placed on the six axis-aligned faces.

- Scintillators
  Top and bottom trigger paddles used for energy-deposition / coincidence triggering.

Default dimensions in the current framework are approximately:

- Core            : 100 x 100 x 100 mm
- Scintillator    : 110 x 10 x 110 mm
- Readout element : 25 x 25 x 1 mm
- Coupling        : 25 x 25 x 1 mm
- Acrylic shell   : 3 mm
- Tyvek           : 0.1 mm

These are defaults, not fixed design constraints.
Most geometry parameters can be changed from macro commands.


2. Readout Configuration
========================

The framework supports all six faces:

    +X
    -X
    +Z
    -Z
    +Y
    -Y

Each face can be configured independently.

Supported layout controls include:

- enabled / disabled
- number of elements in local U direction
- number of elements in local V direction
- pitch U
- pitch V
- offset U
- offset V

Main commands:

    /readout/uniformGrid
    /readout/planeGrid
    /readout/planeAuto

Example:

    +X : 3 x 3
    -X : 2 x 2
    +Z : automatic fill
    -Z : disabled
    +Y : 1 x 2
    -Y : 4 x 4

This can be done without editing C++.

The requested configuration is stored in:

    ReadoutPlaneConfig

The actual generated detector channels are stored in:

    ChannelMap


3. Physics Configuration
========================

Current physics stack:

    FTFP_BERT
    + G4EmStandardPhysics_option4
    + G4OpticalPhysics

Default production cut:

    0.1 mm

Important optical-process defaults:

- Cherenkov      : ON
- Scintillation  : OFF
- Absorption     : ON
- Rayleigh       : ON
- Mie scattering : OFF
- Boundary       : ON
- WLS            : OFF
- WLS2           : OFF

Important charged-particle physics includes:

- ionization
- energy loss
- multiple scattering
- bremsstrahlung for electrons
- hadronic physics from FTFP_BERT
- muon decay
- secondary-particle production

The framework has validation tests for:

- muon dE/dx
- multiple scattering
- electron ionization/bremsstrahlung
- muon lifetime
- Michel-electron endpoint
- particle ancestry
- Cherenkov angle
- Frank-Tamm photon yield
- optical-boundary behavior
- SiPM PDE response


4. Cherenkov Physics
====================

For beta approximately equal to 1:

Water:
    Frank-Tamm yield over 200-900 nm
    approximately 818 photons/cm

Acrylic:
    approximately 1018 photons/cm

For a 3 GeV kinetic-energy muon:
    beta = 0.99942111

Water Cherenkov yield:
    approximately 817 photons/cm

Typical Cherenkov angle in Water:
    300 nm : 42.63 deg
    400 nm : 41.90 deg
    500 nm : 41.57 deg
    600 nm : 41.39 deg


5. Optical Materials
====================

5.1 Water
---------
Canonical profile:
    Water_Default

Includes wavelength-dependent:
- refractive index
- absorption length
- Rayleigh scattering


5.2 Acrylic / PMMA
------------------
Canonical profile:
    PMMA_SNO_A

Important distinction:
    300-500 nm = measured region
    200-300 nm and 500-900 nm = explicit transport-policy extensions

Representative absorption lengths:
    300 nm : approximately 24.4 mm
    400 nm : approximately 1.11 m
    500 nm : approximately 1.43 m


5.3 Optical Coupling
--------------------
Canonical profile:
    PDMS_Generic

Current absorption-length assumption:
    10 m

This is a generic placeholder / engineering assumption rather than a complete
material-specific measured absorption curve.


5.4 Tyvek
---------
The framework includes wavelength-dependent reflectivity.

Representative value:
    R(450 nm) approximately 0.97

The current surface model is a practical generic optical model.
It is not a complete experimentally measured Tyvek angular BRDF model.


6. Sensor Modes
===============

6.1 DirectSiPM
--------------
Simulation chain:

    Optical photon
         |
         v
      Coupling
         |
         v
    Sensor surface
         |
         v
    wavelength-dependent PDE
         |
         v
    DetectedSignal

The current PDE profile is based on a Hamamatsu S1336-3050CS-like response.

Representative value:
    PDE(450 nm) approximately 0.40


6.2 ConcentratorInput
---------------------
This mode does NOT perform full concentrator ray tracing inside Geant4.

Instead:

    Optical photon
         |
         v
    virtual concentrator entrance
         |
         v
    phase-space output
         |
         v
    external optical study / ROBAST


7. Event Generators
===================

7.1 Gun
-------
Useful for deterministic studies and validation.

Configurable:
- particle type
- kinetic energy
- position
- direction
- polarization
- time
- automatic source-Y placement

Recommended first example:
    example_gun_direct_sipm.mac


7.2 CRY
-------
Used for realistic cosmic-ray generation.

Typical configurable quantities:
- particle species
- date
- latitude
- altitude
- subbox length
- time mode

The framework supports events with multiple CRY primaries.


8. Trigger System
=================

Top and bottom scintillator paddles are available for trigger studies.

Event-level information includes:
- EdepTop
- EdepBottom
- timing
- TriggerPassed
- NPrimaryNonOptical
- NNonOpticalTracks
- NOpticalTracks
- NIncident
- NDetected
- EdepCoreTotal


9. CoreBoundaryCrossing Diagnostic
==================================

CoreBoundaryCrossing is an optional diagnostic recorder.

It is NOT a physical detector component.

Purpose:
    record the first successful optical-photon transmission out of the Core

Reflection or total internal reflection:
    not recorded as firstExit

Successful transmission:
    recorded once per photon

Face convention:
    0 = +X
    1 = -X
    2 = +Z
    3 = -Z
    4 = +Y
    5 = -Y

Commands:
    /output/coreBoundaryMode off
    /output/coreBoundaryMode firstExit
    /output/coreBoundaryDetail basic
    /output/coreBoundaryDetail phaseSpace


10. Dependencies
================

Main software dependencies:
- Geant4
- ROOT
- CMake
- C++ compiler
- Python 3

The project can also use:
- CRY
- ROBAST


11. Build Instructions
======================

    cd cerenkov
    mkdir build
    cd build

    cmake ..
    cmake --build . -j8

Main executable:
    ./cherenkov


12. First Simulation
====================

Run:

    ./cherenkov example_gun_direct_sipm.mac 12345

Expected output:
    example_gun_direct_sipm_run0.root


13. Check the First Output
==========================

    root -l -b -q \
    '../analysis/inspect_output.C("example_gun_direct_sipm_run0.root",100,24)'

Healthy run:
    OUTPUT SANITY: PASS


14. Main User Examples
======================

example_gun_direct_sipm.mac
    Deterministic Gun + DirectSiPM

example_direct_sipm.mac
    CRY + DirectSiPM

example_concentrator_input.mac
    CRY + ConcentratorInput


15. Common Configuration Changes
================================

Core size:
    /detector/coreSize ...

Core profile:
    /detector/coreProfile Water_Default

Enable/disable outer shell:
    /detector/useOuterShell true

Outer-shell thickness:
    /detector/outerShellThickness ...

Readout element size:
    /detector/readoutSize ...

Uniform grid:
    /readout/uniformGrid ...

Per-face grid:
    /readout/planeGrid ...

Auto-fill:
    /readout/planeAuto ...

Sensor mode:
    /detector/sensorMode DirectSiPM
or:
    /detector/sensorMode ConcentratorInput

Generator:
    /generator/mode Gun
or:
    /generator/mode CRY

Core-boundary diagnostic:
    /output/coreBoundaryMode off
or:
    /output/coreBoundaryMode firstExit


16. Recommended Workflow for a New Study
========================================

1. Copy an existing public example macro.
2. Change the detector configuration.
3. Change the output file name.
4. Select a random seed.
5. Run the simulation.
6. Check RunConfig.
7. Run inspect_output.C.
8. Analyze physics quantities.
9. Compare configurations.

Only modify C++ when the study requires a capability not available through macros.


17. ROOT Output Structure
=========================

Current schema:
    general-v1.10

Trees:

TriggerEvent
    Event-level trigger and summary information.

ScintDeposit
    Energy deposits in scintillator volumes.

Track
    Non-optical track information.

TrackVolume
    Track/volume traversal information.

Photon
    Optical-photon history and lifetime summary.

SensorIncident
    Optical photons reaching readout/sensor surfaces.

DetectedSignal
    Photons accepted by the sensor PDE model.

ChannelEvent
    Per-channel / per-event summary.

RunConfig
    Simulation configuration and provenance.

ChannelMap
    Actual generated channel geometry and mapping.

VolumeMap
    Volume-ID mapping.

RunStatistics
    Per-worker/run statistics.

BoundaryDebug
    Optional detailed optical-boundary debugging.

ReadoutPlaneConfig
    Requested six-face readout configuration.

PrimaryTruth
    Generator-level primary-particle truth.

CoreBoundaryCrossing
    First successful optical-photon exit from the Core.


18. Data Flow
=============

    PrimaryTruth
         |
         v
       Track
         |
         v
       Photon
         |
         v
    CoreBoundaryCrossing
         |
         v
    SensorIncident
         |
         v
    DetectedSignal


19. Reproducibility
===================

RunConfig stores:
- schema version
- project revision
- Git commit
- Git branch
- dirty/clean source state
- detector geometry
- material profiles
- generator configuration
- random seed
- number of threads
- sensor mode
- readout configuration
- optical-process configuration
- output policy
- CoreBoundary settings

Validated framework release:
    e3c429f8a682d2fe6307a8c0508994d983f8f69e

Freeze state:
    GitDirty = 0


20. Validation
==============

Static validation includes:
- ntuple/schema consistency
- analysis/schema consistency
- macro-command registration
- configuration-to-metadata contract
- macro catalog
- Cherenkov physics sanity
- material optical sanity
- CRY magnetic-latitude check

Runtime validation includes:
- flexible asymmetric readout
- DirectSiPM optical chain
- CoreBoundary off/basic/phaseSpace
- TIR handling
- deterministic six-face first-exit tests
- multithreaded integrated regression
- metadata/provenance validation


21. Example Final Freeze Run
============================

Validated integrated run:
- 1000 events
- 4 threads
- 24 DirectSiPM channels
- six readout faces
- fullAll output
- CoreBoundary OFF

Representative output:
    TriggerEvent        1000
    ScintDeposit       31647
    Track              20665
    TrackVolume        27443
    Photon           9992438
    SensorIncident     214684
    DetectedSignal      58158
    ChannelEvent        23991
    RunConfig               1
    ChannelMap              24
    VolumeMap                8
    RunStatistics            4
    ReadoutPlaneConfig       6
    PrimaryTruth          1000
    CoreBoundaryCrossing     0

Validation:
    OUTPUT SANITY: PASS
    FRAMEWORK METADATA VALIDATION: PASS
    FRAMEWORK FREEZE REGRESSION: PASS


22. Important Limitations
=========================

Current geometry scope:
- cuboid Core
- optional uniform outer shell
- optional Tyvek
- rectangular coupling/readout elements
- six axis-aligned readout faces
- top/bottom scintillator paddles

Current limitations:
- no arbitrary CAD geometry through macros
- PDMS absorption uses a generic assumption
- PMMA outside measured range uses explicit model extensions
- Tyvek uses a generic optical-surface model rather than full measured BRDF
- SiPM PDE uses the currently curated sensor profile
- ConcentratorInput exports phase space rather than doing full Geant4 concentrator ray tracing


23. When Should C++ Be Modified?
===============================

Do NOT modify C++ just to change ordinary detector configuration.

Use macros first.

Modify framework code only when a study requires a general capability the framework does not provide.

New framework features should follow:

    implementation
        |
        v
    configuration interface
        |
        v
    metadata
        |
        v
    validation
        |
        v
    documentation


24. Project Philosophy
======================

The goal is NOT to define one optimal detector.

The goal is to provide a physically validated and reproducible simulation framework
so users can investigate different detector configurations.

Typical studies:
- sensor count
- sensor layout
- active faces
- detector dimensions
- source conditions
- DirectSiPM vs ConcentratorInput
- cosmic-ray angular response
- timing
- wavelength response
- optical collection
- material sensitivity


25. More Documentation
======================

README_GENERAL.md
    General framework and detector information.

README_PROJECT_GUIDE.md
    Project structure and usage guide.

README_TREE_INFO.md
    ROOT tree definitions.

README_VALIDATION.md
    Validation overview.

SCHEMA.md
    Formal schema documentation.

FRAMEWORK_STATUS.md
    Framework status and supported scope.

FRAMEWORK_CONFIG_AUDIT.md
    Configuration and metadata audit.

VALIDATION_RUNBOOK.md
    Detailed validation commands.

macros/README.md
    Macro catalog and user examples.


26. Recommended First Steps
===========================

1. Build:

       mkdir build
       cd build
       cmake ..
       cmake --build . -j8

2. Run:

       ./cherenkov example_gun_direct_sipm.mac 12345

3. Inspect:

       root -l -b -q \
       '../analysis/inspect_output.C("example_gun_direct_sipm_run0.root",100,24)'

4. Confirm:

       OUTPUT SANITY: PASS

5. Copy the example macro.
6. Change one detector parameter at a time.
7. Compare ROOT outputs.
8. Always inspect RunConfig before interpreting results.

End of beginner guide.
