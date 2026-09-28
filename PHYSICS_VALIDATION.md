# Physics validation plan — general schema

The validation strategy is intentionally separated from CRY production. Small deterministic particle-gun tests should establish that geometry, charged-particle physics, optical generation, transport, and sensor response behave as expected before large cosmic-ray runs are interpreted.

## 0. Validation order

Run cheap non-optical tests first, then optical tests:

1. ROOT/MT output-merging sanity, then generator/geometry/path-length sanity;
2. charged-particle dE/dx and secondary production;
3. muon decay proper time / decay products;
4. Cherenkov threshold, angle and Frank–Tamm yield;
5. absorption and Rayleigh free paths;
6. Fresnel/TIR boundary scans;
7. Tyvek reflection/timing;
8. SiPM PDE or concentrator-input phase space;
9. CRY energy/zenith/multiplicity distributions;
10. underground overburden validation.

`validation_output_mt.mac` should pass before physics output is interpreted. `Track` and `TrackVolume` are sufficient for the non-optical physics stages. Optical tables are not needed there and the validation macros disable optical production to save CPU.

## 1. Charged-particle dE/dx

Use `macros/validation_dedx_muon.mac` or the generic gun with another species. Select tracks by PDG rather than by tree name.

For each selected track and material keep two distinct observables:

- local deposition: `TrackVolume.Edep_MeV / TrackVolume.Length_cm`;
- primary kinetic-energy loss: `TrackVolume.DeltaKE_MeV / TrackVolume.Length_cm`.

For a clean through-going primary stopping-power comparison, prefer `DeltaKE/L`; `Edep/L` answers the different question of how much energy was deposited locally in the volume after secondary transport. Compare against an external stopping-power reference appropriate to the same particle, energy and material. Do not combine electrons, muons, protons and pions into one Bethe–Bloch validation curve.

Recommended plots:
- dE/dx distribution by PDG and VolumeID;
- vertex kinetic energy (or beta*gamma offline) vs dE/dx;
- TrackVolume length distribution;
- secondary PDG and creator-process distributions.

## 2. Muon decay

`macros/validation_muon_decay.mac` uses a stopped `mu+` sample. The generic `Track` table stores end proper time and end process, so no muon-specific output tree is required.

Validation plots:
- proper-time distribution for tracks ending by decay;
- exponential fit for an explicitly selected `mu+` sample;
- electron/positron kinetic-energy spectrum from muon-decay daughters (Michel spectrum).

Do not interpret a stopped `mu-` sample in matter as a free-muon lifetime measurement because nuclear capture competes with decay.

## 3. Cherenkov generation

Use a monoenergetic charged-particle gun and a simple radiator geometry.

### Threshold

Verify that photons are absent below the kinematic condition `beta*n(lambda) > 1` and appear above it.

### Angle

For each Cherenkov photon with a finite `Photon.EmissionAngle_deg`, compare with

`cos(theta_C) = 1 / (beta*n(lambda))`.

The preferred diagnostic is a 2-D `theta_C vs wavelength` plot with the analytical curve overlaid, not only a 1-D angle histogram.

### Frank–Tamm yield

For each selected parent track, compare the generated Cherenkov count with the Frank–Tamm integral over the optical wavelength range actually represented by the material RINDEX table. Useful diagnostics are:
- N_gamma vs TrackVolume length;
- N_gamma_MC / N_gamma_FrankTamm per track;
- generated wavelength spectrum.

`validate_physics.py` provides an independent table-based numerical check without Geant4.

## 4. Optical absorption

Launch monochromatic optical photons through a controlled path length. Verify

`N(L)/N(0) = exp(-L/L_abs)`.

Test Water and PMMA separately. The production PMMA absorption preset is the published SNO UV-transmitting Sample-A reference table. It is a reference material, not proof that an arbitrary laboratory PMMA sheet has the same attenuation.


### PDMS_Generic coupling absorption

`PDMS_Generic` uses
`OpticalGrease_PDMS_ABS_Generic10m.txt` as an explicit generic low-loss
transport assumption. The constant `ABSLENGTH = 10 m` is **not** a vendor
measurement and is not assigned an empirical measured wavelength range.

The profile therefore records:

- numerical/tabulated ABS support: 200--900 nm;
- empirical ABS validity: undeclared, stored as `(-1,-1)` nm;
- status: `Generic10m_Placeholder_NotVendorMeasurement`;
- Geant4 table policy outside numerical support: endpoint clamp.

A direct 1-mm coupling transport test at 400 nm validated the runtime unit
conversion and `G4OpAbsorption` behavior. For a 0.9999-mm path with
`L_abs = 10 m`, 1,000,000 photons gave 101 bulk losses versus 99.985
expected (`z = 0.102`), with zero safety kills.

For the nominal 1-mm coupling thickness, Beer--Lambert single-pass
sensitivity is approximately:

| `L_abs` assumption | loss at 0 deg | loss at 45 deg | loss at 70 deg |
|---:|---:|---:|---:|
| 0.1 m | 0.995% | 1.404% | 2.881% |
| 1 m | 0.100% | 0.141% | 0.292% |
| 10 m | 0.010% | 0.014% | 0.029% |
| 100 m | 0.001% | 0.0014% | 0.0029% |
| transparent | 0 | 0 | 0 |

These numbers are a single-pass coupling-transport sensitivity study, not a
full detector-efficiency systematic. The 10-m value remains a generic
design assumption and must be replaced/revalidated when the actual optical
grease is selected.

## 5. Rayleigh scattering

In a sufficiently large water volume, verify that the first Rayleigh-scattering path follows the exponential distribution implied by the configured `RAYLEIGH` mean free path. `Water_Default` now uses the Stage-3E Zhang-Hu PMH pure-water table at 21.5 C; the historical Morel-seawater-like table remains only for regression/comparison. The production `Photon.NRayleigh` counter is an event-history summary; Stage-3D debug instrumentation can tag the first scatter for angular/polarization validation. The remaining `delta=0.039` versus Geant4 ideal-dipole angular-law difference is documented as a small systematic rather than claimed exact.

## 6. Fresnel / TIR interfaces

Use controlled monochromatic optical photons and scan incidence angle through the relevant interfaces. Enable `/output/boundaryDebug true` only for this small sample.

At minimum test the geometries that can actually exist in the detector configuration:
- Core -> Acrylic;
- Acrylic -> Coupling;
- Core -> Coupling (no acrylic shell);
- Coupling -> Core/Acrylic;
- Coupling -> Air side walls;
- Coupling -> Readout surface.

Compare `BoundaryDebug.BoundaryStatus` transitions with Snell/Fresnel expectations calculated from the exact RINDEX files used in that run. Do not encode an interface interpretation in the production schema when `(PreVolumeID, PostVolumeID, BoundaryStatus)` can express the raw fact.


### PMMA <-> PDMS coupling interface

The production coupling profile uses
`OpticalGrease_PDMS_RINDEX_Roychowdhury2018_Cauchy3Anchor.txt`, a
three-anchor Cauchy approximation rather than a raw measured spectral table.

At 400 nm the production tables give approximately:

- `n_PMMA = 1.502601863`;
- `n_PDMS = 1.417962440`.

The actual detector `Acrylic -> Coupling` boundary was validated with
100,000 photons per polarization/angle against analytical Fresnel and Snell
predictions. The scan included the Brewster region and the PMMA->PDMS
critical-angle transition. At 70.67 deg the interface remained Fresnel,
while at 70.68 deg all 100,000 photons underwent total internal reflection.

The reverse PDMS->PMMA direction was also validated against Fresnel/Snell
predictions, including the approximately 46.66-deg Brewster minimum. No
unexpected boundary status or total internal reflection occurred.

The explicit named `RadiatorGelSurface` currently exists for the forward
radiator->coupling direction. The reverse direction falls back to the bulk
dielectric boundary and has been verified to give the same expected optical
physics for the present polished dielectric-dielectric model. If a future
coupling interface adds roughness, coatings, or surface-specific optical
properties, directional border-surface symmetry must be re-audited.

## 7. Tyvek

Baseline model: `dielectric_metal`, UNIFIED, `ground`, with wavelength-dependent measured reflectivity. Validate total reflected fraction versus wavelength first. A precise BRDF/angular validation requires a measured/literature-backed surface model for the actual Tyvek grade and environment.

Useful detector-level QA:
- Tyvek-reflection count per photon;
- photon arrival time vs Tyvek-reflection count;
- late-light timing tail.

## 8. Readout modes

### DirectSiPM

`SensorIncident` records photons arriving at the readout entrance before the effective PDE decision. `DetectedSignal` records Geant4 boundary detections.

At fixed wavelength and normal incidence verify

`N_detected / N_incident = configured PDE`

within binomial uncertainty. If `SiPM_PDE.txt` is a manufacturer PDE curve, do not multiply a separate fill factor into the simulation again.

### ConcentratorInput

No fake SiPM detection is generated. Use `SensorIncident` as the optical-concentrator input. In `phaseSpace` detail it stores local position, local direction, polarization, wavelength and time. Validate coordinate transforms by injecting rays with known positions/directions on each face.

## 9. Trigger

The hardware-like trigger uses separate quantities:
- scintillator pulse integration window;
- deposited-energy threshold;
- top-bottom coincidence window.

It does not require the same Geant4 TrackID in both paddles.

`TriggerEvent` is intentionally small. In `thinAll`, it is stored for all generated events. `ScintDeposit` stores sparse raw non-zero per-step deposits `(paddle, track, time, Edep)`. This allows the scintillator integration window, energy threshold and coincidence window to be reconstructed or scanned offline without rerunning Geant4.

Recommended plots:
- EdepTop vs EdepBottom;
- t_bottom - t_top;
- trigger fraction vs threshold;
- trigger fraction vs coincidence window;
- trigger truth composition using `Track` only as an offline diagnostic, never as a hardware condition.

## 10. CRY production

Before using detector results, check the generated primary distributions:
- PDG/species;
- kinetic energy;
- zenith angle;
- x-z source-plane distribution;
- multiplicity and relative timing.

CRY's `latitude` input is magnetic latitude, not geographic latitude. For the Quy Nhon example, the bundled `tools/cry_magnetic_latitude.py` uses the IGRF-14 degree-1 tilted-dipole approximation and gives about 4.63 deg for 2026-06-18. Recompute this value for the actual run date/location rather than inserting the geographic 13.78 deg directly.

CRY supplies atmospheric secondaries at its supported source elevation; an underground laboratory additionally requires propagation through a validated overburden model or a validated underground spectrum. Surface-generator correctness and overburden correctness are separate validation layers.

## 11. Output-size validation

Physics correctness must also be checked against output policy:
- `thinAll`: all TriggerEvent/ScintDeposit rows, detailed rows only for triggered events;
- `triggered`: only triggered-event data plus RunStatistics denominator;
- `fullAll`: full details for small validation samples.

For a representative production run, record bytes/event and rows/event for `Photon`, `SensorIncident`, and `DetectedSignal`. Full `BoundaryDebug` must remain OFF in production unless a measured size budget justifies it.

## 12. Reproducibility

`RunConfig` records the base RNG seed, generator/CRY configuration, material dataset names, trigger/output settings, actual channel map and per-face readout configuration. A fixed-seed sequential validation run should be exactly repeatable. Multi-threaded event ordering should not be used as a reproducibility identifier; use RunID/EventID/TrackID joins.
