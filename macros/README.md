# Macro catalog

All macros configure the same simulation executable. Run them from the build directory after CMake exposes `macros/*.mac` there as symlinks.

## Ordinary examples

| Macro | Source and purpose |
|---|---|
| `example_gun_direct_sipm.mac` | 100 downward 3 GeV muons, 24 effective readout channels, fullAll output; a deterministic starting point. |
| `example_direct_sipm.mac` | CRY muons with effective DirectSiPM response and thinAll output. |
| `example_concentrator_input.mac` | CRY muons with entrance phase-space output and thinAll output. |

Copy an example for a new study. These configurations illustrate usage; they are not optimized or calibrated detector designs. Readout dimensions differ between examples: geometry must be selected for the study rather than inferred from the sensor-model label.

```bash
./cherenkov example_gun_direct_sipm.mac 12345
```

## Controlled validation

`validation_*.mac` files isolate processes or exercise output bookkeeping. They may disable optical processes, use artificial geometry or inject monochromatic photons. Consult the matching sample/validator table in [Validation](../README_VALIDATION.md).

## Visualization

`init_vis.mac` supplies the interactive configuration when launching without a macro argument. `vis.mac` sets visualization options. Their settings are not production defaults.

## Configuration order

Set detector/readout geometry before `/run/initialize`; configure Gun/CRY afterwards, then call `/run/beamOn N`. CRY input lines need `/CRY/update`. Setting `/generator/position` disables auto-Y; re-enable it afterwards if intended.

For multithreaded geometry or material scans, use a new process per configuration. Give each run a distinct output prefix and retain its macro and seed.

## Output choices

thinAll retains detailed optical records only for triggered events. Use fullAll when studying optical behaviour without trigger selection. incidentDetail auto selects minimal in DirectSiPM and phaseSpace in ConcentratorInput.

CoreBoundary recording is optional and off by default. firstExit records successful Core transmission, not reflection; phaseSpace additionally saves polarization. It follows the detailed-event output policy.
