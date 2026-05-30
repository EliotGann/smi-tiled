# Calibration JSON format

`src/smi_tiled/data/saxs_calibration.json` is read by
`smi_tiled.loader._apply_calibration_override` at import time and
overrides module-level constants.

## Schema

```json
{
  "_doc": "Free-text description.",
  "source_uid": "b900e711-…",
  "energy_ev": 16099.99,
  "wavelength_nm": 0.07700884277,
  "active_beamstop": "pin",
  "scan_fixed_motor_x_mm": 60.009,
  "scan_fixed_motor_y_mm": 10.000,
  "regressions": {
    "bc_col": {
      "intercept": 1099.112,
      "per_motor_z_mm": 0.0002654,
      "per_piezo_z_um": -4.14e-06,
      "rms_px": 3.454,
      "_note": "Fit AT motor_x=60.009, motor_y=10.000. ..."
    },
    "bc_row": {
      "intercept": 1170.957,
      "per_motor_z_mm": -0.0003495,
      "per_piezo_z_um": -5.25e-06,
      "rms_px": 2.833
    },
    "sdd_mm": {
      "intercept": -13.676,
      "per_motor_z_mm": 0.990479,
      "per_piezo_z_um": 0.0009035,
      "rms_mm": 28.491
    }
  },
  "constants": {
    "_SAXS_DEFAULT_DISTANCE_DELTA_MM": -13.676,
    "_SAXS_BEAM_COL_PX_PER_MOTOR_Z_MM": 0.000265,
    "_SAXS_BEAM_ROW_PX_PER_MOTOR_Z_MM": -0.000350,
    "_SAXS_SDD_DELTA_MM_PER_PIEZO_Z_UM": 0.000904
  }
}
```

## Sections

### `regressions` (informational)

Human-readable record of the fits.  Not consumed by the loader; useful
for spot-checking the calibration quality.

| Field | Meaning |
|---|---|
| `intercept` | constant term of the linear fit |
| `per_motor_z_mm` | slope vs `pil2M_motor_z` |
| `per_piezo_z_um` | slope vs `piezo_z` |
| `rms_px` / `rms_mm` | residual RMS of the fit |
| `_note` | free-text caveats (e.g. fixed motor values) |

### `constants` (consumed by the loader)

A dict whose keys must match `_SAXS_*` names in `smi_tiled.loader`.
Recognized keys:

| Key | Default | Type |
|---|---|---|
| `_SAXS_DEFAULT_DISTANCE_DELTA_MM` | `0.0` | `float`, mm |
| `_SAXS_MOTOR_X_REF_MM` | `1.88` | `float`, mm |
| `_SAXS_MOTOR_Y_REF_MM` | `2.45` | `float`, mm |
| `_SAXS_MOTOR_Z_REF_MM` | `0.0` | `float`, mm |
| `_SAXS_PIEZO_Z_REF_UM` | `0.0` | `float`, μm |
| `_SAXS_BEAM_COL_PX_PER_MOTOR_X_MM` | `0.0` | `float`, px/mm (PV already tracks motor_x) |
| `_SAXS_BEAM_ROW_PX_PER_MOTOR_Y_MM` | `0.0` | `float`, px/mm (PV already tracks motor_y) |
| `_SAXS_BEAM_COL_PX_PER_MOTOR_Z_MM` | `0.0` | `float`, px/mm |
| `_SAXS_BEAM_ROW_PX_PER_MOTOR_Z_MM` | `0.0` | `float`, px/mm |
| `_SAXS_SDD_DELTA_MM_PER_PIEZO_Z_UM` | `0.000988` | `float`, mm/μm |
| `_SAXS_DEFAULT_BEAM_DELTA_ROW_PX` | `0.0` | `float`, px |
| `_SAXS_DEFAULT_BEAM_DELTA_COL_PX` | `0.0` | `float`, px |

Unknown keys emit a warning but don't fail the import.

## Application order in resolve_saxs_geometry

```
beam_row = baseline_y_px  (taking abs)
beam_col = baseline_x_px  (taking abs)
dist_mm  = motor_z_mm
```

then the loader applies (in order):

```
beam_col += (motor_x_mm - _SAXS_MOTOR_X_REF_MM) * _SAXS_BEAM_COL_PX_PER_MOTOR_X_MM
beam_row += (motor_y_mm - _SAXS_MOTOR_Y_REF_MM) * _SAXS_BEAM_ROW_PX_PER_MOTOR_Y_MM
beam_col += (motor_z_mm - _SAXS_MOTOR_Z_REF_MM) * _SAXS_BEAM_COL_PX_PER_MOTOR_Z_MM
beam_row += (motor_z_mm - _SAXS_MOTOR_Z_REF_MM) * _SAXS_BEAM_ROW_PX_PER_MOTOR_Z_MM
dist_mm  += (piezo_z_um - _SAXS_PIEZO_Z_REF_UM) * _SAXS_SDD_DELTA_MM_PER_PIEZO_Z_UM
dist_mm  += _SAXS_DEFAULT_DISTANCE_DELTA_MM
beam_row += _SAXS_DEFAULT_BEAM_DELTA_ROW_PX
beam_col += _SAXS_DEFAULT_BEAM_DELTA_COL_PX
```

```{note}
The `motor_x`/`motor_y` beam-center slopes default to `0.0`, so the two
motor_x/y lines above are no-ops by default: the resolved beam center equals the
live EPICS PV, which already tracks the detector translation. They are kept 0.0
to avoid double-counting (a non-zero slope pushed the center ~120 px off on AGB
scans). Do not restore the historical `5.82`/`5.99` values unless a future
detector reports a *static* beam-center PV.
```

## When to regenerate

- After a physical change to the SAXS detector or beamstop mount.
- After re-aligning the beam.
- After a long shutdown or major maintenance.
- When a known sample's q-position drifts from its expected value by
  more than the calibration's stated RMS.

Use `scripts/calibrate_smi_z_scan.py` (see {doc}`../user-guide/calibration`).
