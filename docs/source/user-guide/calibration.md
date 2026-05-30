# Calibration

The SAXS geometry depends on a handful of constants that must be
calibrated against a known sample.  These constants live as
module-level defaults in `smi_tiled.loader`:

| Constant | Default | Role |
|---|---|---|
| `_SAXS_DEFAULT_DISTANCE_DELTA_MM` | from calibration JSON | Additive correction: `SDD = motor_z + delta` |
| `_SAXS_MOTOR_X_REF_MM`, `_SAXS_MOTOR_Y_REF_MM` | `1.88`, `2.45` | Reference motor positions (unused while the x/y slopes are 0.0; kept for the override API) |
| `_SAXS_MOTOR_Z_REF_MM` | `0.0` | Reference motor_z for the small motor_z→BC drift |
| `_SAXS_PIEZO_Z_REF_UM` | `0.0` | Reference piezo_z for the piezo_z→SDD effect |
| `_SAXS_BEAM_COL_PX_PER_MOTOR_X_MM` | `0.0` | px / mm slope, motor_x → beam column. **0.0**: the EPICS PV already tracks motor_x (non-zero double-counts) |
| `_SAXS_BEAM_ROW_PX_PER_MOTOR_Y_MM` | `0.0` | px / mm slope, motor_y → beam row. **0.0**: the EPICS PV already tracks motor_y (non-zero double-counts) |
| `_SAXS_BEAM_COL_PX_PER_MOTOR_Z_MM` | from calibration JSON | px / mm slope, motor_z → beam column |
| `_SAXS_BEAM_ROW_PX_PER_MOTOR_Z_MM` | from calibration JSON | px / mm slope, motor_z → beam row |
| `_SAXS_SDD_DELTA_MM_PER_PIEZO_Z_UM` | from calibration JSON | mm / μm slope, piezo_z → SDD |

These are applied during {func}`~smi_tiled.resolve_saxs_geometry`.

## JSON override

A calibration JSON at `src/smi_tiled/data/saxs_calibration.json` is read
at **import time** and overrides any of the constants in a `"constants"`
block.  The bundled file currently represents the b900e711-… z-grid
scan calibration.

To inspect what's active:

```python
from smi_tiled.loader import (
    _SAXS_CALIBRATION_OVERRIDE,
    _SAXS_DEFAULT_DISTANCE_DELTA_MM,
    _SAXS_BEAM_COL_PX_PER_MOTOR_Z_MM,
)
print(_SAXS_CALIBRATION_OVERRIDE["source_uid"])
print(_SAXS_DEFAULT_DISTANCE_DELTA_MM)        # active value
```

## Re-fitting from an AGB grid scan

Two scripts ship under `scripts/`:

### `calibrate_smi_saxs.py`

Fits the **detector-translation** dependence (motor_x, motor_y, piezo_z)
from a grid scan that varies those motors at fixed motor_z.  Reference
scan: `b0f165c4-…`.

```bash
pixi run python scripts/calibrate_smi_saxs.py b0f165c4-…
```

### `calibrate_smi_z_scan.py`

Fits the **detector-distance** dependence (motor_z) plus the small
motor_z→BC drift, AND extracts the **beamstop centering offset table**
needed by the SMI collection code.  Reference scan: `b900e711-…`.

```bash
pixi run python scripts/calibrate_smi_z_scan.py b900e711-…
```

Both produce a `/tmp/saxs_calibration.json` file with a `"constants"`
block.  Copy it into the package data dir to make the override active:

```bash
cp /tmp/saxs_calibration.json src/smi_tiled/data/saxs_calibration.json
```

And restart Python (the override is applied at import time, so it only
takes effect on a fresh interpreter).

## How the fit works

`calibrate_smi_z_scan.py`:

1. For each frame of an AGB grid scan, locate the pin transmission
   bright spot (initial seed).
2. Find the **pin beamstop shadow centroid** — a stable seed across
   motor_z because the shadow disk is large and consistent.
3. Compute the AGB ring radius via a chi-averaged radial profile.
4. Sample the ring at multiple chi sectors and least-squares fit a
   circle → refined beam center.
5. Convert ring radius to SDD via Bragg: `tan(2θ_AGB) = r_mm / SDD_mm`.

After looping all frames it regresses:

```
SDD(motor_z, piezo_z) = a_s + b_z · motor_z + b_p · piezo_z
BC_col(motor_z)       = a_c + b_zc · motor_z
BC_row(motor_z)       = a_r + b_zr · motor_z
```

The coefficients end up in the `"constants"` block of the output JSON.

## What about motor_x / motor_y / piezo_z?

The b900e711-… scan holds motor_x and motor_y essentially fixed, so
their slopes (`_SAXS_BEAM_*_PX_PER_MOTOR_*_MM`) come from
`calibrate_smi_saxs.py` against a scan like b0f165c4-… that varies
those motors.  Future re-calibrations should ideally use **one big
multi-dimensional grid scan** that varies all five motors so a single
regression captures the full picture.

## Sanity checking your calibration

After loading the new JSON, run a known scan and compare to the
expected ring radius:

```python
from smi_tiled import reduce_smi_combined

result = reduce_smi_combined(uid="b900e711-…", image_cache_path=None)
print(result.merged_iq.q.min(), result.merged_iq.q.max())
# At 16.1 keV, ring 1 of AGB is at q ≈ 1.076 nm⁻¹.
# Look for a peak at that q in result.merged_iq["I"].
```

For a multi-distance check, run a few scans across the motor_z range
and verify the AGB ring 1 q-position is consistent.

## See also

- {doc}`beamstop-centering` — the spline-based offset table also
  produced by `calibrate_smi_z_scan.py`
- {doc}`../reference/calibration-format` — full schema of the JSON
  override file
- {func}`~smi_tiled.resolve_saxs_geometry`
