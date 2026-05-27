# Beamstop offsets JSON format

`src/smi_tiled/data/smi_beamstop_offsets.json` is intended for the SMI
collection code to consume — it carries the offset between the
calibrated beam center and the pin beamstop shadow vs detector z.

## Schema

```json
{
  "_doc": "Free-text description.",
  "source_uid": "b900e711-…",
  "active_beamstop": "pin",
  "current_motor_positions": {
    "saxs_beamstop_x_pin": -227.000,
    "saxs_beamstop_y_pin": 6.800
  },
  "offsets_by_motor_z": [
    {
      "motor_z_mm": 1700.0,
      "d_bsx_mm": 0.220,
      "d_bsy_mm": -1.416,
      "d_bsx_px": 1.28,
      "d_bsy_px": -8.23,
      "n_frames": 11
    },
    ...
  ],
  "spline": {
    "available": true,
    "model": "scipy.interpolate.BSpline (cubic, smoothing)",
    "degree": 3,
    "sigma_dx_mm": 0.126,
    "sigma_dy_mm": 0.149,
    "knots_d_bsx_mm": [t0, t1, ...],
    "coefs_d_bsx_mm": [c0, c1, ...],
    "knots_d_bsy_mm": [t0, t1, ...],
    "coefs_d_bsy_mm": [c0, c1, ...],
    "rms_residual_dx_mm": 0.126,
    "rms_residual_dy_mm": 0.149,
    "_consumer_note": "...",
    "dense_grid": {
      "motor_z_mm": [...],   // 200 points
      "d_bsx_mm":   [...],
      "d_bsy_mm":   [...]
    }
  }
}
```

## Sections

### `current_motor_positions`

The beamstop motor positions during the source calibration scan.  The
offsets are relative to these; the new (centered) positions are
`current + d_bs{x,y}_mm`.

### `offsets_by_motor_z`

Per-`motor_z` averaged measurements.  At each `motor_z` value, the
script averages over the inner `piezo_z` loop frames.  `n_frames` is
the number of frames that contributed to that average.

| Field | Units | Meaning |
|---|---|---|
| `motor_z_mm` | mm | detector z position |
| `d_bsx_mm` | mm | offset in pixel-column direction, in mm |
| `d_bsy_mm` | mm | offset in pixel-row direction, in mm |
| `d_bsx_px` | px | same offset, in pixels |
| `d_bsy_px` | px | same offset, in pixels |

### `spline`

Smoothing cubic B-spline fits to the per-z measurements.  Tracks
non-linear meanders that a linear fit would miss.

| Field | Meaning |
|---|---|
| `degree` | spline degree (always 3 for these fits) |
| `sigma_dx_mm`, `sigma_dy_mm` | per-z scatter estimate used as smoothing factor |
| `knots_*` | full knot vector with boundary multiplicity (degree+1) |
| `coefs_*` | spline coefficients |
| `rms_residual_*` | RMS of raw points vs spline (typically 0.1-0.2 mm) |
| `dense_grid` | 200-point evaluation across the calibrated range |

## Reconstructing the spline

```python
import json
import numpy as np
from scipy.interpolate import BSpline

with open("smi_beamstop_offsets.json") as f:
    bs = json.load(f)

# d_bsx_mm at any motor_z in the calibrated range:
spl_dx = BSpline(
    np.array(bs["spline"]["knots_d_bsx_mm"]),
    np.array(bs["spline"]["coefs_d_bsx_mm"]),
    bs["spline"]["degree"],
)
print(float(spl_dx(5500.0)))   # → e.g. 0.222 mm
```

## Reconstructing without scipy

If the consumer doesn't have scipy, fall back to linear interpolation
through the `dense_grid`:

```python
import numpy as np

z  = np.array(bs["spline"]["dense_grid"]["motor_z_mm"])
dx = np.array(bs["spline"]["dense_grid"]["d_bsx_mm"])
dy = np.array(bs["spline"]["dense_grid"]["d_bsy_mm"])

d_bsx = float(np.interp(5500.0, z, dx))
d_bsy = float(np.interp(5500.0, z, dy))
```

The dense grid has 200 points across the calibrated range; the
interpolation error is well below the spline's residual scatter.

## Sign convention

```
d_bsx > 0:   move the beamstop in the +col direction on the detector
d_bsy > 0:   move the beamstop in the +row direction on the detector
```

The column direction corresponds to the **x** axis of the beamstop
motor stage; row corresponds to **y**.  Verify the sign mapping for
your specific stage before automating corrections.

## See also

- {doc}`../user-guide/beamstop-centering` — narrative + integration
  example
- {doc}`../user-guide/calibration` — the same script produces both
  files
