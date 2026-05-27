# Beamstop centering

A pin beamstop sits in the beam path upstream of the detector.  If
it's not centered on the beam, the SAXS signal is partially blocked or
the beamstop transmission leaks signal into the data.  The
`calibrate_smi_z_scan.py` script extracts the per-`motor_z` offset
needed to keep the beamstop centered.

## What gets exported

`scripts/calibrate_smi_z_scan.py` writes two files:

- `/tmp/saxs_calibration.json` — the loader-side constants (see
  {doc}`calibration`).
- `/tmp/smi_beamstop_offsets.json` — the beamstop offset table.

The package also ships the current best calibration at
`src/smi_tiled/data/smi_beamstop_offsets.json` so the collection code
can read it directly.

## Schema

```json
{
  "_doc": "Pin beamstop centering offsets vs detector motor_z, derived from <UID>.",
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
    …
  ],
  "spline": {
    "available": true,
    "model": "scipy.interpolate.BSpline (cubic, smoothing)",
    "degree": 3,
    "sigma_dx_mm": 0.126,
    "sigma_dy_mm": 0.149,
    "knots_d_bsx_mm": [...],
    "coefs_d_bsx_mm": [...],
    "knots_d_bsy_mm": [...],
    "coefs_d_bsy_mm": [...],
    "rms_residual_dx_mm": 0.126,
    "rms_residual_dy_mm": 0.149,
    "_consumer_note": "To evaluate: from scipy.interpolate import BSpline; ...",
    "dense_grid": {
      "motor_z_mm": [1700.0, 1738.4, …],
      "d_bsx_mm":   [0.220, 0.224, …],
      "d_bsy_mm":   [-1.416, -1.408, …]
    }
  }
}
```

## Why a spline (and not just a linear fit)?

The offset vs motor_z **meanders** — it's not a simple linear function.
A smoothing cubic B-spline tracks the structure while filtering out
per-z noise from averaging over `piezo_z` at each `motor_z`.

The σ-floored smoothing factor (`sigma_dx`, `sigma_dy`) is estimated
from the median absolute deviation between adjacent points and capped
at 0.05 mm so the fit doesn't degenerate to interpolation on
suspiciously-smooth runs.

Typical residuals: **0.1–0.15 mm** in both axes over a 0.5–2 mm range.

## Evaluating the spline

If the consumer has scipy:

```python
import json
import numpy as np
from scipy.interpolate import BSpline

with open("src/smi_tiled/data/smi_beamstop_offsets.json") as f:
    bs = json.load(f)

spl_dx = BSpline(
    np.array(bs["spline"]["knots_d_bsx_mm"]),
    np.array(bs["spline"]["coefs_d_bsx_mm"]),
    bs["spline"]["degree"],
)
spl_dy = BSpline(
    np.array(bs["spline"]["knots_d_bsy_mm"]),
    np.array(bs["spline"]["coefs_d_bsy_mm"]),
    bs["spline"]["degree"],
)

# Evaluate at any motor_z in the calibrated range:
mz = 5500.0
d_bsx, d_bsy = float(spl_dx(mz)), float(spl_dy(mz))
print(f"At motor_z={mz} mm, beamstop offset: ({d_bsx:+.4f}, {d_bsy:+.4f}) mm")
```

If the consumer doesn't have scipy, fall back to linear interpolation
through `spline.dense_grid` (200 points across the calibrated range):

```python
import numpy as np

z = np.array(bs["spline"]["dense_grid"]["motor_z_mm"])
dx = np.array(bs["spline"]["dense_grid"]["d_bsx_mm"])
dy = np.array(bs["spline"]["dense_grid"]["d_bsy_mm"])
mz = 5500.0
d_bsx = float(np.interp(mz, z, dx))
d_bsy = float(np.interp(mz, z, dy))
```

## Sign convention

```
d_bsx > 0:   the beam center is at +col of the beamstop shadow
             → move beamstop in the +col direction to center
d_bsy > 0:   the beam center is at +row of the beamstop shadow
             → move beamstop in the +row direction to center
```

The mm → motor-units conversion depends on the geometry of the
beamstop stage; for a beamstop mounted close to the detector, the
pixel-to-motor ratio is ≈ 1.0.  Consult the beamline mechanical
documentation for the exact mapping.

## Integration with the collection code

The intended consumer is the SMI collection code (e.g. `smi-browser`
or a custom ophyd plan).  It reads the current `motor_z` setpoint,
evaluates the spline at that z, and writes the corrected beamstop
positions back to EPICS.

A minimal example:

```python
def adjust_beamstop_for_motor_z(motor_z_setpoint, current_bsx, current_bsy):
    import json
    import numpy as np
    from scipy.interpolate import BSpline
    from pathlib import Path

    bs_file = Path(__file__).parent / "smi_beamstop_offsets.json"
    bs = json.loads(bs_file.read_text())
    spl_dx = BSpline(
        np.array(bs["spline"]["knots_d_bsx_mm"]),
        np.array(bs["spline"]["coefs_d_bsx_mm"]),
        bs["spline"]["degree"],
    )
    spl_dy = BSpline(
        np.array(bs["spline"]["knots_d_bsy_mm"]),
        np.array(bs["spline"]["coefs_d_bsy_mm"]),
        bs["spline"]["degree"],
    )
    d_bsx = float(spl_dx(motor_z_setpoint))
    d_bsy = float(spl_dy(motor_z_setpoint))
    return current_bsx + d_bsx, current_bsy + d_bsy
```

## Re-running the calibration

When the beamstop is physically realigned, the offset table goes
stale.  Re-run:

```bash
pixi run python scripts/calibrate_smi_z_scan.py <new_AGB_z_scan_UID>
cp /tmp/smi_beamstop_offsets.json src/smi_tiled/data/smi_beamstop_offsets.json
```

The collection code should re-read the file periodically (or on each
plan startup) so the new offsets take effect.

## See also

- {doc}`calibration` — the geometry constants the same script produces
- {doc}`../reference/beamstop-offsets-format` — full schema
- `scripts/calibrate_smi_z_scan.py` — the analysis source
