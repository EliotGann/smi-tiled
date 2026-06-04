# SAXS geometry (Pilatus 2M)

This file covers the SAXS detector layout, beam-center & SDD math, q-from-pixel derivation, and the calibration JSON schema.

## Detector specs

- **Detector**: Pilatus 2M (Dectris).
- **Image shape**: `(rows, cols) = (1679, 1475)`. Row index = vertical (y), column index = horizontal (x).
- **Pixel size**: `0.172 mm × 0.172 mm` (square).
- Constants: `PILATUS_PIXEL_SIZE_M = 0.172e-3` (`loader.py:60`).

## Coordinate convention

- Beam center: `(beam_center_row_px, beam_center_col_px)` — both in pixel units, both positive (loader takes `abs` to recover from negative-stored values, `loader.py:1291-1294`).
- pyFAI poni: `poni1 = bc_row_px * pixel_size`, `poni2 = bc_col_px * pixel_size` (meters).
- Sample-to-detector axis: along `+z`. Beam direction is `+z`.
- Detector pixel `(r, c)` is at lab-frame `(x, y, z)` where:
  - `x = (c - bc_col) * pixel_size` (horizontal, sign convention: + = away from beam in column direction)
  - `y = -(r - bc_row) * pixel_size` (vertical; **negative sign** because row index increases downward, but lab y increases upward)
  - `z = dist_m`

Source: `integrator.py:2960-2962`.

## SAXSGeometry dataclass

`SAXSGeometry` is the resolved summary returned by `resolve_saxs_geometry` (`loader.py:1119-1146`):

```python
@dataclass
class SAXSGeometry:
    dist_m: float                # sample-detector distance (meters)
    poni1_m: float               # beam center row × pixel_size
    poni2_m: float               # beam center col × pixel_size
    pixel1_m: float = PILATUS_PIXEL_SIZE_M
    pixel2_m: float = PILATUS_PIXEL_SIZE_M
    energy_ev: float
    wavelength_m: float
    beam_center_row_px: float
    beam_center_col_px: float
    active_beamstop: str = "rod"            # "rod" | "pin"
    beamstop_pos_mm: dict = ...             # {"rod": {"x": ..., "y": ...}, "pin": {...}}
    motor_x_mm: float | None = None         # reference value (first frame)
    motor_y_mm: float | None = None
    motor_z_mm: float | None = None
    piezo_z_um: float | None = None
```

## Default constants

In `loader.py:79-122` (some overridable from `data/saxs_calibration.json`):

| Constant | Value | Purpose |
|---|---|---|
| `_SAXS_DEFAULT_DISTANCE_MM` | `2000.0` | SDD fallback if all metadata absent |
| `_SAXS_DEFAULT_BEAM_ROW_PX` | `1165.0` | Beam center row fallback |
| `_SAXS_DEFAULT_BEAM_COL_PX` | `746.0` | Beam center col fallback |
| `_SAXS_DEFAULT_DISTANCE_DELTA_MM` | `-13.68` (calibrated) | Additive correction to motor_z readback → true SDD |
| `_SAXS_DEFAULT_BEAM_DELTA_ROW_PX` | `0.0` | Additive correction to beam_center_row |
| `_SAXS_DEFAULT_BEAM_DELTA_COL_PX` | `0.0` | Additive correction to beam_center_col |
| `_SAXS_MOTOR_X_REF_MM` | `1.88` | Reference motor_x (slope is 0; constant kept for override API) |
| `_SAXS_MOTOR_Y_REF_MM` | `2.45` | Reference motor_y |
| `_SAXS_MOTOR_Z_REF_MM` | `0.0` | Reference motor_z |
| `_SAXS_PIEZO_Z_REF_UM` | `0.0` | Reference piezo_z |
| `_SAXS_BEAM_COL_PX_PER_MOTOR_X_MM` | `0.0` | (**by design** — see note below) |
| `_SAXS_BEAM_ROW_PX_PER_MOTOR_Y_MM` | `0.0` | (**by design**) |
| `_SAXS_BEAM_COL_PX_PER_MOTOR_Z_MM` | `+2.65e-4` (calibrated) | Beam center drift with detector translation |
| `_SAXS_BEAM_ROW_PX_PER_MOTOR_Z_MM` | `-3.50e-4` (calibrated) | |
| `_SAXS_SDD_DELTA_MM_PER_PIEZO_Z_UM` | `+9.03e-4` (calibrated) | SDD correction with sample piezo |

> **CRITICAL NOTE on motor_x/y slopes (`loader.py:90-110`)**: The EPICS beam center PVs (`pil2M_beam_center_x_px`, `pil2M_beam_center_y_px`) **already track** detector translations. They drift with `motor_y` at ~6 px/mm naturally. Setting a non-zero `BEAM_COL_PX_PER_MOTOR_X_MM` slope here would **double-count** the correction. The historical regression that found ~5.82/5.99 px/mm is the rate at which the PV itself moves, not an extra correction. Default 0.0 is correct unless the underlying PV behavior changes.

## Motor-driven corrections — the math

`resolve_saxs_geometry` applies these adjustments to the metadata-resolved beam center & distance (`loader.py:1339-1361`):

```python
# Beam center (px)
beam_col += (motor_x_mm - motor_x_ref_mm) * col_per_mx
beam_row += (motor_y_mm - motor_y_ref_mm) * row_per_my
beam_col += (motor_z_mm - motor_z_ref_mm) * col_per_mz
beam_row += (motor_z_mm - motor_z_ref_mm) * row_per_mz

# SDD (mm)
dist_mm  += (piezo_z_um - piezo_z_ref_um) * sdd_per_pz
dist_mm  += distance_delta_mm

# Final additive corrections
beam_row += beam_delta_row_px
beam_col += beam_delta_col_px
```

The last two `beam_delta_*` are intended for fine-tuning at the call site (kwargs to `load_saxs_raw`). The deltas can also come from the calibration JSON.

## Calibration JSON

`src/smi_tiled/data/saxs_calibration.json` (loaded at import time by `_apply_calibration_override`, `loader.py:166-209`).

Schema (all keys optional):
```json
{
  "_doc": "...",
  "source_uid": "<uid of the AGB calibration scan>",
  "energy_ev": 16099.99,
  "wavelength_nm": 0.0770088,
  "active_beamstop": "pin",
  "scan_fixed_motor_x_mm": 60.009,
  "scan_fixed_motor_y_mm": 10.000,
  "regressions": {
    "bc_col": {"intercept": ..., "per_motor_z_mm": ..., "per_piezo_z_um": ..., "rms_px": ...},
    "bc_row": {...},
    "sdd_mm": {...}
  },
  "constants": {
    "_SAXS_DEFAULT_DISTANCE_DELTA_MM": -13.676,
    "_SAXS_BEAM_COL_PX_PER_MOTOR_Z_MM":  0.000265,
    "_SAXS_BEAM_ROW_PX_PER_MOTOR_Z_MM": -0.000350,
    "_SAXS_SDD_DELTA_MM_PER_PIEZO_Z_UM": 0.000903
  }
}
```

The `constants` block is the only part `_apply_calibration_override` honors. Allowed keys (`loader.py:191-201`):

```
_SAXS_DEFAULT_DISTANCE_DELTA_MM
_SAXS_MOTOR_X_REF_MM, _SAXS_MOTOR_Y_REF_MM
_SAXS_MOTOR_Z_REF_MM, _SAXS_PIEZO_Z_REF_UM
_SAXS_BEAM_COL_PX_PER_MOTOR_X_MM
_SAXS_BEAM_ROW_PX_PER_MOTOR_Y_MM
_SAXS_BEAM_COL_PX_PER_MOTOR_Z_MM
_SAXS_BEAM_ROW_PX_PER_MOTOR_Z_MM
_SAXS_SDD_DELTA_MM_PER_PIEZO_Z_UM
_SAXS_DEFAULT_BEAM_DELTA_ROW_PX
_SAXS_DEFAULT_BEAM_DELTA_COL_PX
```

Unknown keys produce `UserWarning` and are ignored.

The current bundled JSON ships with the calibration from AgB scan `b900e711-35a8-4dbc-8afa-2a1e20056608` (run by `scripts/calibrate_smi_z_scan.py`).

## q-from-pixel derivation

The reduction pipeline computes a per-pixel `(q, chi, qh, qv, qz)` map from `(dist_m, poni1, poni2, pixel1, pixel2, wavelength_m, shape)`. Implementation lives in `integrate_saxs` (`integrator.py:2884-2970`):

```python
rr, cc = np.meshgrid(np.arange(ny, dtype=float),
                     np.arange(nx, dtype=float), indexing="ij")
bc_row = poni1_m / pixel1_m
bc_col = poni2_m / pixel2_m

x_m = (cc - bc_col) * pixel2_m            # horizontal, m
y_m = -(rr - bc_row) * pixel1_m           # vertical (sign flip), m
r_m = np.sqrt(x_m**2 + y_m**2 + dist_m**2)  # pixel-to-sample distance

k = 2.0 * np.pi / wavelength_nm           # wavenumber, nm^-1

qh2d = k * x_m / r_m                       # horizontal component (nm^-1)
qv2d = k * y_m / r_m                       # vertical component
qz2d = k * (dist_m / r_m - 1.0)           # along-beam component
q2d  = np.sqrt(qh2d**2 + qv2d**2 + qz2d**2)
chi_deg_2d = np.rad2deg(np.arctan2(qh2d, qv2d))
```

### Derivation sketch

The scattering vector for an elastic interaction is
```
q = k_out - k_in
|k_out| = |k_in| = 2π/λ
k_in    = (0, 0, k)        # along +z (toward detector)
k_out   = k * (x, y, dist) / r    where r = sqrt(x² + y² + dist²)
```

So:
```
q = k * (x/r, y/r, dist/r - 1)
qh = k * x/r       (units: same as k, nm^-1)
qv = k * y/r
qz = k * (dist/r - 1)         # ≤ 0 for forward scattering
|q| = 2k * sin(θ)              where 2θ = scattering angle
```

### Units in attrs vs in integrator

`load_saxs_raw` writes `attrs["wavelength"] = geo.wavelength_m * 1e10` (`loader.py:1985`), i.e. **Ångström**, despite the suggestive `wavelength_m` name on the geometry object. The integrator then converts back: `wavelength_m = float(attrs["wavelength"]) * 1e-10` (`integrator.py:2919`). Net result: `q` in the integrator is in `nm^-1` (because `k = 2π/wavelength_nm`).

> **GOTCHA**: `attrs["wavelength"]` is in **Å**, not meters. The variable name `wavelength_m` in the geometry dataclass refers to meters; the attr conversion happens at the boundary. Don't assume from variable names.

**Conversion**: `q_inv_A = q_inv_nm / 10`.

### Worked example — pin beamstop on scan 1130041

For the SMI defaults at 16.1 keV with SDD = 2.0 m:
- `λ = 12.398 keV·Å / 16.1 keV = 0.770 Å = 0.0770 nm`
- `k = 2π / 0.0770 nm = 81.6 nm^-1 = 8.16 Å^-1`

Measured `q_min` = 0.121 Å^-1 (with pin mask). Inverting:
- `q = (4π/λ) sin(θ)` → `sin(θ) = 0.121 × 0.770 / (4π) = 0.00742`
- `2θ ≈ 0.01484 rad ≈ 0.85°`
- `pixel offset = SDD × tan(2θ) ≈ 2000 mm × 0.01484 = 29.7 mm = 173 px`

So the pin polygon edge sits **~173 px** from the beam center along the column direction (not 17 px as one might naively read off the polygon spec). The corresponding column delta `(beam_col + dc) - bc_col ≈ 173 px` puts the inner mask boundary at `q ≈ 0.121 Å^-1`. Cross-checked against the `make_saxs_mask_from_dict` polygon shift in `integrator.py:565-580`.

If you find a `q_min` that doesn't match the expected mask geometry, trace back through:
1. The actual polygon vertices (in `data/masks/pil2M_mask_polygons.json`).
2. The shift applied by `shift_polygon` (which adds `(dc, dr)` to every vertex).
3. The beam center used by `integrate_saxs` (`bc_col = poni2/pixel2`).
4. The wavelength/SDD used to convert pixel offset → q.

### chi convention

`chi_deg_2d = arctan2(qh, qv)` (`integrator.py:2970`). With:
- `qh > 0, qv = 0` (pixel directly to the right of beam): `chi = +90°`
- `qh = 0, qv > 0` (pixel directly above beam): `chi = 0°`
- `qh = 0, qv < 0` (pixel directly below beam): `chi = ±180°`
- `qh < 0, qv = 0` (pixel directly to the left of beam): `chi = -90°`

So **chi=0 points up** (vertical), **chi=+90 points right** (horizontal). That's `(qv, qh)` order in `arctan2`, not the conventional `(y, x)` of `arctan2`.

## Solid-angle correction

Per-pixel solid angle factor (`integrator.py:2972-2974`):
```
sa[r,c] = pixel_area / r³ * max(dist_m, 0.0)
        = (pixel1·pixel2 · dist) / (x² + y² + dist²)^(3/2)
```

This corrects pixel intensities for the cos(2θ) factor that arises because pixels far from the beam center subtend smaller solid angles (geometric foreshortening). Applied via division (`I_corrected = I_raw / sa`) when `solid_angle_correction=True` (default in `reduce_smi_combined`).

## Geometry cache

`_SAXS_GEOMETRY_CACHE` (`integrator.py:2945`) memoizes `(q2d, qh2d, qv2d, chi_deg_2d, sa_base)` keyed on `_saxs_cache_key(dist_m, poni1_m, poni2_m, pixel1_m, pixel2_m, wavelength_m, shape)`. See `reference/caching.md`.

> **WARNING**: If you edit the geometry build code without bumping any of these inputs, the cache will hand back a stale map. Use `clear_geometry_cache()` after edits.

## Beam center & motor positions in DataArray attrs

After `load_saxs_raw`, the resolved scalars are mirrored into `da.attrs` (`loader.py:1980-2000`):

```python
da.attrs = {
    "dist": geo.dist_m,
    "poni1": geo.poni1_m,
    "poni2": geo.poni2_m,
    "pixel1": pixel_size,
    "pixel2": pixel_size,
    "energy": geo.energy_ev,
    "wavelength": geo.wavelength_m,            # in METERS (verify before use)
    "rot1": 0, "rot2": 0, "rot3": 0,
    "smi_motor_x_mm":  ...,
    "smi_motor_y_mm":  ...,
    "smi_motor_z_mm":  ...,
    "smi_piezo_z_um":  ...,
    "smi_active_beamstop": "rod" | "pin",
    "smi_beamstop_pos_mm": {"rod": {"x":..., "y":...}, "pin": {...}},
}
```

Per-frame variation (e.g. when a scan steps through `pil2M_motor_y`) is stored as `coords` on the `frame` dim, **not** as attrs. Attrs are scalar references.

## Calibrating new scans

To calibrate a new (energy, sample, detector position):
1. Run `pixi run python scripts/calibrate_smi_saxs.py` for a single AGB exposure to fit beam center directly.
2. Run `pixi run python scripts/calibrate_smi_z_scan.py <UID>` for a z-grid scan to fit the per-motor and per-piezo slopes.
3. Inspect the produced `saxs_calibration.json`; commit when satisfied.
4. Restart any long-running process so `_apply_calibration_override` re-runs at import.

See also `reference/validation.md`.
