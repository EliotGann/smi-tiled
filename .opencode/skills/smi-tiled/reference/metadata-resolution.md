# Metadata resolution

Loading SMI scans involves resolving ~30 separate scalars (energy, beam center, SDD, motor positions, beamstop choice, …) from heterogeneous sources. This file documents the resolution helpers, fallback chains, and gotchas for each.

## The canonical fallback chain

For most scalars, the loader walks this priority list (highest → lowest):

```
1. User override (kwarg passed to resolve_*_geometry / load_saxs_raw)
2. Primary stream — per-frame value, first frame in chronological order
3. Baseline stream — snapshot at scan start
4. Primary configuration metadata (run["primary"].metadata["configuration"])
5. Start-doc metadata / sample_name string parsing
6. Hardcoded instrument default in defaults.py
```

This is documented inline at `loader.py:18-19` and is implemented field-by-field in `resolve_saxs_geometry` (`loader.py:1180-1378`) and `resolve_waxs_geometry` (`loader.py:1380-1471`).

There are exceptions — energy uses a different ordering (start metadata wins over baseline), and `resolve_incident_angle` uses a 5-tier chain. They are documented below.

## Low-level helpers

### `_run_uid(run) -> str` (`loader.py:616`)
Extract the run's UID for cache keying. Used by `_BASELINE_CACHE`, `_BASELINE_COLUMNS_CACHE`, and the geometry caches.

### `_as_scalar(value) -> Any` (`loader.py:604`)
Reduce an array-like value to a single Python scalar. Returns `None` if the array is empty. Handles 0-d numpy arrays, 1-element 1-d arrays, and Python scalars uniformly.

### `_apply_primary_sort(arr, run) -> ndarray` (`loader.py:664`)
Reorder a primary-stream array so element 0 is the chronologically-first frame. Necessary because tiled may return tabular data in arbitrary order when frames were saved in parallel.

## Stream readers

### `_read_baseline(run) -> xr.Dataset | None` (`loader.py:681`)
Read the full baseline stream as an `xr.Dataset`, with backwards/forwards-compat for both layouts:
- Old: `run["baseline"].read()` returns an xr.Dataset directly.
- New (`bluesky-tiled-plugins`): `run["baseline"]["internal"]` is a `DataFrameClient`; the helper reads its `pd.DataFrame` and converts via `xr.Dataset.from_dataframe`.

Results are cached per UID in `_BASELINE_CACHE` (process-global).

> **WARNING**: line 699 catches `(KeyError, Exception)`, swallowing **everything** silently. If baseline appears empty, set a breakpoint here and inspect.

### `_baseline_scalar(run, key) -> Any` (`loader.py:720`)
Fast path: pull one column from baseline. Order:
1. If `_BASELINE_CACHE[uid]` exists, read column from cached Dataset.
2. Otherwise: `run["baseline"]["internal"][key].read()` — single-column tiled fetch (avoids reading all 564 columns).
3. Fallback: `_read_baseline(run)` then `_dataset_scalar`.

The single-column path is critical for performance — full baseline can be 564 columns × multiple values, and `xr.Dataset.from_dataframe` is slow.

### `_dataset_scalar(ds, key) -> Any` (`loader.py:750`)
Read a column from an xr.Dataset (returned by `_read_baseline`). Returns `None` if missing.

### `_primary_conf(run, det_key) -> dict` (`loader.py:760`)
Read the per-detector configuration metadata. Returns `run["primary"].metadata["configuration"][det_key]["data"]`. Empty dict on any failure. The `det_key` is the detector handle (`"pil2M"`, `"pil900KW"`).

### `_conf_scalar(conf, key) -> Any` (`loader.py:756`)
Read a scalar from a `_primary_conf`-returned dict. Same semantics as `_dataset_scalar` for plain dicts.

### `_primary_scalar(run, field) -> float | None` (`loader.py:772`)
Read a scalar **only if** the field has a single unique value across all frames. Returns `None` if the field is varying (i.e., is a scan axis) — for varying axes use `_read_scan_axis` instead. Sorts to chronological order before taking element 0.

### `_read_first_scalar(run, field, baseline_keys=(), baseline_ds=None) -> float | None` (`loader.py:800`)
Walk the primary→baseline chain:
1. `_primary_scalar(run, field)`
2. For each `key` in `baseline_keys` (or `(field,)` if not given):
   1. `_baseline_scalar(run, key)`
   2. `_dataset_scalar(baseline_ds, key)`
3. Return `None`.

`baseline_keys` is typically a tuple like `("pil2M_motor_z_user_setpoint", "pil2M_motor_z")` so the setpoint takes priority over the readback.

### `_has_primary_field(run, field)` / `_has_primary_internal_field(run, field)`
Existence checks before attempting to read. `_primary_scalar` doesn't gate on these and returns `None` on any error.

### `_read_primary_internal_array(run, field) -> ndarray | None` (`loader.py:845`)
Read a per-frame array from `run["primary"]["internal"]`. Returns `None` if absent.

### `_read_scan_axis(run, field, cache_path=None) -> ndarray | None` (`loader.py:1799`)
Find a varying scan axis. Walks: cache → primary internal → primary data → primary base/internal → baseline. The result is sorted into chronological order via `_apply_primary_sort` before being returned.

### `_primary_per_frame_motor(run, field, n_frames, cache_path=None) -> ndarray | None` (`loader.py:999`)
Used by `resolve_incident_angle`. Returns a length-`n_frames` array iff the field is *measured in primary* (gates on `_has_primary_field`/`_has_primary_internal_field`). Without that gating, `_read_scan_axis` would silently fall back to baseline, which would fool downstream logic that intends to detect "is this a scanned axis?".

A length-1 result is broadcast to `n_frames`. A mismatched length yields `None`.

## High-level resolvers

### `resolve_saxs_geometry(run, energy_kev=None, **overrides) -> SAXSGeometry` (`loader.py:1180`)

Returns a `SAXSGeometry` dataclass (`loader.py:1119`):

```python
@dataclass
class SAXSGeometry:
    dist_m: float
    poni1_m: float                # = beam_center_row_px * pixel_size
    poni2_m: float                # = beam_center_col_px * pixel_size
    energy_ev: float
    wavelength_m: float
    beam_center_row_px: float
    beam_center_col_px: float
    active_beamstop: str = "rod"  # "rod" | "pin"
    beamstop_pos_mm: dict = ...   # {"rod": {"x": ..., "y": ...}, "pin": {...}}
    motor_x_mm: float | None = None
    motor_y_mm: float | None = None
    motor_z_mm: float | None = None
    piezo_z_um: float | None = None
```

Per-field resolution chains:

| Field | Chain |
|---|---|
| `energy_kev` | override → `start.energy` → `baseline.energy_energy / 1000` → `sample_name._keV` → `DEFAULT_ENERGY_KEV` |
| `beam_center_row_px` | override → `baseline.pil2M_beam_center_y_px` → `baseline_ds.…` → `conf.…` → `_SAXS_DEFAULT_BEAM_ROW_PX` |
| `beam_center_col_px` | override → `baseline.pil2M_beam_center_x_px` → … → default |
| `dist_mm` | override → `primary.pil2M_motor_z` → `baseline.pil2M_motor_z_user_setpoint` → `baseline.pil2M_motor_z` → … → `conf.pil2M_sdd_mm` → `sample_name._sddXm * 1000` → `_SAXS_DEFAULT_DISTANCE_MM` |
| `active_beamstop` | override → `baseline.pil2M_active_beamstop` → … → `conf.pil2M_active_beamstop` → `"rod"` |
| `beamstop_pos_mm` | baseline `saxs_beamstop_{x,y}_{rod,pin}_user_setpoint` → `saxs_beamstop_{x,y}_{rod,pin}` |

> **WARNING — Silent default to "rod"**. If `pil2M_active_beamstop` is missing everywhere, `resolve_saxs_geometry` returns `active_beamstop="rod"` with no warning (`loader.py:1245-1251`). Pin-beamstop scans missing this PV will get the wrong mask.

Then the resolver applies **motor-driven corrections** (`loader.py:1296-1361`):

```python
beam_col += (motor_x_mm - motor_x_ref_mm) * col_per_mx
beam_row += (motor_y_mm - motor_y_ref_mm) * row_per_my
beam_col += (motor_z_mm - motor_z_ref_mm) * col_per_mz
beam_row += (motor_z_mm - motor_z_ref_mm) * row_per_mz
dist_mm  += (piezo_z_um - piezo_z_ref_um) * sdd_per_pz
beam_row += beam_delta_row_px
beam_col += beam_delta_col_px
dist_mm  += distance_delta_mm
```

The slope constants (`col_per_mx`, etc.) all default to `0.0` in `defaults.py`. They become non-zero only when calibrated and committed to `data/saxs_calibration.json` or passed as overrides. See `reference/saxs-geometry.md` for derivation.

> **WARNING — Motor-slope double-count risk**. If a calibration JSON sets `motor_x_ref_mm` and `beam_col_px_per_motor_x_mm` to non-zero values **and** a scan stores a different motor_x reference frame in baseline, the correction can compound. Verify slope sources before applying.

### `resolve_waxs_geometry(run, energy_kev=None, **overrides) -> WAXSGeometry` (`loader.py:1380`)

Returns:

```python
@dataclass
class WAXSPanelGeometry:
    col_start: int
    col_end: int
    offset_deg: float

@dataclass
class WAXSGeometry:
    dist_m: float
    beam_center_row_px: float
    beam_center_col_px: float
    energy_ev: float
    wavelength_m: float
    theta_zero_deg: float = 0.0
    sample_offset_x_mm: float = 0.0
    sample_offset_z_mm: float = 0.0
    panels: tuple[WAXSPanelGeometry, ...] = ()
```

Per-field chains:

| Field | Chain |
|---|---|
| `energy_kev` | (same as SAXS) |
| `dist_mm` | override → `primary.pil900KW_motor_z` → `baseline.pil900KW_motor_z_user_setpoint` → `baseline.pil900KW_motor_z` → … → `conf.pil900KW_sdd_mm` → `_WAXS_DEFAULT_DISTANCE_MM` |
| `beam_center_row_px` | override → `_WAXS_DEFAULT_BEAM_ROW_PX` |
| `beam_center_col_px` | override → `_WAXS_DEFAULT_BEAM_COL_PX` |
| `panel_offsets_deg` | override → `_WAXS_DEFAULT_PANEL_OFFSETS_DEG` (= `(-7.0, 0.0, 7.0)`) |
| `panel_col_ranges` | override → `_WAXS_DEFAULT_PANEL_COL_RANGES` (= `((0,206),(206,413),(413,619))`) |

> **NOTE — WAXS does not read beam center from baseline**. The ophyd config stores WAXS beam center in the raw (pre-rotation) coordinate system, which is incompatible with the rotated `MultiPanelArcDetector` frame. The resolver hard-codes the calibrated default; only `beam_delta_*_px` overrides take effect (`loader.py:1425-1436`).

### `resolve_incident_angle(run, n_frames, *, manual_override=None, theta_offset=0.0, cache_path=None) -> tuple[ndarray | None, str]` (`loader.py:1025`)

5-tier chain returning per-frame `ai` array AND a source description string:

| Tier | Source | Offset applied? |
|---|---|---|
| 0 | `manual_override` | No |
| 1 | Primary `stage_th + piezo_th` (per-frame) | Yes (+`theta_offset`) |
| 2 | Per-frame `target_file_name` parsed `_ai{X}_` | No |
| 3 | start-doc `sample_name` parsed `_ai` (then `_th`) | No |
| 4 | Baseline `stage_th + piezo_th` | Yes (+`theta_offset`) |

The `theta_offset` argument converts a raw motor reading to a true incident angle. It is applied **only** to motor-derived paths (1 and 4), never to angles parsed from filenames or sample names (which already record the commanded angle).

For tier 1, individual motors fall back to baseline if not in primary — so a scan that scans only `piezo_th` while `stage_th` sits in baseline still produces a per-frame array.

### `parse_sample_name_geometry(sample_name) -> dict[str, float]` (`loader.py:949`)

Regex parse of the sample_name string. Recognized patterns:

| Pattern | Returns |
|---|---|
| `_wa<digits>` | `{"waxs_arc_deg": float}` |
| `_sdd<digits>m?` | `{"sdd_m": float}` |
| `_<digits>keV` | `{"energy_kev": float}` |
| `_<digits>eV` | `{"energy_kev": float / 1000}` |
| `_ai<digits>` | `{"incident_angle_deg": float}` |
| `_th<digits>` | `{"theta_deg": float}` |

This is a **last-resort** fallback. Always prefer baseline-derived values when present.

### `_read_target_file_name_geometry(run, cache_path=None)` (`loader.py:918`)

Per-frame variant of `parse_sample_name_geometry`. Reads the `target_file_name` column (if present), decodes from bytes, and runs the regex per frame. Used by `resolve_incident_angle` tier 2.

The column lives in different places across layouts; the helper tries (in order):
1. The HDF5 disk cache.
2. `primary["internal"]`.
3. `primary.base["internal"]` (newer migration catalog).
4. `primary` / `primary["data"]` direct.

## DataArray attrs after `load_*_raw`

After resolution, geometry is attached to the returned DataArray as attrs (`loader.py:1980-2000`):

```python
da.attrs["dist"]                  = geo.dist_m
da.attrs["poni1"]                 = geo.poni1_m
da.attrs["poni2"]                 = geo.poni2_m
da.attrs["pixel1"]                = PILATUS_PIXEL_SIZE_M
da.attrs["pixel2"]                = PILATUS_PIXEL_SIZE_M
da.attrs["energy"]                = geo.energy_ev
da.attrs["wavelength"]            = geo.wavelength_m
da.attrs["smi_active_beamstop"]   = geo.active_beamstop
da.attrs["smi_beamstop_pos_mm"]   = geo.beamstop_pos_mm
da.attrs["smi_motor_x_mm"]        = geo.motor_x_mm   # may be None
... and so on
```

These are **scalar** attrs only. Per-frame variation (motor positions during a scan, ai, waxs_arc) is stored as DataArray **coords** on the appropriate dim (`frame`, `waxs_arc`).

## Resolution-debugging recipe

When you don't trust a resolved value:

```python
from smi_tiled.loader import (
    TiledSMISWAXSLoader, resolve_saxs_geometry, resolve_waxs_geometry,
    _baseline_scalar, _primary_scalar, _primary_conf, _read_baseline,
    parse_sample_name_geometry,
)

loader = TiledSMISWAXSLoader()
run = loader._get_run(uid)

# Inspect each layer of the chain individually:
print("override:    (none)")
print("primary:    ", _primary_scalar(run, "pil2M_motor_z"))
print("baseline:   ", _baseline_scalar(run, "pil2M_motor_z_user_setpoint"))
print("baseline2:  ", _baseline_scalar(run, "pil2M_motor_z"))
print("conf:       ", _primary_conf(run, "pil2M").get("pil2M_sdd_mm"))
print("name_geo:   ", parse_sample_name_geometry(run.metadata["start"].get("sample_name", "")))

# Now compare with what the resolver chooses:
geo = resolve_saxs_geometry(run)
print("resolved dist_m:", geo.dist_m)
```

If the resolver chose a value you don't expect, the chain inspection will tell you which layer it came from.

## What to look at when a number looks wrong

| Symptom | Most likely culprit |
|---|---|
| Energy off by 1000× | `start.energy` is in keV, baseline is in eV — confirm units |
| Beam center swapped X↔Y | `beam_center_row_px` ↔ `beam_center_col_px` confusion (rows = y/vertical, cols = x/horizontal) |
| SDD off by ~50 mm | piezo_z slope (`sdd_per_pz`) using a stale reference |
| Pin-beamstop scan masked as rod | `pil2M_active_beamstop` missing → silent default to `"rod"` |
| Per-frame motor not detected | `_has_primary_field` false → `_primary_per_frame_motor` returns None |
| ai stuck at one value | tier-3 sample_name caught a value before tier-1 motor lookup; check primary fields |
| target_file_name regex fails | bytes vs str — `_coerce_str` handles it, but custom code may not |

Cross-reference with `reference/gotchas.md` for the full list.
