# Data layout & Tiled access

This document covers how SMI data is organized in Tiled, how the loader reaches it, and which fields/streams hold what.

## Tiled connection

SMI data is served at `https://tiled.nsls2.bnl.gov` under the catalog `smi/migration` (legacy + post-bluesky-tiled-plugins data unified).

```python
DEFAULT_TILED_URI = "https://tiled.nsls2.bnl.gov"
DEFAULT_CATALOG   = "smi/migration"
```
(`src/smi_tiled/defaults.py`)

The `TiledSMISWAXSLoader` class wraps the connection (`loader.py:2588-2630`):

```python
loader = TiledSMISWAXSLoader(
    tiled_uri="https://tiled.nsls2.bnl.gov",
    catalog="smi/migration",
    energy_kev=None,    # override; falls back to run metadata
    api_key=None,       # set when authentication is required
)
```

Internally:
- `_get_root_client()` (`loader.py:2651`) creates the tiled root client lazily on first access.
- `_get_catalog()` resolves the slash-separated catalog path against the root.
- `login(...)` (`loader.py:2661`) and `logout()` are exposed for interactive sessions.

## Run lookup: scan_id vs uid

A "run" in Tiled is keyed by **uid** (a UUID), but humans usually reference scans by **scan_id** (an integer counter).

- `loader.searchCatalog(...)` (`loader.py:2860`) — query by any combination of `scan_id`, `plan_name`, `sample_name`, `proposal_id`, `cycle`, `user_name`, `institution`, `detector`, `since`/`until`, plus arbitrary `key=value` kwargs. Returns a `pandas.DataFrame` with columns including `uid`, `scan_id`, `sample_name`, `plan_name`, `start_time`, `n_frames`.
- `loader.peekAtMd(uid)` (`loader.py:2725`) — minimal metadata-only inspection of a run; returns a dict-of-dicts of `start`, `stop`, `streams`, etc. without loading any image.
- `loader.summarizeRun(uid)` (`loader.py:3034`) — slightly richer summary (detectors present, n_frames per detector, sample_name, plan_name, scan_id).
- `loader.loadRun(uid, ...)` (`loader.py:2599`) — load both detectors at once.
- `loader.loadSingleImage(uid, frame_idx, detector="saxs"|"waxs")` (`loader.py:2746`) — single-frame fast path.

> **WARNING.** `searchCatalog` eagerly imports `Key, Regex, TimeRange` from `tiled.queries` (`loader.py:2927`). Some installed `tiled` versions do not export `TimeRange`, breaking the call. Workaround:
> ```python
> from tiled.queries import Key
> node = loader._get_catalog().search(Key("scan_id") == int(scan_id))
> uid, run = next(iter(node.items()))
> ```

## Run streams

Each run is a `BlueskyRun`-shaped node with these streams:

| Stream | Purpose | How to read |
|---|---|---|
| `start` | Scan start metadata (sample_name, plan_name, energy, scan_id, …) | `run.metadata.get("start", {})` |
| `stop` | Scan stop status (success/failure, exit time) | `run.metadata.get("stop", {})` |
| `primary` | Per-frame: images, motors, detectors, scan-axis values | `run["primary"]["data"][field]` or `run["primary"][field]` |
| `primary.config[<det_key>]` | Per-detector configuration scalars (one-shot per scan) | `_primary_conf(run, "pil2M")` |
| `baseline` | Snapshot of all PVs at start (and sometimes end) of scan | `run["baseline"].read()` (legacy) or `run["baseline"]["internal"]` (new layout) |

The "new layout" referenced above is the `bluesky-tiled-plugins` rework where baseline data became a `DataFrameClient` accessed via `run["baseline"]["internal"]` rather than `.read()`. `_read_baseline` (`loader.py:681`) handles both.

> **WARNING.** `_read_baseline` swallows all exceptions silently (`except (KeyError, Exception)` at `loader.py:699`). If baseline is unexpectedly empty, that's the first place to check by hand.

## Image fields and dimensions

The loader uses these field-name constants (re-exported in `__init__.py:29-32`):

```python
SAXS_IMAGE_FIELD = "pil2M_image"      # Pilatus 2M, the SAXS detector
WAXS_IMAGE_FIELD = "pil900KW_image"   # Pilatus 900KW, the WAXS arc detector
WAXS_ARC_FIELD   = "waxs_arc"         # WAXS arc-rotation angle, degrees
WAXS_BSX_FIELD   = "waxs_bsx"         # WAXS beamstop X motor, mm
```

### SAXS image dims and shape

`load_saxs_raw(...)` (`loader.py:1899`) returns an `xarray.DataArray`:

| Scan type | Dims | Shape |
|---|---|---|
| Single-frame `count` | `("pix_y", "pix_x")` | `(1679, 1475)` |
| Multi-frame plan | `("frame", "pix_y", "pix_x")` | `(N, 1679, 1475)` |
| WAXS-arc-stepped scan | `("waxs_arc", "pix_y", "pix_x")` | `(N_arc, 1679, 1475)` |

Per-pixel size: `0.172 mm` (`PILATUS_PIXEL_SIZE_M = 0.172e-3`, `loader.py`).

The DataArray carries pyFAI-compatible geometry attrs (mirroring PyHyperScattering's `template_xr` convention):
```
dist, poni1, poni2, rot1, rot2, rot3, pixel1, pixel2, energy, wavelength
```

Where:
- `poni1 = beam_center_row_px * pixel1` (meters)
- `poni2 = beam_center_col_px * pixel2` (meters)
- `dist`, `wavelength` in meters; `energy` in eV.

SMI-specific extras live under the `smi_` prefix (panel geometry, arc angle, beamstop, motor positions).

### WAXS image dims and shape

`load_waxs_raw(...)` (`loader.py:2047`) returns:

| Scan type | Dims | Shape |
|---|---|---|
| Multi-arc (typical) | `("waxs_arc", "pix_y", "pix_x")` | `(N_arc, 619, 1475)` after rotation, OR `(N_arc, 195, 619)` per panel pre-rotation |

The "raw" image from the detector is `(195, 619)` — the 3 panels concatenated horizontally. The loader applies a `np.fliplr(np.rot90(image, k=3))` rotation at integration time (`integrator.py:488 rotate_image_and_mask`) so that the WAXS image matches the SAXS coordinate convention.

> **WARNING.** WAXS DataArrays carry pyFAI geometry attrs but the underlying detector is a 3-panel folded arc. The single-flat-panel pyFAI assumption is violated. Loader emits a warning at `loader.py:2615-2619`. Use `MultiPanelArcDetector` / `integrate_waxs` directly; do NOT pipe WAXS through `PFGeneralIntegrator`.

## Per-frame motor positions

Per-frame motor values live on the primary stream and can be read either via the loader's DataArray coords or through the helper `_primary_scalar(run, field)` (`loader.py:772`).

Most commonly used:

| Field | Type | Purpose |
|---|---|---|
| `pil2M_motor_x` | float, mm | SAXS detector X translation |
| `pil2M_motor_y` | float, mm | SAXS detector Y translation |
| `pil2M_motor_z` | float, mm | SAXS detector Z (= sample-detector distance) |
| `piezo_z` | float, μm | Sample piezo Z (fine SDD adjustment) |
| `pil900KW_motor_z` | float, mm | WAXS detector Z |
| `waxs_arc` | float, deg | WAXS rotation arc angle |
| `waxs_bsx` | float, mm | WAXS beamstop X position |
| `bsx`, `bsy` | float, mm | SAXS beamstop X/Y position (per active beamstop) |
| `pil2M_active_beamstop` | string | `"rod"` or `"pin"` — which SAXS beamstop is in beam |
| `energy_energy` | float, eV | Photon energy |

The `_primary_scalar(run, field)` helper enforces "this must have a single value" — if the field varies across frames (i.e., it's a scan axis), it returns `None`. For varying scan axes, use `_read_scan_axis` instead.

For chronological-order semantics (because Tiled may return tables in arbitrary order for parallel-saved data), `_apply_primary_sort` (`loader.py`, near `_primary_scalar`) reorders to seq-num-ascending before indexing.

## Baseline snapshot fields

Baseline contains all monitored PVs at start (and often end) of scan. Common keys you'll resolve from baseline:

| Key | Purpose |
|---|---|
| `pil2M_beam_center_x_px` | SAXS beam center, column px |
| `pil2M_beam_center_y_px` | SAXS beam center, row px |
| `pil2M_motor_z_user_setpoint` | Reference SDD (mm) |
| `pil2M_active_beamstop` | `"rod"` or `"pin"` |
| `saxs_beamstop_x_rod`, `saxs_beamstop_y_rod` | Rod beamstop position (mm) |
| `saxs_beamstop_x_pin`, `saxs_beamstop_y_pin` | Pin beamstop position (mm) |
| `energy_energy` | Photon energy (eV) |
| `beam_center_x`, `beam_center_y` | Beam center alternative names |
| `target_file_name` | Last-saved sample target filename (sometimes encodes geometry like `_ai0.5_th30`) |

Many of these have "_user_setpoint" companions; the loader's `_read_first_scalar` (`loader.py:800`) tries both forms.

## Encoded sample-name geometry

Some scans encode geometry into the `sample_name` start-doc field, e.g. `EG_AGB_16.10keV_wa14.5_sdd2.0m`. `parse_sample_name_geometry(name)` (`loader.py`) extracts:

- `energy_kev` — from a `<digits>.<digits>keV` token
- `sdd_m` — from `sdd<digits>.<digits>m`
- `wa<digits>` — WAXS arc reference angle (used by some calibration scans)
- `_ai<digits>` — incident angle (alpha_i, used for grazing incidence)
- `_th<digits>` — sample theta

This is a **lowest-priority** fallback in the metadata chain; baseline / primary always win.

## What the loader DataArrays carry as attrs

After `load_saxs_raw(uid)`:

```python
da.attrs == {
    "dist":       <SDD_m>,
    "poni1":      <bc_row_px * pixel1>,
    "poni2":      <bc_col_px * pixel2>,
    "rot1": 0, "rot2": 0, "rot3": 0,
    "pixel1":     PILATUS_PIXEL_SIZE_M,   # 0.172e-3
    "pixel2":     PILATUS_PIXEL_SIZE_M,
    "energy":     <eV>,
    "wavelength": <m>,
    # SMI extras:
    "smi_motor_x_mm":  <ref_value>,
    "smi_motor_y_mm":  <ref_value>,
    "smi_motor_z_mm":  <ref_value>,
    "smi_piezo_z_um":  <ref_value>,
    "smi_active_beamstop": "rod" | "pin",
    "smi_beamstop_pos_mm": {"rod": {"x": ..., "y": ...}, "pin": {"x": ..., "y": ...}},
}
```

Per-frame variation is exposed as **coords** on the `frame` (or `waxs_arc`) dim, not as attrs — attrs are scalars.

For WAXS, additional attrs include the panel layout (`smi_panel_offsets_deg`, `smi_panel_col_ranges`) and `smi_waxs_arc` (the per-frame arc array as a coord).

## Image data are lazy

`load_saxs_raw` and `load_waxs_raw` return `xr.DataArray` whose `.data` is a `dask.array` backed by tiled. Reading is deferred until `.values` or `.compute()` is called, or until `np.asarray(...)` is applied. The integrator pulls data in blocks of `IMAGE_BLOCK_FRAMES = 50` (`integrator.py`, near top) to bound peak memory.

`load_saxs_raw` and `load_waxs_raw` accept these args (`loader.py:1899-2046`):

```python
load_saxs_raw(run, energy_kev=None, **geometry_overrides) -> xr.DataArray
load_waxs_raw(run, energy_kev=None, **geometry_overrides) -> xr.DataArray
```

Geometry overrides flow through to `resolve_*_geometry(...)`.

## Disk-cache layout

`populate_cache(uid, run=None, include_images=True, path=None)` (`loader.py:423-480`) writes a per-UID HDF5 file with:

```
<uid>.h5
├── primary/
│   ├── pil2M_motor_x        (1-d array, n_frames)
│   ├── pil2M_motor_y
│   ├── ... (all primary scalar fields)
│   ├── pil2M_image          (3-d array, n × 1679 × 1475)
│   └── pil900KW_image       (3-d array, n × 195 × 619)
├── baseline/                (Dataset of all baseline scalars)
└── attrs                    (start metadata mirrored as group attrs)
```

Path resolution (`loader.py:231-258`):
1. `$SMI_BROWSER_CACHE_DIR/<uid>.h5` if env var set.
2. `$TMPDIR/smi_browser_cache/<uid>.h5` (or `/tmp/smi_browser_cache/<uid>.h5`).
3. Last fallback: `~/.cache/smi_browser/<uid>.h5`.

When you pass `image_cache_path=...` to `reduce_smi_combined`, that explicit path is used directly.

The cache is best-effort: write failures (read-only filesystem, disk full) are caught and ignored. The reduction pipeline doesn't fail because the cache write failed — it just refetches from tiled next time.

See `reference/caching.md` for the full cache architecture (4 cache layers).

## Discovery cheat sheet

```python
# What runs were taken on a sample, in a date range?
df = loader.searchCatalog(
    sample_name="AgB",
    since="2024-01-01", until="2024-12-31",
    detector="saxs",
)

# Just check what a single run looks like, no images:
md = loader.peekAtMd(uid)
# md is a dict; md["start"] / md["stop"] / md["streams"]

# Count frames per detector, plan name, sample:
summary = loader.summarizeRun(uid)

# Get the actual run object for low-level access:
run = loader._get_run(uid)
```

For deep diagnostics (motor values, baseline content), use the helpers in `reference/metadata-resolution.md`.
