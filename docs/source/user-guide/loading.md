# Loading raw images

Raw SMI data lives in a Tiled catalog (production at
`https://tiled.nsls2.bnl.gov`, in `smi/migration`).  Each scan is a
*run* identified by a UID.  Two detectors may be present:

- **Pilatus 2M** — SAXS, ~1679 × 1475 pixels, single flat panel.
- **Pilatus 900KW** — WAXS, ~619 × 1475 pixels split into three folded
  panels (offsets typically −7°, 0°, +7°).

This page covers how {class}`smi_tiled.TiledSMISWAXSLoader` resolves
geometry, reads images, and exposes per-frame motor data.

## The loader

```python
from smi_tiled import TiledSMISWAXSLoader

loader = TiledSMISWAXSLoader(
    tiled_uri="https://tiled.nsls2.bnl.gov",   # default
    catalog="smi/migration",                   # default
    energy_kev=None,                           # auto-resolve from run
    api_key=None,                              # or pass / use TILED_API_KEY env
)
```

`TiledSMISWAXSLoader` is **not a `FileLoader` subclass** — SMI's data
source is a remote HTTP catalog indexed by UID, not a directory of
files indexed by path.  The closest pattern in the broader Python
scattering ecosystem is `PyHyperScattering.SST1RSoXSDB`.

## Loading one detector

```python
saxs = loader.loadSingleImage(uid="6e61b977-…", detector="saxs")
print(saxs)
# <xarray.DataArray (frame: 14, pix_y: 1679, pix_x: 1475)>
# Coordinates: frame, plus per-frame motor positions
# Attributes: dist, poni1, poni2, pixel1, pixel2, energy, wavelength,
#             smi_detector, smi_active_beamstop, uid, scan_id, sample_name,
#             ...
```

For a single-frame scan, the leading axis is squeezed away
(`['pix_y', 'pix_x']` only).  When the scan axis is `waxs_arc`, the
leading dim is named `waxs_arc` rather than `frame`.

## Loading both detectors

```python
both = loader.loadRun(uid="…")
saxs_raw = both["saxs"]
waxs_raw = both["waxs"]
```

`None` is returned for detectors absent from the run.

## Geometry resolution

The loader's
{func}`smi_tiled.resolve_saxs_geometry` /
{func}`smi_tiled.resolve_waxs_geometry` use a strict fallback chain:

```
1. User override (kwarg, e.g. beam_center_row_px=…)
2. Primary stream  (per-frame, when motor is scanned)
3. Baseline stream (start-of-scan snapshot — always present)
4. Primary configuration metadata
5. Start metadata / sample_name encoding (e.g. _wa20.0_, _sdd9m, _16.10keV)
6. Hardcoded instrument defaults
```

The "trust motors over metadata" rule applies: when a motor is being
*scanned* in the primary stream, its first-frame value (not the
baseline snapshot, which is from before the scan started) is the
reference.

## Motor-driven corrections

The default beam center and SDD come from the EPICS baseline PVs.  The
beam-center PVs (`pil2M_beam_center_x/y_px`) **already track the detector
translation motors** (`pil2M_motor_x/y`), so by default the loader applies
**no** additional motor_x/y correction — the `_SAXS_BEAM_*_PX_PER_MOTOR_X/Y_MM`
slopes default to `0.0`.  Only the small motor_z drift and the `piezo_z`→SDD
effect remain active (calibrated from AGB scans — see {doc}`calibration`):

```python
beam_col += (motor_x_mm - _SAXS_MOTOR_X_REF_MM) * _SAXS_BEAM_COL_PX_PER_MOTOR_X_MM   # 0.0 by default
beam_row += (motor_y_mm - _SAXS_MOTOR_Y_REF_MM) * _SAXS_BEAM_ROW_PX_PER_MOTOR_Y_MM   # 0.0 by default
beam_col += (motor_z_mm - _SAXS_MOTOR_Z_REF_MM) * _SAXS_BEAM_COL_PX_PER_MOTOR_Z_MM
beam_row += (motor_z_mm - _SAXS_MOTOR_Z_REF_MM) * _SAXS_BEAM_ROW_PX_PER_MOTOR_Z_MM
dist_mm  += (piezo_z_um - _SAXS_PIEZO_Z_REF_UM) * _SAXS_SDD_DELTA_MM_PER_PIEZO_Z_UM
dist_mm  += _SAXS_DEFAULT_DISTANCE_DELTA_MM
```

```{warning}
A non-zero motor_x/y slope **double-counts** the detector translation (the PV
already encodes it) and pushes the beam center ~120 px off — which also
misplaces the beamstop sub-mask.  Keep these slopes at 0.0 unless a future
detector reports a static, position-independent beam-center PV.
```

These constants live in `smi_tiled.loader` and can be overridden
**at runtime** (`reduce_smi_combined(saxs_beam_delta_px=(d_row, d_col), …)`)
or **at import time** via `src/smi_tiled/data/saxs_calibration.json`
(see {doc}`calibration`).

## Browsing the catalog

```python
df = loader.searchCatalog(
    sample="AGB",            # case-insensitive regex on start.sample_name
    plan="scan",
    since="2026-05-01",
    until="2026-05-31",
    detector="saxs",         # filter on start.detectors
    limit=50,
)
```

Returns a `pandas.DataFrame` keyed by `scan_id`, `sample_name`,
`plan_name`, `detectors`, `num_points`, `uid`, etc.  Pass
`outputType="all"` to include `cycle`, `user_name`, `institution`,
`proposal_id`.

For a Jupyter interactive grid:

```python
loader.browseCatalog(sample="AGB")   # requires `pip install ipyaggrid`
```

## HDF5 cache (optional)

Tiled round-trips are slow; the loader can read from a pre-populated
HDF5 cache:

```python
result = reduce_smi_combined(
    uid=uid,
    image_cache_path="/tmp/smi_browser_cache/" + uid + ".h5",
)
```

The cache layout matches what `smi-browser` (or
{func}`smi_tiled.populate_cache`) writes:

```
/primary/<field>    — 1-D arrays of per-frame scalars
/baseline/<field>   — 1-D arrays (length 1–2)
/images/<field>     — 3-D detector stacks (N, H, W)
```

Missing fields transparently fall back to Tiled.  Set
`image_cache_path="auto"` (the default in `reduce_smi_combined`) to
auto-detect the cache file in `$SMI_BROWSER_CACHE_DIR` or
`$TMPDIR/smi_browser_cache/`.

## Robustness: chunked-read fallback

The Tiled server occasionally returns HTTP 500 for bulk multi-frame
detector reads (the SMI `pil2M_image` field at ~99 MB per chunk is
prone to this).  The loader detects oversize chunks and pre-emptively
reads frame-by-frame in parallel:

```python
# Internal: _BULK_READ_MAX_CHUNK_BYTES = 64 MiB
# When a node's chunks exceed this, the loader skips the bulk read
# and uses a ThreadPoolExecutor to read N frames in parallel,
# retrying transient failures.
```

This is automatic — no user action needed.

## Inspecting metadata without reading the image

```python
md = loader.peekAtMd(uid="…", detector="saxs")
print(md["dist_m"], md["beam_center_row_px"], md["beam_center_col_px"])
```

`peekAtMd` resolves geometry without downloading the image stack.

## See also

- {doc}`reduction` — turning raw images into `I(q)`
- {doc}`calibration` — re-fitting the geometry constants
- {class}`smi_tiled.TiledSMISWAXSLoader` (API reference)
