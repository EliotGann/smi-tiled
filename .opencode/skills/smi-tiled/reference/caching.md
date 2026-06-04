# Caching architecture

`smi-tiled` uses **four independent cache layers**, each addressing a different bottleneck. They have different keys, lifetimes, and invalidation rules — getting them wrong causes silently-stale results.

| Layer | Where | Key | Stored | Lifetime |
|---|---|---|---|---|
| 1. Baseline cache | Process memory | run UID | `xr.Dataset` of all baseline columns | Process lifetime |
| 2. Geometry cache (SAXS) | Process memory | `(dist, poni1, poni2, pixel1, pixel2, wavelength, shape)` | `(q2d, qh2d, qv2d, chi_deg_2d, sa_base)` | Process lifetime |
| 3. Geometry cache (WAXS qmap) | Process memory | `(theta_deg, wavelength_nm)` rounded | `MultiPanelArcDetector` instance | Process lifetime |
| 4. Disk cache | HDF5 file `~/<uid>.h5` | UID-named file | All primary scalars + baseline + images | Until manually deleted |

## Layer 1: Baseline cache (`loader.py:38-39`)

```python
_BASELINE_CACHE: dict[str, xr.Dataset | None] = {}     # keyed by run UID
_BASELINE_COLUMNS_CACHE: dict[str, list[str]] = {}     # keyed by run UID
```

Why: parsing baseline takes seconds when the new layout's `xr.Dataset.from_dataframe` runs over 564 columns. Caching avoids that overhead on subsequent calls.

When populated:
- First call to `_read_baseline(run)` → cached in `_BASELINE_CACHE`.
- First call to `_baseline_scalar(run, key)` → column list cached in `_BASELINE_COLUMNS_CACHE` (so `key in columns` is fast).
- `_prepopulate_caches_from_h5(run, cache_path)` (`loader.py:581`) — pre-fills from disk cache.

When invalidated:
- `clear_baseline_cache()` (loader.py:50-51).
- Process exit.

> **GOTCHA**: `_BASELINE_CACHE[uid] = None` is a valid value (cached "this run has no baseline"). The check is `if _rid in _BASELINE_CACHE`, not `if _BASELINE_CACHE.get(_rid)` — so a `None` value is sticky.

## Layer 2: SAXS geometry cache (`integrator.py:169`)

```python
_SAXS_GEOMETRY_CACHE: dict[tuple, tuple] = {}
```

Key: `_saxs_cache_key(dist_m, poni1_m, poni2_m, pixel1_m, pixel2_m, wavelength_m, shape)` (`integrator.py:252`).

Stored: `(q2d, qh2d, qv2d, chi_deg_2d, sa_base)` — the per-pixel q maps and the un-toggled solid angle. The `sa` actually used is `sa_base if solid_angle_correction else None`, so the cache is independent of that toggle.

Memory footprint: `~5 × 1679 × 1475 × 8 bytes ≈ 95 MB per cached entry`. Two entries (e.g. one with motor scan, one without) ≈ 190 MB.

When populated:
- `integrate_saxs(...)` with `cache_geometry=True` (default).

When invalidated:
- `clear_geometry_cache()` (`integrator.py:172`).
- Manual: `_SAXS_GEOMETRY_CACHE.clear()`.
- **Not** invalidated when:
  - `solid_angle_correction` toggle changes (the cache stores `sa_base`, applied conditionally on read).
  - Mask changes.
  - `pixel_splitting` changes.

> **CRITICAL — stale cache after code edits**: If you modify the q-map computation (e.g. fix a sign convention), the cache hands back the old result. Geometry cache is keyed on **inputs**, not on code version. Call `clear_geometry_cache()` after edits.

## Layer 3: WAXS qmap cache

`_WAXS_QMAP_CACHE` (in `integrate_waxs`, `integrator.py:3200-3300`).

Key: `(theta_rounded_to_3dp, wavelength_rounded_to_15dp)` — rounding ensures floating-point drift doesn't fragment the cache.

Stored: `MultiPanelArcDetector` instance with its computed `qmap` dataset.

Why: at fixed wavelength, frames at the same arc angle share a qmap. For an N-frame scan stepping arc by 0.5°, that's 1 qmap per arc value (often a few unique values).

When invalidated: only by process exit. There is **no** explicit `clear_waxs_qmap_cache()` function — `clear_geometry_cache()` does NOT clear it (despite the name).

> **GOTCHA**: `clear_geometry_cache()` (`integrator.py:172`) only clears `_SAXS_GEOMETRY_CACHE`. The WAXS qmap cache and the split-bin-plan cache live elsewhere and persist.

## Layer 4: Disk cache (`loader.py:231-532`)

The most user-visible cache. HDF5 file at `<cache_dir>/<uid>.h5`:

```
<uid>.h5
├── primary/
│   ├── pil2M_motor_x        (n_frames,)
│   ├── pil2M_motor_y        (n_frames,)
│   ├── ... (all primary scalar fields: motors, signals)
│   ├── target_file_name     (n_frames,)  variable-length UTF-8
│   ├── pil2M_image          (n_frames, 1679, 1475)  chunked, gzip-compressed
│   └── pil900KW_image       (n_frames, 195, 619)
├── baseline/
│   └── (all baseline scalars, one Dataset per key)
└── attrs (mirrored start metadata)
```

### Path resolution

`_auto_cache_path(uid)` (`loader.py:231`):
1. `$SMI_BROWSER_CACHE_DIR/<uid>.h5` if env var set.
2. Otherwise `$TMPDIR/smi_browser_cache/<uid>.h5` (or `/tmp/smi_browser_cache/<uid>.h5` if no `$TMPDIR`).

`_cache_dir()` (`loader.py:248`) is the same logic but creates the directory if missing and returns the directory path.

User-supplied `image_cache_path` overrides both.

### Population

`populate_cache(uid, run, cache_path=None, include_images=True)` (`loader.py:423`):
1. Open HDF5 in append mode.
2. Read all primary scalar fields, write each as a 1-D dataset.
3. Read baseline as `xr.Dataset`, write each numeric column as a dataset.
4. (Optional) Read `pil2M_image` and `pil900KW_image` arrays from primary, write chunked + gzip-compressed (`compression_opts=2`, chunk shape `(1, ny, nx)` so individual frames are addressable).
5. Mirror start metadata into root group's attrs.

Failure modes are caught: `try ... except Exception: pass`. Failure to populate the cache **does not** prevent reduction; the next call just re-fetches from tiled.

### Read paths

| Function | Reads |
|---|---|
| `_read_cached_images(cache_path, field)` | Full image stack as ndarray |
| `_read_cached_images_lazy(cache_path, field)` | Returns dask array; HDF5 file kept open |
| `_read_cached_primary_field(cache_path, field)` | Single 1-D array (e.g. motor_z values) |
| `_read_cached_baseline(cache_path)` | All baseline as `dict[str, ndarray]` |
| `_read_cached_baseline_field(cache_path, field)` | Single baseline scalar |
| `_prepopulate_caches_from_h5(run, cache_path)` | Fills `_BASELINE_CACHE` from disk |

### When the disk cache is consulted

When you pass `image_cache_path=...` to `load_saxs_raw`, `load_waxs_raw`, `_read_scan_axis`, etc., the helpers try the cache first. If the file or the requested field is missing, they fall through to tiled. So a partial cache is harmless.

When `reduce_smi_combined` is called with `image_cache_path="auto"`:
- The auto-resolved path is consulted before tiled.
- `populate_disk_cache=True` (default) populates the cache after a successful tiled fetch.

### Invalidation

The disk cache has **no automatic invalidation**. To rebuild:
```bash
rm $TMPDIR/smi_browser_cache/<uid>.h5
# or
rm $SMI_BROWSER_CACHE_DIR/<uid>.h5
```

Or pass `image_cache_path=None` to bypass disk entirely.

> **GOTCHA — silent corruption**: If `populate_cache` is interrupted mid-write, the HDF5 file may be truncated. Subsequent reads will return partial data. The fix is `rm` and re-run.

## SplitBinPlan cache

There is also a **lazy SplitBinPlan cache** keyed on `(geometry_hash, q_edges_id, chi_edges_id, pixel_splitting)`. Currently inline-only (not module-global). Each call to `integrate_saxs` builds one SplitBinPlan and reuses it across all frames.

If the geometry varies per-frame (which is rare in transmission scans but common in z-grid scans), each unique geometry triggers a new plan. Memory cost: `~tens of MB` per plan (sparse matrix).

## Cache stats helper

`integrator.py:209-216`:
```python
{
    "saxs_entries":  len(_SAXS_GEOMETRY_CACHE),
    "saxs_est_mb":   <estimated MB>,
    ...
}
```

Useful for diagnosing memory pressure.

## When to clear caches

| Symptom | Likely fix |
|---|---|
| Edited `integrate_saxs` q-map code, results unchanged | `clear_geometry_cache()` |
| Edited mask JSON, mask still missing polygons | Pass `mask_path=...` explicitly; static masks aren't cached |
| Edited `populate_cache` image format | `rm <cache_dir>/<uid>.h5`, re-run |
| Edited `_read_baseline` logic | `_BASELINE_CACHE.clear(); _BASELINE_COLUMNS_CACHE.clear()` |
| Process memory growing per scan | `clear_geometry_cache()` after each scan; or batch scans by similar geometry |
| WAXS results stale despite geometry edit | Restart Python — there is no `clear_waxs_qmap_cache` |

## Recipe: full cache reset

```python
from smi_tiled.integrator import clear_geometry_cache
from smi_tiled.loader import _BASELINE_CACHE, _BASELINE_COLUMNS_CACHE

clear_geometry_cache()        # Layer 2
_BASELINE_CACHE.clear()       # Layer 1 part A
_BASELINE_COLUMNS_CACHE.clear()  # Layer 1 part B
# Layer 3 (WAXS qmap) requires process restart
# Layer 4 (disk): rm files manually
```

For Layer 3, the simplest workaround is to start a fresh Python process.
