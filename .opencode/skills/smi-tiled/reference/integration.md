# Integration pipeline

This document covers the per-frame reduction from `(image, mask, qmap)` → `(q, chi)` histogram → `I(q)`. It also covers cross-detector merging (SAXS + WAXS) and pixel splitting.

## Top-level entry: `reduce_smi_combined`

`integrator.py:3523-3700+`. The end-to-end pipeline:

```python
result = reduce_smi_combined(
    uid: str,
    tiled_uri: str = "https://tiled.nsls2.bnl.gov",
    catalog: str = "smi/migration",
    n_q: int = 2000,
    n_chi: int = 360,
    solid_angle_correction: bool = True,
    saxs_mask: str | Path | dict | None = None,
    waxs_mask: str | Path | dict | None = None,
    saxs_kwargs: dict | None = None,
    waxs_kwargs: dict | None = None,
    geometry: str = "transmission",        # or "grazing_incidence"
    incident_angle_deg: float = 0.0,
    saxs_beam_delta_px: tuple | None = None,
    waxs_beam_delta_px: tuple | None = None,
    saxs_distance_delta_mm: float | None = None,
    saxs_q_cutoff: float | None = None,
    saxs_agbh_ring_order: int = 5,
    saxs_q_margin_fraction: float = 0.01,
    dezinger_threshold: float | None = 3000.0,
    dezinger_kernel: int = 5,
    waxs_beam_col_per_arc_deg: float = 0.08,
    cache_geometry: bool = True,
    pixel_splitting: int = 1,
    build_detector_ds: bool | None = None,
    build_frame_qchi: bool = True,
    frame_qchi_store: str | Path | None = "auto",   # zarr path for per-frame
    image_cache_path: str | Path | None = "auto",   # HDF5 path
    populate_disk_cache: bool = True,
    progress: ProgressCallback | None = None,
    # Derived stages
    virtual_axes=None,
    line_cuts=None,
    peak_fits=None,
) -> CombinedReductionResult
```

Stages:
1. **Connect to tiled** and load the run.
2. **Resolve geometry** for both detectors.
3. **Load lazy DataArrays** (`load_saxs_raw`, `load_waxs_raw`).
4. **Populate disk cache** if `populate_disk_cache=True` and `image_cache_path` is set.
5. **Integrate SAXS** via `integrate_saxs` (skips if no SAXS frames).
6. **Integrate WAXS** via `integrate_waxs` (skips if no WAXS frames).
7. **Merge** to common `(q, chi)` grid via `merge_q_chi_weighted`.
8. **Compute merged I(q)** via `merge_iq_profiles`.
9. **Build per-frame I(q)** via `_build_per_frame_iq`.
10. **Optional derived stages** (line cuts, peak fits, virtual axes).
11. Return `CombinedReductionResult`.

## SAXS integration: `integrate_saxs`

`integrator.py:2884-3200`. Direct pixel-space q-map and histogram binning. Steps:

### 1. Build geometry (`integrator.py:2913-2978`)

Read `dist, poni1, poni2, pixel1, pixel2, wavelength` from `da.attrs` (units: `wavelength` is in **Å**, the rest in **m**; the integrator multiplies by `1e-10` then `1e9` to get nm).

Compute `q2d, qh2d, qv2d, chi_deg_2d, sa_base` — see `reference/saxs-geometry.md` for the math.

Result keyed by `_saxs_cache_key(...)` and stored in `_SAXS_GEOMETRY_CACHE` if `cache_geometry=True`.

### 2. Apply masks (`integrator.py:2984-3050`)

```python
base_valid = np.isfinite(q2d) & np.isfinite(chi_deg_2d)
if mask_use is not None:
    base_valid &= mask_use
# Per-frame dynamic mask combines shadow + aperture (see reference/masks.md)
```

### 3. Build q/chi grids and SplitBinPlan

```python
q_edges   = np.linspace(q_min, q_max, n_q + 1)
chi_edges = np.linspace(chi_min, chi_max, n_chi + 1)
plan = _SplitBinPlan(q2d, chi_deg_2d, q_edges, chi_edges, pixel_splitting)
```

The `_SplitBinPlan` (`integrator.py:1487-1581`) precomputes a sparse `(n_bins, n_pixels)` matrix `M` so that `I_hist = M @ (img.ravel() * valid.ravel())`. For pixel_splitting > 1, multiple sub-pixel-shifted matrices are accumulated. See "Pixel splitting" below.

### 4. Integrate frames in blocks

`IMAGE_BLOCK_FRAMES = 32` (`loader.py:73`). For each block:
1. Compute `imgs[s:e] = np.asarray(da.data[s:e], dtype=float)` (lazy dask compute).
2. Apply solid-angle correction: `imgs /= sa[None, :, :]` if enabled.
3. Optionally dezinger: `_dezinger_per_frame(imgs, kernel, threshold)`.
4. For each frame in block:
   - `valid_f = base_valid & dynamic_mask[i] & np.isfinite(imgs[i])`
   - `I_hist, N_hist = plan.integrate_frame(imgs[i], valid_f)`
   - Accumulate into `accum_I, accum_N`.
   - If `frame_qchi_store` is set, write `(I_hist, N_hist)` to zarr.

### 5. Reduce to outputs

```python
mean_I = np.where(accum_N > 0, accum_I / accum_N, np.nan)
qchi = xr.Dataset({"intensity": ..., "counts": ...}, coords={"q": ..., "chi": ...})
total_N = accum_N.sum(axis=1)
I_1d = np.where(total_N > 0, (mean_I_no_nan * accum_N).sum(axis=1) / total_N, np.nan)
iq = xr.Dataset({"I": ..., "counts": ...}, coords={"q": ...})
```

This is implemented in `_qchi_and_iq` (`integrator.py:1583-1609`).

### Returns

```python
{
    "q_chi":        xr.Dataset (q, chi) — accumulated mean intensity + counts
    "iq":           xr.Dataset (q,)     — azimuthally averaged
    "qchi_frames":  xr.Dataset (frame, q, chi) | lazy zarr | None
    "iq_frames":    xr.Dataset (frame, q)
    "qmap":         xr.Dataset (row, col) — q2d, qh, qv, chi_deg, valid_mask
}
```

## WAXS integration: `integrate_waxs`

`integrator.py:3200-3500`. Per-frame `MultiPanelArcDetector` because the panel geometry depends on `waxs_arc`.

### Per-frame geometry building

For each frame `f`:
1. Read `θ_f = waxs_arc[f]` and `λ_f = wavelength_per_frame[f]`.
2. Look up cached `MultiPanelArcDetector` keyed on `(θ_rounded, λ_rounded)`. Two frames sharing `(θ, λ)` reuse the same detector.
3. Call `detector.qmap(θ_f)` to get `(qx, qy, qz, qabs, solid_angle)`.
4. Compute per-frame mask via `make_mask_for_angle(θ_f, ...)`.
5. Run `_qchi_and_iq` on the masked image.

### Returns

Same shape as SAXS:
```python
{
    "q_chi":  xr.Dataset (q, chi),
    "iq":     xr.Dataset (q,),
    "qchi_frames": xr.Dataset (frame, q, chi) | None,
    "iq_frames":   xr.Dataset (frame, q),
    "qmap":   xr.Dataset (row, col),
}
```

## Pixel splitting (`integrator.py:1487-1581`)

When `pixel_splitting > 1`, each pixel is subdivided into an `n × n` grid. The `(q, chi)` of each sub-pixel is interpolated using gradients:

```python
dq_dr = np.gradient(q2d, axis=0)
dq_dc = np.gradient(q2d, axis=1)
dchi_dr = np.gradient(chi2d, axis=0)
dchi_dc = np.gradient(chi2d, axis=1)
offsets = np.linspace(-0.5 + 0.5/n, 0.5 - 0.5/n, n)   # sub-pixel centers
for dr in offsets:
    for dc in offsets:
        q_sub = q2d + dr * dq_dr + dc * dq_dc
        chi_sub = chi2d + dr * dchi_dr + dc * dchi_dc
```

Each sub-pixel gets weight `1 / (n*n)`. Per-frame integration becomes a sum over `n²` sparse matrix-vector products.

`pixel_splitting=1` (default) skips this and uses single-point binning (`q_sub = q2d`).

## Solid-angle correction

Applied **before** binning (`integrator.py` near 3050):
```python
imgs /= sa[None, :, :]    # corrects for cos(2θ) angular factor
```

This is critical for absolute intensity work — without it, peak heights vs q are systematically biased low at high q.

The correction is `pixel_area · max(dist, 0) / r³` (see `reference/saxs-geometry.md`). When `solid_angle_correction=False`, divide by 1 (no-op).

## Histogram + binning math

For each frame:
```
For each pixel (r, c):
    q   = q2d[r, c]
    chi = chi2d[r, c]
    q_bin   = searchsorted(q_edges, q)   - 1
    chi_bin = searchsorted(chi_edges, chi) - 1
    if valid[r, c] and 0 <= q_bin < n_q and 0 <= chi_bin < n_chi:
        I_hist[q_bin, chi_bin] += img[r, c]
        N_hist[q_bin, chi_bin] += 1
```

The `_SplitBinPlan` precomputes the bin assignments **once** (per geometry), so per-frame integration is just a sparse mat-vec — no `searchsorted` per frame.

## Cross-detector merge: `merge_q_chi_weighted`

`integrator.py:1742-1813`. Inputs are two q-chi datasets (SAXS, WAXS), each on its own grid.

Steps:
1. Build common `q_grid = linspace(q_min_overall, q_max_overall, n_q)`.
2. Build common `chi_grid = linspace(chi_min_overall, chi_max_overall, n_chi)`.
3. **Regrid** both inputs onto the common grid via `RegularGridInterpolator(method="nearest")`. Counts use `fill_value=0`; intensities use `fill_value=nan`.
4. **Count-weighted average**:
   ```python
   total_N = saxs_N + waxs_N
   merged_I = (saxs_I * saxs_N + waxs_I * waxs_N) / total_N
   ```
   where `nan` intensities are replaced with 0 in the numerator (so that bins with only one detector's data still get that detector's intensity).
5. Return Dataset with vars `intensity, counts, saxs_intensity, saxs_counts, waxs_intensity, waxs_counts`.

> **NOTE — nearest-neighbor regridding**: This is a lossy operation. For high-precision work consider raising `n_q` or interpolating linearly. Default `n_q=2000` is dense enough that nearest-neighbor is rarely a problem.

> **NOTE — chi grids may not overlap**: SAXS chi spans roughly ±180°; WAXS chi spans a much narrower band (because WAXS is below the equator). The merged chi range is the union of both, so most chi-bins will only have data from one detector.

## Merged I(q): `merge_iq_profiles`

`integrator.py:1816-1858`. Two-step:
1. Re-azimuthally-average the **merged** q-chi map (count-weighted).
2. Also interpolate the original SAXS-only and WAXS-only I(q) onto the merged q grid (carried as `saxs_I` and `waxs_I` data vars).

Result variables: `I, counts, saxs_I, waxs_I`.

Note: `merge_iq_profiles` does **not** itself do count-weighted merging — that comes from step 1 (the q-chi merge). The `saxs_I` and `waxs_I` carried alongside are unmerged interpolated views.

## Per-frame I(q): `_build_per_frame_iq`

`integrator.py:1861-1971`. For per-frame analysis (each frame's profile separately):
1. Take `saxs_iq_frames` and `waxs_iq_frames` (each `(frame, q)`).
2. Interpolate each onto the merged q grid.
3. Count-weighted merge per frame.
4. Attach any per-frame primary-stream scalars (motor positions, etc.) from `scan_info["step_candidates"]` as data vars on the `frame` dim.

Returned shape: `(frame, q)` Dataset with vars `I, saxs_I, waxs_I` plus per-frame motor positions etc.

## Detector-stack DataSet

`build_detector_ds=True` returns a `(frame, row, col)` xarray Dataset with the **raw** detector frames (post-mask, pre-bin) on disk. Useful for QC and diagnosis but expensive — for large scans (`n_frames > 50`) the auto-default switches to `False`. See `_DETECTOR_DS_AUTO_MAX_FRAMES = 50` (`integrator.py:2785`).

> **GOTCHA**: For 51+ frame scans, the auto-default silently drops the detector DS. If you need it, pass `build_detector_ds=True` explicitly. Memory cost is ~16 GB for 800 SAXS frames.

## Per-frame qchi store

`build_frame_qchi=True` and `frame_qchi_store="auto"` (or a path) writes per-frame `(q, chi)` data to a zarr store on disk:
- Layout: chunks `(1, n_q, n_chi)` so each frame is its own chunk.
- Backed by `_ZarrFrameWriter` (`integrator.py:1653-1695`).
- Returned via `result.saxs["qchi_frames"]` as a lazy dask-backed Dataset.
- Consumers materialize one frame at a time via indexing.

For "auto", the store goes to `{TMPDIR or /tmp}/smi_qchi_{uid}.zarr/{saxs|waxs}`.

## Geometry caches

| Cache | Key | Stores |
|---|---|---|
| `_SAXS_GEOMETRY_CACHE` | `_saxs_cache_key(dist, poni1, poni2, pixel1, pixel2, wavelength, shape)` | `(q2d, qh2d, qv2d, chi_deg_2d, sa_base)` |
| `_WAXS_QMAP_CACHE` | `(theta_rounded_to_3dp, wavelength_rounded_to_15dp)` | `MultiPanelArcDetector + qmap` |
| `_SPLIT_PLAN_CACHE` | hash of geometry + grid + pixel_splitting | `_SplitBinPlan` instance |
| `_BASELINE_CACHE`, `_BASELINE_COLUMNS_CACHE` | run UID | parsed baseline content |

`clear_geometry_cache()` (`integrator.py:172`) wipes all of these.

> **WARNING — stale cache after edits**: If you modify geometry build code without bumping any input field, the cache hands back the old (stale) result. Call `clear_geometry_cache()` after edits.

## Performance notes

| Bottleneck | Mitigation |
|---|---|
| Per-frame `searchsorted` in histograms | `_SplitBinPlan` precomputes bin indices |
| Repeated dask compute on already-loaded data | Block-level `np.asarray` materialization in `IMAGE_BLOCK_FRAMES = 32` chunks |
| Holding `(n_frames, ny, nx)` mask | `_LargeAreaMaskLazy` for n_frames > 50 |
| Holding per-frame qchi arrays | `_ZarrFrameWriter` streaming |
| Re-reading from tiled in repeated runs | `populate_cache` writes HDF5 disk cache |

## Key file locations

| Function | Location |
|---|---|
| `polygons_to_mask`, `shift_polygon`, `make_mask_for_angle` | `integrator.py:546-591` |
| `make_saxs_mask_from_dict`, `make_saxs_mask_from_spec` | `integrator.py:612-790` |
| `_make_waxs_shadow_mask`, `_LargeAreaMaskLazy` | `integrator.py:1219-1287` |
| `_make_aperture_mask`, `make_saxs_large_area_masks`, `_LargeAreaMaskCombined` | `integrator.py:1290-1370` |
| `_histogram2d_pixel_split` | `integrator.py:1378-1485` |
| `_SplitBinPlan` | `integrator.py:1487-1581` |
| `_qchi_and_iq` | `integrator.py:1583-1609` |
| `_stack_qchi_frames`, `_stack_iq_frames` | `integrator.py:1612-1650` |
| `_ZarrFrameWriter` | `integrator.py:1653-1698` |
| `merge_q_chi_weighted` | `integrator.py:1742-1813` |
| `merge_iq_profiles` | `integrator.py:1816-1858` |
| `_build_per_frame_iq` | `integrator.py:1861-1971` |
| `merge_multiple_qchi`, `merge_multiple_iq` | `integrator.py:1978-2080` |
| `CombinedReductionResult.to_dataarray` | `integrator.py:2086-2200` |
| `integrate_saxs` | `integrator.py:2884-3200` |
| `integrate_waxs` | `integrator.py:3200-3520` |
| `reduce_smi_combined` | `integrator.py:3523-onwards` |
