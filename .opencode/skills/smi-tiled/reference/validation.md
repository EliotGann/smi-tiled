# Validation & sanity checks

How to verify a `smi-tiled` reduction is healthy without ground-truth comparison: inspect, run the bundled scripts, check sentinel values.

## Quick inspection (no image read)

### `loader.peekAtMd(uid, detector="saxs")` (`loader.py:2725-2744`)

Returns geometry metadata only — no images, no primary read:

```python
{
    "energy_kev":         16.1,
    "dist_m":             5.013,
    "beam_center_row_px": 1107.0,
    "beam_center_col_px": 743.5,
    "active_beamstop":    "rod",        # or "pin"
}
```

For `detector="waxs"`:
```python
{
    "energy_kev":         16.1,
    "dist_m":             0.273,
    "beam_center_row_px": 80.0,
    "beam_center_col_px": 310.0,
    "n_panels":           3,
}
```

Use this first when diagnosing a run: a sanity-check dict is one HTTP round trip.

### `loader.summarizeRun(uid)` (`loader.py:3034-3061`)

Headline metadata, no metadata-fallback chain:

```python
{
    "uid":              "c849dd2c-…",
    "scan_id":          1130041,
    "sample_name":      "AGB_scan_x_y_9m",
    "plan_name":        "count",
    "detectors":        ["pil2M_image", "pil900KW_image"],
    "detector_kinds":   ["saxs", "waxs"],     # via classify_detector_field
    "start_time":       1717023456.789,
    "num_points":       11,
}
```

Returns `start.get(...)` directly, so missing keys come back as `None` rather than an exception. Use this to confirm the right scan was selected before paying for a full reduction.

### `loader.searchCatalog(...)` (`loader.py:2860-3009`)

Discover runs without enumerating the catalog. Examples:

```python
df = loader.searchCatalog(scan_id=1130041)
df = loader.searchCatalog(sample="AgB", limit=10)
df = loader.searchCatalog(plan="grid_scan", since="2024-05-01")
df = loader.searchCatalog(detector="saxs", outputType="all")
```

Returns a `pandas.DataFrame`. `outputType="scans"` returns a 1-column DataFrame of scan_ids.

> **GOTCHA**: `searchCatalog` eagerly imports `Key, Regex, TimeRange` from `tiled.queries` (`loader.py:2927`). Some `tiled` versions miss `TimeRange`, breaking the call. Bypass:
> ```python
> from tiled.queries import Key
> node = loader._get_catalog().search(Key("scan_id") == 1130041)
> uid, run = next(iter(node.items()))
> ```

### `infer_detectors_and_steps(run)` (`loader.py:2339`)

The introspection used by `reduce_smi_combined` to discover scan dimensions. Avoids `.read()` on detector image fields (which can return multi-GB and trigger HTTP 500). Returns:

```python
{
    "n_frames": 33,
    "sample_name": "...",
    "scan_id": 1130041,
    "detectors": ["saxs", "waxs"],
    "step_candidates": [
        {"name": "pil2M_motor_x", "n_unique": 3, "min": 0.0, "max": 60.0,
         "values": np.array([0., 30., 60., 0., 30., 60., …])},
        {"name": "piezo_z",       "n_unique": 11, "min": -10000., "max": 10000.,
         "values": …},
    ],
}
```

`step_candidates` lists every 1-D primary field with > 1 unique value. Use this to find the scanned axes when you don't know them ahead of time.

## Sentinel values to verify

Calling a reduction "looks right" depends on knowing the expected ranges:

| Quantity | Healthy range (transmission) | Source |
|---|---|---|
| `geo.dist_m` (SAXS) | 1.7 – 9.3 m | `pil2M_motor_z` motor |
| `geo.dist_m` (WAXS) | 0.273 m (calibrated, `_DEFAULT_CAL`) | `integrator.py:3967-3968` |
| `geo.energy_ev` | 14000 – 17500 eV | beamline cycle |
| `geo.wavelength_m` | 7e-11 – 9e-11 m | E → λ via E·λ = hc |
| `geo.beam_center_row_px` (SAXS) | ~1107 ± 50 | varies with motor_y |
| `geo.beam_center_col_px` (SAXS) | ~744 ± 50 | varies with motor_x |
| `geo.active_beamstop` | `"rod"` or `"pin"` | from baseline |
| AGB ring 1 q | 1.076 nm⁻¹ | `D_AGB = 5.838 nm`, `q = 2π/D` |
| AGB ring 1 in I(q) plot | sharp peak at q = 1.076 nm⁻¹ | scripts/calibrate_smi_saxs.py:11 |
| Pin disk masked fraction | ≥ 95% in disk of r=22 px around BC + (5, 5) | gotcha #2 in SKILL.md |
| Solid-angle correction range | 1.0 – 1.5 (faster falloff at large 2θ) | `integrator.py:1465-1480` |
| Number of valid pixels per scan | 1679×1475 - mask area ≈ 2.0 × 10⁶ for SAXS | `_smi_run_raw_shape` |
| Per-frame I(q) `counts` at small q | < 100 (small annulus) | bin geometry |
| Per-frame I(q) `counts` at q=1 nm⁻¹ | 10⁴ – 10⁵ | bin geometry |

If any of these are off by orders of magnitude, the reduction has a problem.

### Default loader calibration deltas

`LOADER_DEFAULTS` (`defaults.py:234`):
```python
LoaderCalibration(
    saxs_row_delta_px=…,
    saxs_col_delta_px=…,
    waxs_row_delta_px=…,
    waxs_col_delta_px=…,
    saxs_distance_delta_mm=…,
)
```

These are the *static offsets* applied to the live EPICS PVs. If a calibration script regresses the live PV against a motor, its slope `b_x` is **not** something to add — the PV already moves; doing so double-counts. See gotcha #1 in `gotchas.md`.

### `BSX_PER_ARC_DEG` (`defaults.py:195`)

`-4.39 mm/deg` — mechanical linkage between `waxs_arc` and `waxs_bsx`. When `waxs_arc` rotates by Δθ, `waxs_bsx` moves by `-4.39 × Δθ` mm because they share a kinematic mount. Used in WAXS qmap construction to compute the per-frame beam center column (`waxs-geometry.md`).

If your WAXS reduction shows the beamstop at the wrong column for non-zero `waxs_arc`, this constant is the first place to check.

## Bundled scripts

All under `scripts/`. Run via `pixi run python scripts/<name>.py [args]`.

### `validate_large_scan.py` (123 lines)

Stress test: runs `reduce_smi_combined` on a known large scan with progress callbacks, prints peak memory (via `tracemalloc`), reports timing per stage and progress-bar granularity.

```bash
pixi run python scripts/validate_large_scan.py
```

Hard-coded UID `ce94c000-369d-444a-8078-a9ed3c36b872`. Reports:
- Number of progress callbacks per stage
- Min/max/mean callback intervals
- Per-stage timing (load, mask, saxs, waxs, merge)
- Peak memory usage
- Output dataset shapes

Use after large refactors to confirm no perf regression.

### `calibrate_smi_saxs.py` (457 lines)

AgB grid-scan calibration of SAXS beam center and SDD as functions of `pil2M_motor_x/y/z` and `piezo_z`. Assumes silver behenate (D = 5.838 nm) is in the beam. Steps:

1. **Find bright spot**: locate largest connected bright region (filters hot pixels), find its peak, take ±15 px intensity-weighted centroid.
2. **Radial profile**: chi-averaged radial profile around that point.
3. **Find ring 1**: peak in the radial profile gives ring radius (px).
4. **Refine BC**: sample the ring at multiple chi angles, least-squares fit a circle.
5. **Bragg → SDD**: `tan(2θ_AGB) = r_mm / SDD_mm`.

Then regress:
```
beam_col_px = a_x + b_x · motor_x_mm
beam_row_px = a_y + b_y · motor_y_mm
SDD_mm      = a_z + b_z · piezo_z_um
```

> **CRITICAL — see warning at `scripts/calibrate_smi_saxs.py:36-44`**: The fitted `b_x ≈ 5.82 px/mm` is how the live EPICS BC PV moves with motor_x — the PV already tracks the detector. Pasting it into `_SAXS_BEAM_COL_PX_PER_MOTOR_X_MM` double-counts and throws the BC ~120 px off. Default = 0.0 and **stay 0.0** unless a future detector reports a static, position-independent BC PV.

### `calibrate_smi_z_scan.py` (657 lines)

Multi-dimensional calibration: BC, SDD, and pin shadow centroid as functions of `pil2M_motor_z` and `piezo_z`. Outputs:

- `/tmp/agb_z_calibration_results.npz` — per-frame measurements
- `/tmp/saxs_calibration.json` — JSON override file the loader can read
- `/tmp/smi_beamstop_offsets.json` — `motor_z` → `(d_bsx, d_bsy)` table for the data-collection code to keep the pin beamstop centered on the beam

Reference scan: `b900e711-…` (`AGB_scan_z`), 79 z-values × 11 piezo-z values. Active beamstop: pin.

### `benchmark_gpu_histogram.py` (184 lines)

Compares CPU `_SplitBinPlan` vs GPU `TorchBinPlan` on synthetic data mimicking pil2M geometry. Verifies correctness (must match within float tolerance) and reports throughput (frames/sec). Requires PyTorch + GPU (CUDA or MPS).

Use this before promoting GPU integration to default.

## AgB ring sanity check

The fastest "is the geometry right?" test:

```python
import numpy as np
result = reduce_smi_combined(uid, …)
da = result.to_dataarray("merged_iq")
q = np.asarray(da.q)
I = np.asarray(da.values)
finite = np.isfinite(I)
peak_q = q[finite][np.argmax(I[finite][(q[finite] > 0.9) & (q[finite] < 1.2)])]
print(f"AgB ring 1 at q = {peak_q:.3f} nm⁻¹  (expected 1.076)")
```

If the run had AgB:
- Off by < 1% → calibration good.
- Off by 5-10% → SDD wrong; check `geo.dist_m`.
- Off by > 20% or no peak → BC wrong, mask wrong, or no AgB.

Same idea works for WAXS using a known-d-spacing standard (LaB6, Si).

## Healthy reduction signature

A successful `reduce_smi_combined` should produce all of:
- `result.timing["total"] < 60s` for typical 11-frame scans (assuming cache populated).
- `result.merged_iq["I"]` has finite values for a contiguous q range.
- `result.merged_iq["counts"]` is monotone-increasing from low q to mid q (more bins fit in a thicker annulus), then drops as you exit the detector.
- `result.merged_qchi["intensity"]` shows roughly uniform intensity in chi at fixed q (apart from anisotropy of the sample).
- `result.saxs["iq"]` and `result.waxs["iq"]` overlap reasonably in q ≈ 1-2 nm⁻¹ where their q-ranges meet.

If `merged_qchi.saxs_intensity` and `merged_qchi.waxs_intensity` disagree by > 2× in the overlap region, suspect:
- WAXS panel rotation wrong (`MultiPanelArcDetector`).
- WAXS BC col wrong (`_BSX_PER_ARC_DEG` mismatch).
- Solid-angle correction toggle differs between detectors (it's applied via the same `solid_angle_correction` flag — not per-detector).

## Mask coverage check

```python
from smi_tiled import resolve_saxs_geometry, make_saxs_mask_from_spec
from smi_tiled.defaults import resolve_mask_path

run = loader._get_run(uid)
geo = resolve_saxs_geometry(run)
mask = make_saxs_mask_from_spec(
    image_shape=(1679, 1475),
    mask_path=resolve_mask_path(None, "saxs"),
    active_beamstop=geo.active_beamstop,
    beamstop_pos_mm=geo.beamstop_pos_mm,
    beam_center_px=(geo.beam_center_row_px, geo.beam_center_col_px),
)
print(f"valid fraction: {mask.mean():.2%}")  # expect ~95% for rod, ~98% for pin
```

If valid fraction < 80%, the mask is over-masking — likely wrong polygon vertices in the JSON or wrong beam center.

## When in doubt, dump intermediates

The `result.saxs["q_chi"]` and `result.waxs["q_chi"]` Datasets are tiny — saving them as NetCDF or pickle for later comparison is cheap:

```python
result.saxs["q_chi"].to_netcdf(f"/tmp/saxs_qchi_{uid[:8]}.nc")
result.waxs["q_chi"].to_netcdf(f"/tmp/waxs_qchi_{uid[:8]}.nc")
```

Compare against a known-good run's output. Diffs in `intensity` / `counts` at the same q-chi pixel point at: mask change, calibration change, or solid-angle toggle change.

## Healthy detector-space frame

If you set `build_detector_ds=True` (note: auto-disabled for > 50 frames; pass explicit `True` to override):

```python
ds = result.saxs["ds"]
img = ds["intensity"].isel(frame=0).values
mask = ds["mask"].isel(frame=0).values
print(f"raw max:    {img.max():.0e}")          # ~1e6 typical
print(f"raw min:    {img.min()}")               # 0 or -1 (dead pixels)
print(f"mask sum:   {mask.sum()}")              # ~2e6 valid pixels
print(f"q range:    {ds['q_abs'].min().values:.3f} – {ds['q_abs'].max().values:.3f} nm⁻¹")
```

Pilatus pixel saturation is `~1e6 cps` (counts per second integrated over exposure). Counts above that mean dead time / overflow.

## Summary: validation toolkit

| Step | Command | Detects |
|---|---|---|
| 1. Headlines | `loader.summarizeRun(uid)` | Wrong scan, missing detectors. |
| 2. Geometry | `loader.peekAtMd(uid)` | Wrong BC, wrong SDD, wrong beamstop. |
| 3. Steps | `infer_detectors_and_steps(run)` | Wrong scan dimensions, missing motor. |
| 4. AgB ring | Reduction + peak finder at q ≈ 1.076 | Wrong calibration. |
| 5. Mask coverage | Manual mask + `valid_frac` | Bad polygons, wrong BC. |
| 6. Detector-space | `build_detector_ds=True` | Hot pixels, dezinger problems, panel-edge artifacts. |
| 7. Stress test | `validate_large_scan.py` | Memory leaks, perf regressions. |
