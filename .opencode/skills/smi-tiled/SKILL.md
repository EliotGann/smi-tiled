---
name: smi-tiled
description: Use when working with the smi-tiled package, SMI beamline (NSLS-II 12-ID) data, Pilatus 2M (pil2M / SAXS) or Pilatus 900KW (pil900KW / WAXS arc) detectors, beamstop masks (rod / pin), q-chi reduction, reduce_smi_combined, TiledSMISWAXSLoader, scan_id lookups against tiled.nsls2.bnl.gov / smi/migration, or debugging missing or incorrect masking, geometry, or calibration. Covers the metadata fallback chain (override → baseline → primary → defaults), polygon mask formats, the q/chi math used by integrate_saxs and integrate_waxs, the SAXS+WAXS merge, the 4 cache layers, and the "where to check first" decision tree for diagnosing reductions.
---

# smi-tiled

The `smi-tiled` package is a Tiled-native loader, integrator, and uploader for SMI beamline data (NSLS-II 12-ID). This skill captures the operational knowledge needed to read its outputs correctly, diagnose problems, and extend it safely. Source of truth is the code under `src/smi_tiled/`; this skill points at canonical lines.

## Trigger

Use this skill when:
- A user message mentions any of: `smi-tiled`, `SMI`, `12-ID`, `pil2M`, `pil900KW`, `Pilatus 2M`, `Pilatus 900KW`, `beamstop`, `rod` / `pin`, `tiled.nsls2.bnl.gov`, `smi/migration`, `reduce_smi_combined`, `TiledSMISWAXSLoader`, `SAXSGeometry`, `WAXSCalibration`, `MultiPanelArcDetector`, `waxs_arc`, `bsx`, `merged_qchi`, `merged_iq`, `AgB` / `silver behenate` calibration of SAXS, or scan IDs in the smi catalog.
- The task is reading, editing, or debugging code in `/nsls2/users/egann/git/smi/smi-tiled` (or a checkout of it).
- The user observes "missing masking", "wrong q range", "wrong beam center", "missing pin / rod beamstop", or geometry/intensity that disagrees with expectations.

Do NOT use this skill for unrelated x-ray analysis, generic pyFAI work, or PyHyperScattering questions that don't touch SMI data.

## Repo orientation

Single Python package under `src/smi_tiled/`, ~5k lines. Public API re-exported from `__init__.py`.

| Module | Lines | Purpose |
|---|---|---|
| `loader.py` | ~3000 | `TiledSMISWAXSLoader`, `SAXSGeometry`, `resolve_*_geometry`, raw image loading, baseline/primary metadata helpers, baseline cache. |
| `integrator.py` | ~4000 | `reduce_smi_combined`, `integrate_saxs`, `integrate_waxs`, `MultiPanelArcDetector`, mask construction, q/chi math, sparse pixel→bin plan, merging. |
| `defaults.py` | ~600 | Bundled mask paths, calibration JSONs, `resolve_mask_path`, `load_mask_polygons`, `save_mask_polygons`. |
| `data/masks/{pil2M,900KW}_mask_polygons.json` | — | Static masks. Pin / rod beamstops as `polygon_offsets_from_beam`. |
| `data/saxs_calibration.json` | — | Default SAXS calibration constants (slopes, references). |
| `derived/` | ~600 | Optional: virtual axes, line cuts, peak fits (added in commit `26a46e2`). |
| `upload/` | — | Optional Tiled write-back for derived products. |
| `scripts/` | — | `validate_large_scan.py`, `calibrate_smi_saxs.py`, `calibrate_smi_z_scan.py`, `benchmark_gpu_histogram.py`. |

Public API (from `src/smi_tiled/__init__.py:67-94`):
- Loader: `TiledSMISWAXSLoader`, `SAXSGeometry`, `WAXSGeometry`, `WAXSPanelGeometry`, `resolve_saxs_geometry`, `resolve_waxs_geometry`, `parse_sample_name_geometry`, `infer_detectors_and_steps`, `load_saxs_raw`, `load_waxs_raw`, `populate_cache`, `clear_baseline_cache`.
- Integrator: `reduce_smi_combined`, `reduce_smi_gi`, `CombinedReductionResult`, `GIReductionResult`, `integrate_saxs`, `integrate_waxs`, `MultiPanelArcDetector`, `WAXSCalibration`, `PanelSpec`, `ProgressCallback`.
- Masks: `make_saxs_mask_from_spec`, `make_saxs_mask_from_dict`, `make_waxs_mask_callable`, `make_waxs_mask_callable_from_dict`, `polygons_to_mask`, `mask_for_frame`.
- Caches: `clear_geometry_cache`, `geometry_cache_info`.
- Merge: `merge_q_chi_weighted`, `merge_iq_profiles`, `merge_multiple_qchi`, `merge_reduction_results`.

## Always-loaded facts

These constants and conventions are referenced everywhere; assume them unless proven otherwise.

### Detector inventory

| Detector | Code key | Image field | Shape (rows, cols) | Pixel size | Layout |
|---|---|---|---|---|---|
| Pilatus 2M (SAXS) | `pil2M` | `pil2M_image` | `(1679, 1475)` | 0.172 mm | Single flat panel |
| Pilatus 900KW (WAXS) | `pil900KW` | `pil900KW_image` | `(195, 487)` per panel × 3 | 0.172 mm | 3-panel folded arc, panels at ≈ −7°, 0°, +7° relative to `waxs_arc` |

Field constants (`src/smi_tiled/defaults.py:47-50`, re-exported from `__init__.py:29-32`):
```
SAXS_IMAGE_FIELD = "pil2M_image"
WAXS_IMAGE_FIELD = "pil900KW_image"
WAXS_ARC_FIELD   = "waxs_arc"   # arc angle, deg
WAXS_BSX_FIELD   = "waxs_bsx"   # bsx motor, mm (links to WAXS beamstop position)
```

Tiled defaults (`src/smi_tiled/defaults.py`):
```
DEFAULT_TILED_URI = "https://tiled.nsls2.bnl.gov"
DEFAULT_CATALOG   = "smi/migration"
```

### Streams on a run

A run typically has:
- `primary` — per-frame data, including images and per-frame motor positions.
- `primary.config[<detector>]` — per-detector configuration (one-shot scalars from the bluesky configuration dict).
- `baseline` — start- and end-of-run snapshots of all PVs.

### Per-frame fields you'll commonly need (from `loader.py:2023-2029`)

Per-frame motor positions live on the primary stream:
- `pil2M_motor_x`, `pil2M_motor_y`, `pil2M_motor_z` — SAXS detector translation (mm)
- `piezo_z` — sample piezo Z (mm)
- `waxs_arc` (deg), `waxs_bsx` (mm) — WAXS arc and beamstop X
- `bsx`, `bsy` — SAXS beamstop position (mm); also accessible as the per-beamstop entries in `beamstop_pos_mm`
- `pil2M_active_beamstop` — string, `"rod"` or `"pin"`. Critical for which beamstop polygon is used.

### Metadata fallback chain (essential)

Almost every scalar (SDD, beam center, energy, beamstop, incident angle) resolves through this priority order, documented at `loader.py:18-19`:

```
override (caller-supplied)
  → baseline (cached one-shot scalar)
  → primary stream (per-frame, first value)
  → primary.conf[detector] (configuration scalar)
  → start doc
  → bundled default
```

Helper functions (`loader.py`):
- `_baseline_scalar(run, key)` — `loader.py:720`. Cached read of baseline.
- `_dataset_scalar(ds, key)` — `loader.py:750`. Pull a scalar out of a baseline-shaped Dataset.
- `_primary_scalar(run, field)` — `loader.py:772`. First per-frame value.
- `_primary_conf(run, det_key)` — `loader.py:760`. Returns `primary.config[det_key]` as dict.
- `_read_first_scalar` — first-of-multiple sources (also in `loader.py`).

> **WARNING — silent fallthrough.** If neither baseline nor primary nor conf has a key, the caller's default is used silently. The default for `active_beamstop` is `"rod"` (`loader.py:1140`, `integrator.py:615, 770`). If `pil2M_active_beamstop` is not recorded, the pin polygon is never applied — without warning. See `reference/gotchas.md`.

> **WARNING — bare exception swallow.** `_read_baseline` in `loader.py:681-717` catches `(KeyError, Exception)`; backend errors are silenced. If baseline reads "look empty", verify with `peekAtMd(uid)` directly.

### Public functions you call most

```python
from smi_tiled import (
    TiledSMISWAXSLoader,                # discover + load runs
    reduce_smi_combined,                # full SAXS+WAXS reduction
    resolve_saxs_geometry,              # SAXSGeometry from a run
    resolve_waxs_geometry,              # WAXSGeometry from a run
    make_saxs_mask_from_spec,           # SAXS mask from path/dict/JSON
    make_waxs_mask_callable_from_dict,  # WAXS mask (callable per arc angle)
    populate_cache,                     # write per-UID HDF5 cache
    clear_geometry_cache,               # invalidate _SplitBinPlan cache
)
```

## "Where to check first" — decision tree

This is the diagnostic ladder for "result looks wrong". Walk top-to-bottom; stop at the first deviation.

### Symptom: missing or wrong beamstop masking

```
1. active_beamstop resolution
   - _baseline_scalar(run, "pil2M_active_beamstop")           expect: "pin" or "rod"
   - _dataset_scalar(_read_baseline(run), "pil2M_active_beamstop")
   - _conf_scalar(_primary_conf(run, "pil2M"), "pil2M_active_beamstop")
   - resolve_saxs_geometry(run).active_beamstop               expect: matches above
   FAIL MODE: returns "rod" by silent fallback when key absent.

2. Mask file
   - resolve_mask_path(None, detector="saxs")                 expect: bundled JSON path
   - JSON keys: ["image_shape", "static_regions", "beamstops"]
   - beamstops[active].polygon_offsets_from_beam present?     expect: True
   FAIL MODE: round-trip through load_mask_polygons → save_mask_polygons drops
              polygon_offsets_from_beam; only legacy `polygon` (a stale square)
              remains. See defaults.py:303,313,393.

3. Geometry
   - resolve_saxs_geometry(run): beam_center_row_px, beam_center_col_px,
     dist_m, energy_ev, beamstop_pos_mm.{rod,pin}.{x,y}
   FAIL MODE: setting _SAXS_BEAM_COL_PX_PER_MOTOR_X_MM (etc.) non-zero on top
              of a tracked beam-center PV double-counts the motor (~120 px
              error). See scripts/calibrate_smi_saxs.py:36-44.

4. Built mask
   - make_saxs_mask_from_spec(image_shape, mask_path=None,
                              active_beamstop=geo.active_beamstop,
                              beamstop_pos_mm=geo.beamstop_pos_mm,
                              beam_center_px=(geo.beam_center_row_px,
                                              geo.beam_center_col_px))
   - For pin: 98%+ of pixels in the disk centered at
     (bc_row, bc_col + 5), radius 22 px, should be masked.
   - For rod: vertical strip from beam center upward.

5. Integration result
   - For 16.1 keV / SDD ≈ 2 m / pin: q_min should land near 0.12 Å⁻¹
     (pin polygon edge at ~173 px from beam center along +col direction;
     `q = (4π/λ) sin(½ × atan(173·0.172/2000))` ≈ 0.121 Å⁻¹). Use
     `reference/saxs-geometry.md` for the full derivation.
   - Compare WITH-mask vs NO-mask reductions: the chi sector at
     chi ≈ ±90° (pin location, due to chi = arctan2(qh, qv) and pin
     offset along +column → +90°) should be empty in WITH-mask.
   - reduce_smi_combined returns a CombinedReductionResult; the SAXS-only
     iq lives at `result.saxs["iq"]` (NOT result.iq).
```

### Symptom: wrong q range / wrong beam center

1. Verify `beam_center_x`, `beam_center_y`, `detector_distance` in baseline.
2. If a calibration delta is being used: confirm `saxs_beam_delta_px` and `saxs_distance_delta_mm` arguments to `reduce_smi_combined` reflect the offset from the live PV (not the absolute value).
3. If motor slopes (`_SAXS_BEAM_COL_PX_PER_MOTOR_X_MM` etc.) are nonzero, confirm the beam-center PV is *not* already tracking detector translation. Default slopes are 0.0 because the PV does track; setting them adds a second correction.
4. AgB ring positions: with bundled calibration, the 1st AgB ring should land at q = 2π/5.838 nm = 1.076 Å⁻¹ (= 0.1076 nm⁻¹). Drift > 0.5 % means recalibrate.

### Symptom: TimeRange / searchCatalog ImportError

Pre-existing. `searchCatalog` imports `TimeRange` from `tiled.queries` eagerly. Some installed `tiled` versions don't have it. Workaround: query directly with `Key`:
```python
from tiled.queries import Key
node = loader._get_catalog().search(Key("scan_id") == int(scan_id))
items = list(node.items())
```

For the full per-symptom diagnostic patterns and runnable recipe, see `reference/debug-workflows.md`.

## Must-know gotchas (top 10)

These are the silent-fail traps. Treat as load-bearing.

1. **`active_beamstop` defaults to `"rod"` silently** when no source recorded it. Pin polygon is never applied. (`loader.py:1140`, `integrator.py:615,770`.)
2. **Mask round-trip data loss.** `defaults._coerce_polygon` (line 303) preserves only `polygon`, not `polygon_offsets_from_beam`. Calling `save_mask_polygons` on a loaded mask destroys beamstop offset semantics. Never re-save bundled masks.
3. **Empirical fudge `1.088`** in `make_mask_for_angle` (`integrator.py:587`) when converting `bsx` mm → pixel offset. Undocumented in defaults; comes from instrument calibration.
4. **`_BSX_PER_ARC_DEG = -4.39`** at `integrator.py:3925` shadows `defaults.BSX_PER_ARC_DEG`. Recalibration must update both sites.
5. **Motor-slope double-count.** `_SAXS_BEAM_COL_PX_PER_MOTOR_X_MM` (and y, z) default to 0.0 because the EPICS beam-center PV already tracks detector translation. Setting them non-zero is an error unless the PV is broken. (`scripts/calibrate_smi_saxs.py:36-44`.)
6. **Bare exception swallow** in `_read_baseline` (`loader.py:699`): `except (KeyError, Exception)` silently hides backend errors. If baseline reads look empty, log directly.
7. **Lazy tiled import** in `searchCatalog` (`loader.py:2927`). Eagerly imports `Key, Regex, TimeRange`. Some `tiled` versions miss `TimeRange`, breaking the call. Bypass with direct `Key` queries.
8. **`_DETECTOR_DS_AUTO_MAX_FRAMES = 50`** (`integrator.py:2785`). Per-frame detector ds (raw images + per-frame masks) is auto-disabled silently for scans > 50 frames. To force on for big scans, pass `build_detector_ds=True`.
9. **Geometry cache key incomplete.** `_SAXS_GEOMETRY_CACHE` keys do not include all calibration constants. Editing slopes mid-session can leave a stale cached `_SplitBinPlan`. Call `clear_geometry_cache()` after any calibration edit.
10. **WAXS does not satisfy pyFAI single-panel assumption.** WAXS DataArrays carry pyFAI geometry attrs but the underlying detector is a 3-panel arc. Warning at `loader.py:2615-2619`. Use `MultiPanelArcDetector` / `integrate_waxs` directly; do NOT pipe WAXS through `PFGeneralIntegrator`.

## Reference catalog

When the task needs depth on a topic, read the matching reference file under `reference/`.

| File | Topic |
|---|---|
| `reference/data-layout.md` | Tiled access, run/stream layout, field-name conventions, `searchCatalog`/`peekAtMd`/`summarizeRun`, image array dims (`pix_y`, `pix_x`, `frame`, `waxs_arc`). |
| `reference/metadata-resolution.md` | Full fallback chain, every `_*_scalar` helper, `resolve_incident_angle` 5-tier example, where each canonical scalar is sourced. |
| `reference/saxs-geometry.md` | `SAXSGeometry` dataclass fields, `resolve_saxs_geometry`, motor delta math, calibration JSON schema, the q-from-pixel derivation. |
| `reference/waxs-geometry.md` | `MultiPanelArcDetector`, `PanelSpec`, `WAXSCalibration`, `beam_center_at_angle`, `qmap`, panel rotation math. |
| `reference/masks.md` | Static + beamstop polygon JSON schema, `make_saxs_mask_from_dict`, `polygons_to_mask`, `shift_polygon`, the WAXS callable mask, `_make_waxs_shadow_mask`, `_make_aperture_mask`, `_LargeAreaMaskLazy`, `mask_for_frame`. |
| `reference/integration.md` | q/chi map construction, solid-angle correction, `_SplitBinPlan` (sparse pixel→bin matrix, sub-pixel splitting), `integrate_saxs` per-frame loop, `integrate_waxs` per-panel multi-arc, count-weighted merge formulas. |
| `reference/caching.md` | The 4 cache layers (baseline, SAXS geometry, WAXS qmap, HDF5 disk), env-var paths, invalidation, plus the per-frame zarr streaming store. |
| `reference/outputs.md` | `CombinedReductionResult`, per-detector dict keys, `to_dataarray`, line-cut helpers. |
| `reference/validation.md` | Existing scripts, peek/summarize, AgB ring sanity checks, expected sentinel values. |
| `reference/gotchas.md` | Full failure-mode catalog with file:line, symptoms, and fixes. |
| `reference/debug-workflows.md` | Diagnostic decision trees by symptom, with runnable recipes. |

## Common workflows (quick recipes)

### 1. Look up a run by `scan_id`

```python
from smi_tiled import TiledSMISWAXSLoader

loader = TiledSMISWAXSLoader()  # defaults to tiled.nsls2.bnl.gov / smi/migration
df = loader.searchCatalog(scan_id=1130041)
uid = df.iloc[0]["uid"]
```

If `searchCatalog` errors on `TimeRange`, query directly:
```python
from tiled.queries import Key
node = loader._get_catalog().search(Key("scan_id") == 1130041)
uid, run = next(iter(node.items()))
```

### 2. Reduce a run

```python
from smi_tiled import reduce_smi_combined

result = reduce_smi_combined(
    uid="c849dd2c-…",
    solid_angle_correction=True,
    geometry="transmission",   # or "grazing_incidence"
    n_q=3000,
    pixel_splitting=3,
    saxs_mask_path=None,       # None → bundled default
    waxs_mask_path=None,
)

iq_saxs = result.saxs["iq"]            # xarray.Dataset, dims (q,)
qchi_saxs = result.saxs["q_chi"]       # xarray.Dataset, dims (q, chi)
merged = result.merged_iq              # xarray.Dataset, dims (q,)
da = result.to_dataarray("merged_iq")  # for PyHyperScattering accessors
```

### 3. Build a SAXS mask for a specific run

```python
from smi_tiled import resolve_saxs_geometry, make_saxs_mask_from_spec

geo = resolve_saxs_geometry(run)
mask = make_saxs_mask_from_spec(
    image_shape=(1679, 1475),
    mask_path=None,                     # None → bundled JSON
    active_beamstop=geo.active_beamstop,
    beamstop_pos_mm=geo.beamstop_pos_mm,
    beam_center_px=(geo.beam_center_row_px, geo.beam_center_col_px),
)
# mask is bool, shape (1679, 1475); True = valid pixel.
```

### 4. Inspect a run without loading images

```python
md = loader.peekAtMd(uid)        # tiled metadata only (loader.py:2725)
summary = loader.summarizeRun(uid)  # detectors, n_frames, sample, plan (loader.py:3034)
```

### 5. Diagnose missing pin masking

See `reference/debug-workflows.md` for the full decision-tree. The shortest version:

```python
from smi_tiled.loader import _baseline_scalar, resolve_saxs_geometry
from smi_tiled.integrator import make_saxs_mask_from_spec, _smi_run_raw_shape
from smi_tiled.defaults import resolve_mask_path, SAXS_IMAGE_FIELD

run = loader._get_run(uid)
print("baseline:", _baseline_scalar(run, "pil2M_active_beamstop"))   # "pin"?
geo = resolve_saxs_geometry(run)
print("resolved:", geo.active_beamstop)                              # "pin"?
mask_path = resolve_mask_path(None, detector="saxs")                 # bundled JSON
shape = _smi_run_raw_shape(run, SAXS_IMAGE_FIELD)
mask = make_saxs_mask_from_spec(
    image_shape=shape, mask_path=mask_path,
    active_beamstop=geo.active_beamstop,
    beamstop_pos_mm=geo.beamstop_pos_mm,
    beam_center_px=(geo.beam_center_row_px, geo.beam_center_col_px),
)
# Spot-check the pin disk around (bc_row, bc_col + 5), r=22:
import numpy as np
rr, cc = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), indexing="ij")
disk = (rr - geo.beam_center_row_px)**2 + (cc - geo.beam_center_col_px - 5)**2 <= 22**2
masked_frac = (~mask[disk]).sum() / disk.sum()
print(f"pin disk masked: {masked_frac:.2%}")  # expect >= 95% for pin
```

If `geo.active_beamstop` is `"rod"` but you expect `"pin"`, the resolve chain fell through. If the pin disk is < 95 % masked, the bundled JSON is broken (see gotcha #2).

## House style

When editing `smi-tiled`:
- Follow the existing import style: `from smi_tiled.<module> import …` for cross-module references inside the package; absolute paths for CLI scripts.
- The codebase already uses `numpy`, `xarray`, `dask`, `scipy.sparse`, `zarr`, `h5py`, `tiled`. Don't add new optional dependencies without a clear reason.
- Comments tend to be paragraph-style and explain *why*, not *what*. Match this. The integrator in particular has long comment blocks above tricky code; preserve them.
- Tests live in `tests/`. Run via `pixi run pytest` (env: `default`).
- Pixi project; `pixi run python -c '...'` to invoke.

## Repo conventions

- Branch: `master` is the integration branch.
- Lockfile commits are made with `pixi update`; don't hand-edit `pixi.lock`.
- The bundled mask JSONs and `saxs_calibration.json` are tracked. Regenerate masks via the calibration scripts under `scripts/`, never by `save_mask_polygons` (gotcha #2).

## When unsure, escalate to source

This skill summarizes; the code is authoritative. The integrator and loader are heavily commented — when behavior surprises you, read the function above the call site. For mask schema questions, the JSON itself is the most reliable spec.
