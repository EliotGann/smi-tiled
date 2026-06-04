# Debug workflows

Diagnostic decision trees for common failure modes. Each section: symptom → checks (in order, cheapest first) → likely fix.

The workflows assume you have:
- `loader = TiledSMISWAXSLoader()` and a `uid`,
- `run = loader._get_run(uid)`,
- a working Python session that can import `smi_tiled`.

> No standalone scripts are shipped — these are workflows you run interactively or paste into a notebook.

---

## Symptom: missing or wrong beamstop masking

(SAXS pin or rod beamstop not appearing in the integrated I(q,χ); high counts in chi sector where beamstop should be.)

### Decision tree

```
1. ┌─ Did the run record `pil2M_active_beamstop`?
   │
   │  _baseline_scalar(run, "pil2M_active_beamstop")
   │       expect: "pin" or "rod"
   │       FAIL: returns None  ← run never recorded it
   │
2. └─→ Does resolve_saxs_geometry return the right value?
   │
   │  geo = resolve_saxs_geometry(run)
   │  geo.active_beamstop
   │       expect: matches step 1
   │       FAIL: returns "rod" (silent default)
   │
3. └─→ Does the bundled mask JSON have the polygon?
   │
   │  resolve_mask_path(None, "saxs")
   │       expect: bundled JSON path
   │  load_mask_polygons(path)["beamstops"][geo.active_beamstop]
   │       expect: non-empty list of vertices
   │       FAIL: empty list ← polygon was lost (gotcha #2)
   │
4. └─→ Does the built mask actually mask the beamstop region?
   │
   │  mask = make_saxs_mask_from_spec(image_shape, mask_path=None,
   │            active_beamstop=geo.active_beamstop,
   │            beamstop_pos_mm=geo.beamstop_pos_mm,
   │            beam_center_px=(geo.beam_center_row_px, geo.beam_center_col_px))
   │  Test: pin disk at (bc_row, bc_col + 5), r=22 → ≥95% masked
   │        FAIL: mask polygons placed wrong ← bad calibration
   │
5. └─→ Does the integrated result still show the beamstop?
   │
   │  Compare merged_qchi.intensity at chi ≈ +90° (pin) or chi ≈ +90°/−90° (rod)
   │       Should be NaN / very low counts in the masked region
   │       FAIL: mask not applied to integration ← caller bypassed mask
```

### Recipe

```python
from smi_tiled import (
    TiledSMISWAXSLoader, resolve_saxs_geometry,
    make_saxs_mask_from_spec,
)
from smi_tiled.loader import _baseline_scalar
from smi_tiled.defaults import resolve_mask_path, load_mask_polygons
from smi_tiled.integrator import _smi_run_raw_shape
from smi_tiled.defaults import SAXS_IMAGE_FIELD
import numpy as np

loader = TiledSMISWAXSLoader()
uid = "your-uid-here"
run = loader._get_run(uid)

# Step 1
print("baseline:", _baseline_scalar(run, "pil2M_active_beamstop"))

# Step 2
geo = resolve_saxs_geometry(run)
print(f"resolved: bc=({geo.beam_center_row_px}, {geo.beam_center_col_px}) "
      f"dist={geo.dist_m} bs={geo.active_beamstop}")

# Step 3
mask_path = resolve_mask_path(None, detector="saxs")
mask_dict = load_mask_polygons(mask_path)
print(f"polygons for {geo.active_beamstop}: "
      f"{len(mask_dict['beamstops'].get(geo.active_beamstop, []))} vertices")

# Step 4
shape = _smi_run_raw_shape(run, SAXS_IMAGE_FIELD)
mask = make_saxs_mask_from_spec(
    image_shape=shape,
    mask_path=mask_path,
    active_beamstop=geo.active_beamstop,
    beamstop_pos_mm=geo.beamstop_pos_mm,
    beam_center_px=(geo.beam_center_row_px, geo.beam_center_col_px),
)
rr, cc = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), indexing="ij")
disk = (rr - geo.beam_center_row_px)**2 + (cc - geo.beam_center_col_px - 5)**2 <= 22**2
masked_frac = (~mask[disk]).sum() / disk.sum()
print(f"pin-disk masked fraction: {masked_frac:.2%}  (expect >= 95% for pin)")

# Step 5: full reduction + check chi sector
from smi_tiled import reduce_smi_combined
result = reduce_smi_combined(uid=uid)
qchi = result.merged_qchi
chi = qchi.chi.values
chi_sector = (chi > 80) & (chi < 100)  # pin sector for SMI
mean_intensity_in_sector = np.nanmean(qchi.intensity.values[:, chi_sector])
print(f"mean intensity in pin sector (should be near 0): {mean_intensity_in_sector:.2f}")
```

### Most common causes

- `geo.active_beamstop == "rod"` but the run used pin (gotcha #1) → pass `active_beamstop="pin"` explicitly.
- `mask_dict["beamstops"]["pin"]` is empty (gotcha #2) → bundled mask was overwritten; restore from git.
- Beam center off by ~100 px (gotcha #5) → check the `_SAXS_BEAM_*_PER_MOTOR_*_MM` constants are 0.0.

---

## Symptom: wrong q range / shifted AgB ring

### Decision tree

```
1. ┌─ Energy correct?
   │  geo.energy_ev → expect 14000–17500 eV
   │  Convert: λ = 1239.84193 / energy_ev nm
   │
2. └─→ Distance correct?
   │  geo.dist_m → expect 1.7–9.3 m for SAXS
   │  This is `pil2M_motor_z` + offset. Check baseline.
   │
3. └─→ Beam center correct?
   │  geo.beam_center_row_px, geo.beam_center_col_px
   │  For SAXS at default motor positions: expect (~1107, ~744) ± 50
   │
4. └─→ Bragg check
   │  AgB d = 5.838 nm → q1 = 2π/d = 1.076 nm⁻¹
   │  In I(q): peak should be at 1.076 ± 0.01 nm⁻¹
```

### Recipe

```python
from smi_tiled import resolve_saxs_geometry
import numpy as np

geo = resolve_saxs_geometry(run)
hc_ev_nm = 1239.84193
wavelength_nm = hc_ev_nm / geo.energy_ev
print(f"E = {geo.energy_ev} eV, λ = {wavelength_nm:.5f} nm")
print(f"SDD = {geo.dist_m:.4f} m, BC = ({geo.beam_center_row_px}, {geo.beam_center_col_px})")
print(f"Active beamstop: {geo.active_beamstop}")

# AgB ring 1 expected pixel radius:
two_theta = 2 * np.arcsin(wavelength_nm / (2 * 5.838))
r_mm = geo.dist_m * 1000.0 * np.tan(two_theta)
r_px = r_mm / 0.172
print(f"AgB ring 1 expected at radius = {r_px:.1f} px from beam center")

# After reduction, find where AgB peak actually lands
result = reduce_smi_combined(uid=uid)
q = result.merged_iq["q"].values
I = result.merged_iq["I"].values
finite = np.isfinite(I) & (q > 0.9) & (q < 1.2)
peak_q = q[finite][np.argmax(I[finite])]
print(f"AgB peak found at q = {peak_q:.4f} nm⁻¹  (expect 1.076)")
print(f"  drift: {(peak_q - 1.076) / 1.076 * 100:+.2f}%")
```

### Likely causes by drift magnitude

| Drift | Check |
|---|---|
| < 1% | Within calibration tolerance — accept. |
| 1-5% | Likely energy mis-set; check `geo.energy_ev` vs storage ring. |
| 5-20% | Distance wrong; check `pil2M_motor_z` calibration. |
| > 20% | Beam center wrong, or wrong detector. |
| No peak found | Mask masked the ring (over-masked), or the run had no AgB. |

---

## Symptom: WAXS reduction produces wrong q

### Decision tree

```
1. ┌─ All 3 panels appearing?
   │  result.waxs["ds"] should have 3 panel-like regions per frame.
   │  If only 1 visible: panel mask too aggressive, or panel rotation wrong.
   │
2. └─→ Per-frame beam center correct?
   │  Check WAXSCalibration.beam_center_at_angle(arc) for each frame
   │  Expect: bc_col = bc_col_arc0 + (-_BSX_PER_ARC_DEG * arc) / pixel_size * fudge
   │
3. └─→ Panel rotations correct?
   │  PanelSpec angles: [-7°, 0°, +7°] relative to waxs_arc
   │  Wrong sign → panels mirror about the central panel
   │
4. └─→ Bsx reference correct?
   │  See gotcha #14: waxs_bsx_ref must reflect arc-0 position, not first frame
```

### Recipe

```python
from smi_tiled import resolve_waxs_geometry, reduce_smi_combined
from smi_tiled.integrator import WAXSCalibration, MultiPanelArcDetector
import numpy as np

geo = resolve_waxs_geometry(run)
print(f"WAXS geo: dist={geo.dist_m:.4f} m, "
      f"bc=({geo.beam_center_row_px}, {geo.beam_center_col_px}), "
      f"E={geo.energy_ev} eV, n_panels={len(geo.panels)}")

# Inspect panel specs
for i, panel in enumerate(geo.panels):
    print(f"  panel {i}: angle={panel.angle_deg}°, "
          f"row_offset={panel.row_offset_px}, col_offset={panel.col_offset_px}")

# Run reduction with detector_ds for inspection
result = reduce_smi_combined(uid=uid, build_detector_ds=True)

if result.waxs is None:
    print("No WAXS data in this run.")
else:
    ds = result.waxs["ds"]
    print(f"detector_ds shape: {ds['intensity'].shape}")
    print(f"q range per frame:")
    for fi in range(ds.sizes["frame"]):
        q = ds["q_abs"].isel(frame=fi).values
        finite = np.isfinite(q) & (q > 0)
        print(f"  frame {fi}: arc={ds['waxs_arc'].isel(frame=fi).values:.2f}°, "
              f"q ∈ [{q[finite].min():.3f}, {q[finite].max():.3f}] nm⁻¹")
```

### Common WAXS bugs

- Panels mirrored: `flip_horizontal` toggle. Try `waxs_kwargs={"flip_horizontal": True}`.
- bsx reference off: pass explicit `waxs_kwargs={"waxs_bsx_ref": <known_arc0_bsx>}`.
- Constant `_BSX_PER_ARC_DEG` wrong: edit at `integrator.py:3925` AND `defaults.py:195` (gotcha #4).
- Stale qmap cache after edit: restart Python (gotcha #11).

---

## Symptom: SAXS+WAXS merge has visible seam

(Discontinuity in `merged_iq["I"]` around q ≈ 1-2 nm⁻¹ where SAXS ends and WAXS begins.)

### Decision tree

```
1. ┌─ Do per-detector intensities agree in the overlap region?
   │  result.saxs["iq"] vs result.waxs["iq"] at q in overlap
   │  > 2× disagreement: scale factor or geometry mismatch
   │
2. └─→ Solid-angle correction applied to both?
   │  Single bool, applied to both detectors (gotcha #16)
   │  Default True; toggle via solid_angle_correction kwarg
   │
3. └─→ Are mask coverage fractions similar?
   │  SAXS valid_frac and WAXS valid_frac should both be ~95%+
   │  Big difference → one detector is over-masked, biasing intensity
   │
4. └─→ Diagnostic: split the merged data
   │  merged_qchi.saxs_intensity vs merged_qchi.waxs_intensity
   │  Plot both at fixed chi → should overlap in the seam region
```

### Recipe

```python
import numpy as np
import matplotlib.pyplot as plt

result = reduce_smi_combined(uid=uid)
saxs_q = result.saxs["iq"]["q"].values
saxs_I = result.saxs["iq"]["I"].values
waxs_q = result.waxs["iq"]["q"].values
waxs_I = result.waxs["iq"]["I"].values
merged_q = result.merged_iq["q"].values
merged_I = result.merged_iq["I"].values

# Overlap region
q_lo = max(saxs_q.min(), waxs_q.min())
q_hi = min(saxs_q.max(), waxs_q.max())
print(f"Overlap region: [{q_lo:.3f}, {q_hi:.3f}] nm⁻¹")

# Sample both detectors at the centre of the overlap
mid = (q_lo + q_hi) / 2
saxs_at_mid = np.interp(mid, saxs_q, saxs_I)
waxs_at_mid = np.interp(mid, waxs_q, waxs_I)
print(f"At q={mid:.3f}: SAXS I={saxs_at_mid:.4g}, WAXS I={waxs_at_mid:.4g}")
print(f"Ratio: {saxs_at_mid / waxs_at_mid:.3f}  (expect 1.0 ± 20%)")

plt.loglog(saxs_q, saxs_I, "b", label="SAXS")
plt.loglog(waxs_q, waxs_I, "r", label="WAXS")
plt.loglog(merged_q, merged_I, "k--", label="merged")
plt.axvspan(q_lo, q_hi, alpha=0.2, color="gray", label="overlap")
plt.xlabel("q (nm⁻¹)"); plt.ylabel("I"); plt.legend()
plt.show()
```

### Common causes

- Solid-angle correction differs (it shouldn't — gotcha #16).
- WAXS panel orientation wrong → q-values shifted.
- One detector saturated → intensity capped.
- Mask over-masking one detector → counts lost.

---

## Symptom: reduction returns None / silent failure

### Decision tree

```
1. ┌─ Did the run have the detector?
   │  loader.summarizeRun(uid)["detectors"]
   │  expect: list containing "pil2M_image" and/or "pil900KW_image"
   │
2. └─→ Does the primary stream have the field?
   │  _has_primary_field(run, "pil2M_image")  # → True/False
   │
3. └─→ Is the run reachable from tiled?
   │  loader._get_run(uid)  ← raises if unreachable
   │
4. └─→ Cache poisoning?
   │  Try with image_cache_path=None to bypass disk cache
```

### Recipe

```python
from smi_tiled.loader import _has_primary_field
from smi_tiled.defaults import SAXS_IMAGE_FIELD, WAXS_IMAGE_FIELD

print(loader.summarizeRun(uid))
print(f"Has SAXS field: {_has_primary_field(run, SAXS_IMAGE_FIELD)}")
print(f"Has WAXS field: {_has_primary_field(run, WAXS_IMAGE_FIELD)}")

# Force re-fetch from tiled, bypassing cache
result = reduce_smi_combined(uid=uid, image_cache_path=None)
```

If `image_cache_path=None` succeeds where `"auto"` failed: cache is corrupted (gotcha #13), `rm` it.

---

## Symptom: `searchCatalog` ImportError

```
ImportError: searchCatalog requires `tiled.queries`
```

(But `tiled` IS installed.)

### Cause

Eager import of `Key, Regex, TimeRange` at `loader.py:2927`. Some `tiled` versions don't export `TimeRange` (gotcha #7).

### Workaround

```python
from tiled.queries import Key
node = loader._get_catalog().search(Key("scan_id") == 1130041)
items = list(node.items())
uid, run = items[0]
```

### Permanent fix

Modify `loader.py:2927` to import lazily, or check `hasattr(queries, "TimeRange")` before constructing the `TimeRange` query.

---

## Symptom: cache returns stale results after code edit

(Reduction returns the OLD result for the same UID after editing q-map / calibration code.)

### Cache layers (see `caching.md` for full detail)

| Layer | Cleared by |
|---|---|
| Baseline | `clear_baseline_cache()` |
| SAXS geometry | `clear_geometry_cache()` |
| WAXS qmap | Process restart only |
| Disk (HDF5) | `rm $TMPDIR/smi_browser_cache/<uid>.h5` |
| SplitBinPlan | Per-call (no global cache) |

### Recipe

```python
from smi_tiled.integrator import clear_geometry_cache
from smi_tiled.loader import _BASELINE_CACHE, _BASELINE_COLUMNS_CACHE

clear_geometry_cache()           # SAXS qmap
_BASELINE_CACHE.clear()          # baseline xr.Dataset
_BASELINE_COLUMNS_CACHE.clear()  # baseline column list

# WAXS qmap requires Python restart — do this last:
# exit() / restart kernel
```

For disk cache:

```bash
rm -f $TMPDIR/smi_browser_cache/*.h5
# or with explicit env var:
rm -f $SMI_BROWSER_CACHE_DIR/*.h5
```

Or per-uid:
```bash
rm -f $TMPDIR/smi_browser_cache/<uid>.h5
```

---

## Symptom: out-of-memory on a large scan

### Trigger checks

- `result.saxs["ds"]` and `result.waxs["ds"]` populated for > 50 frames? (gotcha #8) — these are the biggest memory hogs.
- `result.per_frame_qchi["saxs"]` / `["waxs"]` are full in-memory `(n_frames, n_q, n_chi)` Datasets.
- Geometry cache holding multiple entries (~95 MB each for SAXS).

### Recipe

```python
from smi_tiled.integrator import geometry_cache_info, clear_geometry_cache

print(geometry_cache_info())
# {"saxs_entries": 3, "saxs_est_mb": 285, ...}

# Reduce memory pressure for large scans:
result = reduce_smi_combined(
    uid=uid,
    build_detector_ds=False,                          # don't store per-frame raw images
    build_frame_qchi=True,
    frame_qchi_store=f"/tmp/qchi_{uid}.zarr",         # stream to disk
    image_cache_path="auto",                          # use disk cache
)
```

The zarr-backed `q_chi_frames` keeps peak RAM at ~1 frame instead of N frames. See `outputs.md` for the lazy zarr mode.

For very large scans, use `validate_large_scan.py` as a benchmark:
```bash
pixi run python scripts/validate_large_scan.py
```
It uses `tracemalloc` to report peak memory per stage.

---

## Symptom: incident angle wrong

(GI reduction shows wrong qz, or transmission scan accidentally treated as GI.)

### Decision tree

`resolve_incident_angle` walks 5 sources (`metadata-resolution.md`):
1. Explicit override
2. `target_file_name` (per-frame string parse)
3. `sample_name` (start metadata string parse)
4. `stage_th + piezo_th` (motor sum)
5. Default 0°

### Recipe

```python
from smi_tiled.loader import resolve_incident_angle

ai_deg, source = resolve_incident_angle(run, n_frames=11)
print(f"αᵢ = {ai_deg} (resolved via: {source})")
```

If `source` is `"default 0°"` for a GI scan, none of the upstream sources had usable info. Pass an explicit override:
```python
result = reduce_smi_gi(uid=uid, incident_angle_deg=0.5)
```

---

## Symptom: structured per-frame strings not parsed (`fn:*` columns missing)

(You expect `result.per_frame_iq["fn:eV"]` but it's not there.)

### Checks

1. Source field exists?
   ```python
   "target_file_name" in result.per_frame_iq.data_vars
   ```
2. Source field is non-empty per frame?
   ```python
   names = result.per_frame_iq["target_file_name"].values
   print(names[:3])  # bytes or str?
   ```
3. Did `apply_virtual_axes` run? (Default yes; check `virtual_axes` kwarg.)
4. Was the column filtered by `min_fill`?

### Recipe

```python
from smi_tiled.derived import VirtualAxesConfig, apply_virtual_axes

# Force re-run with permissive filter
cfg = VirtualAxesConfig(min_fill=0.01)
apply_virtual_axes(result, cfg)

print([v for v in result.per_frame_iq.data_vars if v.startswith("fn:")])
```

If the column appears with `min_fill=0.01` but not `0.5`, your scan has < 50% frames with the structured string in `target_file_name`. Either fix the source data or accept the lower fill.

---

## Symptom: fitting peaks gives all-NaN or all-False success

### Checks (in order)

1. `result.per_frame_iq["I"]` has finite values?
2. The peak window `[q_min, q_max]` actually contains the peak?
3. The peak is above noise (SNR > 3)?
4. Width hasn't hit upper bound (gives `success=False`)?

### Recipe

```python
import numpy as np
from smi_tiled.derived import PeakDef, apply_peak_fits

peak = PeakDef(name="agb_1", q_min=0.9, q_max=1.2, model="gaussian", baseline="linear")

# Sanity-check the data first
q = result.per_frame_iq["q"].values
I = result.per_frame_iq["I"].values
in_window = (q >= peak.q_min) & (q <= peak.q_max)
print(f"Window: {in_window.sum()} q-bins, "
      f"frame 0 finite: {np.isfinite(I[0, in_window]).sum()}")
print(f"Frame 0 max in window: {np.nanmax(I[0, in_window]):.4g}")

apply_peak_fits(result, [peak])
print(result.peak_fits)

# If all NaN, lower thresholds (debug-only):
# patch: smi_tiled.derived.peakfit.MIN_SNR = 1.0; MIN_R2 = 0.0
```

---

## Symptom: `populate_cache` slow / hangs

`populate_cache` reads detector image stacks from tiled. For a 1679×1475 detector with 200 frames, that's ~2 GB over HTTP. If the network is slow, expect minutes.

### Checks

```bash
# Test tiled bandwidth:
pixi run python -c "
from smi_tiled.loader import TiledSMISWAXSLoader
import time
t0 = time.perf_counter()
loader = TiledSMISWAXSLoader()
run = loader._get_run('your-uid')
img = run['primary']['data']['pil2M_image'].read()
print(f'shape={img.shape}, time={time.perf_counter() - t0:.1f}s')
"
```

If > 30s for a small scan, the network is the bottleneck. Options:
- Run `populate_cache` once during a low-traffic period; subsequent reductions are fast.
- Use `include_images=False` for cache-only-metadata mode.
- Move to a host with better network access to `tiled.nsls2.bnl.gov`.

---

## General debugging principles

1. **Start cheap**: `summarizeRun` → `peekAtMd` → `resolve_*_geometry`. These don't read images.
2. **Bypass caches when in doubt**: pass `image_cache_path=None`, call `clear_geometry_cache()`, restart Python. Stale caches are #1 source of "but I edited the code!" confusion.
3. **Diagnose in detector space**: `build_detector_ds=True` lets you visualize raw frames + masks + q-maps. Confirms whether the bug is in input data or reduction.
4. **Compare against known-good**: save reduction outputs as NetCDF for known-good UIDs; diff against current results.
5. **Read the source comment**: the code has paragraph-style explanatory comments above tricky logic. Check `loader.py` and `integrator.py` near the line of interest.

## Symptom-to-file map

| Symptom | First check | Reference |
|---|---|---|
| Wrong beamstop masking | `_baseline_scalar`, `geo.active_beamstop` | `masks.md`, `metadata-resolution.md` |
| Wrong q range | `geo.dist_m`, `geo.beam_center_*`, AgB peak | `saxs-geometry.md`, `validation.md` |
| WAXS q wrong | Panel angles, `_BSX_PER_ARC_DEG`, qmap cache | `waxs-geometry.md`, `gotchas.md` #4, #11 |
| Merge seam | Per-detector I(q) at overlap | `integration.md`, `gotchas.md` #16 |
| Reduction returns None | `summarizeRun`, `_has_primary_field` | `data-layout.md` |
| `searchCatalog` import error | tiled version | `gotchas.md` #7 |
| Stale results after edit | Caches | `caching.md`, `gotchas.md` #9, #11 |
| OOM on large scan | `build_detector_ds`, `frame_qchi_store` | `outputs.md`, `gotchas.md` #8 |
| Wrong αᵢ | `resolve_incident_angle` chain | `metadata-resolution.md` |
| `fn:*` columns missing | `target_file_name`, `min_fill` | `outputs.md` |
| Peak fits fail | I(q) values, SNR | `outputs.md` peak fits section |
| populate_cache slow | Network bandwidth | `caching.md` |
