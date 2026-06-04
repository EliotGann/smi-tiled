# Gotchas

Silent-fail traps in `smi-tiled`. Each entry: location, symptom, why, fix.

These are load-bearing — most reduction bugs trace to one of these.

## 1. `active_beamstop` defaults to `"rod"` silently

**Location**: `loader.py:1140`, `integrator.py:615`, `integrator.py:770`, `integrator.py:1134`.

```python
# loader.py:1140 (SAXSGeometry default)
active_beamstop: str = "rod"
```

**Symptom**: Pin disk masking missing in reductions, even though the run used the pin beamstop. Visible as a bright spot at (bc_row, bc_col + 5) in `merged_qchi`, and as anomalously high counts in the chi ≈ 90° sector of I(q,χ).

**Why**: `resolve_saxs_geometry` walks the metadata fallback chain looking for `pil2M_active_beamstop` (override → primary config → baseline → primary stream → defaults). If none of those have the key, it returns `"rod"` without warning. The bundled mask JSON keys polygons by beamstop name (`"rod"` vs `"pin"`), so the wrong polygon set is applied.

**Fix**:
1. Verify `_baseline_scalar(run, "pil2M_active_beamstop")` returns `"pin"` (or `"rod"`).
2. If it returns `None`, the run never recorded which beamstop was used. Pass an explicit override:
   ```python
   geo = resolve_saxs_geometry(run, active_beamstop="pin")
   ```
3. Or rebuild the mask explicitly:
   ```python
   mask = make_saxs_mask_from_spec(
       image_shape=(1679, 1475),
       mask_path=None,
       active_beamstop="pin",        # <- explicit
       beamstop_pos_mm=...,
       beam_center_px=(...),
   )
   ```

## 2. Mask round-trip data loss

**Location**: `defaults.py:303` (`_coerce_polygon`), `defaults.py:393` (`save_mask_polygons`).

```python
# defaults.py:303
def _coerce_polygon(verts: Any) -> list:
    out = []
    for v in verts:
        c, r = v[0], v[1]
        out.append([float(c), float(r)])
    return out
```

**Symptom**: After calling `save_mask_polygons` on a loaded mask, the bundled JSON loses its `polygon_offsets_from_beam` key. Subsequent reductions silently skip the beamstop polygon (because the offsets were how the polygon got positioned relative to the beam center).

**Why**: `_coerce_polygon` only preserves the bare 2-D vertex list. Wrapper dicts like `{"polygon": [...], "polygon_offsets_from_beam": [...]}` are unwrapped on load (`load_mask_polygons` at `defaults.py:368-376`), and on save only the `polygon` list is written back (`defaults.py:406-407`). The `polygon_offsets_from_beam` field is **lost**.

**Fix**:
- **Never** `save_mask_polygons` on the bundled mask. The bundled JSON is a reference and shouldn't be re-saved.
- If you need to edit a mask, edit the JSON file directly.
- If you need to programmatically generate a new mask, do so in the JSON schema directly using a text editor or `json.dump(...)`, preserving `polygon_offsets_from_beam`.
- Use the calibration scripts under `scripts/` to (re)generate masks; they write the JSON in the correct schema.

## 3. Empirical fudge factor `1.088`

**Location**: `integrator.py:587`, `integrator.py:891`, `integrator.py:51` (docstring).

```python
# integrator.py:587
bs_shift_px = (bs_shift_mm / pixel_size_mm) * 1.088  # empirical fudge factor
```

**Symptom**: WAXS beamstop polygon mis-positioned by ~9% in pixels.

**Why**: The naive `bsx_mm / pixel_size_mm` conversion under-counts the actual pixel shift by 8.8%. The factor was determined empirically from instrument calibration; the cause is likely a small geometric correction (panel tilt, off-axis projection) not captured in the simple linear model. **Not documented in `defaults.py`** — only lives in `integrator.py`.

**Fix**:
- Don't change this number unless you've recalibrated against a known beamstop position.
- If you write new mask code, mirror the same `* 1.088` factor to stay consistent with the existing reductions.

## 4. `_BSX_PER_ARC_DEG = -4.39` shadows public constant

**Location**: `integrator.py:3925` (private), `defaults.py:195` (public).

```python
# integrator.py:3925
_BSX_PER_ARC_DEG = -4.39  # mm/deg, SMI mechanical linkage

# defaults.py:195
BSX_PER_ARC_DEG: float = -4.39
```

**Symptom**: After recalibrating the bsx-vs-arc linkage, only one site is updated. Beam center column is wrong for non-zero `waxs_arc` frames.

**Why**: The constant is duplicated. Worse, in the function at `integrator.py:1098` there's:
```python
from smi_tiled.defaults import BSX_PER_ARC_DEG as _BSX_PER_ARC_DEG_PUBLIC
```
— a third name to keep track of.

**Fix**:
- If you recalibrate, update **all three**:
  - `defaults.py:195` (`BSX_PER_ARC_DEG`)
  - `integrator.py:3925` (`_BSX_PER_ARC_DEG` local)
  - Any test fixtures that pin a value.
- Search regex: `BSX_PER_ARC_DEG`.

## 5. Motor-slope double-count

**Location**: `loader.py:_SAXS_BEAM_COL_PX_PER_MOTOR_X_MM`, `_SAXS_BEAM_ROW_PX_PER_MOTOR_Y_MM`, etc. Default value 0.0.

Documented at `scripts/calibrate_smi_saxs.py:36-44`:

```
.. warning::
    The fitted ``b_x`` / ``b_y`` (~5.82 / 5.99 px/mm) describe how the EPICS
    beam-center PV *itself* moves with motor_x/y — the PV already tracks the
    detector translation.  Do NOT paste them into
    ``_SAXS_BEAM_COL_PX_PER_MOTOR_X_MM`` / ``_SAXS_BEAM_ROW_PX_PER_MOTOR_Y_MM``:
    the loader reads that same live PV, so adding ``(motor - ref) * slope`` on
    top double-counts and throws the beam center ~120 px off.
```

**Symptom**: Beam center off by ~100 px after running `calibrate_smi_saxs.py`. AgB ring centred at the wrong q.

**Why**: The EPICS BC PV (`pil2M_beam_center_x_PV`) **already tracks** detector translation — it is updated by the IOC in response to motor moves. The calibration script regresses that PV against the motor and finds slope ≈ 5.82 px/mm. Pasting that slope into `_SAXS_BEAM_COL_PX_PER_MOTOR_X_MM` makes the loader compute:
```
bc_col = bc_pv + (motor - ref) * 5.82
```
which double-counts the motor's effect.

**Fix**:
- These slopes default to 0.0. **Keep them at 0.0**.
- Set non-zero only if a future detector reports a *static, position-independent* BC PV.
- If your reduction is off by ~100 px and you recently edited `loader.py`, check these constants first.

## 6. Bare exception swallow in `_read_baseline`

**Location**: `loader.py:699`.

```python
try:
    result = run["baseline"].read()
except (KeyError, Exception):
    pass
```

**Symptom**: Baseline data appears empty for a run that should have it. No exception, no log.

**Why**: `except (KeyError, Exception)` is equivalent to `except Exception` (the `KeyError` is redundant because `Exception` already catches it). All errors silently fall through to the next code path. If the tiled backend is misbehaving, you get an empty Dataset instead of a connection error.

**Fix**:
- When debugging "missing baseline":
  ```python
  result = run["baseline"].read()  # let it raise
  ```
- Or temporarily change the except clause to log:
  ```python
  except (KeyError, Exception) as exc:
      print(f"baseline read failed: {exc!r}")
  ```
- The sloppy `except` clause is a code-style artifact; consider tightening when you touch this function next.

## 7. Lazy tiled import / `TimeRange` missing

**Location**: `loader.py:2927`.

```python
try:
    from tiled.queries import Key, Regex, TimeRange  # type: ignore
except Exception as exc:
    raise ImportError(...)
```

**Symptom**: `searchCatalog(...)` raises `ImportError` even when `tiled` is installed.

**Why**: Some `tiled` versions (older or partial installs) don't expose `TimeRange` from `tiled.queries`. The import is eager — even if you don't use `since`/`until`, the import fails at module load.

**Fix**:
- **Workaround** without modifying smi-tiled:
  ```python
  from tiled.queries import Key
  node = loader._get_catalog().search(Key("scan_id") == scan_id)
  uid, run = next(iter(node.items()))
  ```
- **Permanent fix**: change the import at `loader.py:2927` to lazy `try/except` per name, or guard `TimeRange` behind a check that only imports it when `since`/`until` is provided.

## 8. `_DETECTOR_DS_AUTO_MAX_FRAMES = 50`

**Location**: `integrator.py:2785`.

```python
_DETECTOR_DS_AUTO_MAX_FRAMES = 50

# Used at integrator.py:2798:
if n_frames <= _DETECTOR_DS_AUTO_MAX_FRAMES:
    build_detector_ds = True   # silent default
else:
    print(f"... {n_frames} frames (> {_DETECTOR_DS_AUTO_MAX_FRAMES}) to bound memory; ...")
    build_detector_ds = False
```

**Symptom**: `result.saxs["ds"]` is `None` for a 60-frame scan, even though you wanted detector-space output.

**Why**: For a 1679 × 1475 detector, each frame's `(intensity, q_abs, q_horizontal, q_vertical, mask)` tuple is ~50 MB float64. At 50 frames that's already 2.5 GB. The auto-disable bounds memory.

**Fix**: Force on with explicit kwarg:
```python
result = reduce_smi_combined(
    uid=uid,
    build_detector_ds=True,   # explicit, overrides auto-disable
    ...
)
```
Be sure your machine has enough memory. For a 200-frame scan with both detectors, peak RAM can hit 30+ GB.

## 9. Geometry cache key incomplete

**Location**: `integrator.py:252` (`_saxs_cache_key`), `integrator.py:169` (`_SAXS_GEOMETRY_CACHE`).

```python
_SAXS_GEOMETRY_CACHE: dict[tuple, tuple] = {}
# Key: (dist_m, poni1_m, poni2_m, pixel1_m, pixel2_m, wavelength_m, shape)
```

**Symptom**: After editing q-map computation code (e.g. fixing a sign convention or a calibration constant), reductions return the OLD result for the same UID.

**Why**: The cache key is the geometry **inputs**, not the **code version**. The cached `(q2d, qh2d, qv2d, chi_deg_2d, sa_base)` tuple is reused for any reduction with the same key. Edits to `_compute_saxs_qmaps` don't invalidate.

**Fix**:
```python
from smi_tiled.integrator import clear_geometry_cache
clear_geometry_cache()
```
Or restart the Python process. Make this part of your workflow when editing q-map / calibration code.

> See also `caching.md` for the **WAXS qmap cache** which has no clear function — it requires process restart. And there's a SplitBinPlan cache built per-call (not module-global, so it's fine).

## 10. WAXS does not satisfy pyFAI single-panel assumption

**Location**: `loader.py:2615-2619` (warning), `integrator.py` `MultiPanelArcDetector`.

```
WAXS raw DataArrays (dims ['waxs_arc', 'pix_y', 'pix_x']) describe a
3-panel folded arc detector whose geometry cannot be modeled with the
single flat-panel pyFAI assumption used by PFGeneralIntegrator.
Always use this package's own reduce_smi_combined for WAXS reduction.
```

**Symptom**: Piping `waxs_raw` into `PyHyperScattering.PFGeneralIntegrator(geomethod='template_xr')` produces wrong q values. Rings are mis-located.

**Why**: WAXS is 3 panels at angles [-7°, 0°, +7°] relative to `waxs_arc`. Each panel has its own beam center, distance, and orientation in the lab frame. pyFAI assumes one flat panel.

**Fix**: Use `integrate_waxs(...)` directly, or `reduce_smi_combined(...)` which calls it. Never pipe WAXS DataArrays through pyFAI.

## 11. WAXS qmap cache has no clear function

**Location**: `_WAXS_QMAP_CACHE` inside `integrate_waxs` (`integrator.py:3200-3300`).

**Symptom**: Same as #9 but for WAXS — edits to WAXS q-map computation are stale.

**Why**: `clear_geometry_cache()` (`integrator.py:172`) only clears `_SAXS_GEOMETRY_CACHE`. The WAXS qmap cache is keyed on `(theta_rounded, wavelength_rounded)` and is never explicitly cleared.

**Fix**: Restart the Python process. There is no `clear_waxs_qmap_cache()` function. (TODO: add one.)

## 12. `_BASELINE_CACHE` sticky `None` value

**Location**: `loader.py:38`, `loader.py:693-694`.

```python
_BASELINE_CACHE: dict[str, xr.Dataset | None] = {}
# ...
if _rid in _BASELINE_CACHE:
    return _BASELINE_CACHE[_rid]   # may be None!
```

**Symptom**: First call to `_read_baseline(run)` returns `None` (transient backend error). All subsequent calls also return `None` for the same UID — the cache hit returns the cached `None`.

**Why**: The check is `if _rid in _BASELINE_CACHE`, which is True for `_BASELINE_CACHE[_rid] = None`. The cache treats "no baseline" and "we don't know yet" identically.

**Fix**:
- Call `clear_baseline_cache()` after a transient error.
- Or directly: `del _BASELINE_CACHE[uid]` for that run.

## 13. Population mid-write corruption

**Location**: `loader.py:423` (`populate_cache`).

**Symptom**: HDF5 disk cache returns partial data; subsequent reductions are inconsistent across runs of the same scan.

**Why**: `populate_cache` writes to HDF5 in append mode (`h5py.File(cache_path, "a")`). If interrupted mid-write (Ctrl-C, crash, OOM), the file is left in a partial state. h5py doesn't atomically rename — it edits the file in-place.

**Fix**:
```bash
rm $TMPDIR/smi_browser_cache/<uid>.h5
# or
rm $SMI_BROWSER_CACHE_DIR/<uid>.h5
```
Then re-run the reduction with `image_cache_path="auto"`; it will re-populate.

## 14. WAXS bsx reference inferred from first frame

**Location**: `integrator.py:3919-3941`.

```python
# waxs_bsx_ref is the bsx position where the mask polygon was
# drawn (typically arc ≈ 0°).  If the scan started at a different
# arc angle the first-frame bsx will be offset and we must NOT
# use it as the reference.
```

**Symptom**: WAXS beamstop polygon mis-positioned for arc-scanned WAXS data, even when `_BSX_PER_ARC_DEG` is correct.

**Why**: The reference position for the mask polygon is taken from the first scan frame's `waxs_bsx`. If the scan starts at non-zero arc, the first-frame bsx is already offset by `bsx_per_arc_deg * arc[0]`. Using it as the reference would shift every subsequent frame's mask by that constant.

**Fix** (already in code, but note): the integrator either fits the bsx-vs-arc slope (when arc is scanned) or extrapolates back to arc=0 using the known `_BSX_PER_ARC_DEG`. If you need to override:
```python
result = reduce_smi_combined(
    ...,
    waxs_kwargs={"waxs_bsx_ref": <known_arc0_bsx_mm>},
)
```

## 15. Frozen dataclass mutation

**Location**: `integrator.py:2085` (`@dataclass(frozen=True)`).

**Symptom**: `result.line_cuts = {...}` raises `FrozenInstanceError`.

**Why**: `CombinedReductionResult` is intentionally frozen so its top-level fields can't be reassigned by accident. But the optional derived stages (`apply_line_cuts`, `apply_peak_fits`, `apply_virtual_axes`) need to attach new products.

**Fix**: They use `object.__setattr__` to bypass the frozen check:
```python
# linecuts.py:291
object.__setattr__(result, "line_cuts", out)
```
- If you write a new derived stage, follow the same pattern.
- Don't reassign top-level fields (`merged_iq`, `saxs`, etc.) post-construction.

## 16. Solid-angle correction is global, not per-detector

**Location**: `integrator.py:3987` (passed to integrate_waxs), `integrator.py:3877` (passed to integrate_saxs).

**Symptom**: SAXS and WAXS show inconsistent intensity scaling in the merged output. The seam at q ≈ 1-2 nm⁻¹ has a discontinuity.

**Why**: `solid_angle_correction` is a single bool that's applied to **both** detectors with the same toggle. There's no `solid_angle_correction_saxs` / `solid_angle_correction_waxs`.

**Fix**: For physically-correct merged data, leave `solid_angle_correction=True` (default). If you must disable it for one detector, you'll need to call `integrate_saxs` and `integrate_waxs` directly with different toggles, then merge manually with `merge_q_chi_weighted`.

## 17. SAXS `q_abs` 2-D, WAXS `q_abs` 3-D

**Location**: `integrator.py:3177-3186` (SAXS), `integrator.py:3490` (WAXS).

**Symptom**: Code that consumes `result.saxs["ds"]["q_abs"]` and tries the same on `result.waxs["ds"]["q_abs"]` breaks (`q_abs` doesn't have a `frame` dim for SAXS).

**Why**: SAXS geometry is fixed across frames in a transmission scan, so q-maps are stored once as 2-D `(row, col)`. WAXS arc rotates, so q-maps are 3-D `(frame, row, col)`.

**Fix**: Branch on detector:
```python
ds = result.saxs["ds"]
q_abs = ds["q_abs"].values  # shape (1679, 1475), per-frame share
# vs
ds = result.waxs["ds"]
q_abs = ds["q_abs"].values  # shape (n_frames, 195, 619)
q_first = ds["q_abs"].isel(frame=0).values
```

## 18. Single-detector merged_qchi passthrough

**Location**: `integrator.py:1755-1759`.

```python
if saxs_qchi is None:
    saxs_qchi = _empty_qchi_like(waxs_qchi)
if waxs_qchi is None:
    waxs_qchi = _empty_qchi_like(saxs_qchi)
```

**Symptom**: Code that expects `result.merged_qchi` to be `None` for SAXS-only or WAXS-only scans gets a fully-populated Dataset instead.

**Why**: `merge_q_chi_weighted` substitutes a zero-counts Dataset for the missing detector. The merge becomes effectively a regrid of the present detector. This is intentional — downstream code can always use `merged_qchi` regardless of which detectors fired.

**Fix**: Either branch on `result.saxs is None or result.waxs is None`, or just use `merged_qchi` directly (it's correct).

## 19. `target_file_name` cached as bytes, not str

**Location**: `loader.py:484-486`, `loader.py:2090` (decode at parse time).

```python
primary_grp.create_dataset(
    "target_file_name",
    data=np.asarray(_names, dtype=object),
    dtype=h5py.string_dtype(encoding="utf-8"),
)
```

**Symptom**: After cache hit, `target_file_name` comes back as `bytes` (`b"sample_2450eV"`), not `str`. String operations (`.startswith`, regex) raise `TypeError`.

**Why**: HDF5 variable-length string columns return bytes in h5py by default. Tiled's primary stream returns str. The `_coerce_str` helper at `loader.py:~480` handles both, but only at specific call sites. If you read from cache directly, you must decode.

**Fix**: When reading `target_file_name` from cache yourself:
```python
names = _read_cached_primary_field(cache_path, "target_file_name")
names = [n.decode("utf-8", "replace") if isinstance(n, bytes) else str(n) for n in names]
```

## 20. Energy in eV vs keV

**Location**: throughout. `SAXSGeometry.energy_ev` is in eV; `peekAtMd` returns `energy_kev`.

**Symptom**: λ computed from energy is off by 10³.

**Why**: The codebase is mostly in eV (SI-adjacent), but the user-facing helpers (`peekAtMd`, `WAXSCalibration`) use keV. Cross-module conversions are silent.

**Fix**: Always convert at the boundary. The standard formula:
```python
hc_ev_nm = 1239.84193  # eV·nm
wavelength_nm = hc_ev_nm / energy_ev
wavelength_m = wavelength_nm * 1e-9
```
Don't mix `energy_ev / 1000.0` with `wavelength_m` casually.

## 21. Auto-disabled detector_ds is silent without explicit print

**Location**: `integrator.py:2798-2802`, `integrator.py:2846-2850`.

```python
if n_frames <= _DETECTOR_DS_AUTO_MAX_FRAMES:
    build_detector_ds = True
elif build_detector_ds is None:
    print(f"... {n_frames} frames (> {_DETECTOR_DS_AUTO_MAX_FRAMES}) to bound memory; ...")
    build_detector_ds = False
```

**Symptom**: `result.saxs["ds"]` is None despite passing nothing to `build_detector_ds`. The auto-disable only prints when `build_detector_ds is None`, but the user expected on-by-default.

**Why**: "Off when n > 50" is the new contract; pre-50-frame scans get on-by-default. Tests and docs may pre-date this change.

**Fix**: Pass explicit `build_detector_ds=True` when you need it. See gotcha #8 for memory implications.

## Summary table

| # | File:line | Class |
|---|---|---|
| 1 | `loader.py:1140` | Silent default |
| 2 | `defaults.py:303` | Round-trip data loss |
| 3 | `integrator.py:587` | Undocumented constant |
| 4 | `integrator.py:3925` | Duplicate constant |
| 5 | `loader.py:_SAXS_BEAM_*` | Double-counting |
| 6 | `loader.py:699` | Exception swallow |
| 7 | `loader.py:2927` | Eager import |
| 8 | `integrator.py:2785` | Silent auto-disable |
| 9 | `integrator.py:169` | Stale cache |
| 10 | `loader.py:2615` | Wrong abstraction |
| 11 | `integrator.py:_WAXS_QMAP_CACHE` | Stale cache, no clear |
| 12 | `loader.py:38` | Sticky None |
| 13 | `loader.py:423` | Mid-write corruption |
| 14 | `integrator.py:3919` | Reference inference |
| 15 | `integrator.py:2085` | Frozen dataclass |
| 16 | `integrator.py:3987` | Global toggle |
| 17 | `integrator.py:3177/3490` | Inconsistent shape |
| 18 | `integrator.py:1755` | Hidden passthrough |
| 19 | `loader.py:484` | bytes vs str |
| 20 | mixed | eV vs keV units |
| 21 | `integrator.py:2798` | Silent default |
