# Reduction outputs

This file documents the structure of `CombinedReductionResult` and `GIReductionResult` (the public outputs of `reduce_smi_combined` and `reduce_smi_gi`), plus the per-detector dicts and merged datasets they wrap.

## `CombinedReductionResult` (`integrator.py:2085-2199`)

Frozen dataclass returned from `reduce_smi_combined` for transmission scans.

```python
@dataclass(frozen=True)
class CombinedReductionResult:
    uid: str
    scan_info: dict[str, Any]
    saxs: dict[str, Any] | None
    waxs: dict[str, Any] | None
    merged_qchi: xr.Dataset | None
    merged_iq: xr.Dataset | None
    per_frame_iq: xr.Dataset | None = None
    timing: dict[str, float] | None = None
    geometry: str = "transmission"
    incident_angle_deg: float = 0.0
    # Optional, attached by smi_tiled.derived stages:
    per_frame_qchi: dict[str, xr.Dataset] | None = None
    line_cuts: dict[str, xr.Dataset] | None = None
    peak_fits: xr.Dataset | None = None
```

| Field | Set when | Purpose |
|---|---|---|
| `uid` | always | Tiled run UID. |
| `scan_info` | always | `infer_detectors_and_steps` output: `n_frames`, `sample_name`, `step_candidates`, scan_id, etc. |
| `saxs` | SAXS detector present | Per-detector intermediate dict (see below). |
| `waxs` | WAXS detector present | Per-detector intermediate dict (see below). |
| `merged_qchi` | both detectors OR single-detector passthrough | Count-weighted (q,χ) merge. |
| `merged_iq` | merged_qchi exists | Azimuthally-averaged 1-D profile. |
| `per_frame_iq` | per-frame data exists | Per-step I(q) on the merged grid. |
| `timing` | always | Per-stage seconds (load, mask, saxs, waxs, merge, total). |
| `geometry` | always | `"transmission"` or `"grazing_incidence"`. |
| `incident_angle_deg` | always | Resolved αᵢ. 0 for transmission. |
| `per_frame_qchi` | `build_frame_qchi=True` (default) | `{"saxs": Dataset, "waxs": Dataset}`, each with dims `(frame, q, chi)`. |
| `line_cuts` | passed `line_cuts=[...]` to reduce_smi_combined | `{name: Dataset}` of 1-D cuts. |
| `peak_fits` | passed `peak_fits=[...]` | `Dataset` with dims `(peak, frame)`. |

## Per-detector dict: `saxs` / `waxs`

Built by `integrate_saxs` (`integrator.py:2884`, returns at line 3193) and `integrate_waxs` (line 3200, returns at 3516). Both have the same shape:

```python
{
    "q_chi": xr.Dataset,           # always set, scan-summed
    "iq": xr.Dataset,              # always set, scan-summed
    "iq_frames": xr.Dataset,       # always set, per-frame
    "q_chi_frames": xr.Dataset | None,  # set unless build_frame_qchi=False
    "ds": xr.Dataset | None,       # set if build_detector_ds=True
}
```

### `q_chi` (scan-summed)

`integrator.py:1591-1597`. Built by `_qchi_and_iq` after pixel-summing into `accum_I, accum_N` over all frames.

```python
xr.Dataset(
    {"intensity": (("q", "chi"), mean_I),
     "counts":    (("q", "chi"), accum_N)},
    coords={"q": q_grid, "chi": chi_grid},
)
```

| Variable | Dims | Meaning |
|---|---|---|
| `intensity` | `(q, chi)` | `accum_I / accum_N` where `accum_N > 0`, else NaN. |
| `counts` | `(q, chi)` | Sum of valid (mask & finite) pixels in each (q, χ) bin. |

`q` units: nm⁻¹ for SAXS and WAXS (radial wavevector). `chi`: degrees, range depends on detector.

### `iq` (scan-summed)

`integrator.py:1605-1607`. Azimuthally collapsed.

```python
{"I":      ("q",), "counts": ("q",)}
# I_1d = sum(intensity * counts, axis=chi) / sum(counts, axis=chi)
```

### `iq_frames` (per-frame)

Stacked from per-frame `_qchi_and_iq` results. Dims: `(frame, q)`. Variables: `I, counts`.

### `q_chi_frames` (per-frame)

Three modes depending on `build_frame_qchi` and `frame_qchi_store`:

| Mode | Trigger | Storage |
|---|---|---|
| In-memory stack | `build_frame_qchi=True` (default), no `frame_qchi_store` | `xr.Dataset` with dims `(frame, q, chi)`, fully eager. |
| Lazy zarr | `frame_qchi_store="path/to/qchi.zarr"` | `xr.Dataset` backed by dask + zarr — peak memory ~ 1 frame. |
| Skipped | `build_frame_qchi=False` | `None` |

The lazy path uses `_ZarrFrameWriter`. Use it for scans > 50 frames where in-memory `(N, n_q, n_chi)` would blow up RAM (e.g. 200 frames × 2000 × 360 × 8 B ≈ 1.2 GB).

### `ds` (detector-space)

Only built if `build_detector_ds=True`. Auto-disabled silently for scans > `_DETECTOR_DS_AUTO_MAX_FRAMES = 50` (`integrator.py:2785`).

SAXS shape (`integrator.py:3177-3186`):
```python
xr.Dataset(
    {"intensity":   (("frame", "row", "col"), images),
     "q_abs":       (("row", "col"),          q2d),       # 2-D, shared
     "q_horizontal":(("row", "col"),          qh2d),
     "q_vertical":  (("row", "col"),          qv2d),
     "mask":        (("frame", "row", "col"), masks)},
    coords={"frame": np.arange(n_frames)},
)
```

WAXS shape (`integrator.py:3484-3509`):
```python
xr.Dataset(
    {"intensity":   (("frame", "row", "col"), images),
     "q_abs":       (("frame", "row", "col"), qabs_frames),  # 3-D, varies per arc
     "q_horizontal":(("frame", "row", "col"), qh_frames),
     "q_vertical":  (("frame", "row", "col"), qv_frames),
     "mask":        (("frame", "row", "col"), masks),
     "waxs_arc":    (("frame",),              arc_angles)},
    coords={"frame": np.arange(len(arc_angles))},
)
```

Difference: SAXS q-maps are 2-D (geometry is fixed), WAXS q-maps are 3-D (arc angle changes per frame).

## `merged_qchi` (`integrator.py:1742-1813`)

Built by `merge_q_chi_weighted(saxs_qchi, waxs_qchi)`. Steps:

1. Determine common q grid: `np.linspace(min(saxs_q, waxs_q), max(saxs_q, waxs_q), n_q)`.
2. Same for `chi_grid`.
3. Nearest-neighbor regrid each detector onto the common grid (using `scipy.interpolate.RegularGridInterpolator`).
4. Count-weighted merge:
   ```
   merged_I = (s_I * s_N + w_I * w_N) / (s_N + w_N)   where (s_N + w_N) > 0
            = NaN                                       elsewhere
   ```

```python
xr.Dataset(
    {"intensity":      (("q", "chi"), merged_I),       # the merged answer
     "counts":         (("q", "chi"), s_N + w_N),
     "saxs_intensity": (("q", "chi"), s_I_interp),      # diagnostic
     "saxs_counts":    (("q", "chi"), s_N_interp),
     "waxs_intensity": (("q", "chi"), w_I_interp),
     "waxs_counts":    (("q", "chi"), w_N_interp)},
    coords={"q": q_grid, "chi": chi_grid},
)
```

The `saxs_*` / `waxs_*` variables let you split out which detector contributed where — useful for diagnosing seam artifacts in the overlap region (~ q = 1-2 nm⁻¹ at typical SMI geometry).

### Single-detector passthrough

If only one detector is present, `merge_q_chi_weighted` calls `_empty_qchi_like` to build a zero-counts dataset for the missing detector, then runs the same merge. Result: `merged_qchi` is essentially the present detector's data, regridded.

## `merged_iq` (`integrator.py:1816-1858`)

Built by `merge_iq_profiles(merged_qchi, saxs_iq, waxs_iq)`.

Procedure: azimuthally collapse `merged_qchi`, then attach the per-detector 1-D I(q) interpolated onto the merged q grid (for diagnostics).

```python
xr.Dataset(
    {"I":      ("q",), counts: ("q",),     # the merged answer
     "saxs_I": ("q",),                      # interpolated, diagnostic
     "waxs_I": ("q",)},                     # interpolated, diagnostic
    coords={"q": q_grid},
)
```

`I = sum(merged_intensity * merged_counts, axis=chi) / sum(merged_counts, axis=chi)`.

## `per_frame_iq` (`integrator.py:1861-1976`)

Per-step I(q) on the merged q grid. Dims `(frame, q)`.

```python
{
    "I":       (("frame", "q"), merged_I_2d),    # count-weighted merge per frame
    "saxs_I":  (("frame", "q"), saxs_I_2d),      # interpolated
    "waxs_I":  (("frame", "q"), waxs_I_2d),
    # plus any per-frame primary scalars:
    "pil2M_motor_x": (("frame",), motor_values),
    "energy_energy": (("frame",), energy_values),
    "fn:eV":         (("frame",), parsed_from_target_file_name),
    ...
}
```

Per-frame primary scalars (motor positions, energy, ring current, etc.) are attached as 1-D data variables on `frame`. The `fn:*` columns come from `apply_virtual_axes` parsing structured strings like `target_file_name` (see Virtual axes below).

## `to_dataarray` (`integrator.py:2131-2199`)

Helper to extract a single-variable `xr.DataArray` view of a merged product, for compatibility with PyHyperScattering accessors (`da.rsoxs.*`, `da.fit.*`).

```python
da_iq = result.to_dataarray("merged_iq")           # I on (q,)
da_qchi = result.to_dataarray("merged_qchi")       # intensity on (q, chi)
da_pf = result.to_dataarray("per_frame_iq", "saxs_I")  # saxs_I on (frame, q)
```

Default variable: `"I"` for 1-D, `"intensity"` for 2-D.

Attached attrs: `uid, scan_id, sample_name, geometry, incident_angle_deg, source` plus the original Dataset's attrs.

Raises `ValueError` if:
- key is not one of `merged_iq | merged_qchi | per_frame_iq`,
- the dataset is `None` (e.g. requesting `merged_qchi` for a SAXS-only scan that wasn't merged),
- `variable` not in dataset.

## `GIReductionResult` (`integrator.py:2202-2271`)

Returned from `reduce_smi_gi` for grazing-incidence WAXS scans.

```python
@dataclass(frozen=True)
class GIReductionResult:
    uid: str
    sample_name: str
    scan_motor: str                # e.g. "piezo_th"
    scan_motor_values: np.ndarray
    alpha_i_deg: np.ndarray         # per-frame
    alpha_i_source: str             # description, e.g. "stage_th + piezo_th"
    qxy_grid: np.ndarray            # 1-D bin centres (nm⁻¹)
    qz_grid: np.ndarray             # 1-D bin centres (nm⁻¹)
    frames: list[np.ndarray]        # per-frame I(qxy, qz)
    summed: np.ndarray              # mean over frames
    q_chi_frames: xr.Dataset | None  # dims (frame, qxy, qz)
    summed_ds: xr.Dataset | None     # dims (qxy, qz), vars intensity, counts
    timing: dict[str, float] | None
    line_cuts: dict[str, xr.Dataset] | None = None
    peak_fits: xr.Dataset | None = None
```

GI uses sample-frame (qxy, qz) instead of (q, chi). The lab-to-sample rotation is `lab_to_sample_frame` (`integrator.py:2278-2303`):

```
qx_s = qx_lab
qy_s = qy_lab * sin(αᵢ) - qz_lab * cos(αᵢ)
qz_s = qy_lab * cos(αᵢ) + qz_lab * sin(αᵢ)
```

### Line cut helpers (built-in, no `apply_line_cuts` needed)

```python
qxy, I_qxy = result.line_cut_qxy(qz_center=0.5, qz_width=0.05, frame=None)
qz,  I_qz  = result.line_cut_qz(qxy_center=1.0, qxy_width=0.1,  frame=3)
```

`frame=None` uses the summed image; an integer indexes `self.frames`. The cut averages over the band; if the band is too narrow to contain any bin, returns NaN.

## Optional derived products

These are **opt-in** — they're attached only if you pass the relevant kwarg to `reduce_smi_combined`.

### `per_frame_qchi` (`integrator.py:4035-4042`)

Dict `{"saxs": ds, "waxs": ds}` promoted from `saxs_result["q_chi_frames"]` and `waxs_result["q_chi_frames"]`. Always populated when `build_frame_qchi=True` (default) and the relevant detector has frames.

This is **separate from** `result.saxs["q_chi_frames"]` only by reference layout — same Dataset, exposed at the top level for downstream consumers (line cuts on `saxs_qchi` / `waxs_qchi` targets, upload schemas, GI pipeline).

### `line_cuts` (`smi_tiled/derived/linecuts.py:247`)

Triggered by passing `line_cuts=[LineCutSpec(...), ...]` to `reduce_smi_combined`. Attached as `result.line_cuts: dict[name, xr.Dataset]`.

`LineCutSpec` (`linecuts.py:84`):
```python
LineCutSpec(
    kind="h" | "v",                  # h: I(x) over y-band; v: I(y) over x-band
    center=0.5,                       # band centre
    width=0.05,                       # band full width
    target="merged_qchi" | "saxs_qchi" | "waxs_qchi" | "qxy_qz",
    name="optional_stable_id",
)
```

Per-frame is the default — `apply_line_cuts(result, cuts, per_frame=True)` returns a Dataset with `frame` dim. Pass `per_frame=False` to collapse via `nanmean(axis=frame)` first.

Targets:
- `merged_qchi` — single 2-D image (`merged_qchi.intensity`)
- `saxs_qchi` / `waxs_qchi` — per-frame stacks (`per_frame_qchi[detector].intensity`)
- `qxy_qz` — GI per-frame stack (`q_chi_frames` of GI result)

### `peak_fits` (`smi_tiled/derived/peakfit.py:434`)

Triggered by passing `peak_fits=[PeakDef(...), ...]` to `reduce_smi_combined`. Attached as `result.peak_fits: xr.Dataset`.

`PeakDef` (`peakfit.py:47`):
```python
PeakDef(
    name="primary_peak",
    q_min=0.4, q_max=0.6,
    model="gaussian" | "lorentzian",
    baseline="none" | "linear",
    link="independent" | "linked" | "tracked",
    bg_factor=2.0,           # widens fit window for baseline anchoring
)
```

Output:
```python
xr.Dataset(
    {"amplitude": (("peak", "frame"), …),
     "center":    (("peak", "frame"), …),
     "fwhm":      (("peak", "frame"), …),
     "area":      (("peak", "frame"), …),
     "success":   (("peak", "frame"), bool),
     "peak_key":  (("peak",),         …)},  # SHA1 hash for cache-staleness
    coords={"peak": [name1, name2, …]},
)
```

Quality gates (`peakfit.py:42-44`): `MIN_SNR=3.0`, `MIN_R2=0.2`, `WIDTH_BOUND_TOL=0.97`. A fit fails (`success=False`) if it doesn't pass.

### Virtual axes (`smi_tiled/derived/virtual_axes.py`)

`apply_virtual_axes` runs **by default** (with `VirtualAxesConfig()`) at the end of `reduce_smi_combined` (`integrator.py:4062-4065`). It parses structured per-frame strings (e.g. `target_file_name`) into numeric `fn:*` columns attached to `result.per_frame_iq`.

Example:
- Source: `target_file_name = "Lucas_sample2_pos1_2450.00eV_ai0.50_degC100.0"`
- Parsed (per-frame): `{"sample": 2, "pos": 1, "eV": 2450.0, "ai": 0.5, "degC": 100.0}`
- Result: `per_frame_iq["fn:eV"]`, `per_frame_iq["fn:ai"]`, etc.

Filter: a column is kept only if its non-NaN fraction ≥ `min_fill` (default 0.5).

To disable: pass `virtual_axes=VirtualAxesConfig(enabled=False)`.

## Quick lookup table

| Want… | Read |
|---|---|
| Final 1-D I(q) over full q range | `result.merged_iq["I"]` |
| Final 2-D I(q, χ) | `result.merged_qchi["intensity"]` |
| SAXS-only 1-D I(q) | `result.saxs["iq"]["I"]` |
| WAXS-only 1-D I(q) | `result.waxs["iq"]["I"]` |
| Per-frame I(q) | `result.per_frame_iq["I"]` (dims `(frame, q)`) |
| Per-frame SAXS (q,χ) | `result.per_frame_qchi["saxs"]["intensity"]` |
| Frame motor positions | `result.per_frame_iq["pil2M_motor_x"]` etc. |
| Energy axis (parsed) | `result.per_frame_iq["fn:eV"]` |
| Detector-space images | `result.saxs["ds"]["intensity"]` (only if `build_detector_ds=True`) |
| Reduction timing | `result.timing` |
| Run scan_id | `result.scan_info["scan_id"]` |
| For PyHyperScattering | `result.to_dataarray("merged_iq")` |
