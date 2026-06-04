# Masks (SAXS + WAXS)

SMI uses **polygon-based masks** for both detectors, stored as JSON. This file documents the schema, the build pipeline, the dynamic (per-frame) masks, and the round-trip pitfalls.

## Schema (nested form, used by SAXS)

`src/smi_tiled/data/masks/pil2M_mask_polygons.json` — bundled SAXS mask:

```json
{
  "image_shape": [1679, 1475],
  "static_regions": {
    "gap_col_1": [[col, row], [col, row], ...],
    "gap_col_2": [...],
    "gap_row_1": [...],
    ...
  },
  "beamstops": {
    "rod": {
      "polygon_offsets_from_beam": [
        [-14, -459], [14, -459], [14, 575], [-14, 575]
      ],
      "polygon": [...],            // legacy; ignored when offsets present
      "x_motor_key": "saxs_beamstop_x_rod_user_setpoint",
      "y_motor_key": "saxs_beamstop_y_rod_user_setpoint",
      "reference_mm": {"x": 6.8, "y": 3.44},
      "pixels_per_mm": {"x": 5.81, "y": 5.81}
    },
    "pin": {
      "polygon_offsets_from_beam": [
        [[27.0, 0.0], [25.9, 6.8], ..., [25.9, -6.8]],   // pin disc (20 verts)
        [[-12, 18], [20, 18], [11, 676], [-20, 676]],    // pin rod (extends downward)
        [[-155, 572], ..., [-155, 676]],                 // ancillary mask
        [[-95, 560], ..., [-95, 676]]                    // ancillary mask
      ],
      "polygon": [...],
      ...
    }
  }
}
```

Coordinates are in **raw detector indexing** `(col, row)` (not (row, col)!).

`polygon_offsets_from_beam` is the **preferred** form — offsets are relative to the resolved beam center, so motor-driven beam-center drift is automatically tracked. Each beamstop entry can carry **either**:
- A single polygon (`[[col, row], ...]`).
- A list of polygons (`[[[c, r], ...], [[c, r], ...], ...]`).

The pin entry uses the list-of-polygons form because the pin beamstop has a disc plus an L-shaped support structure plus auxiliary masks.

## Schema (flat form, used by WAXS)

`src/smi_tiled/data/masks/900KW_mask_polygons.json` — bundled WAXS mask is **flat**:

```json
{
  "image_shape": [195, 619],
  "<region_name>": [[col, row], ...],
  "<region_name>": [...],
  "beamstop": [...]                  // any key containing "beamstop"
}
```

`load_mask_polygons` (`defaults.py:313-390`) routes any key containing `"beamstop"` (case-insensitive) into `beamstops`, all other keys into `static_regions`. The result is normalized to the nested schema regardless of input.

## load_mask_polygons (`defaults.py:313`)

Always returns:
```python
{
    "image_shape": [rows, cols] | None,
    "static_regions": {name: [[col, row], ...], ...},
    "beamstops":      {name: [[col, row], ...], ...},
}
```

For nested-schema beamstops with `{"polygon": [...], ...}` wrappers, the wrapper is **unwrapped** — only the polygon survives.

> **CRITICAL ROUND-TRIP GOTCHA**: `load_mask_polygons` calls `_coerce_polygon` (`defaults.py:303`), which **drops** any non-polygon keys in a beamstop entry. That includes `polygon_offsets_from_beam`, `x_motor_key`, `reference_mm`, `pixels_per_mm`, etc. So `save_mask_polygons(load_mask_polygons(path), path2)` is **not lossless** for the SAXS schema — it produces a flat-polygon-only file that cannot be motor-corrected. Fix: edit JSON directly, or build a custom serializer.

## save_mask_polygons (`defaults.py:393`)

Writes back using the nested schema, but only `static_regions` + `beamstops` (each as bare polygons). Does not preserve `polygon_offsets_from_beam` etc. on round-trip.

## resolve_mask_path (`defaults.py:91`)

Resolves a user-supplied mask argument:
- `None` → bundled default for `detector` (`"saxs"` → pil2M, `"waxs"` → 900KW).
- Absolute path / existing file → returned as-is.
- Bare filename matching the bundled basename → returns the bundled path.
- Anything else → returned unchanged (caller will see FileNotFoundError later).

## Mask polarity

**Throughout SMI tiled code: `True = valid pixel, False = masked-out.`**

Inside `polygons_to_mask` (`integrator.py:546`):
```python
mask = np.ones(shape, dtype=bool)   # start: everything valid
for poly in polygons:
    if not poly:
        continue
    cols = np.array([p[0] for p in poly], dtype=float)
    rows = np.array([p[1] for p in poly], dtype=float)
    rr, cc = skpoly(rows, cols, shape=shape)
    mask[rr, cc] = False             # mark polygon interiors as INVALID
return mask
```

The `skpoly` import is `from skimage.draw import polygon as skpoly`.

## SAXS mask build: `make_saxs_mask_from_dict` (`integrator.py:612`)

```python
def make_saxs_mask_from_dict(
    image_shape: tuple[int, int],
    mask_spec: dict,
    active_beamstop: str = "rod",
    beamstop_pos_mm: dict | None = None,
    beam_center_px: tuple[float, float] | None = None,
) -> np.ndarray:
```

Returns a static (non-per-frame) mask. Operation order:

1. Collect all `static_regions[*]` polygons.
2. Look up `beamstops[active_beamstop]` (typically `"rod"` or `"pin"`).
3. **Beamstop dict cases**:
   - **`polygon_offsets_from_beam`** present → for each offset polygon, shift by `(beam_center_col + dc, beam_center_row + dr)`. Requires `beam_center_px=(row, col)` arg.
     - List-of-lists handling: if first element is a 2-vertex (e.g. `[27, 0]`), it's one polygon. If first element is a list-of-vertices, it's multiple polygons.
   - **`polygon`** present (legacy) → use absolute coords. If `beamstop_pos_mm` provided, shift by `(cur_motor - reference_mm) * pixels_per_mm`.
4. **Bare polygon**:
   - `bs` is a list — could be a single polygon `[[c,r], ...]` or a list-of-polygons. Detected by checking if `bs[0][0]` is a list/tuple.
   - Used as-is (no beam-center anchoring or motor shift).
5. **Empty / missing** → no beamstop polygon added (warning if `polygon_offsets_from_beam` would have been used but `beam_center_px` missing).
6. Return `polygons_to_mask(image_shape, polys)`.

> **WARNING — silent default to `"rod"`**: `active_beamstop` defaults to `"rod"` everywhere — `loader.py:1140`, `integrator.py:615`, `integrator.py:770`. If `pil2M_active_beamstop` is missing in baseline, `resolve_saxs_geometry` returns `"rod"` silently. A pin-beamstop scan integrates with the wrong mask and gets a wrong q_min.

## SAXS mask build: `make_saxs_mask_from_spec` (`integrator.py:767`)

Thin wrapper around `make_saxs_mask_from_dict`. `mask_path` may be either a file path or a parsed dict; `_resolve_mask_spec` (`integrator.py:598`) handles both.

## WAXS mask build: `make_mask_for_angle` (`integrator.py:573`)

WAXS masks are **per-arc-angle** because the beamstop and the open-edge boundary move with `waxs_arc`.

```python
def make_mask_for_angle(
    image_shape_raw: tuple[int, int],   # (195, 619), raw orientation
    static_regions: dict,
    beamstop_region: list,              # single polygon
    waxs_bsx: float,                    # current waxs_bsx motor (mm)
    waxs_bsx_ref: float,                # reference waxs_bsx (mm)
    pixel_size_mm: float = 0.172,
    rotation_k: int = 3,
    include_beamstop: bool = True,
) -> np.ndarray:
```

Math:
```python
bs_shift_mm = waxs_bsx - waxs_bsx_ref
bs_shift_px = (bs_shift_mm / pixel_size_mm) * 1.088   # *** EMPIRICAL FUDGE FACTOR ***
shifted = shift_polygon(beamstop_region, dx_px=0.0, dy_px=bs_shift_px)
raw_mask = polygons_to_mask(image_shape_raw, all_polys)
mask_rot = rotate_image_and_mask(raw_mask, k=rotation_k)
return mask_rot
```

> **CRITICAL — empirical 1.088 fudge factor at integrator.py:587**: This factor was determined empirically from earlier mask alignment work (likely needed because of differential WAXS pixel-size mismatch with raw vs rotated orientation). It has no derivation, no calibration scan documented; if mask alignment changes after a detector swap or remount, this factor must be re-determined.

## Dynamic SAXS shadow mask: `_make_waxs_shadow_mask` (`integrator.py:1219`)

When the WAXS arc detector is in the SAXS beam path (small SDD scans), it physically blocks part of the SAXS detector. The blocked region depends on `waxs_arc`.

The shadow boundary is a single column index per frame:
```python
boundary_col = beam_center_col_px + beam_visible_offset_px +
               (waxs_arc - beam_visible_deg) / (clear_edge_deg - beam_visible_deg) *
               (clear_col - beam_center_col_px - beam_visible_offset_px)
```

where:
- `beam_visible_deg = 14.5` — at this arc angle, the WAXS edge sits at the beam.
- `clear_edge_deg = 18.0` — at this arc angle, the WAXS edge has fully cleared the SAXS detector.
- `clear_col = nx - 1 - edge_margin_px` — the rightmost detector column.

Pixels with `col <= boundary_col[frame]` are kept; the shadow blocks `col > boundary_col[frame]`. Frames where `waxs_arc >= clear_edge_deg` are fully unmasked.

For large scans (n_frames > 50), returns a `_LargeAreaMaskLazy` (`integrator.py:1264`) that materializes a single-frame mask on `__getitem__`. Avoids holding `(n_frames, ny, nx)` boolean array in memory.

> **Threshold to remember**: 50 frames. Above that, the lazy path kicks in and `mask[i]` returns a freshly-built ndarray each call.

## Aperture mask: `_make_aperture_mask` (`integrator.py:1290`)

Models the WAXS arc detector blocking part of the **SAXS angular range**. Uses a q-cutoff (which is angle for fixed wavelength; q being wavelength-independent is convenient for cross-energy comparison).

Default q_cutoff is determined from AGB d-spacing (`D_nm = 5.838`):
```python
q_cutoff = (2π / D_nm) * agbh_ring_order * (1 + q_margin_fraction)
```
e.g. for `agbh_ring_order=5` and `q_margin_fraction=0.08`: `q_cutoff = (2π/5.838) × 5 × 1.08 ≈ 5.81 nm^-1 ≈ 0.581 Å^-1`.

If `q_cutoff > detector_max_q`, no occlusion is applied (the cutoff falls beyond the detector corner — typical for long-SDD scans where the WAXS arc is far below the SAXS detector's q-range). Otherwise:
```python
return (q_abs <= q_cutoff)
```

## Combined dynamic mask: `make_saxs_large_area_masks` (`integrator.py:1329`)

Returns `(combined, shadow_only, aperture_only)` where:
- `combined = shadow & aperture` (per-frame).
- For lazy mode (large scans), returns `_LargeAreaMaskCombined` (`integrator.py:1355`).

## Mask combination at integration time

Inside `integrate_saxs` (`integrator.py:2984-3000`):
```python
base_valid = np.isfinite(q2d) & np.isfinite(chi_deg_2d)
if mask_use is not None:                        # static mask from JSON
    base_valid &= mask_use
if dynamic_mask is not None:                    # large-area shadow + aperture
    valid_per_frame = base_valid & dynamic_mask[frame_idx]
```

For WAXS (`integrator.py:3300-3400`), the per-frame `make_mask_for_angle` result is intersected with the per-pixel `qmap`-finite mask.

## Mask debug recipe

```python
import numpy as np
import matplotlib.pyplot as plt
from smi_tiled.defaults import default_saxs_mask_path, load_mask_polygons
from smi_tiled.integrator import make_saxs_mask_from_dict
from smi_tiled.loader import TiledSMISWAXSLoader, resolve_saxs_geometry

loader = TiledSMISWAXSLoader()
run = loader._get_run("<uid>")
geo = resolve_saxs_geometry(run)

# Use the raw (un-normalized) JSON for full beamstop info
import json
with open(default_saxs_mask_path()) as f:
    spec = json.load(f)

mask = make_saxs_mask_from_dict(
    image_shape=(1679, 1475),
    mask_spec=spec,
    active_beamstop=geo.active_beamstop,
    beam_center_px=(geo.beam_center_row_px, geo.beam_center_col_px),
)
plt.imshow(mask, origin="lower")  # SAXS display orientation: flipud
plt.axhline(geo.beam_center_row_px, color="r", linewidth=0.5)
plt.axvline(geo.beam_center_col_px, color="r", linewidth=0.5)
plt.title(f"active={geo.active_beamstop}, masked={(~mask).sum()} px")
plt.show()
```

Verify:
- Number of masked pixels for rod (~14 × 1034 ≈ 14k) vs pin (~22²π + rod ≈ 18k+).
- All masked pixels lie within the polygon boundaries.
- The chi sector at +90° (right of beam) is empty if pin polygon offset is +col-direction.

## Display orientation gotcha

`orient_frame_for_display(arr, "saxs")` returns `np.flipud(arr)` (`defaults.py:248-258`). When overlaying a mask polygon on a displayed image:
- The mask is in raw `(row, col)` coords.
- The displayed image has rows flipped: `display_row = nrows - 1 - raw_row`.
- A polygon vertex `(col_raw, row_raw)` displays at `(x, y) = (col_raw, raw_h - row_raw)` — see `orient_polygon_xy` (`defaults.py:261`).

WAXS uses `np.fliplr(np.rot90(arr, k=3))` ≈ transpose, so the polygon orient becomes `(x, y) = (row_raw, col_raw)`.

## Files of interest

| File | Purpose |
|---|---|
| `src/smi_tiled/data/masks/pil2M_mask_polygons.json` | Bundled SAXS static + beamstop mask (nested schema) |
| `src/smi_tiled/data/masks/900KW_mask_polygons.json` | Bundled WAXS static + beamstop mask (flat schema) |
| `src/smi_tiled/defaults.py:303-410` | `_coerce_polygon`, `load_mask_polygons`, `save_mask_polygons` |
| `src/smi_tiled/integrator.py:546-714` | `polygons_to_mask`, `shift_polygon`, `make_mask_for_angle`, `make_saxs_mask_from_dict`, `make_saxs_mask_from_spec` |
| `src/smi_tiled/integrator.py:1219-1370` | `_make_waxs_shadow_mask`, `_LargeAreaMaskLazy`, `_make_aperture_mask`, `make_saxs_large_area_masks`, `_LargeAreaMaskCombined` |
