# WAXS geometry (Pilatus 900KW arc detector)

The WAXS detector at SMI is a **3-panel folded arc** built from a single Pilatus 900KW rotated 90° and split logically into three panels separated by 7° kinks. The geometry is custom to SMI; pyFAI's flat-panel assumptions do **not** apply.

## Detector specs

- **Detector**: Pilatus 900KW (1 Pilatus module, 195 × 1475 raw).
- **Raw image shape**: `(195, 619)` in pre-rotation coords. After `np.fliplr(np.rot90(image, k=3))` (`integrator.py:493`), the working shape is `(619, 195)`. The integrator operates on the **rotated** image internally, so:
  - "rows" (`ny=619`) span the **arc curvature** (panel-to-panel direction).
  - "cols" (`nx=195`) span the panel **width** (radial direction; vertical on detector).
  - That's the opposite of the SAXS convention where rows are vertical.
- **Pixel size**: `0.172 mm × 0.172 mm` (same Pilatus module).
- **Panels**: 3, indexed 0/1/2.

## Default panel layout

(`loader.py:130-131`, `integrator.py:297-298`)

| Panel | Column range | Offset (deg) |
|---|---|---|
| 0 | `[0, 206)` | `-7.0°` |
| 1 | `[206, 413)` | `0.0°` (reference, central) |
| 2 | `[413, 619)` | `+7.0°` |

The **offset_deg** is the in-plane (about y, the vertical axis) tilt of each panel relative to the central panel. Adjacent panels meet at fold lines `col = 205.5` and `col = 412.5` (`integrator.py:409, 423`).

## WAXSCalibration dataclass

Defined in `integrator.py:290-340`. Defaults:

```python
@dataclass
class WAXSCalibration:
    energy_kev: float = 16.1
    sample_distance_mm: float = 273.0          # 2-3 cm WAXS, calibrated
    pixel_size_mm: float = 0.172
    beam_center_row: float = 217.0             # row in rotated coords
    beam_center_col: float = 319.0             # col in rotated coords
    panel_col_ranges: Tuple = ((0,206), (206,413), (413,619))
    panel_offsets_deg: Tuple = (-7.0, 0.0, 7.0)
    panel_row_shifts: Tuple   = (0.0, 0.0, 0.0)
    panel_col_shifts: Tuple   = (0.0, 0.0, 0.0)
    panel_delta_deg: Tuple    = (0.0, 0.0, 0.0)  # per-panel angle correction
    theta_zero_deg: float = 0.0                  # arc-angle zero offset
    sample_offset_x_mm: float = 0.0              # sample horizontal offset
    sample_offset_z_mm: float = 0.0              # sample along-beam offset
    beam_col_per_arc_deg: float = 0.0            # BC drift with arc rotation
    q_horizontal_sign: float = -1.0              # convention switch
    q_vertical_sign: float = -1.0
    rotation_k: int = 3
```

The cal is a stand-alone, tunable calibration object for `MultiPanelArcDetector`. It is constructed from a `WAXSGeometry` (and arc angle) inside `integrate_waxs` (`integrator.py:2432`).

## WAXSGeometry dataclass (loader-side)

`loader.py:1147-1180`:

```python
@dataclass
class WAXSGeometry:
    dist_m: float
    beam_center_row_px: float
    beam_center_col_px: float
    pixel_m: float = PILATUS_PIXEL_SIZE_M
    energy_ev: float = ...
    wavelength_m: float = ...
    theta_zero_deg: float = 0.0
    sample_offset_x_mm: float = 0.0
    sample_offset_z_mm: float = 0.0
    panels: tuple[WAXSPanelGeometry, ...] = ()
```

The loader's `WAXSGeometry` is the wire format from metadata resolution. Inside `integrate_waxs` it gets translated into a per-frame `MultiPanelArcDetector` with the panel-arc machinery.

> **NOTE — beam center is hardcoded** (`loader.py:1425-1436`): The ophyd configuration stores WAXS beam center in the **raw, pre-rotation** coordinate system. That's incompatible with the rotated-frame `MultiPanelArcDetector`. So `resolve_waxs_geometry` does **not** read beam center from baseline; it uses `_WAXS_DEFAULT_BEAM_ROW_PX = 217.0`, `_WAXS_DEFAULT_BEAM_COL_PX = 319.0`. Override via `beam_delta_row_px` / `beam_delta_col_px` only.

## MultiPanelArcDetector — the model

`integrator.py:355-486`. Constructor signature:
```python
MultiPanelArcDetector(
    image_shape: tuple[int, int],         # (ny, nx) of rotated image
    panel_specs: Sequence[PanelSpec],     # 3-element list
    wavelength_nm: float,
    pixel_size_mm: float = 0.172,
    sample_distance_mm: float = 300.0,
    beam_center_px: tuple[float, float] = (0.0, 0.0),   # (row, col) in rotated frame
    theta_zero_deg: float = 0.0,
    sample_offset_x_mm: float = 0.0,
    sample_offset_z_mm: float = 0.0,
)
```

A `PanelSpec` carries:
```python
@dataclass
class PanelSpec:
    image_cols: slice         # range of columns this panel occupies
    offset_deg: float         # nominal angle of the panel relative to central
    row_shift_px: float = 0.0 # per-panel row offset (calibration)
    col_shift_px: float = 0.0 # per-panel col offset (calibration)
    delta_deg: float = 0.0    # per-panel angle correction (calibration)
```

The fold geometry is reconstructed in `qmap()` (`integrator.py:381-485`):

1. Pick the **reference panel** as the one whose `offset_deg` is closest to 0 (typically the central panel).
2. Place its center using the beam-center offset: `bc_u = -(bc_col - panel_mid_ref) * pixel_mm`, then
   ```
   center_x[ref] = -bc_u * cos(α_ref)
   center_z[ref] =  R - bc_u * sin(α_ref)
   ```
   where `R = sample_distance_mm`, `α = panel_offset + delta_deg`.
3. Walk **outward** from the reference panel (both directions). Two adjacent panels share a fold line at the boundary column (`c1` of the inner or `c0` of the outer). The fold-line position must be the **same physical point** seen from both panels:
   ```
   u_fold_inner = -(fold_col - panel_mid_inner) * pixel_mm
   u_fold_outer = -(fold_col - panel_mid_outer) * pixel_mm
   center[outer] = center[inner] + u_inner * (cos α_inner, sin α_inner)
                                 - u_outer * (cos α_outer, sin α_outer)
   ```
   This propagates the geometry across panels while keeping the seams continuous.
4. For each panel, compute per-pixel `(x_det, y_det, z_det)` in detector frame:
   ```
   u    = -(col - panel_mid[idx]) * pixel_mm                # along the panel surface
   y_det = -(row - (bc_row + row_shift[idx])) * pixel_mm    # vertical (sign flip — see SAXS)
   x_det =  center_x[idx] + u * cos(α)
   z_det =  center_z[idx] + u * sin(α)
   ```
5. Apply the **arc rotation** (the WAXS detector physically rotates about y by `theta_deg`):
   ```
   x_lab = x_det * cos(θ) - z_det * sin(θ)
   z_lab = x_det * sin(θ) + z_det * cos(θ)
   y_lab = y_det
   ```
   where `θ = θ_arc + θ_zero`.
6. Apply sample offsets:
   ```
   p_x = x_lab - sample_offset_x_mm
   p_z = z_lab - sample_offset_z_mm
   p_y = y_det
   ```
7. Compute q from same formula as SAXS:
   ```
   r = sqrt(p_x² + p_y² + p_z²)
   k = 2π / wavelength_nm
   qx = k * p_x / r
   qy = k * p_y / r
   qz = k * (p_z / r - 1)
   qabs = sqrt(qx² + qy² + qz²)
   solid_angle = pixel_area * max(p_z, 0) / r³
   ```
8. Return as an `xr.Dataset` with vars `qx, qy, qz, qabs, solid_angle` and dims `(row, col)`.

## beam_center_at_angle (`WAXSCalibration`, integrator.py:330)

The WAXS detector **rotates about y** as the arc moves. The beam center on the detector face shifts with arc rotation in two ways:
1. **Sample offset projection**: if the sample isn't at the rotation axis, the beam center drifts as the detector rotates.
2. **Direct linear drift** (mechanical wobble), captured by `beam_col_per_arc_deg`.

The combined formula (`integrator.py:330-340`):

```python
def beam_center_at_angle(self, theta_deg: float) -> tuple[float, float]:
    row = float(self.beam_center_row)
    col = float(self.beam_center_col)
    if self.sample_offset_x_mm != 0 or self.sample_offset_z_mm != 0:
        th = np.deg2rad(theta_deg + self.theta_zero_deg)
        dx_mm = (self.sample_offset_x_mm * (np.cos(th) - 1.0)
                 - self.sample_offset_z_mm * np.sin(th))
        col += dx_mm / self.pixel_size_mm
    if self.beam_col_per_arc_deg != 0:
        col += self.beam_col_per_arc_deg * theta_deg
    return (row, col)
```

The `(cos θ - 1)` term is small for small `θ` and comes from rotating the sample-offset displacement vector by `-θ` (the detector frame's view of the sample motion).

> **NOTE**: beam_col_per_arc_deg is a per-frame correction with sign convention "+" = beam center moves +col with +arc. Calibration scripts in `scripts/` measure this from AGB scans.

## waxs_arc and waxs_bsx coupling

The arc rotation also drives a beamstop-y motor (`waxs_bsx`) so the WAXS beamstop tracks the direct beam:

```python
BSX_PER_ARC_DEG = -4.39       # (defaults.py:195)
_BSX_PER_ARC_DEG = -4.39      # (integrator.py:3925)
```

The two constants ought to be identical and currently are. **GOTCHA**: the integrator-internal `_BSX_PER_ARC_DEG` shadows the public `defaults.BSX_PER_ARC_DEG`. Only the integrator one is used by `_waxs_bsx_for_arc`, so editing `defaults.BSX_PER_ARC_DEG` will not affect the integrator. See `reference/gotchas.md`.

The relationship is:
```
waxs_bsx_expected_mm = waxs_bsx_at_arc_zero + BSX_PER_ARC_DEG * waxs_arc_deg
```

This is used inside `_make_waxs_shadow_mask` to compute where the dynamic shadow mask should sit.

## Per-panel calibration deltas

The 6 calibration deltas (`panel_row_shifts`, `panel_col_shifts`, `panel_delta_deg`, all length-3 tuples) capture the small departures of the real detector from its nominal layout. They are typically determined by:
1. AGB diffraction ring fitting at multiple arc angles.
2. The fold lines between adjacent panels are set so that AGB rings are **continuous** across the seams.

Currently all defaults are 0 (uncalibrated). When set, they should typically be in the range `±0.5 px` and `±0.05°`.

## How `integrate_waxs` builds the detector per frame

`integrator.py:3200-3500` (slice).

For each frame `f`:
1. Pull this frame's arc angle `θ_f` from the per-frame coord on the DataArray.
2. Pull the per-frame wavelength `λ_f` (from per-frame energy, `integrator.py:3287`).
3. Build a `WAXSCalibration` from the resolved geometry + θ + λ.
4. Look up the cached `MultiPanelArcDetector` keyed on `(θ_f, λ_f)`. Two frames with the same arc angle and energy share a `qmap`.
5. Call `qmap(θ_f)` to get `(qx, qy, qz, qabs, solid_angle)` per pixel.
6. Construct the static + dynamic mask (`_make_waxs_shadow_mask`) and run `_qchi_and_iq` on the masked image.

## Frame rotation: critical to remember

The image arrays come off tiled in **raw** orientation `(195, 619)`. The integrator immediately calls:

```python
rotate_image_and_mask(img, mask, k=3)
# = (np.fliplr(np.rot90(img, k=3)), np.fliplr(np.rot90(mask, k=3)))
```

(`integrator.py:488-495`). This is equivalent to a transpose for typical k=3 input, so the working shape becomes `(619, 195)`.

The rotation puts:
- the **arc axis** along the row direction (so panels are side-by-side in column index).
- the **panel-width axis** (radial / vertical-on-detector) along the column direction.

This is why `panel_col_ranges` are in column space `[0, 619)` even though the raw image only had 619 columns.

Masks in `data/masks/900KW_mask_polygons.json` are stored in raw coordinates and must be rotated before use. `make_saxs_mask_from_dict` and the WAXS mask helpers handle this internally.

## Cross-detector merging (SAXS + WAXS)

When `reduce_smi_combined` produces both SAXS and WAXS reductions:
- Both produce `(q, chi, intensity)` tables on independent `q` grids.
- `merge_q_chi_weighted` (`integrator.py:1742`) interleaves them into a single common `q` grid using count-weighted averaging.
- The overlap region is typically `q ≈ 0.3 — 1.0 Å^-1` (SAXS reaches to ~1 Å^-1, WAXS starts from ~0.3 Å^-1).
- `merge_iq_profiles` does the analogous merge on the chi-averaged 1D `I(q)`.

The full SAXS+WAXS merging logic uses linear interpolation onto a common log-spaced grid, then count-weighted averaging in the overlap region. See `reference/integration.md`.
