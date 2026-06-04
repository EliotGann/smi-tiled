"""
SMISWAXSIntegrator
===================
pyFAI-backed azimuthal integration for the SMI WAXS+SAXS instrument,
designed to work with :class:`SMISWAXSLoader.TiledSMISWAXSLoader`.

The module is self-contained: all geometry (multi-panel arc detector,
single-panel SAXS), masking, histogram binning, and merging utilities
are included so the only external runtime dependencies are numpy, xarray,
and (optionally) pyFAI.

Architecture
------------
1.  Raw images arrive as xr.DataArray from the loader with calibration in attrs.
2.  SAXS uses a flat-panel pixel-space q-map (matching SinglePanelSAXSDetector).
3.  WAXS uses per-frame q-maps from :class:`MultiPanelArcDetector` which
    computes exact 3D pixel positions for the folded 3-panel geometry.
4.  Both detectors bin into (q, chi) grids via histogram2d.
5.  Merged I(q, chi) and I(q) are produced via weighted overlap merge.

Key classes / functions
-----------------------
- ``MultiPanelArcDetector``  – WAXS 3-panel geometry model
- ``integrate_saxs``         – SAXS integration from raw DataArray
- ``integrate_waxs``         – WAXS integration from raw DataArray
- ``reduce_smi_combined``    – Full SAXS + WAXS pipeline returning merged result

Mask architecture
-----------------
SMI masking is built up in three layers, applied as a logical AND in this
order:

**Layer 1 — Fixed instrument geometry.**  Inter-module gaps and bad-pixel
regions on each detector.  These are physically wired into the hardware
and never move.  Shipped in the bundled JSON mask files
(``pil2M_mask_polygons.json``, ``900KW_mask_polygons.json``) under the
``static_regions`` block (SAXS) or as flat top-level polygons (WAXS).

**Layer 2 — Beamstop polygon, position-corrected from motor readings.**
SAXS has two beamstop variants (``rod`` and ``pin``) under
``beamstops``; the integrator reads ``pil2M_active_beamstop`` from
baseline to choose.  Polygon position uses one of:

* ``polygon_offsets_from_beam`` *(preferred)* — anchored to the
  per-frame beam center, which already tracks ``pil2M_motor_x/y`` and
  ``piezo_z``.
* ``polygon`` + ``reference_mm`` + ``pixels_per_mm`` *(legacy)* —
  shifted by ``(motor − reference) × px_per_mm``.

WAXS has a single beamstop polygon that shifts vertically by
``(waxs_bsx − waxs_bsx_ref) / pixel_size × 1.088``.  The reference is
derived from the SMI mechanical linkage
``waxs_bsx_ref = waxs_bsx − BSX_PER_ARC_DEG × waxs_arc``.  Auto-disabled
when ``|waxs_arc| > beamstop_max_abs_arc_deg`` (default 15°) — the
beamstop has cleared the active area.

**Layer 3 — Dynamic per-frame occlusion (SAXS only).**  Two extra masks
computed during integration:

* **WAXS-shadow mask** — the WAXS detector physically blocks part of the
  SAXS detector's view; the boundary column moves with ``waxs_arc``
  (``_make_waxs_shadow_mask``).
* **Aperture mask** — q-cutoff anchored to an AgBh ring order (default
  5) at the current SDD; auto-disables when the cutoff falls beyond the
  detector edge (long SDD).  (``_make_aperture_mask``)

Override knobs
~~~~~~~~~~~~~~
Layers 2 and 3 are independently togglable via :func:`reduce_smi_combined`::

    reduce_smi_combined(
        uid,
        saxs_mask=my_dict_or_path_or_None,    # Layer 1+2 source
        waxs_mask=my_dict_or_path_or_None,    # Layer 1+2 source (WAXS)
        saxs_kwargs={
            "dynamic_saxs_kwargs": {
                "aperture":    {"enabled": False},    # Layer 3a off
                "waxs_shadow": {"enabled": False},    # Layer 3b off
            }
        },
        saxs_q_cutoff=0.6,           # force aperture cutoff (nm⁻¹)
        saxs_agbh_ring_order=8,      # change anchor ring
    )

Mask inputs (Path *or* dict)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
All public mask entry points accept either a path to a JSON file or an
already-parsed dict with the same schema.  The dict form lets notebook
users compose and edit masks in memory:

.. code-block:: python

    import json
    spec = json.load(open(default_saxs_mask_path()))
    spec["static_regions"]["my_extra_blob"] = [[100, 200], [150, 200], ...]
    result = reduce_smi_combined(uid, saxs_mask=spec)

Functions taking either form:

* :func:`make_saxs_mask_from_spec` / :func:`make_saxs_mask_from_dict`
* :func:`make_waxs_mask_callable` / :func:`make_waxs_mask_callable_from_dict`
* :func:`mask_for_frame` (via ``mask_path=`` accepting dict)
* :func:`reduce_smi_combined` (``saxs_mask=`` / ``waxs_mask=`` kwargs)
* :func:`reduce_smi_gi` (``waxs_mask=`` kwarg)
"""
from __future__ import annotations

import json
import time as _time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence, Tuple

#: Callback signature for progress reporting from reduction pipelines.
#: ``(stage, current, total)`` — *stage* is a short string identifier
#: (e.g. ``"saxs_integrate"``), *current* is the 1-based step index,
#: and *total* is the number of steps in that stage.
ProgressCallback = Callable[[str, int, int], None]

import numpy as np
import xarray as xr

from smi_tiled.loader import (
    IMAGE_BLOCK_FRAMES,
    SAXSGeometry,
    WAXSGeometry,
    resolve_saxs_geometry,
)


# ===================================================================
# Geometry cache – persists across calls within the same Python process
# ===================================================================
#
# GUI integration
# ---------------
# The geometry cache stores precomputed per-pixel q-maps so that repeated
# reductions with the same detector geometry (energy, distance, beam center,
# panel offsets, masks, etc.) skip the expensive trigonometry.
#
# Usage from a GUI or batch script:
#
#     from smi_tiled.integrator import (
#         reduce_smi_combined,
#         clear_geometry_cache,
#         geometry_cache_info,
#     )
#
#     # Process many scans — geometry is computed once, then reused:
#     for uid in uid_list:
#         result = reduce_smi_combined(uid, cache_geometry=True, ...)
#
#     # When the user changes calibration (beam center, distance, energy,
#     # panel offsets, etc.), clear the stale cache:
#     clear_geometry_cache()
#
#     # Or inspect current cache size:
#     info = geometry_cache_info()
#     print(f"Cache holds {info['waxs_entries']} WAXS geometries, "
#           f"~{info['estimated_mb']:.1f} MB")
#
# The cache is keyed on the full set of geometry parameters, so if you
# change *any* calibration value the old entry simply won't match and a
# new one will be computed (no stale-data risk).  Call
# ``clear_geometry_cache()`` only to free memory.

_WAXS_GEOMETRY_CACHE: dict[tuple, dict[float, tuple]] = {}
_SAXS_GEOMETRY_CACHE: dict[tuple, tuple] = {}


def clear_geometry_cache() -> None:
    """Clear the module-level geometry cache to free memory.

    Call this when you want to reclaim memory in a long-running process
    (e.g. a GUI).  It is *not* necessary to call this when calibration
    parameters change — the cache is keyed on the full parameter set, so
    a changed parameter simply produces a new cache entry.

    Example (in a GUI callback when user clicks "Reset Calibration")::

        from smi_tiled.integrator import clear_geometry_cache
        clear_geometry_cache()
    """
    _WAXS_GEOMETRY_CACHE.clear()
    _SAXS_GEOMETRY_CACHE.clear()


def geometry_cache_info() -> dict[str, Any]:
    """Return a summary of the current geometry cache state.

    Returns
    -------
    dict with keys:
        waxs_entries : int – number of distinct WAXS calibration keys cached
        waxs_angles_total : int – total number of cached arc-angle q-maps
        saxs_entries : int – number of distinct SAXS geometry keys cached
        estimated_mb : float – rough memory estimate in megabytes
    """
    waxs_angles = sum(len(v) for v in _WAXS_GEOMETRY_CACHE.values())
    # Each cached angle holds ~5 arrays of shape (ny, nx); estimate 619×487
    # ≈ 300k pixels × 8 bytes × 5 arrays ≈ 12 MB per angle
    est_per_angle_mb = 12.0
    # Each SAXS entry holds ~5 arrays of shape (ny, nx); estimate 1475×1679
    # ≈ 2.5M pixels × 8 bytes × 5 ≈ 100 MB per entry
    est_per_saxs_mb = 100.0
    estimated_mb = (
        waxs_angles * est_per_angle_mb
        + len(_SAXS_GEOMETRY_CACHE) * est_per_saxs_mb
    )
    return {
        "waxs_entries": len(_WAXS_GEOMETRY_CACHE),
        "waxs_angles_total": waxs_angles,
        "saxs_entries": len(_SAXS_GEOMETRY_CACHE),
        "estimated_mb": estimated_mb,
    }


def _waxs_cache_key(
    cal: "WAXSCalibration",
    image_shape: tuple[int, int],
    flip_horizontal: bool,
    qx_shift_nm: float,
    qy_shift_nm: float,
) -> tuple:
    """Build a hashable key from all parameters that affect WAXS q-maps."""
    return (
        cal.energy_kev,
        cal.sample_distance_mm,
        cal.pixel_size_mm,
        cal.beam_center_row,
        cal.beam_center_col,
        tuple(cal.panel_col_ranges),
        tuple(cal.panel_offsets_deg),
        tuple(cal.panel_row_shifts),
        tuple(cal.panel_col_shifts),
        tuple(cal.panel_delta_deg),
        cal.theta_zero_deg,
        cal.sample_offset_x_mm,
        cal.sample_offset_z_mm,
        cal.beam_col_per_arc_deg,
        cal.q_horizontal_sign,
        cal.q_vertical_sign,
        cal.rotation_k,
        image_shape,
        flip_horizontal,
        round(qx_shift_nm, 10),
        round(qy_shift_nm, 10),
    )


def _saxs_cache_key(
    dist_m: float,
    poni1_m: float,
    poni2_m: float,
    pixel1_m: float,
    pixel2_m: float,
    wavelength_m: float,
    image_shape: tuple[int, int],
) -> tuple:
    """Build a hashable key from all parameters that affect SAXS q-maps."""
    return (
        round(dist_m, 12),
        round(poni1_m, 12),
        round(poni2_m, 12),
        round(pixel1_m, 12),
        round(pixel2_m, 12),
        round(wavelength_m, 15),
        image_shape,
    )


# ===================================================================
# WAXS detector geometry – ported from waxs_reduce.py
# ===================================================================

def wavelength_nm_from_energy_kev(energy_kev: float) -> float:
    return 1.23984198 / float(energy_kev)


@dataclass(frozen=True)
class PanelSpec:
    image_cols: slice
    offset_deg: float
    row_shift_px: float = 0.0
    col_shift_px: float = 0.0
    delta_deg: float = 0.0


@dataclass
class WAXSCalibration:
    energy_kev: float = 16.1
    sample_distance_mm: float = 273.0
    pixel_size_mm: float = 0.172
    beam_center_row: float = 217.0
    beam_center_col: float = 319.0
    panel_col_ranges: Tuple = ((0, 206), (206, 413), (413, 619))
    panel_offsets_deg: Tuple = (-7.0, 0.0, 7.0)
    panel_row_shifts: Tuple = (0.0, 0.0, 0.0)
    panel_col_shifts: Tuple = (0.0, 0.0, 0.0)
    panel_delta_deg: Tuple = (0.0, 0.0, 0.0)
    theta_zero_deg: float = 0.0
    sample_offset_x_mm: float = 0.0
    sample_offset_z_mm: float = 0.0
    beam_col_per_arc_deg: float = 0.0
    q_horizontal_sign: float = -1.0
    q_vertical_sign: float = -1.0
    rotation_k: int = 3

    @property
    def wavelength_nm(self) -> float:
        return wavelength_nm_from_energy_kev(self.energy_kev)

    def make_panel_specs(self) -> list[PanelSpec]:
        specs = []
        for i, ((c0, c1), off) in enumerate(
            zip(self.panel_col_ranges, self.panel_offsets_deg)
        ):
            specs.append(
                PanelSpec(
                    image_cols=slice(int(c0), int(c1)),
                    offset_deg=float(off),
                    row_shift_px=float(self.panel_row_shifts[i]),
                    col_shift_px=float(self.panel_col_shifts[i]),
                    delta_deg=float(self.panel_delta_deg[i]),
                )
            )
        return specs

    def beam_center_at_angle(self, theta_deg: float) -> Tuple[float, float]:
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


# Default calibration matching legacy waxs_reduce._DEFAULT_CAL
_DEFAULT_CAL = dict(
    energy_kev=16.1,
    sample_distance_mm=273,
    beam_center_row=217.0,
    beam_center_col=319.0,
    panel_offsets_deg=(-7.0, 0.0, 7.0),
    theta_zero_deg=0,
    sample_offset_z_mm=0.0,
)


class MultiPanelArcDetector:
    """3-panel folded arc WAXS detector geometry model."""

    def __init__(
        self,
        image_shape: Tuple[int, int],
        panel_specs: Sequence[PanelSpec],
        wavelength_nm: float,
        pixel_size_mm: float = 0.172,
        sample_distance_mm: float = 300.0,
        beam_center_px: Tuple[float, float] = (0.0, 0.0),
        theta_zero_deg: float = 0.0,
        sample_offset_x_mm: float = 0.0,
        sample_offset_z_mm: float = 0.0,
    ) -> None:
        self.ny, self.nx = image_shape
        self.panel_specs = list(panel_specs)
        self.wavelength_nm = float(wavelength_nm)
        self.pixel_size_mm = float(pixel_size_mm)
        self.sample_distance_mm = float(sample_distance_mm)
        self.beam_center_row_px = float(beam_center_px[0])
        self.beam_center_col_px = float(beam_center_px[1])
        self.theta_zero_deg = float(theta_zero_deg)
        self.sample_offset_x_mm = float(sample_offset_x_mm)
        self.sample_offset_z_mm = float(sample_offset_z_mm)

    def qmap(self, theta_deg: float) -> xr.Dataset:
        """Compute per-pixel q-vectors and solid angle for a given arc angle."""
        p_mm = self.pixel_size_mm
        R    = self.sample_distance_mm
        n_panels = len(self.panel_specs)
        rows = np.arange(self.ny, dtype=float)

        panel_c0s, panel_c1s, panel_mids, alphas_det = [], [], [], []
        for ps in self.panel_specs:
            c0 = 0 if ps.image_cols.start is None else ps.image_cols.start
            c1 = self.nx if ps.image_cols.stop is None else ps.image_cols.stop
            panel_c0s.append(c0)
            panel_c1s.append(c1)
            panel_mids.append(0.5 * (c0 + c1 - 1) + ps.col_shift_px)
            alphas_det.append(np.deg2rad(ps.offset_deg + ps.delta_deg))

        ref_idx = min(
            range(n_panels),
            key=lambda i: abs(self.panel_specs[i].offset_deg),
        )
        bc_u = -(self.beam_center_col_px - panel_mids[ref_idx]) * p_mm
        alpha_r = alphas_det[ref_idx]
        centers_x = [None] * n_panels
        centers_z = [None] * n_panels
        centers_x[ref_idx] = -bc_u * np.cos(alpha_r)
        centers_z[ref_idx] = R - bc_u * np.sin(alpha_r)

        for i in range(ref_idx - 1, -1, -1):
            fold = panel_c1s[i] - 0.5
            u1 = -(fold - panel_mids[i + 1]) * p_mm
            u0 = -(fold - panel_mids[i]) * p_mm
            centers_x[i] = (
                centers_x[i + 1]
                + u1 * np.cos(alphas_det[i + 1])
                - u0 * np.cos(alphas_det[i])
            )
            centers_z[i] = (
                centers_z[i + 1]
                + u1 * np.sin(alphas_det[i + 1])
                - u0 * np.sin(alphas_det[i])
            )
        for i in range(ref_idx + 1, n_panels):
            fold = panel_c0s[i] - 0.5
            u1 = -(fold - panel_mids[i - 1]) * p_mm
            u0 = -(fold - panel_mids[i]) * p_mm
            centers_x[i] = (
                centers_x[i - 1]
                + u1 * np.cos(alphas_det[i - 1])
                - u0 * np.cos(alphas_det[i])
            )
            centers_z[i] = (
                centers_z[i - 1]
                + u1 * np.sin(alphas_det[i - 1])
                - u0 * np.sin(alphas_det[i])
            )

        theta_rad = np.deg2rad(float(theta_deg) + self.theta_zero_deg)
        cos_th, sin_th = np.cos(theta_rad), np.sin(theta_rad)

        px_mm = np.full((self.ny, self.nx), np.nan)
        py_mm = np.full((self.ny, self.nx), np.nan)
        pz_mm = np.full((self.ny, self.nx), np.nan)

        for idx, ps in enumerate(self.panel_specs):
            c0, c1 = panel_c0s[idx], panel_c1s[idx]
            cols_p = np.arange(c0, c1, dtype=float)
            _rr, _cc = np.meshgrid(rows, cols_p, indexing="ij")
            u = -(_cc - panel_mids[idx]) * p_mm
            y_det = -(
                _rr - (self.beam_center_row_px + ps.row_shift_px)
            ) * p_mm
            alpha = alphas_det[idx]
            x_det = centers_x[idx] + u * np.cos(alpha)
            z_det = centers_z[idx] + u * np.sin(alpha)
            x_lab = x_det * cos_th - z_det * sin_th
            z_lab = x_det * sin_th + z_det * cos_th
            px_mm[:, c0:c1] = x_lab - self.sample_offset_x_mm
            py_mm[:, c0:c1] = y_det
            pz_mm[:, c0:c1] = z_lab - self.sample_offset_z_mm

        r = np.sqrt(px_mm**2 + py_mm**2 + pz_mm**2)
        k = 2.0 * np.pi / self.wavelength_nm
        with np.errstate(invalid="ignore", divide="ignore"):
            qx = k * px_mm / r
            qy = k * py_mm / r
            qz = k * (pz_mm / r - 1.0)
        qabs = np.sqrt(qx**2 + qy**2 + qz**2)

        pixel_area_mm2 = p_mm * p_mm
        with np.errstate(invalid="ignore", divide="ignore"):
            solid_angle = pixel_area_mm2 * np.maximum(pz_mm, 0.0) / (r**3)

        return xr.Dataset(
            {
                "qx": (("row", "col"), qx),
                "qy": (("row", "col"), qy),
                "qz": (("row", "col"), qz),
                "qabs": (("row", "col"), qabs),
                "solid_angle": (("row", "col"), solid_angle),
            },
            coords={
                "row": np.arange(self.ny),
                "col": np.arange(self.nx),
            },
        )


def rotate_image_and_mask(
    image: np.ndarray,
    mask: np.ndarray | None = None,
    k: int = 3,
) -> tuple[np.ndarray, np.ndarray | None]:
    img_rot = np.fliplr(np.rot90(np.asarray(image), k=k))
    mask_rot = np.fliplr(np.rot90(np.asarray(mask), k=k)) if mask is not None else None
    return img_rot, mask_rot


def dezinger(
    image: np.ndarray,
    kernel_size: int = 5,
    threshold: float = 5.0,
) -> np.ndarray:
    """Detect hot/dead pixels via local-mean outlier rejection.

    Uses a fast separable uniform (mean) filter instead of median_filter.
    For photon-counting detectors, zingers are extreme outliers (10-100×
    above background) so mean-based detection is effective and ~10× faster
    than median-based for large images.

    Pixels deviating by more than ``threshold`` × σ from the local mean
    are flagged.

    Parameters
    ----------
    image : ndarray
        2-D detector image.
    kernel_size : int
        Side length of the square filter kernel (default 5).
    threshold : float
        Number of σ above the local mean to flag as bad (default 5).

    Returns
    -------
    ndarray[bool]
        True for *valid* pixels, False for outliers (matches mask convention).
    """
    from scipy.ndimage import uniform_filter

    img = np.asarray(image, dtype=float)
    local_mean = uniform_filter(img, size=kernel_size, mode="reflect")
    # Local variance via E[X²] - E[X]²
    local_sq = uniform_filter(img * img, size=kernel_size, mode="reflect")
    local_var = np.maximum(local_sq - local_mean**2, 0.0)
    local_std = np.sqrt(local_var)
    # Floor the std to avoid flagging pixels in truly uniform regions
    global_std = np.nanstd(img[np.isfinite(img)]) if np.any(np.isfinite(img)) else 1.0
    local_std = np.maximum(local_std, 0.01 * global_std)
    diff = img - local_mean
    return ~(np.abs(diff) > threshold * local_std)


# ===================================================================
# Mask utilities – ported from waxs_reduce.py / saxs_reduce.py
# ===================================================================

def polygons_to_mask(
    shape: tuple[int, int],
    polygons: list,
) -> np.ndarray:
    """Build a boolean mask (True = valid) from a list of polygon regions."""
    from skimage.draw import polygon as skpoly

    mask = np.ones(shape, dtype=bool)
    ny, nx = shape
    for poly in polygons:
        if not poly:
            continue
        cols = np.array([p[0] for p in poly], dtype=float)
        rows = np.array([p[1] for p in poly], dtype=float)
        rr, cc = skpoly(rows, cols, shape=shape)
        mask[rr, cc] = False
    return mask


def shift_polygon(
    polygon: list[list[float]],
    dx_px: float = 0.0,
    dy_px: float = 0.0,
) -> list[list[float]]:
    return [[col + dx_px, row + dy_px] for col, row in polygon]


def make_mask_for_angle(
    image_shape_raw: tuple[int, int],
    static_regions: dict,
    beamstop_region: list,
    waxs_bsx: float,
    waxs_bsx_ref: float,
    pixel_size_mm: float = 0.172,
    rotation_k: int = 3,
    include_beamstop: bool = True,
) -> np.ndarray:
    """Build a per-angle WAXS mask accounting for beamstop motor position."""
    polys = list(static_regions.values())
    if include_beamstop and beamstop_region:
        bs_shift_mm = waxs_bsx - waxs_bsx_ref
        bs_shift_px = (bs_shift_mm / pixel_size_mm) * 1.088  # empirical fudge factor
        polys.append(shift_polygon(beamstop_region, dx_px=0.0, dy_px=bs_shift_px))
    raw_mask = polygons_to_mask(image_shape_raw, polys)
    mask_rot, _ = rotate_image_and_mask(raw_mask, k=rotation_k)
    return mask_rot


# ===================================================================
# SAXS mask builders
# ===================================================================

def _resolve_mask_spec(spec: "str | Path | dict") -> dict:
    """Resolve a mask spec input to a parsed dict.

    Accepts either a file path (str or pathlib.Path) pointing to a JSON
    mask file, or an already-parsed dict in the same schema.  The dict
    form lets notebook users compose and edit masks in memory without
    writing temp files.
    """
    if isinstance(spec, dict):
        return spec
    with open(spec) as f:
        return json.load(f)


def make_saxs_mask_from_dict(
    image_shape: tuple[int, int],
    mask_spec: dict,
    active_beamstop: str = "rod",
    beamstop_pos_mm: dict | None = None,
    beam_center_px: tuple[float, float] | None = None,
) -> np.ndarray:
    """Build a static SAXS mask (True = valid) from a parsed mask dict.

    This is the in-memory variant of :func:`make_saxs_mask_from_spec` —
    use it when you have the mask polygons in a Python dict (perhaps
    built up from a notebook UI or assembled programmatically) and don't
    want to serialize to a temp file first.

    The expected schema matches the bundled JSON mask files:

    .. code-block:: python

        {
            "image_shape": [rows, cols],          # optional, informational
            "static_regions": {
                "gap_1": [[col, row], [col, row], ...],   # one polygon per key
                "gap_2": [...],
            },
            "beamstops": {
                "rod": {                          # one entry per beamstop variant
                    "polygon_offsets_from_beam": [[d_col, d_row], ...],
                    # OR (legacy form):
                    "polygon": [[col, row], ...],
                    "reference_mm": {"x": 0.0, "y": 0.0},
                    "pixels_per_mm": {"x": 5.81, "y": 5.81},
                },
                "pin": { ... },
            },
        }

    A beamstop entry may be a ``dict`` (specified in one of two ways) or a
    bare polygon:

    * ``polygon_offsets_from_beam`` — list of ``[d_col, d_row]`` offsets
      relative to the resolved beam center.  This is the preferred form
      because the beam center already accounts for detector motor
      positions (pil2M_motor_x/y) and sample-z offsets (piezo_z).  Requires
      *beam_center_px* to be supplied.  Also accepts a list of polygons.

    * ``polygon`` (legacy) — absolute pixel coordinates.  The polygon is
      then shifted by ``(cur_motor - reference_mm) * pixels_per_mm`` if
      *beamstop_pos_mm* is provided.  Retained for backwards compatibility
      with mask files written before beam-center anchoring was supported.

    * **bare polygon** — the beamstop maps directly to a ``[[col, row], ...]``
      polygon (or a list of such polygons) in absolute detector coordinates,
      with no enclosing dict.  This is the normalized form emitted by
      :func:`smi_tiled.defaults.load_mask_polygons` and produced by hand-edited
      / exported mask files.  Used as-is (no beam-center anchoring or motor
      shift).

    Parameters
    ----------
    image_shape : (rows, cols)
        Raw detector image shape.
    mask_spec : dict
        Parsed polygon spec (see schema above).
    active_beamstop : {'rod', 'pin'} or other key present in mask spec
        Which beamstop is currently in the beam.  Selects the matching
        polygon (or polygons) from the ``beamstops`` block.
    beamstop_pos_mm : dict | None
        Per-beamstop motor positions (legacy ``polygon`` only).
    beam_center_px : (row, col) | None
        Resolved beam center in raw-detector pixel coordinates.  Required
        when a beamstop entry uses ``polygon_offsets_from_beam``.

    Returns
    -------
    np.ndarray[bool]
        Boolean mask; ``True`` marks a valid pixel.
    """
    static_regions = mask_spec.get("static_regions", {})
    beamstops_spec = mask_spec.get("beamstops", {})
    polys = list(static_regions.values())

    bs = beamstops_spec.get(active_beamstop, {})

    if isinstance(bs, (list, tuple)):
        # Bare polygon(s) in absolute detector pixel coords — the normalized
        # form emitted by :func:`smi_tiled.defaults.load_mask_polygons` and used
        # by hand-edited / exported mask files, where a beamstop maps directly
        # to a ``[[col, row], ...]`` polygon (or a list of such polygons) rather
        # than to a ``{"polygon_offsets_from_beam": ...}`` dict.  Appended as-is:
        # the coordinates are already absolute, so no beam-center anchoring or
        # motor shift is applied.
        if bs:
            first = bs[0]
            is_list_of_polys = (
                isinstance(first, (list, tuple))
                and len(first) > 0
                and isinstance(first[0], (list, tuple))
            )
            bare_polys = bs if is_list_of_polys else [bs]
            for poly in bare_polys:
                if poly:
                    polys.append([[float(c), float(r)] for c, r in poly])
        return polygons_to_mask(image_shape, polys)

    bs_offsets = bs.get("polygon_offsets_from_beam")
    if bs_offsets is not None:
        if beam_center_px is None:
            warnings.warn(
                f"Beamstop {active_beamstop!r} uses polygon_offsets_from_beam "
                "but no beam_center_px was provided — skipping beamstop mask.",
                stacklevel=2,
            )
        else:
            bc_row, bc_col = float(beam_center_px[0]), float(beam_center_px[1])
            # Normalize to a list of polygons (single polygon → list-of-one).
            if bs_offsets and isinstance(bs_offsets[0][0], (list, tuple)):
                offset_polys = bs_offsets
            else:
                offset_polys = [bs_offsets]
            for offs in offset_polys:
                if not offs:
                    continue
                shifted = [[bc_col + float(dc), bc_row + float(dr)]
                           for dc, dr in offs]
                polys.append(shifted)
    else:
        bs_poly = bs.get("polygon")
        if bs_poly is not None:
            bs_ref = bs.get("reference_mm", {})
            bs_ref_x = bs_ref.get("x", 0.0)
            bs_ref_y = bs_ref.get("y", 0.0)
            px_per_mm_map = bs.get("pixels_per_mm", {})
            px_per_mm_x = (
                px_per_mm_map.get("x", 1.0 / 0.172)
                if isinstance(px_per_mm_map, dict)
                else float(px_per_mm_map)
            )
            px_per_mm_y = (
                px_per_mm_map.get("y", 1.0 / 0.172)
                if isinstance(px_per_mm_map, dict)
                else float(px_per_mm_map)
            )
            if beamstop_pos_mm:
                cur = beamstop_pos_mm.get(active_beamstop) or {}
                cur_x = cur.get("x") or bs_ref_x
                cur_y = cur.get("y") or bs_ref_y
                dx = (cur_x - bs_ref_x) * px_per_mm_x
                dy = (cur_y - bs_ref_y) * px_per_mm_y
                polys.append(shift_polygon(bs_poly, dx_px=dx, dy_px=dy))
            else:
                polys.append(bs_poly)

    return polygons_to_mask(image_shape, polys)


def make_saxs_mask_from_spec(
    image_shape: tuple[int, int],
    mask_path: "str | Path | dict",
    active_beamstop: str = "rod",
    beamstop_pos_mm: dict | None = None,
    beam_center_px: tuple[float, float] | None = None,
) -> np.ndarray:
    """Build a static SAXS mask from a JSON mask file *or* parsed dict.

    Thin wrapper around :func:`make_saxs_mask_from_dict`: when
    *mask_path* is a string or :class:`~pathlib.Path`, the JSON is loaded
    from disk; when it's a ``dict`` the parsed contents are used
    directly (no file I/O).

    All other arguments and the return value are documented in
    :func:`make_saxs_mask_from_dict`.
    """
    mask_spec = _resolve_mask_spec(mask_path)
    return make_saxs_mask_from_dict(
        image_shape=image_shape,
        mask_spec=mask_spec,
        active_beamstop=active_beamstop,
        beamstop_pos_mm=beamstop_pos_mm,
        beam_center_px=beam_center_px,
    )


def make_waxs_mask_callable_from_dict(
    mask_data: dict,
    waxs_bsx_ref: float = 0.0,
    beamstop_max_abs_arc_deg: float | None = 15.0,
):
    """Return a per-frame WAXS mask function built from a parsed dict.

    Two schemas are accepted:

    1. **Nested** (matches the bundled SAXS schema)::

           {
               "static_regions": {"gap_upper": [...], "bad_module": [...]},
               "beamstops": {"beamstop": <polygon-or-wrapper>},
           }

       where ``<polygon-or-wrapper>`` is either a bare polygon
       (``[[col, row], ...]``) or a wrapper ``{"polygons": [[...]]}``.

    2. **Flat** (legacy, matches the bundled WAXS file)::

           {
               "gap_upper":  [[col, row], ...],
               "bad_module": [[col, row], ...],
               "beamstop":   [[col, row], ...],
           }

       The ``beamstop`` key is the moving region; everything else is
       treated as static.

    Returns
    -------
    callable
        ``mask_fn(image_shape_raw, theta_deg, waxs_bsx) -> ndarray[bool]``,
        returning a mask already rotated to display orientation
        (``np.fliplr(np.rot90(..., k=3))`` applied internally).

    Parameters
    ----------
    mask_data : dict
        Parsed polygon spec.
    waxs_bsx_ref : float
        Reference ``waxs_bsx`` value used to compute beamstop polygon
        shift.  Per-frame ``(waxs_bsx − waxs_bsx_ref)`` drives a vertical
        polygon shift via the SMI mechanical linkage.
    beamstop_max_abs_arc_deg : float or None
        Skip the beamstop polygon when ``|theta_deg|`` exceeds this
        value (the beamstop has cleared the active area).  ``None``
        keeps the beamstop active at every angle.
    """
    if "static_regions" in mask_data or "beamstops" in mask_data:
        static_regions = mask_data.get("static_regions", {})
        bs_entry = mask_data.get("beamstops", {}).get("beamstop", {})
        if isinstance(bs_entry, list):
            # beamstop is stored directly as a polygon (list of [x, y] pairs)
            beamstop_region = bs_entry
        else:
            beamstop_region = (
                bs_entry.get("polygons", [[]])[0]
                if bs_entry.get("polygons")
                else []
            )
    else:
        beamstop_region = mask_data.get("beamstop", [])
        static_regions = {
            key: value
            for key, value in mask_data.items()
            if key != "beamstop"
        }

    # Pre-compute the static mask once (static regions never change per frame).
    # Only the beamstop shifts per-frame via waxs_bsx.
    _static_mask_cache: dict[tuple[int, int], np.ndarray] = {}
    _rotated_cache: dict[tuple[tuple[int, int], bool, float], np.ndarray] = {}

    def _get_static_mask(image_shape_raw):
        if image_shape_raw not in _static_mask_cache:
            _static_mask_cache[image_shape_raw] = polygons_to_mask(
                image_shape_raw, list(static_regions.values())
            )
        return _static_mask_cache[image_shape_raw]

    def mask_fn(image_shape_raw, theta_deg, waxs_bsx):
        include_beamstop = (
            beamstop_max_abs_arc_deg is None
            or abs(float(theta_deg)) <= float(beamstop_max_abs_arc_deg)
        )
        # Cache key: (shape, beamstop_included, rounded bsx)
        cache_key = (image_shape_raw, include_beamstop,
                     round(float(waxs_bsx), 4) if include_beamstop else 0.0)
        cached = _rotated_cache.get(cache_key)
        if cached is not None:
            return cached

        static_mask = _get_static_mask(image_shape_raw)
        if include_beamstop and beamstop_region:
            bs_shift_mm = float(waxs_bsx) - float(waxs_bsx_ref)
            bs_shift_px = (bs_shift_mm / 0.172) * 1.088
            shifted_bs = shift_polygon(beamstop_region, dx_px=0.0, dy_px=bs_shift_px)
            bs_mask = polygons_to_mask(image_shape_raw, [shifted_bs])
            raw_mask = static_mask & bs_mask
        else:
            raw_mask = static_mask
        mask_rot, _ = rotate_image_and_mask(raw_mask, k=3)
        _rotated_cache[cache_key] = mask_rot
        return mask_rot

    return mask_fn


def make_waxs_mask_callable(
    mask_path: "str | Path | dict",
    waxs_bsx_ref: float = 0.0,
    beamstop_max_abs_arc_deg: float | None = 15.0,
):
    """Return a per-frame WAXS mask function built from a path or dict.

    Thin wrapper around :func:`make_waxs_mask_callable_from_dict`: loads
    the JSON when *mask_path* is a string or :class:`~pathlib.Path`,
    passes through when it's a ``dict``.  All other parameters are
    documented in :func:`make_waxs_mask_callable_from_dict`.

    Returns
    -------
    callable
        ``mask_fn(image_shape_raw, theta_deg, waxs_bsx) -> bool mask``
    """
    mask_data = _resolve_mask_spec(mask_path)
    return make_waxs_mask_callable_from_dict(
        mask_data=mask_data,
        waxs_bsx_ref=waxs_bsx_ref,
        beamstop_max_abs_arc_deg=beamstop_max_abs_arc_deg,
    )


# ===================================================================
# Single-frame mask convenience (browser / notebook helper)
# ===================================================================

def _smi_run_field_at(run: Any, field: str, frame_idx: int) -> float:
    """Return ``run.primary[field][frame_idx]`` as a Python float.

    Tolerates several access patterns so it works against bluesky/tiled
    runs *and* dict-like fakes in tests::

        run["primary"]["data"][field]              # tiled, bluesky-tiled
        run["primary"][field]                      # legacy tiled
        run.primary[field]                         # attribute-style
        run["primary"].read()[field]               # full read fallback
        run[field]                                 # bare dict fallback (tests)
    """
    primary = None
    try:
        primary = run["primary"]
    except Exception:
        primary = getattr(run, "primary", None)

    candidates = []
    if primary is not None:
        # bluesky-tiled: primary/data/<field>
        try:
            candidates.append(primary["data"][field])
        except Exception:
            pass
        try:
            candidates.append(primary[field])
        except Exception:
            pass
        try:
            candidates.append(getattr(primary, field))
        except Exception:
            pass
    # Bare dict-like at top level (used by simple test fakes)
    try:
        candidates.append(run[field])
    except Exception:
        pass

    for node in candidates:
        if node is None:
            continue
        try:
            values = node.read() if hasattr(node, "read") else node
            arr = np.asarray(values).reshape(-1)
            if arr.size == 0:
                continue
            idx = int(frame_idx) if arr.size > 1 else 0
            return float(arr[idx])
        except Exception:
            continue

    # Fall back to a full primary.read() — slower, but always correct.
    if primary is not None:
        try:
            ds = primary.read()
            arr = np.asarray(ds[field].values).reshape(-1)
            idx = int(frame_idx) if arr.size > 1 else 0
            return float(arr[idx])
        except Exception:
            pass

    raise KeyError(field)


def _smi_run_raw_shape(run: Any, image_field: str) -> tuple[int, int]:
    """Return ``(rows, cols)`` of one frame without downloading the data."""
    primary = None
    try:
        primary = run["primary"]
    except Exception:
        primary = getattr(run, "primary", None)

    nodes = []
    if primary is not None:
        try:
            nodes.append(primary["data"][image_field])
        except Exception:
            pass
        try:
            nodes.append(primary[image_field])
        except Exception:
            pass
    try:
        nodes.append(run[image_field])
    except Exception:
        pass

    for node in nodes:
        try:
            shp = tuple(getattr(node, "shape", ()))
            if len(shp) >= 2:
                return (int(shp[-2]), int(shp[-1]))
        except Exception:
            continue
    raise KeyError(f"could not determine shape of {image_field!r} on run")


def mask_for_frame(
    run_or_uid: Any,
    frame_idx: int,
    detector: str,
    *,
    mask_path: "str | Path | dict | None" = None,
    orient_for_display: bool = False,
    tiled_uri: str | None = None,
    catalog: str | None = None,
    raw_shape: tuple[int, int] | None = None,
    beamstop_max_abs_arc_deg: float | None = 15.0,
) -> np.ndarray:
    """Return the boolean validity mask (True = valid) for one frame.

    Thin wrapper over :func:`make_saxs_mask_from_spec` /
    :func:`make_waxs_mask_callable` that pulls the per-frame motor
    positions a browser/notebook would otherwise have to fetch by hand.

    Parameters
    ----------
    run_or_uid
        A bluesky/tiled run object, **or** a uid string.  If a string is
        given, ``tiled_uri`` and ``catalog`` are used to resolve it
        (defaults from :mod:`smi_tiled.defaults`).
    frame_idx : int
        Frame index along the scan axis (e.g. ``waxs_arc``).  Ignored for
        SAXS in current SMI configuration but accepted for symmetry.
    detector : {'saxs', 'waxs'}
        Which detector's mask to build.
    mask_path : str, Path, dict, or None
        Polygon-mask JSON path *or* already-parsed polygon dict.
        ``None`` selects the bundled default from
        :func:`smi_defaults.resolve_mask_path`.
    orient_for_display : bool, optional
        If True, return the mask already aligned with
        :func:`smi_defaults.orient_frame_for_display`.  For WAXS the
        underlying mask builder *already* returns a display-oriented
        array (it applies ``np.fliplr(np.rot90(..., k=3))`` internally),
        so no extra orientation pass is performed in that case.
    tiled_uri, catalog : str | None
        Used only when ``run_or_uid`` is a uid string.
    raw_shape : tuple[int, int] | None
        Override for the raw detector shape ``(rows, cols)``.  Useful for
        testing; normally read from the run's primary stream metadata.
    beamstop_max_abs_arc_deg : float | None
        Forwarded to :func:`make_waxs_mask_callable`.

    Returns
    -------
    np.ndarray[bool]
        Boolean mask, ``True`` where the pixel is valid for integration.

    Raises
    ------
    KeyError
        If WAXS is requested and ``waxs_arc`` / ``waxs_bsx`` cannot be
        located on the run's primary stream.
    ValueError
        If ``detector`` is not ``'saxs'`` or ``'waxs'``.
    """
    from smi_tiled.defaults import (
        DEFAULT_TILED_URI as _DEFAULT_TILED_URI,
        DEFAULT_CATALOG as _DEFAULT_CATALOG,
        SAXS_IMAGE_FIELD as _SAXS_IMAGE_FIELD,
        WAXS_IMAGE_FIELD as _WAXS_IMAGE_FIELD,
        WAXS_ARC_FIELD as _WAXS_ARC_FIELD,
        WAXS_BSX_FIELD as _WAXS_BSX_FIELD,
        BSX_PER_ARC_DEG as _BSX_PER_ARC_DEG_PUBLIC,
        orient_frame_for_display as _orient_frame_for_display,
        resolve_mask_path as _resolve_mask_path,
    )

    det = str(detector).lower()
    if det not in ("saxs", "waxs"):
        raise ValueError(f"detector must be 'saxs' or 'waxs', got {detector!r}")

    # Resolve uid → run if necessary.
    if isinstance(run_or_uid, str):
        from tiled.client import from_uri  # heavy, lazy
        cat_path = catalog or _DEFAULT_CATALOG
        uri = tiled_uri or _DEFAULT_TILED_URI
        run = from_uri(uri)[cat_path][run_or_uid]
    else:
        run = run_or_uid

    image_field = _SAXS_IMAGE_FIELD if det == "saxs" else _WAXS_IMAGE_FIELD
    if raw_shape is None:
        raw_shape = _smi_run_raw_shape(run, image_field)
    raw_shape = (int(raw_shape[0]), int(raw_shape[1]))

    # Dict inputs bypass the path-based resolver; otherwise let the
    # resolver apply bundled-default fallback.
    if isinstance(mask_path, dict):
        resolved_mask_path = mask_path
    else:
        resolved_mask_path = _resolve_mask_path(mask_path, detector=det)

    if det == "saxs":
        # Pull the actual active beamstop + per-run motor positions from
        # the run's baseline / configuration so the dynamic mask reflects
        # this scan's geometry, not just the polygon file's reference
        # positions.  Falls back gracefully if the resolver fails (e.g.
        # mocked test runs without a baseline).
        active_bs = "rod"
        bs_pos: dict | None = None
        saxs_geo = None
        try:
            from smi_tiled.loader import resolve_saxs_geometry
            saxs_geo = resolve_saxs_geometry(run)
            active_bs = saxs_geo.active_beamstop or "rod"
            bs_pos = saxs_geo.beamstop_pos_mm
        except Exception:
            pass

        bc_px = None
        if saxs_geo is not None:
            bc_px = (saxs_geo.beam_center_row_px, saxs_geo.beam_center_col_px)
        mask = make_saxs_mask_from_spec(
            image_shape=raw_shape,
            mask_path=resolved_mask_path,
            active_beamstop=active_bs,
            beamstop_pos_mm=bs_pos,
            beam_center_px=bc_px,
        )

        # AND in the per-frame WAXS-shadow occlusion (depends on
        # ``waxs_arc``).  This is the same shadow that
        # ``_integrate_saxs_batch`` applies during full reduction, so the
        # dynamic-mask overlay matches what the integrator actually sees.
        # Skipped silently if waxs_arc / beam_center are unavailable
        # (e.g. SAXS-only runs without the WAXS arc motor).
        try:
            waxs_arc = _smi_run_field_at(run, _WAXS_ARC_FIELD, frame_idx)
            beam_col = (
                saxs_geo.beam_center_col_px
                if saxs_geo is not None else None
            )
            if waxs_arc is not None and beam_col is not None:
                shadow = _make_waxs_shadow_mask(
                    raw_shape, [float(waxs_arc)], float(beam_col),
                )[0]
                mask = mask & shadow
        except Exception:
            pass

        if orient_for_display:
            mask = _orient_frame_for_display(mask, "saxs")
        return mask

    # WAXS: pull arc + bsx, derive bsx_ref via the SMI mechanical linkage.
    waxs_arc = _smi_run_field_at(run, _WAXS_ARC_FIELD, frame_idx)
    waxs_bsx = _smi_run_field_at(run, _WAXS_BSX_FIELD, frame_idx)
    waxs_bsx_ref = waxs_bsx - _BSX_PER_ARC_DEG_PUBLIC * waxs_arc

    mask_fn = make_waxs_mask_callable(
        resolved_mask_path,
        waxs_bsx_ref=waxs_bsx_ref,
        beamstop_max_abs_arc_deg=beamstop_max_abs_arc_deg,
    )
    # ``make_waxs_mask_callable`` returns a mask already in the display
    # orientation (rot90+fliplr applied inside ``make_mask_for_angle``).
    # Therefore for WAXS we never reapply ``orient_frame_for_display``.
    return mask_fn(raw_shape, theta_deg=waxs_arc, waxs_bsx=waxs_bsx)


# SAXS large-area masks (WAXS shadow + aperture)
_DEFAULT_SAXS_WAXS_SHADOW = {
    "enabled": True,
    "beam_visible_deg": 14.15,
    "clear_edge_deg": 18.0,
    "beam_visible_offset_px": -3.0,
    "edge_margin_px": 0.0,
}
_DEFAULT_SAXS_APERTURE = {
    "enabled": True,
    "agbh_ring_order": 5,
    "q_margin_fraction": 0.01,
    "q_cutoff": None,
}


def _silver_behenate_q_rings(max_q: float, max_order: int = 20) -> np.ndarray:
    D = 5.838  # nm
    orders = np.arange(1, max_order + 1)
    q_rings = 2.0 * np.pi / D * orders
    return q_rings[q_rings <= max_q]


def _make_waxs_shadow_mask(
    image_shape, waxs_arc, beam_center_col_px, **kwargs
) -> np.ndarray:
    """Return per-frame shadow boundary columns (1-D) for lazy evaluation.

    Instead of materializing an (n_frames, ny, nx) boolean array, returns a
    tuple (boundary_cols, always_clear_mask) where boundary_cols is shape
    (n_frames,) and the actual mask for frame i can be computed as:
        cols <= boundary_cols[i]  (broadcast over ny, nx)
    Frames where always_clear_mask[i] is True are fully unmasked.

    For backward compat when n_frames is small (<= 50), still returns the
    materialized array.
    """
    ny, nx = image_shape
    enabled = kwargs.get("enabled", True)
    if not enabled or waxs_arc is None:
        return np.ones((1, ny, nx), dtype=bool)
    waxs_arc = np.asarray(waxs_arc, dtype=float).reshape(-1)
    beam_visible_deg = float(kwargs.get("beam_visible_deg", 14.5))
    clear_edge_deg = float(kwargs.get("clear_edge_deg", 18.0))
    clear_span = clear_edge_deg - beam_visible_deg
    if abs(clear_span) < 1e-6:
        return np.ones((waxs_arc.size, ny, nx), dtype=bool)
    start_col = float(beam_center_col_px) + float(
        kwargs.get("beam_visible_offset_px", 0.0)
    )
    clear_col = float(nx - 1) - float(kwargs.get("edge_margin_px", 0.0))
    boundary_col = start_col + (
        (waxs_arc - beam_visible_deg) / clear_span
    ) * (clear_col - start_col)
    boundary_col = np.clip(boundary_col, -1.0, float(nx))

    # For small scans, materialize directly (cheap in memory).
    if waxs_arc.size <= 50:
        cols = np.arange(nx, dtype=float)[np.newaxis, np.newaxis, :]
        keep = cols <= boundary_col[:, np.newaxis, np.newaxis]
        keep[waxs_arc >= clear_edge_deg] = True
        return np.broadcast_to(keep, (waxs_arc.size, ny, nx)).copy()

    # For large scans, return a _LargeAreaMaskLazy object that computes
    # single-frame slices on demand without materializing the full array.
    return _LargeAreaMaskLazy(nx, ny, boundary_col, waxs_arc, clear_edge_deg)


class _LargeAreaMaskLazy:
    """Lazy per-frame shadow mask that mimics ndarray indexing [idx]."""

    __slots__ = ("nx", "ny", "boundary_col", "clear_mask", "_shape", "_cols")

    def __init__(self, nx, ny, boundary_col, waxs_arc, clear_edge_deg):
        self.nx = nx
        self.ny = ny
        self.boundary_col = boundary_col  # (n_frames,)
        self.clear_mask = waxs_arc >= clear_edge_deg  # (n_frames,) bool
        self._shape = (len(boundary_col), ny, nx)
        self._cols = np.arange(nx, dtype=float)

    @property
    def shape(self):
        return self._shape

    def __getitem__(self, idx):
        if self.clear_mask[idx]:
            return np.ones((self.ny, self.nx), dtype=bool)
        return np.broadcast_to(
            (self._cols <= self.boundary_col[idx])[np.newaxis, :],
            (self.ny, self.nx),
        )


def _make_aperture_mask(q_abs, **kwargs) -> np.ndarray:
    """Build an angular (in q) aperture mask: True = keep, False = occlude.

    The aperture models a physical occlusion (typically the WAXS detector arc
    blocking part of the SAXS detector's view).  Occlusion is defined in
    *angle* (here parameterised as a q-cutoff, which is angle for fixed
    wavelength).  Two key behaviours:

    * The cutoff is computed from AGB d-spacing alone, never from "the
      highest visible ring on this detector" — silently clamping to the
      highest visible ring caused over-occlusion at large SDD where the
      requested ring order wasn't even visible.
    * If the cutoff exceeds the detector's maximum q (i.e., the angular
      cutoff falls beyond the corner of the detector), no occlusion is
      applied.  This naturally turns the aperture off at long SDD (e.g.,
      9 m), and turns it on at short SDD (e.g., 2 m) where it physically
      represents what the WAXS detector blocks.
    """
    enabled = kwargs.get("enabled", True)
    if not enabled:
        return np.ones_like(q_abs, dtype=bool)
    q_cutoff = kwargs.get("q_cutoff")
    if q_cutoff is None:
        ring_order = max(int(kwargs.get("agbh_ring_order", 5)), 1)
        # AGB ring n in q-space (depends on d-spacing only — angle by Bragg
        # is wavelength-dependent, but q is the wavelength-independent
        # reciprocal-space coordinate, so q_n = 2π n / D works for any λ).
        D_nm = 5.838
        q_cutoff = (2.0 * np.pi / D_nm) * ring_order
        q_margin = float(kwargs.get("q_margin_fraction", 0.08))
        q_cutoff *= (1.0 + q_margin)
    # Off-detector cutoff → no occlusion (long-SDD short-q-range case).
    finite = np.isfinite(q_abs)
    detector_max_q = float(np.nanmax(q_abs)) if finite.any() else 0.0
    if float(q_cutoff) >= detector_max_q:
        return finite
    return finite & (q_abs <= float(q_cutoff))


def make_saxs_large_area_masks(
    image_shape, q_abs, waxs_arc, *, beam_center_col_px, waxs_shadow=None, aperture=None
):
    shadow_cfg = dict(_DEFAULT_SAXS_WAXS_SHADOW)
    if waxs_shadow:
        shadow_cfg.update(waxs_shadow)
    aperture_cfg = dict(_DEFAULT_SAXS_APERTURE)
    if aperture:
        aperture_cfg.update(aperture)
    shadow_mask = _make_waxs_shadow_mask(
        image_shape, waxs_arc, beam_center_col_px, **shadow_cfg
    )
    aperture_2d = _make_aperture_mask(q_abs, **aperture_cfg)

    # If shadow_mask is lazy (large scan), return a lazy composite.
    if isinstance(shadow_mask, _LargeAreaMaskLazy):
        combined = _LargeAreaMaskCombined(shadow_mask, aperture_2d)
        return combined, shadow_mask, aperture_2d

    # Small scan: materialize as before.
    aperture_mask = np.broadcast_to(
        aperture_2d, shadow_mask.shape
    ).copy()
    return shadow_mask & aperture_mask, shadow_mask, aperture_mask


class _LargeAreaMaskCombined:
    """Lazy combined shadow + aperture mask for large scans."""

    __slots__ = ("_shadow", "_aperture_2d", "_shape")

    def __init__(self, shadow: "_LargeAreaMaskLazy", aperture_2d: np.ndarray):
        self._shadow = shadow
        self._aperture_2d = aperture_2d
        self._shape = shadow.shape

    @property
    def shape(self):
        return self._shape

    def __getitem__(self, idx):
        return self._shadow[idx] & self._aperture_2d


# ===================================================================
# Histogram binning helpers
# ===================================================================


def _histogram2d_pixel_split(
    q2d: np.ndarray,
    chi2d: np.ndarray,
    img: np.ndarray,
    valid: np.ndarray,
    q_edges: np.ndarray,
    chi_edges: np.ndarray,
    pixel_splitting: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Histogram with optional pyFAI-style pixel splitting.

    When *pixel_splitting* > 1 each pixel is subdivided into an
    NxN grid of sub-pixels.  The q/chi position of each sub-pixel is
    estimated via gradient-based interpolation of the full q/chi maps,
    and the pixel intensity is fractionally distributed across the bins
    the sub-pixels fall into.

    Parameters
    ----------
    q2d, chi2d : 2-D arrays
        Per-pixel q and chi maps (same shape as *img*).
    img : 2-D array
        Intensity image (may contain NaN for invalid pixels).
    valid : 2-D bool array
        Mask of pixels to include in the histogram.
    q_edges, chi_edges : 1-D arrays
        Bin edges for the output histogram.
    pixel_splitting : int
        Number of sub-pixel divisions per axis.  1 (default) disables
        splitting and falls back to the standard single-point histogram.

    Returns
    -------
    I_hist, N_hist : 2-D arrays shaped (n_q, n_chi)
    """
    if pixel_splitting <= 1:
        q_sel = q2d[valid].ravel()
        chi_sel = chi2d[valid].ravel()
        I_sel = img[valid].ravel()
        I_hist, _, _ = np.histogram2d(
            q_sel, chi_sel, bins=[q_edges, chi_edges], weights=I_sel,
        )
        N_hist, _, _ = np.histogram2d(
            q_sel, chi_sel, bins=[q_edges, chi_edges],
        )
        return I_hist, N_hist

    # Gradient-based sub-pixel interpolation
    dq_dr = np.gradient(q2d, axis=0)
    dq_dc = np.gradient(q2d, axis=1)
    dchi_dr = np.gradient(chi2d, axis=0)
    dchi_dc = np.gradient(chi2d, axis=1)

    n = pixel_splitting
    offsets = np.linspace(-0.5 + 0.5 / n, 0.5 - 0.5 / n, n)
    weight = 1.0 / (n * n)

    n_q = len(q_edges) - 1
    n_chi = len(chi_edges) - 1
    I_hist = np.zeros((n_q, n_chi), dtype=float)
    N_hist = np.zeros((n_q, n_chi), dtype=float)

    for dr in offsets:
        for dc in offsets:
            q_sub = q2d + dr * dq_dr + dc * dq_dc
            chi_sub = chi2d + dr * dchi_dr + dc * dchi_dc

            sub_valid = valid & np.isfinite(q_sub) & np.isfinite(chi_sub)
            q_sel = q_sub[sub_valid].ravel()
            chi_sel = chi_sub[sub_valid].ravel()
            I_sel = img[sub_valid].ravel() * weight

            I_h, _, _ = np.histogram2d(
                q_sel, chi_sel, bins=[q_edges, chi_edges], weights=I_sel,
            )
            N_h, _, _ = np.histogram2d(
                q_sel, chi_sel, bins=[q_edges, chi_edges],
            )
            I_hist += I_h
            N_hist += N_h * weight

    return I_hist, N_hist


def _bin_indices(vals: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce ``np.histogram``'s bin assignment for explicit edges.

    For an array of edges, ``np.histogramdd`` computes
    ``searchsorted(edges, x, side='right')`` and decrements the entries that
    sit exactly on the rightmost edge (so the last bin is closed).  Values
    outside ``[edges[0], edges[-1]]`` — and NaNs, which sort to the end — are
    flagged invalid.  Reproducing that logic lets us precompute the bin index
    of every pixel **once** and reuse it across all frames (the geometry, and
    therefore the q/chi maps and bin edges, are identical for every frame).

    Returns
    -------
    bin_idx : int64 array
        Zero-based bin index (meaningful only where ``in_range``).
    in_range : bool array
        True where the value falls inside the histogram range.
    """
    nb = len(edges) - 1
    idx = np.searchsorted(edges, vals, side="right")
    idx[vals == edges[-1]] -= 1
    in_range = (idx >= 1) & (idx <= nb)
    return idx - 1, in_range


class _SplitBinPlan:
    """Precomputed pixel->(q,chi) bin mapping for fast repeated integration.

    ``_histogram2d_pixel_split`` recomputes, for every frame, which (q, chi)
    bin each pixel falls into — an expensive ``searchsorted`` over millions of
    pixels.  But for a stack of frames sharing one geometry the mapping never
    changes; only the per-pixel intensities and the per-frame validity mask do.

    This class precomputes, for each sub-pixel offset, a sparse
    ``(n_bins, n_pixels)`` matrix ``M`` whose row is the flattened output bin
    and whose column is the source pixel.  Integrating a frame is then a sparse
    mat-vec: ``I = M @ (img * valid)`` and ``N = M @ valid`` — folding the
    per-frame mask into the input vector.  Results are numerically equivalent
    to :func:`_histogram2d_pixel_split` (to within float summation order).

    The construction mirrors the gradient-based sub-pixel splitting of
    :func:`_histogram2d_pixel_split` exactly so outputs match for any
    ``pixel_splitting``.
    """

    def __init__(
        self,
        q2d: np.ndarray,
        chi2d: np.ndarray,
        q_edges: np.ndarray,
        chi_edges: np.ndarray,
        pixel_splitting: int = 1,
    ) -> None:
        from scipy import sparse

        self.n_q = len(q_edges) - 1
        self.n_chi = len(chi_edges) - 1
        self.n_bins = self.n_q * self.n_chi
        self.n_pixels = q2d.size
        self.shape = q2d.shape

        n = max(int(pixel_splitting), 1)
        self.weight = 1.0 / (n * n)

        if n == 1:
            sub_positions = [(q2d, chi2d)]
        else:
            dq_dr = np.gradient(q2d, axis=0)
            dq_dc = np.gradient(q2d, axis=1)
            dchi_dr = np.gradient(chi2d, axis=0)
            dchi_dc = np.gradient(chi2d, axis=1)
            offsets = np.linspace(-0.5 + 0.5 / n, 0.5 - 0.5 / n, n)
            sub_positions = []
            for dr in offsets:
                for dc in offsets:
                    q_sub = q2d + dr * dq_dr + dc * dq_dc
                    chi_sub = chi2d + dr * dchi_dr + dc * dchi_dc
                    sub_positions.append((q_sub, chi_sub))

        self._mats: list[Any] = []
        pix_idx = np.arange(self.n_pixels, dtype=np.int64)
        for q_sub, chi_sub in sub_positions:
            qb, q_ok = _bin_indices(q_sub.ravel(), q_edges)
            cb, c_ok = _bin_indices(chi_sub.ravel(), chi_edges)
            ok = q_ok & c_ok & np.isfinite(q_sub.ravel()) & np.isfinite(chi_sub.ravel())
            rows = (qb[ok] * self.n_chi + cb[ok]).astype(np.int64)
            cols = pix_idx[ok]
            data = np.ones(rows.size, dtype=np.float64)
            self._mats.append(
                sparse.csr_matrix(
                    (data, (rows, cols)), shape=(self.n_bins, self.n_pixels)
                )
            )

    def integrate_frame(
        self, img: np.ndarray, valid: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Integrate one frame. ``img``/``valid`` are 2-D, same shape as geometry.

        ``valid`` is the boolean per-frame mask (mask AND finite-image AND any
        dezinger flags).  Invalid pixels are zeroed in the weight vector so
        NaNs never reach the histogram.
        """
        valid_f = valid.ravel().astype(np.float64)
        wI = np.where(valid.ravel(), img.ravel().astype(np.float64), 0.0)
        I_hist = np.zeros(self.n_bins, dtype=np.float64)
        N_hist = np.zeros(self.n_bins, dtype=np.float64)
        for M in self._mats:
            I_hist += M.dot(wI)
            N_hist += M.dot(valid_f)
        # ``self.weight`` == 1.0 for pixel_splitting == 1 (single sub-pixel),
        # matching the original unweighted single-point histogram; for >1 it is
        # 1/(n*n), distributing each pixel's intensity/count across sub-pixels.
        I_hist *= self.weight
        N_hist *= self.weight
        return (
            I_hist.reshape(self.n_q, self.n_chi),
            N_hist.reshape(self.n_q, self.n_chi),
        )


def _qchi_and_iq(
    accum_I: np.ndarray,
    accum_N: np.ndarray,
    q_grid: np.ndarray,
    chi_grid: np.ndarray,
) -> dict[str, Any]:
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_I = np.where(accum_N > 0, accum_I / accum_N, np.nan)
    qchi = xr.Dataset(
        {
            "intensity": (("q", "chi"), mean_I),
            "counts": (("q", "chi"), accum_N),
        },
        coords={"q": q_grid, "chi": chi_grid},
    )
    total_N = accum_N.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        I_1d = np.where(
            total_N > 0,
            (np.nan_to_num(mean_I, nan=0.0) * accum_N).sum(axis=1) / total_N,
            np.nan,
        )
    iq = xr.Dataset(
        {"I": ("q", I_1d), "counts": ("q", total_N)},
        coords={"q": q_grid},
    )
    return {"q_chi": qchi, "iq": iq}


def _stack_qchi_frames(frame_qchi: list[xr.Dataset]) -> xr.Dataset:
    first = frame_qchi[0]
    return xr.Dataset(
        {
            "intensity": (
                ("frame", "q", "chi"),
                np.stack([ds["intensity"].values for ds in frame_qchi], axis=0),
            ),
            "counts": (
                ("frame", "q", "chi"),
                np.stack([ds["counts"].values for ds in frame_qchi], axis=0),
            ),
        },
        coords={
            "frame": np.arange(len(frame_qchi), dtype=int),
            "q": np.asarray(first["q"].values, dtype=float),
            "chi": np.asarray(first["chi"].values, dtype=float),
        },
    )


def _stack_iq_frames(frame_iq: list[xr.Dataset]) -> xr.Dataset:
    first = frame_iq[0]
    return xr.Dataset(
        {
            "I": (
                ("frame", "q"),
                np.stack([ds["I"].values for ds in frame_iq], axis=0),
            ),
            "counts": (
                ("frame", "q"),
                np.stack([ds["counts"].values for ds in frame_iq], axis=0),
            ),
        },
        coords={
            "frame": np.arange(len(frame_iq), dtype=int),
            "q": np.asarray(first["q"].values, dtype=float),
        },
    )


class _ZarrFrameWriter:
    """Stream per-frame ``(q, chi)`` maps to a zarr store on disk.

    For large parallel scans the full ``(frame, q, chi)`` stack is tens of GB
    and must not be held in RAM.  This writer creates the on-disk array layout
    up front (via an xarray template written with ``compute=False``, which
    writes coordinates + metadata but no bulk data), then fills one frame at a
    time with a direct zarr slice write — so peak memory is a single frame
    regardless of frame count.

    :meth:`dataset` returns the store opened with :func:`xarray.open_zarr`, i.e.
    a **lazy, dask-backed** ``xr.Dataset`` chunked one-frame-per-chunk.
    Consumers index/compute frames on demand without loading the whole stack.
    """

    def __init__(self, path, n_frames: int, q_grid: np.ndarray, chi_grid: np.ndarray):
        import dask.array as da
        import zarr

        self.path = str(path)
        nq, nchi = len(q_grid), len(chi_grid)
        chunks = (1, nq, nchi)
        template = xr.Dataset(
            {
                "intensity": (("frame", "q", "chi"),
                              da.zeros((n_frames, nq, nchi), chunks=chunks, dtype=float)),
                "counts": (("frame", "q", "chi"),
                           da.zeros((n_frames, nq, nchi), chunks=chunks, dtype=float)),
            },
            coords={
                "frame": np.arange(n_frames, dtype=int),
                "q": np.asarray(q_grid, dtype=float),
                "chi": np.asarray(chi_grid, dtype=float),
            },
        )
        # Writes coords + array metadata now; defers (skips) the lazy zeros.
        template.to_zarr(self.path, mode="w", compute=False)
        self._z = zarr.open(self.path, mode="r+")

    def write(self, idx: int, intensity_2d: np.ndarray, counts_2d: np.ndarray) -> None:
        self._z["intensity"][idx] = np.asarray(intensity_2d, dtype=float)
        self._z["counts"][idx] = np.asarray(counts_2d, dtype=float)

    def dataset(self) -> xr.Dataset:
        return xr.open_zarr(self.path)


# ===================================================================
# Merge utilities – ported from combined_reduce.py
# ===================================================================

def _interp_axis(source_axis, values, target_axis, fill_value):
    out = np.full(target_axis.shape, fill_value, dtype=float)
    finite = np.isfinite(values) & np.isfinite(source_axis)
    if np.count_nonzero(finite) == 0:
        return out
    x = np.asarray(source_axis[finite], dtype=float)
    y = np.asarray(values[finite], dtype=float)
    order = np.argsort(x)
    out = np.interp(target_axis, x[order], y[order], left=fill_value, right=fill_value)
    return out


def _empty_qchi_like(ref: xr.Dataset) -> xr.Dataset:
    """Return a zero-count q-chi dataset with the same grid as *ref*."""
    q = ref["q"].values
    chi = ref["chi"].values
    nq, nc = len(q), len(chi)
    return xr.Dataset(
        {
            "intensity": (("q", "chi"), np.full((nq, nc), np.nan)),
            "counts":    (("q", "chi"), np.zeros((nq, nc))),
        },
        coords={"q": q, "chi": chi},
    )


def _empty_iq_like(ref: xr.Dataset) -> xr.Dataset:
    """Return a NaN I(q) dataset with the same q grid as *ref*."""
    q = ref["q"].values
    return xr.Dataset(
        {
            "I":      ("q", np.full(len(q), np.nan)),
            "counts": ("q", np.zeros(len(q))),
        },
        coords={"q": q},
    )


def merge_q_chi_weighted(
    saxs_qchi: xr.Dataset | None,
    waxs_qchi: xr.Dataset | None,
    n_q: int = 1000,
    n_chi: int = 360,
) -> xr.Dataset | None:
    """Merge SAXS and WAXS q-chi maps on a common grid with count-weighting.

    Returns None if both inputs are None. Returns the single detector's data
    (re-gridded) if only one is present.
    """
    if saxs_qchi is None and waxs_qchi is None:
        return None
    # Single-detector passthrough: use the available one for both slots
    if saxs_qchi is None:
        saxs_qchi = _empty_qchi_like(waxs_qchi)
    if waxs_qchi is None:
        waxs_qchi = _empty_qchi_like(saxs_qchi)

    saxs_q = np.asarray(saxs_qchi["q"].values, dtype=float)
    waxs_q = np.asarray(waxs_qchi["q"].values, dtype=float)
    q_min = min(float(np.nanmin(saxs_q)), float(np.nanmin(waxs_q)))
    q_max = max(float(np.nanmax(saxs_q)), float(np.nanmax(waxs_q)))
    q_grid = np.linspace(q_min, q_max, n_q)

    saxs_chi = np.asarray(saxs_qchi["chi"].values, dtype=float)
    waxs_chi = np.asarray(waxs_qchi["chi"].values, dtype=float)
    chi_min = min(float(np.nanmin(saxs_chi)), float(np.nanmin(waxs_chi)))
    chi_max = max(float(np.nanmax(saxs_chi)), float(np.nanmax(waxs_chi)))
    chi_grid = np.linspace(chi_min, chi_max, n_chi)

    saxs_I = np.asarray(saxs_qchi["intensity"].values, dtype=float)
    saxs_N = np.asarray(saxs_qchi["counts"].values, dtype=float)
    waxs_I = np.asarray(waxs_qchi["intensity"].values, dtype=float)
    waxs_N = np.asarray(waxs_qchi["counts"].values, dtype=float)

    def _regrid_2d(src_q, src_chi, data, target_q, target_chi, fill):
        """Regrid a 2D (q, chi) array onto a new grid via nearest-neighbor."""
        from scipy.interpolate import RegularGridInterpolator
        finite_data = np.where(np.isfinite(data), data, fill)
        interp = RegularGridInterpolator(
            (src_q, src_chi), finite_data,
            method="nearest", bounds_error=False, fill_value=fill,
        )
        tq, tc = np.meshgrid(target_q, target_chi, indexing="ij")
        return interp((tq, tc))

    s_I_interp = _regrid_2d(saxs_q, saxs_chi, saxs_I, q_grid, chi_grid, np.nan)
    s_N_interp = _regrid_2d(saxs_q, saxs_chi, saxs_N, q_grid, chi_grid, 0.0)
    w_I_interp = _regrid_2d(waxs_q, waxs_chi, waxs_I, q_grid, chi_grid, np.nan)
    w_N_interp = _regrid_2d(waxs_q, waxs_chi, waxs_N, q_grid, chi_grid, 0.0)

    total_N = s_N_interp + w_N_interp
    with np.errstate(divide="ignore", invalid="ignore"):
        merged_I = np.where(
            total_N > 0,
            (np.nan_to_num(s_I_interp, nan=0.0) * s_N_interp
             + np.nan_to_num(w_I_interp, nan=0.0) * w_N_interp) / total_N,
            np.nan,
        )

    return xr.Dataset(
        {
            "intensity": (("q", "chi"), merged_I),
            "counts": (("q", "chi"), total_N),
            "saxs_intensity": (("q", "chi"), s_I_interp),
            "saxs_counts": (("q", "chi"), s_N_interp),
            "waxs_intensity": (("q", "chi"), w_I_interp),
            "waxs_counts": (("q", "chi"), w_N_interp),
        },
        coords={"q": q_grid, "chi": chi_grid},
    )


def merge_iq_profiles(
    merged_qchi: xr.Dataset | None,
    saxs_iq: xr.Dataset | None,
    waxs_iq: xr.Dataset | None,
) -> xr.Dataset | None:
    """Produce merged I(q) by azimuthal integration of the merged q-chi map."""
    if merged_qchi is None:
        return None
    q_grid = np.asarray(merged_qchi["q"].values, dtype=float)
    merged_I = np.asarray(merged_qchi["intensity"].values, dtype=float)
    merged_N = np.asarray(merged_qchi["counts"].values, dtype=float)

    total_N = merged_N.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        I_1d = np.where(
            total_N > 0,
            (np.nan_to_num(merged_I, nan=0.0) * merged_N).sum(axis=1) / total_N,
            np.nan,
        )

    if saxs_iq is not None:
        saxs_q_src = np.asarray(saxs_iq["q"].values, dtype=float)
        saxs_I_src = np.asarray(saxs_iq["I"].values, dtype=float)
        saxs_I_interp = _interp_axis(saxs_q_src, saxs_I_src, q_grid, np.nan)
    else:
        saxs_I_interp = np.full(len(q_grid), np.nan)

    if waxs_iq is not None:
        waxs_q_src = np.asarray(waxs_iq["q"].values, dtype=float)
        waxs_I_src = np.asarray(waxs_iq["I"].values, dtype=float)
        waxs_I_interp = _interp_axis(waxs_q_src, waxs_I_src, q_grid, np.nan)
    else:
        waxs_I_interp = np.full(len(q_grid), np.nan)

    return xr.Dataset(
        {
            "I": ("q", I_1d),
            "counts": ("q", total_N),
            "saxs_I": ("q", saxs_I_interp),
            "waxs_I": ("q", waxs_I_interp),
        },
        coords={"q": q_grid},
    )


def _build_per_frame_iq(
    merged_iq: xr.Dataset | None,
    saxs_result: dict[str, Any] | None,
    waxs_result: dict[str, Any] | None,
    scan_info: dict[str, Any] | None = None,
) -> xr.Dataset | None:
    """Build per-frame I(q) Dataset on the same q grid as merged_iq.

    Combines per-frame SAXS and WAXS I(q) via interpolation onto the
    merged q grid, then produces a count-weighted merge per frame.

    If *scan_info* is provided and contains per-frame primary-stream
    scalars (step_candidates with a 'values' key), they are attached
    as data variables on the (frame,) dimension.

    Returns
    -------
    xr.Dataset with dims (frame, q) and variables I, saxs_I, waxs_I,
    plus any per-frame primary scalars, or None if no per-frame data
    is available.
    """
    if merged_iq is None:
        return None

    q_grid = np.asarray(merged_iq["q"].values, dtype=float)
    n_q = len(q_grid)

    saxs_iq_frames = saxs_result["iq_frames"] if saxs_result else None
    waxs_iq_frames = waxs_result["iq_frames"] if waxs_result else None

    if saxs_iq_frames is None and waxs_iq_frames is None:
        return None

    # Determine number of frames from whichever detector is present
    if saxs_iq_frames is not None and waxs_iq_frames is not None:
        n_frames = max(
            len(saxs_iq_frames["frame"]),
            len(waxs_iq_frames["frame"]),
        )
    elif saxs_iq_frames is not None:
        n_frames = len(saxs_iq_frames["frame"])
    else:
        n_frames = len(waxs_iq_frames["frame"])

    saxs_I_2d = np.full((n_frames, n_q), np.nan)
    waxs_I_2d = np.full((n_frames, n_q), np.nan)

    if saxs_iq_frames is not None:
        saxs_q_src = np.asarray(saxs_iq_frames["q"].values, dtype=float)
        saxs_I_src = np.asarray(saxs_iq_frames["I"].values, dtype=float)
        n_saxs = saxs_I_src.shape[0]
        for fi in range(min(n_saxs, n_frames)):
            saxs_I_2d[fi] = _interp_axis(saxs_q_src, saxs_I_src[fi], q_grid, np.nan)

    if waxs_iq_frames is not None:
        waxs_q_src = np.asarray(waxs_iq_frames["q"].values, dtype=float)
        waxs_I_src = np.asarray(waxs_iq_frames["I"].values, dtype=float)
        n_waxs = waxs_I_src.shape[0]
        for fi in range(min(n_waxs, n_frames)):
            waxs_I_2d[fi] = _interp_axis(waxs_q_src, waxs_I_src[fi], q_grid, np.nan)

    # Count-weighted merge per frame (same logic as merge_iq_profiles)
    saxs_N_2d = np.zeros((n_frames, n_q), dtype=float)
    waxs_N_2d = np.zeros((n_frames, n_q), dtype=float)

    if saxs_iq_frames is not None:
        saxs_counts_src = np.asarray(saxs_iq_frames["counts"].values, dtype=float)
        saxs_q_src = np.asarray(saxs_iq_frames["q"].values, dtype=float)
        n_saxs = saxs_counts_src.shape[0]
        for fi in range(min(n_saxs, n_frames)):
            saxs_N_2d[fi] = _interp_axis(saxs_q_src, saxs_counts_src[fi], q_grid, 0.0)

    if waxs_iq_frames is not None:
        waxs_counts_src = np.asarray(waxs_iq_frames["counts"].values, dtype=float)
        waxs_q_src = np.asarray(waxs_iq_frames["q"].values, dtype=float)
        n_waxs = waxs_counts_src.shape[0]
        for fi in range(min(n_waxs, n_frames)):
            waxs_N_2d[fi] = _interp_axis(waxs_q_src, waxs_counts_src[fi], q_grid, 0.0)

    total_N = saxs_N_2d + waxs_N_2d
    with np.errstate(divide="ignore", invalid="ignore"):
        merged_I_2d = np.where(
            total_N > 0,
            (np.nan_to_num(saxs_I_2d, nan=0.0) * saxs_N_2d
             + np.nan_to_num(waxs_I_2d, nan=0.0) * waxs_N_2d) / total_N,
            np.nan,
        )

    data_vars: dict[str, Any] = {
        "I": (("frame", "q"), merged_I_2d),
        "saxs_I": (("frame", "q"), saxs_I_2d),
        "waxs_I": (("frame", "q"), waxs_I_2d),
    }

    # Attach per-frame primary-stream scalars as data variables
    if scan_info is not None:
        for cand in scan_info.get("step_candidates", []):
            vals = cand.get("values")
            if vals is None:
                continue
            vals = np.asarray(vals, dtype=float)
            if vals.shape[0] == n_frames:
                data_vars[cand["name"]] = ("frame", vals)

    return xr.Dataset(
        data_vars,
        coords={
            "q": q_grid,
            "frame": np.arange(n_frames, dtype=int),
        },
    )


# -------------------------------------------------------------------
# Multi-scan merging
# -------------------------------------------------------------------

def merge_multiple_qchi(
    datasets: list[xr.Dataset],
    n_q: int = 2000,
    n_chi: int = 360,
) -> xr.Dataset:
    """Count-weighted merge of N ``(q, chi)`` datasets onto a common grid.

    Each input must have ``intensity`` and ``counts`` variables with
    dimensions ``(q, chi)``.  The output grid spans the union of all
    input q/chi ranges.
    """
    from scipy.interpolate import RegularGridInterpolator

    if not datasets:
        raise ValueError("Need at least one dataset to merge")
    if len(datasets) == 1:
        return datasets[0]

    all_q = [np.asarray(ds["q"].values, dtype=float) for ds in datasets]
    all_chi = [np.asarray(ds["chi"].values, dtype=float) for ds in datasets]
    q_min = min(float(q.min()) for q in all_q)
    q_max = max(float(q.max()) for q in all_q)
    chi_min = min(float(c.min()) for c in all_chi)
    chi_max = max(float(c.max()) for c in all_chi)
    q_grid = np.linspace(q_min, q_max, n_q)
    chi_grid = np.linspace(chi_min, chi_max, n_chi)
    tq, tc = np.meshgrid(q_grid, chi_grid, indexing="ij")
    pts = (tq, tc)

    accum_IN = np.zeros((n_q, n_chi), dtype=float)
    accum_N = np.zeros((n_q, n_chi), dtype=float)

    for ds, src_q, src_chi in zip(datasets, all_q, all_chi):
        I_src = np.asarray(ds["intensity"].values, dtype=float)
        N_src = np.asarray(ds["counts"].values, dtype=float)

        N_fill = np.where(np.isfinite(N_src), N_src, 0.0)
        IN_src = np.where(np.isfinite(I_src), I_src * N_fill, 0.0)

        interp_N = RegularGridInterpolator(
            (src_q, src_chi), N_fill,
            method="nearest", bounds_error=False, fill_value=0.0,
        )
        interp_IN = RegularGridInterpolator(
            (src_q, src_chi), IN_src,
            method="nearest", bounds_error=False, fill_value=0.0,
        )
        accum_N += interp_N(pts)
        accum_IN += interp_IN(pts)

    with np.errstate(divide="ignore", invalid="ignore"):
        merged_I = np.where(accum_N > 0, accum_IN / accum_N, np.nan)

    return xr.Dataset(
        {
            "intensity": (("q", "chi"), merged_I),
            "counts": (("q", "chi"), accum_N),
        },
        coords={"q": q_grid, "chi": chi_grid},
    )


def merge_multiple_iq(
    merged_qchi: xr.Dataset,
) -> xr.Dataset:
    """Azimuthally average a merged q-chi map into I(q)."""
    q_grid = np.asarray(merged_qchi["q"].values, dtype=float)
    merged_I = np.asarray(merged_qchi["intensity"].values, dtype=float)
    merged_N = np.asarray(merged_qchi["counts"].values, dtype=float)

    total_N = merged_N.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        I_1d = np.where(
            total_N > 0,
            (np.nan_to_num(merged_I, nan=0.0) * merged_N).sum(axis=1) / total_N,
            np.nan,
        )
    return xr.Dataset(
        {"I": ("q", I_1d), "counts": ("q", total_N)},
        coords={"q": q_grid},
    )


def merge_reduction_results(
    *results: "CombinedReductionResult",
    n_q: int = 2000,
    n_chi: int = 360,
) -> Tuple[xr.Dataset, xr.Dataset]:
    """Merge multiple :class:`CombinedReductionResult` objects.

    Collects the per-scan merged q-chi maps and combines them with
    count-weighting.

    Returns ``(merged_qchi, merged_iq)``.
    """
    qchi_list = [r.merged_qchi for r in results if r.merged_qchi is not None]
    if not qchi_list:
        raise ValueError("No merged q-chi data in any of the results")
    merged_qchi = merge_multiple_qchi(qchi_list, n_q=n_q, n_chi=n_chi)
    merged_iq = merge_multiple_iq(merged_qchi)
    return merged_qchi, merged_iq


# ===================================================================
# Result dataclasses
# ===================================================================

@dataclass(frozen=True)
class CombinedReductionResult:
    """Reduced SAXS+WAXS data for one SMI run.

    Attributes
    ----------
    uid : str
        Source tiled run UID.
    scan_info : dict
        Output of :func:`SMISWAXSLoader.infer_detectors_and_steps` —
        contains step_candidates, n_frames, sample_name, etc.
    saxs, waxs : dict or None
        Per-detector intermediate products (q-maps, per-frame qchi, iq).
        ``None`` if that detector wasn't present in the scan.
    merged_qchi : xr.Dataset or None
        Count-weighted SAXS+WAXS merge on a common ``(q, chi)`` grid.
        Variables: ``intensity, counts, saxs_intensity, saxs_counts,
        waxs_intensity, waxs_counts``.
    merged_iq : xr.Dataset or None
        Azimuthally-averaged 1-D profile on the merged q grid.
        Variables: ``I, counts, saxs_I, waxs_I``.
    per_frame_iq : xr.Dataset or None
        Per-scan-step I(q).  Dims ``(frame, q)``.  Variables: ``I, saxs_I,
        waxs_I`` plus any per-frame primary scalars attached as data vars.

    Notes
    -----
    For compatibility with PyHyperScattering accessors (``da.rsoxs.*``,
    ``da.fit.*``), use :meth:`to_dataarray` to extract a single-variable
    ``xr.DataArray`` view of one of the merged datasets.
    """
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
    #: Optional derived products attached by the smi_tiled.derived stages.
    per_frame_qchi: dict[str, xr.Dataset] | None = None
    line_cuts: dict[str, xr.Dataset] | None = None
    peak_fits: xr.Dataset | None = None

    def to_dataarray(
        self,
        key: str = "merged_iq",
        variable: str | None = None,
    ) -> xr.DataArray:
        """Extract a single ``xr.DataArray`` view of one merged product.

        The merged outputs are stored as ``xr.Dataset`` (so intensity and
        counts ride together).  PyHyperScattering's xarray accessors
        (``da.rsoxs.slice_chi``, ``da.fit.apply``, ...) operate on
        ``xr.DataArray``, so this helper extracts the intensity variable
        and attaches reduction provenance as attrs.

        Parameters
        ----------
        key : {'merged_iq', 'merged_qchi', 'per_frame_iq'}
            Which reduction product to extract.
        variable : str, optional
            Variable name within the dataset to return.  Defaults to
            ``'I'`` for the 1-D products and ``'intensity'`` for q-chi.

        Returns
        -------
        xr.DataArray
            With dims appropriate to the requested product and attrs
            containing ``uid, scan_id, sample_name, geometry,
            incident_angle_deg`` plus the original dataset's attrs.

        Raises
        ------
        ValueError
            If ``key`` does not name a reduction product, or that product
            is ``None`` (e.g. requesting merged_qchi for a SAXS-only scan
            that wasn't merged).
        """
        sources = {
            "merged_iq":    (self.merged_iq, "I"),
            "merged_qchi":  (self.merged_qchi, "intensity"),
            "per_frame_iq": (self.per_frame_iq, "I"),
        }
        if key not in sources:
            raise ValueError(
                f"key must be one of {list(sources)}, got {key!r}"
            )
        ds, default_var = sources[key]
        if ds is None:
            raise ValueError(
                f"{key!r} is None — that reduction product wasn't produced "
                f"(likely because the relevant detector was absent from the scan)."
            )
        var = variable or default_var
        if var not in ds:
            raise ValueError(
                f"Variable {var!r} not in {key!r} dataset.  "
                f"Available: {list(ds.data_vars)}"
            )
        da = ds[var]
        sample_name = (self.scan_info or {}).get("sample_name", "")
        scan_id = (self.scan_info or {}).get("scan_id")
        new_attrs = dict(da.attrs)
        new_attrs.update({
            "uid":                self.uid,
            "scan_id":            scan_id,
            "sample_name":        sample_name,
            "geometry":           self.geometry,
            "incident_angle_deg": self.incident_angle_deg,
            "source":             key,
        })
        return da.assign_attrs(new_attrs)


@dataclass(frozen=True)
class GIReductionResult:
    """Result of a grazing-incidence WAXS reduction.

    Attributes
    ----------
    uid : str
        Tiled run UID.
    sample_name : str
        Sample name from the start document.
    scan_motor : str
        Name of the scanned motor (e.g. ``'piezo_th'``).
    scan_motor_values : np.ndarray
        Per-frame values of the scanned motor.
    alpha_i_deg : np.ndarray
        Per-frame incident angle (degrees).
    alpha_i_source : str
        Description of how alpha_i was determined.
    qxy_grid : np.ndarray
        1-D q_xy bin centres (nm\ :sup:`-1`).
    qz_grid : np.ndarray
        1-D q_z bin centres (nm\ :sup:`-1`).
    frames : list[np.ndarray]
        Per-frame I(qxy, qz) images (shape ``(n_qxy, n_qz)``).
    summed : np.ndarray
        Averaged I(qxy, qz) over all frames.
    q_chi_frames : xr.Dataset or None
        Per-frame I(qxy, qz) as xr.Dataset with dims (frame, qxy, qz).
    summed_ds : xr.Dataset or None
        Summed I(qxy, qz) as xr.Dataset with dims (qxy, qz) and
        variables ``intensity`` and ``counts``.
    timing : dict[str, float] | None
        Timing breakdown.
    """
    uid: str
    sample_name: str
    scan_motor: str
    scan_motor_values: np.ndarray
    alpha_i_deg: np.ndarray
    alpha_i_source: str
    qxy_grid: np.ndarray
    qz_grid: np.ndarray
    frames: list  # list[np.ndarray]
    summed: np.ndarray
    q_chi_frames: xr.Dataset | None = None
    summed_ds: xr.Dataset | None = None
    timing: dict[str, float] | None = None
    #: Optional derived products attached by the smi_tiled.derived stages.
    line_cuts: dict[str, xr.Dataset] | None = None
    peak_fits: xr.Dataset | None = None

    # -- Line cut helpers --------------------------------------------------

    def line_cut_qxy(self, qz_center: float, qz_width: float = 0.05,
                     frame: int | None = None) -> tuple[np.ndarray, np.ndarray]:
        """I(qxy) at constant qz +/- width.  *frame*=None uses the sum."""
        img = self.summed if frame is None else self.frames[frame]
        mask = (self.qz_grid >= qz_center - qz_width) & (self.qz_grid <= qz_center + qz_width)
        if not mask.any():
            return self.qxy_grid, np.full_like(self.qxy_grid, np.nan)
        return self.qxy_grid, np.nanmean(img[:, mask], axis=1)

    def line_cut_qz(self, qxy_center: float, qxy_width: float = 0.1,
                    frame: int | None = None) -> tuple[np.ndarray, np.ndarray]:
        """I(qz) at constant qxy +/- width.  *frame*=None uses the sum."""
        img = self.summed if frame is None else self.frames[frame]
        mask = (self.qxy_grid >= qxy_center - qxy_width) & (self.qxy_grid <= qxy_center + qxy_width)
        if not mask.any():
            return self.qz_grid, np.full_like(self.qz_grid, np.nan)
        return self.qz_grid, np.nanmean(img[mask, :], axis=0)


# ===================================================================
# Grazing-incidence helpers
# ===================================================================

def lab_to_sample_frame(
    qx_lab: np.ndarray,
    qy_lab: np.ndarray,
    qz_lab: np.ndarray,
    alpha_i_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rotate lab-frame q into the sample frame for grazing incidence.

    Parameters
    ----------
    qx_lab, qy_lab, qz_lab : ndarray
        Lab-frame q components (horizontal, vertical-up, along-beam).
    alpha_i_deg : float
        Incident angle in degrees.

    Returns
    -------
    (qx_s, qy_s, qz_s) where qz_s is along the surface normal.
    """
    ai = np.deg2rad(alpha_i_deg)
    cos_ai = np.cos(ai)
    sin_ai = np.sin(ai)
    qx_s = qx_lab
    qy_s = qy_lab * sin_ai - qz_lab * cos_ai
    qz_s = qy_lab * cos_ai + qz_lab * sin_ai
    return qx_s, qy_s, qz_s


import re as _re

_AI_PATTERNS = [
    _re.compile(r'(?:^|[_\-])ai[_\-]?([0-9]+\.?[0-9]*)', _re.IGNORECASE),
    _re.compile(r'alpha[_]?i?[_\-]?([0-9]+\.?[0-9]*)', _re.IGNORECASE),
    _re.compile(r'incident[_\-]?([0-9]+\.?[0-9]*)', _re.IGNORECASE),
]


def parse_incident_angle_from_string(s: str) -> float | None:
    """Extract incident angle (degrees) from a sample-name string.

    Recognized patterns (case-insensitive):
    ``ai0.12``, ``ai_0.12``, ``alpha0.12``, ``alpha_i0.12``, ``incident0.12``.

    Returns ``float`` or ``None``.
    """
    for pat in _AI_PATTERNS:
        m = pat.search(s)
        if m:
            val = float(m.group(1))
            if 0 < val < 90:
                return val
    return None


def find_incident_angle(
    run,
    n_frames: int,
    manual_override: float | None = None,
    theta_offset: float = 0.0,
) -> tuple[np.ndarray, str]:
    """Determine per-frame incident angle (degrees).

    Thin wrapper over :func:`smi_tiled.loader.resolve_incident_angle`, which
    resolves the angle per-frame with the priority:

    0. ``manual_override`` if not None.
    1. **Primary measured motor** ``stage_th + piezo_th`` ``+ theta_offset``.
    2. **target_file_name** parsed ``_ai`` (per-frame, no offset).
    3. **sample_name** parsed ``_ai``/``_th`` (no offset).
    4. **Baseline** ``stage_th + piezo_th`` ``+ theta_offset`` (last resort).

    The zero offset is applied only to motor-derived paths (1 and 4).

    Returns ``(alpha_i_array, source_description)``.  Raises ``RuntimeError``
    when no source resolves.
    """
    from smi_tiled.loader import resolve_incident_angle

    ai, source = resolve_incident_angle(
        run, n_frames,
        manual_override=manual_override,
        theta_offset=theta_offset,
    )
    if ai is None:
        raise RuntimeError(
            "Cannot determine incident angle. Pass incident_angle_deg manually."
        )
    return ai, source


def integrate_waxs_gi(
    waxs_raw: xr.DataArray,
    mask_fn,
    alpha_i_deg: np.ndarray,
    n_qxy: int = 500,
    n_qz: int = 500,
    cal: WAXSCalibration | None = None,
    dezinger_threshold: float | None = None,
    dezinger_kernel: int = 5,
    pixel_splitting: int = 1,
    progress: ProgressCallback | None = None,
) -> dict[str, Any]:
    """WAXS GI reduction: bin each frame into (qxy, qz) in the sample frame.

    Parameters
    ----------
    waxs_raw : xr.DataArray
        Raw WAXS images from ``TiledSMISWAXSLoader.loadSingleImage``.
    mask_fn : callable or None
        ``mask_fn(image_shape_raw, theta_deg, waxs_bsx) -> bool mask``.
    alpha_i_deg : array-like
        Per-frame incident angle (degrees).
    n_qxy, n_qz : int
        Grid dimensions.
    cal : WAXSCalibration or None
        Detector calibration.  None uses ``_DEFAULT_CAL``.
    dezinger_threshold, dezinger_kernel
        Hot-pixel rejection parameters.
    pixel_splitting : int
        Number of sub-pixel divisions per axis for fractional pixel
        splitting during histogram binning.  1 (default) disables splitting.

    Returns
    -------
    dict with keys ``qxy_grid``, ``qz_grid``, ``frames``, ``summed``,
    ``q_chi_frames``, ``iq_frames``.
    """
    if cal is None:
        cal = WAXSCalibration(**_DEFAULT_CAL)

    images = np.asarray(waxs_raw.values, dtype=float)
    if images.ndim == 2:
        images = images[np.newaxis, :, :]
    n_frames = images.shape[0]

    arc_angles = np.asarray(
        waxs_raw.coords[waxs_raw.dims[0]].values, dtype=float,
    )
    bsx_per_frame = np.asarray(
        waxs_raw.attrs.get("smi_waxs_bsx_per_frame", [0.0] * n_frames),
        dtype=float,
    )
    alpha_arr = np.asarray(alpha_i_deg, dtype=float)
    if alpha_arr.size == 1:
        alpha_arr = np.full(n_frames, float(alpha_arr))

    # Build detector geometry (constant arc for GI scans)
    img_0_rot, _ = rotate_image_and_mask(images[0], k=cal.rotation_k)
    rot_shape = img_0_rot.shape
    theta_arc = float(arc_angles[0])
    bc = cal.beam_center_at_angle(theta_arc)
    det = MultiPanelArcDetector(
        image_shape=rot_shape,
        panel_specs=cal.make_panel_specs(),
        wavelength_nm=cal.wavelength_nm,
        pixel_size_mm=cal.pixel_size_mm,
        sample_distance_mm=cal.sample_distance_mm,
        beam_center_px=bc,
        theta_zero_deg=cal.theta_zero_deg,
        sample_offset_x_mm=cal.sample_offset_x_mm,
        sample_offset_z_mm=cal.sample_offset_z_mm,
    )
    qds = det.qmap(theta_arc)
    # Apply the same sign conventions as the transmission code so that
    # "up" in the lab frame is positive qy.
    qx_lab = cal.q_horizontal_sign * np.asarray(qds["qx"].values, dtype=float)
    qy_lab = cal.q_vertical_sign * np.asarray(qds["qy"].values, dtype=float)
    qz_lab_arr = np.asarray(qds["qz"].values, dtype=float)

    # --- First pass: global qxy/qz range ---
    _qxy_mn, _qxy_mx, _qz_mn, _qz_mx = [], [], [], []
    for fi in range(n_frames):
        qx_s, qy_s, qz_s = lab_to_sample_frame(
            qx_lab, qy_lab, qz_lab_arr, alpha_arr[fi],
        )
        qxy = np.sign(qx_s) * np.sqrt(qx_s ** 2 + qy_s ** 2)
        ok = np.isfinite(qxy) & np.isfinite(qz_s)
        if ok.any():
            _qxy_mn.append(float(np.nanmin(qxy[ok])))
            _qxy_mx.append(float(np.nanmax(qxy[ok])))
            _qz_mn.append(float(np.nanmin(qz_s[ok])))
            _qz_mx.append(float(np.nanmax(qz_s[ok])))

    qxy_edges = np.linspace(min(_qxy_mn), max(_qxy_mx), n_qxy + 1)
    qz_edges = np.linspace(min(_qz_mn), max(_qz_mx), n_qz + 1)
    qxy_grid = 0.5 * (qxy_edges[:-1] + qxy_edges[1:])
    qz_grid = 0.5 * (qz_edges[:-1] + qz_edges[1:])

    # --- Second pass: histogram each frame ---
    accum_I = np.zeros((n_qxy, n_qz), dtype=float)
    accum_N = np.zeros((n_qxy, n_qz), dtype=float)
    frame_maps: list[np.ndarray] = []

    for fi in range(n_frames):
        theta_f = float(arc_angles[fi])
        ai = alpha_arr[fi]
        bsx = float(bsx_per_frame[fi]) if fi < len(bsx_per_frame) else 0.0

        img_rot, _ = rotate_image_and_mask(images[fi], k=cal.rotation_k)

        mask_rot = None
        if mask_fn is not None:
            try:
                mask_rot = mask_fn(images[fi].shape, theta_f, bsx)
            except Exception as exc:
                warnings.warn(
                    f"mask_fn failed for frame {fi}: {exc}", stacklevel=2,
                )

        if dezinger_threshold is not None:
            dz = dezinger(
                img_rot, kernel_size=dezinger_kernel,
                threshold=dezinger_threshold,
            )
            mask_rot = (mask_rot & dz) if mask_rot is not None else dz

        qx_s, qy_s, qz_s = lab_to_sample_frame(
            qx_lab, qy_lab, qz_lab_arr, ai,
        )
        qxy = np.sign(qx_s) * np.sqrt(qx_s ** 2 + qy_s ** 2)

        valid = np.isfinite(qxy) & np.isfinite(qz_s) & np.isfinite(img_rot)
        if mask_rot is not None:
            valid &= mask_rot

        I_hist, N_hist = _histogram2d_pixel_split(
            qxy, qz_s, img_rot, valid, qxy_edges, qz_edges,
            pixel_splitting=pixel_splitting,
        )
        accum_I += I_hist
        accum_N += N_hist

        with np.errstate(invalid="ignore", divide="ignore"):
            frame_maps.append(np.where(N_hist > 0, I_hist / N_hist, np.nan))

        if progress is not None:
            progress("gi_integrate", fi + 1, n_frames)

    with np.errstate(invalid="ignore", divide="ignore"):
        summed = np.where(accum_N > 0, accum_I / accum_N, np.nan)

    # Build xr.Dataset outputs paralleling the transmission path
    # Per-frame qxy-qz maps as xr.Dataset
    frame_qchi: list[xr.Dataset] = []
    for fmap in frame_maps:
        frame_qchi.append(xr.Dataset(
            {
                "intensity": (("qxy", "qz"), fmap),
            },
            coords={"qxy": qxy_grid, "qz": qz_grid},
        ))

    # Stacked per-frame dataset
    q_chi_frames = xr.Dataset(
        {
            "intensity": (
                ("frame", "qxy", "qz"),
                np.stack(frame_maps, axis=0),
            ),
        },
        coords={
            "frame": np.arange(n_frames, dtype=int),
            "qxy": qxy_grid,
            "qz": qz_grid,
        },
    )

    # Summed (averaged) dataset
    summed_ds = xr.Dataset(
        {
            "intensity": (("qxy", "qz"), summed),
            "counts": (("qxy", "qz"), accum_N),
        },
        coords={"qxy": qxy_grid, "qz": qz_grid},
    )

    return {
        "qxy_grid": qxy_grid,
        "qz_grid": qz_grid,
        "frames": frame_maps,
        "summed": summed,
        "q_chi_frames": q_chi_frames,
        "summed_ds": summed_ds,
    }


# ===================================================================
# Grazing-incidence reduction entry point
# ===================================================================

def reduce_smi_gi(
    uid: str,
    tiled_uri: str = "https://tiled.nsls2.bnl.gov",
    catalog: str = "smi/migration",
    waxs_mask: "str | Path | dict | None" = None,
    waxs_mask_path: "str | Path | None" = None,
    n_qxy: int = 500,
    n_qz: int = 500,
    incident_angle_deg: float | None = None,
    theta_offset: float = -0.5,
    waxs_beam_col_per_arc_deg: float = 0.08,
    beamstop_max_abs_arc_deg: float = 15.0,
    dezinger_threshold: float | None = 30000.0,
    dezinger_kernel: int = 5,
    waxs_cal_overrides: dict[str, Any] | None = None,
    image_cache_path: str | Path | None = "auto",
    populate_disk_cache: bool = True,
    pixel_splitting: int = 1,
    progress: ProgressCallback | None = None,
    # ---- Optional derived-analysis stages (smi_tiled.derived) ----
    line_cuts: "Sequence[Any] | None" = None,
) -> GIReductionResult:
    """Full grazing-incidence WAXS reduction pipeline.

    Parameters
    ----------
    uid : str
        Tiled run UID.
    tiled_uri, catalog : str
        Tiled connection parameters.
    waxs_mask : str, Path, dict, or None
        WAXS mask spec.  Accepts a JSON file path, a Path, an in-memory
        dict (same schema as the bundled JSON; see
        :func:`make_waxs_mask_callable_from_dict`), or ``None`` to use
        the bundled SMI default mask shipped with smi-tiled
        (``smi_tiled.defaults.default_waxs_mask_path``).
    waxs_mask_path : str or Path or None
        Deprecated alias for ``waxs_mask`` (Path only).  When both are
        supplied, ``waxs_mask`` wins.
    n_qxy, n_qz : int
        Output grid dimensions.
    incident_angle_deg : float or None
        Manual incident-angle override.  None = auto-detect from
        sample_name or motor positions.
    theta_offset : float
        Added to ``stage_th + piezo_th`` when auto-detecting the
        incident angle.
    waxs_beam_col_per_arc_deg : float
        Beam-centre drift per degree of waxs_arc.
    beamstop_max_abs_arc_deg : float
        Mask beamstop only for ``|arc| <= this``.
    dezinger_threshold, dezinger_kernel
        Hot-pixel rejection parameters.
    waxs_cal_overrides : dict or None
        Extra overrides for ``WAXSCalibration`` fields.
    image_cache_path : str, Path, None, or "auto"
        Path to a pre-cached HDF5 image file (SMI Browser disk cache).
        If ``"auto"`` (default), automatically checks
        ``$SMI_BROWSER_CACHE_DIR/<uid>.h5``.  If a cache file is found,
        images are read from it instead of tiled; missing fields fall back
        to tiled transparently.  Pass ``None`` to disable cache lookup.
    populate_disk_cache : bool
        If True (default) and no cache file existed, write fetched data
        to the cache after loading from tiled.
    pixel_splitting : int
        Number of sub-pixel divisions per axis for fractional pixel
        splitting during histogram binning.  1 (default) disables splitting.
    progress : callable or None
        Optional callback ``(stage: str, current: int, total: int) -> None``
        invoked to report progress.  Stages: ``"load"``, ``"gi_setup"``,
        ``"gi_integrate"``.  *current* is 1-based; *total* is the number of
        steps in that stage.

    Returns
    -------
    GIReductionResult
    """
    import time as _time
    from smi_tiled.loader import (
        TiledSMISWAXSLoader,
        _auto_cache_path,
        populate_cache,
        resolve_waxs_geometry,
    )

    t0 = _time.perf_counter()

    # Resolve image cache path
    _cache_was_missing = False
    if image_cache_path == "auto":
        image_cache_path = _auto_cache_path(uid)
        if image_cache_path is None:
            _cache_was_missing = True
    elif image_cache_path is not None:
        image_cache_path = Path(image_cache_path) if not isinstance(image_cache_path, Path) else image_cache_path
        if not image_cache_path.exists():
            _cache_was_missing = True
            image_cache_path = None

    # --- Connect & get metadata ---
    from tiled.client import from_uri
    client = from_uri(tiled_uri)
    run = client[catalog + "/" + uid]
    start = run.metadata.get("start", {})
    sample_name = start.get("sample_name", "")
    n_frames = start.get("num_points", 1)
    scan_motor = (start.get("motors") or ["unknown"])[0]

    # --- Incident angle ---
    alpha_i, ai_source = find_incident_angle(
        run, n_frames,
        manual_override=incident_angle_deg,
        theta_offset=theta_offset,
    )

    # --- Scan motor values (for labelling) ---
    scan_motor_values = alpha_i.copy()  # default: use alpha_i
    try:
        primary_ds = run["primary"]
        if scan_motor in primary_ds:
            scan_motor_values = np.asarray(
                primary_ds[scan_motor].read(), dtype=float,
            )
    except Exception:
        pass

    # --- Load WAXS images ---
    loader = TiledSMISWAXSLoader(tiled_uri=tiled_uri, catalog=catalog)
    waxs_raw = loader.loadSingleImage(uid, detector="waxs", image_cache_path=image_cache_path)
    if waxs_raw is None:
        raise RuntimeError(f"No WAXS data in scan {uid}")
    t_load = _time.perf_counter()

    if progress is not None:
        progress("load", 1, 1)

    # Populate disk cache for future runs if data was fetched from tiled
    if populate_disk_cache and _cache_was_missing:
        try:
            populate_cache(uid, run, include_images=True)
        except Exception:
            pass  # cache write is best-effort

    # --- Mask ---
    # waxs_mask (new, accepts dict) wins over waxs_mask_path (legacy)
    if waxs_mask is None:
        waxs_mask = waxs_mask_path
    from smi_tiled.defaults import resolve_mask_path
    if not isinstance(waxs_mask, dict):
        waxs_mask = resolve_mask_path(waxs_mask, detector="waxs")
    waxs_mask_fn = None
    if waxs_mask is not None:
        waxs_mask_fn = make_waxs_mask_callable(
            waxs_mask,                                  # accepts dict or Path
            beamstop_max_abs_arc_deg=beamstop_max_abs_arc_deg,
        )

    # --- WAXS calibration ---
    cal_dict: dict[str, Any] = dict(_DEFAULT_CAL)
    cal_dict["beam_col_per_arc_deg"] = waxs_beam_col_per_arc_deg
    if waxs_cal_overrides:
        cal_dict.update(waxs_cal_overrides)
    waxs_cal = WAXSCalibration(**cal_dict)

    if progress is not None:
        progress("gi_setup", 1, 1)

    # --- Integrate ---
    t_int = _time.perf_counter()
    gi_out = integrate_waxs_gi(
        waxs_raw=waxs_raw,
        mask_fn=waxs_mask_fn,
        alpha_i_deg=alpha_i,
        n_qxy=n_qxy,
        n_qz=n_qz,
        cal=waxs_cal,
        dezinger_threshold=dezinger_threshold,
        dezinger_kernel=dezinger_kernel,
        pixel_splitting=pixel_splitting,
        progress=progress,
    )
    t_done = _time.perf_counter()

    result = GIReductionResult(
        uid=uid,
        sample_name=sample_name,
        scan_motor=scan_motor,
        scan_motor_values=scan_motor_values,
        alpha_i_deg=alpha_i,
        alpha_i_source=ai_source,
        qxy_grid=gi_out["qxy_grid"],
        qz_grid=gi_out["qz_grid"],
        frames=gi_out["frames"],
        summed=gi_out["summed"],
        q_chi_frames=gi_out["q_chi_frames"],
        summed_ds=gi_out["summed_ds"],
        timing={
            "total": t_done - t0,
            "tiled_load": t_load - t0,
            "integrate": t_done - t_int,
        },
    )
    if line_cuts:
        from .derived import apply_line_cuts
        apply_line_cuts(result, list(line_cuts))
    return result


# ===================================================================
# SAXS integration
# ===================================================================

#: Above this frame count, the per-frame detector-space ``ds`` output is
#: skipped by default.  ``ds`` carries a full float64 copy of the image stack
#: plus per-frame geometry/mask arrays — tens of GB on large parallel scans —
#: and is never consumed by the reduction itself (merged products and per-frame
#: I(q) do not depend on it).  Pass ``build_detector_ds=True`` to force it.
_DETECTOR_DS_AUTO_MAX_FRAMES = 50


def _resolve_build_detector_ds(
    build_detector_ds: bool | None, n_frames: int, detector: str
) -> bool:
    """Resolve the tri-state ``build_detector_ds`` flag.

    ``None`` (default) means *auto*: build the detector-space ``ds`` only for
    small scans, skip it (with a warning) for large ones so they don't OOM.
    """
    if build_detector_ds is not None:
        return bool(build_detector_ds)
    if n_frames <= _DETECTOR_DS_AUTO_MAX_FRAMES:
        return True
    warnings.warn(
        f"[integrate_{detector}] skipping the detector-space 'ds' output for "
        f"{n_frames} frames (> {_DETECTOR_DS_AUTO_MAX_FRAMES}) to bound memory; "
        f"pass build_detector_ds=True to force it.",
        stacklevel=3,
    )
    return False


def _dir_is_writable(path: Path) -> bool:
    """Create *path* if needed and confirm we can write a file inside it."""
    try:
        path.mkdir(parents=True, exist_ok=True)
        probe = path / ".smi_write_test"
        probe.touch()
        probe.unlink()
        return True
    except OSError:
        return False


def _resolve_frame_store_dir(
    uid: str,
    frame_qchi_store: "str | Path | None",
    image_cache_path: "str | Path | None",
    n_frames: int,
) -> "Path | None":
    """Resolve a *writable* directory for the streamed per-frame q-chi store.

    Returns ``None`` to mean "keep per-frame q-chi in memory".

    ``"auto"`` streams only for large scans, preferring the directory of the
    provided image cache (known writable in the SMI browser flow), then the
    standard cache dir, then a private temp dir.  Each candidate is probed for
    writability — ``_cache_dir()`` uses ``mkdir(exist_ok=True)`` and so happily
    returns a shared dir owned by another user, which then fails only when we
    try to create a sub-store.  An explicit path is honoured with the same
    fallback.  Only if nothing is writable do we fall back to in-memory.
    """
    import tempfile

    if frame_qchi_store is None:
        return None

    candidates: list[Path] = []
    if frame_qchi_store == "auto":
        if n_frames <= _DETECTOR_DS_AUTO_MAX_FRAMES:
            return None
        if isinstance(image_cache_path, (str, Path)):
            try:
                candidates.append(
                    Path(image_cache_path).expanduser().parent / f"{uid}_qchi"
                )
            except Exception:
                pass
        try:
            from smi_tiled.loader import _cache_dir
            candidates.append(_cache_dir() / f"{uid}_qchi")
        except Exception:
            pass
    else:
        candidates.append(Path(frame_qchi_store))

    for cand in candidates:
        if _dir_is_writable(cand):
            return cand

    try:
        tmp = Path(tempfile.mkdtemp(prefix=f"smi_qchi_{uid}_"))
        warnings.warn(
            f"frame_qchi_store: preferred cache locations were not writable; "
            f"streaming per-frame q-chi to temp dir {tmp} instead.",
            stacklevel=3,
        )
        return tmp
    except OSError:
        warnings.warn(
            "frame_qchi_store: no writable location found; keeping per-frame "
            "q-chi in memory (may use significant RAM on large scans).",
            stacklevel=3,
        )
        return None


def integrate_saxs(
    saxs_raw: xr.DataArray,
    mask: np.ndarray | None,
    n_q: int = 1000,
    n_chi: int = 360,
    solid_angle_correction: bool = False,
    rotate_cw_90: bool = False,
    waxs_arc: np.ndarray | None = None,
    beam_center_col_px: float | None = None,
    dynamic_saxs_mask: bool = False,
    dynamic_saxs_kwargs: dict[str, Any] | None = None,
    dezinger_threshold: float | None = None,
    dezinger_kernel: int = 5,
    cache_geometry: bool = True,
    pixel_splitting: int = 1,
    build_detector_ds: bool | None = None,
    build_frame_qchi: bool = True,
    frame_qchi_store: "str | Path | None" = None,
    progress: ProgressCallback | None = None,
) -> dict[str, Any]:
    """SAXS reduction via direct pixel-space q-map and histogram binning.

    When *frame_qchi_store* is given, the per-frame ``(q, chi)`` maps are
    streamed to that zarr path one frame at a time and ``out['q_chi_frames']``
    is returned as a lazy, dask-backed dataset (see :class:`_ZarrFrameWriter`),
    so peak memory stays at a single frame regardless of frame count.
    """
    _t0_saxs = _time.perf_counter()
    attrs = saxs_raw.attrs
    # Build geometry arrays from attrs
    dist_m = float(attrs["dist"])
    poni1_m = float(attrs["poni1"])
    poni2_m = float(attrs["poni2"])
    pixel1_m = float(attrs["pixel1"])
    pixel2_m = float(attrs["pixel2"])
    wavelength_m = float(attrs["wavelength"]) * 1e-10
    wavelength_nm = wavelength_m * 1e9

    # Lazy image source.  ``saxs_raw.data`` may be a dask array (streamed from
    # tiled, or a chunk-backed HDF5 cache) or a plain ndarray.  We pull frames
    # in blocks via ``_load_block`` so the full stack (~16 GB at 800 frames) is
    # never materialized at once — peak raw memory is one block.
    _img_src = saxs_raw.data
    _is_2d = saxs_raw.ndim == 2
    if _is_2d:
        n_frames = 1
        ny, nx = saxs_raw.shape
    else:
        n_frames = int(saxs_raw.shape[0])
        ny, nx = saxs_raw.shape[-2:]
    shape = (ny, nx)

    def _load_block(s: int, e: int) -> np.ndarray:
        """Frames [s:e) as a float64 ndarray (m, ny, nx); computes if lazy."""
        if _is_2d:
            return np.asarray(_img_src, dtype=float)[np.newaxis, :, :]
        return np.asarray(_img_src[s:e], dtype=float)

    mask_use = None if mask is None else np.asarray(mask, dtype=bool)

    # Check persistent geometry cache
    _saxs_key = _saxs_cache_key(dist_m, poni1_m, poni2_m, pixel1_m, pixel2_m, wavelength_m, shape)
    _cached = _SAXS_GEOMETRY_CACHE.get(_saxs_key) if cache_geometry else None

    if _cached is not None:
        q2d, qh2d, qv2d, chi_deg_2d, sa_base = _cached
        # Recompute sa based on current solid_angle_correction setting
        sa = sa_base if solid_angle_correction else None
    else:
        rr, cc = np.meshgrid(
            np.arange(ny, dtype=float),
            np.arange(nx, dtype=float),
            indexing="ij",
        )
        bc_row = poni1_m / pixel1_m
        bc_col = poni2_m / pixel2_m

        x_m = (cc - bc_col) * pixel2_m
        y_m = -(rr - bc_row) * pixel1_m
        r_m = np.sqrt(x_m**2 + y_m**2 + dist_m**2)
        k = 2.0 * np.pi / wavelength_nm

        qh2d = k * x_m / r_m
        qv2d = k * y_m / r_m
        qz2d = k * (dist_m / r_m - 1.0)
        q2d = np.sqrt(qh2d**2 + qv2d**2 + qz2d**2)
        chi_deg_2d = np.rad2deg(np.arctan2(qh2d, qv2d))

        pixel_area_m2 = pixel1_m * pixel2_m
        with np.errstate(invalid="ignore", divide="ignore"):
            sa_base = pixel_area_m2 * np.maximum(dist_m, 0.0) / (r_m**3)
        sa = sa_base if solid_angle_correction else None

        if cache_geometry:
            _SAXS_GEOMETRY_CACHE[_saxs_key] = (q2d, qh2d, qv2d, chi_deg_2d, sa_base)

    bc_col = poni2_m / pixel2_m
    print(f"  [integrate_saxs] geometry build/cache: "
          f"{_time.perf_counter() - _t0_saxs:.3f}s")

    _t_mask = _time.perf_counter()
    base_valid = np.isfinite(q2d) & np.isfinite(chi_deg_2d)
    if mask_use is not None:
        base_valid &= mask_use

    # Dynamic large-area masks
    if waxs_arc is None:
        dim0 = saxs_raw.dims[0] if saxs_raw.ndim > 2 else None
        if dim0 == "waxs_arc":
            waxs_arc = np.asarray(
                saxs_raw.coords["waxs_arc"].values, dtype=float
            )
    if beam_center_col_px is None:
        beam_center_col_px = bc_col

    mask_options = dict(dynamic_saxs_kwargs or {})
    ws_kw = dict(mask_options.pop("waxs_shadow", {}) or {})
    ap_kw = dict(mask_options.pop("aperture", {}) or {})

    large_area_mask, _, _ = make_saxs_large_area_masks(
        shape,
        q2d,
        waxs_arc,
        beam_center_col_px=float(beam_center_col_px),
        waxs_shadow=ws_kw,
        aperture=ap_kw,
    )

    # Per-frame static validity = base_valid AND the large-area (beamstop
    # shadow / aperture) mask for that frame.  Computed per frame inside the
    # block loop rather than materialized as a full (frame, row, col) stack
    # (~2 GB of bool at 800 frames).  Dezinger (on the raw frame) is folded in
    # there too.
    def _large_slice(idx: int):
        if large_area_mask.shape[0] == 1:
            return large_area_mask[0]
        if large_area_mask.shape[0] >= n_frames:
            return large_area_mask[idx]
        return None
    print(f"  [integrate_saxs] mask+large_area: "
          f"{_time.perf_counter() - _t_mask:.3f}s")

    # Bin edges
    _t_bins = _time.perf_counter()
    q_vals = q2d[base_valid]
    chi_vals = chi_deg_2d[base_valid]
    # Use the full physical range of the unmasked detector for q and chi
    # bin edges.  An earlier version used np.percentile(q_vals, [0.5, 99.5])
    # to clip outliers, but this silently dropped the lowest q values —
    # the most important region for USAXS scans (long-SDD, small q_min).
    # The mask already excludes bad/beamstop pixels, so the literal
    # min/max of surviving q is the right physical range.
    if q_vals.size > 0:
        q_min = float(np.nanmin(q_vals))
        q_max = float(np.nanmax(q_vals))
        q_edges = np.linspace(q_min, q_max, n_q + 1)
    else:
        q_edges = np.linspace(0.0, 10.0, n_q + 1)
    q_grid = 0.5 * (q_edges[:-1] + q_edges[1:])

    if chi_vals.size > 0:
        chi_min = float(np.nanmin(chi_vals))
        chi_max = float(np.nanmax(chi_vals))
        chi_edges = np.linspace(chi_min, chi_max, n_chi + 1)
    else:
        chi_edges = np.linspace(-180.0, 180.0, n_chi + 1)
    chi_grid = 0.5 * (chi_edges[:-1] + chi_edges[1:])
    print(f"  [integrate_saxs] bin edges: "
          f"{_time.perf_counter() - _t_bins:.3f}s")

    accum_I = np.zeros((n_q, n_chi), dtype=float)
    accum_N = np.zeros((n_q, n_chi), dtype=float)
    frame_qchi: list[xr.Dataset] = []
    frame_iq: list[xr.Dataset] = []

    # Precompute the pixel->bin mapping ONCE.  The q/chi maps and bin edges are
    # identical for every frame, so we never need to recompute which bin each
    # pixel falls into — only the per-frame intensities and validity mask vary.
    # ``_SplitBinPlan`` turns per-frame integration into a sparse mat-vec.
    _t_plan = _time.perf_counter()
    plan = _SplitBinPlan(q2d, chi_deg_2d, q_edges, chi_edges,
                         pixel_splitting=pixel_splitting)
    print(f"  [integrate_saxs] bin-plan precompute: "
          f"{_time.perf_counter() - _t_plan:.3f}s")

    qchi_writer = (
        _ZarrFrameWriter(frame_qchi_store, n_frames, q_grid, chi_grid)
        if frame_qchi_store is not None else None
    )

    # Detector-space ds is large and never consumed by the reduction; only
    # accumulate the raw frames/masks for it when actually building it (auto-
    # skipped for big scans).
    make_ds = _resolve_build_detector_ds(build_detector_ds, n_frames, "saxs")
    ds_images = np.empty((n_frames, ny, nx), dtype=float) if make_ds else None
    ds_masks = np.empty((n_frames, ny, nx), dtype=bool) if make_ds else None

    # Solid-angle acceptance mask is frame-independent; precompute once.
    sa_ok = None
    if sa is not None:
        sa_ok = (mask_use if mask_use is not None else np.ones(shape, bool)) & (sa > 0)

    _t_loop = _time.perf_counter()
    _t_io = _t_dez = _t_hist_total = _t_qchi_total = 0.0

    # Double-buffered I/O: prefetch the next block while processing the current one.
    from concurrent.futures import ThreadPoolExecutor
    _block_starts = list(range(0, n_frames, IMAGE_BLOCK_FRAMES))
    _prefetch_exec = ThreadPoolExecutor(max_workers=1)
    _prefetch_future = _prefetch_exec.submit(_load_block, _block_starts[0],
                                             min(_block_starts[0] + IMAGE_BLOCK_FRAMES, n_frames))
    for _bi, bstart in enumerate(_block_starts):
        bend = min(bstart + IMAGE_BLOCK_FRAMES, n_frames)
        _tio = _time.perf_counter()
        block = _prefetch_future.result()
        # Submit prefetch for next block while we process this one.
        if _bi + 1 < len(_block_starts):
            _next_start = _block_starts[_bi + 1]
            _next_end = min(_next_start + IMAGE_BLOCK_FRAMES, n_frames)
            _prefetch_future = _prefetch_exec.submit(_load_block, _next_start, _next_end)
        _t_io += _time.perf_counter() - _tio
        for j in range(bend - bstart):
            idx = bstart + j
            img_raw = block[j]
            # Dezinger flags hot pixels on the RAW frame (before SA scaling).
            _td = _time.perf_counter()
            dz = (dezinger(img_raw, kernel_size=dezinger_kernel,
                           threshold=dezinger_threshold)
                  if dezinger_threshold is not None else None)
            _t_dez += _time.perf_counter() - _td
            # Solid-angle correction.
            if sa_ok is not None:
                with np.errstate(divide="ignore", invalid="ignore"):
                    img = np.where(sa_ok, img_raw / sa, np.nan)
            else:
                img = img_raw
            # Static per-frame validity = base & large-area & dezinger.
            pfv = base_valid
            ls = _large_slice(idx)
            if ls is not None:
                pfv = pfv & ls
            if dz is not None:
                pfv = pfv & dz
            valid = pfv & np.isfinite(img)

            _th = _time.perf_counter()
            i_hist, n_hist = plan.integrate_frame(img, valid)
            _t_hist_total += _time.perf_counter() - _th
            accum_I += i_hist
            accum_N += n_hist
            _tq = _time.perf_counter()
            frame_out = _qchi_and_iq(i_hist, n_hist, q_grid, chi_grid)
            if qchi_writer is not None:
                qchi_writer.write(idx, frame_out["q_chi"]["intensity"].values,
                                  frame_out["q_chi"]["counts"].values)
            elif build_frame_qchi:
                frame_qchi.append(frame_out["q_chi"])
            frame_iq.append(frame_out["iq"])
            _t_qchi_total += _time.perf_counter() - _tq

            if make_ds:
                ds_images[idx] = img_raw
                ds_masks[idx] = pfv

            if progress is not None:
                progress("saxs_integrate", idx + 1, n_frames)

    _prefetch_exec.shutdown(wait=False)
    print(f"  [integrate_saxs] per-frame loop ({n_frames} frames): "
          f"{_time.perf_counter() - _t_loop:.3f}s "
          f"(io={_t_io:.3f}s, dez={_t_dez:.3f}s, hist={_t_hist_total:.3f}s, "
          f"qchi={_t_qchi_total:.3f}s)")
    _t_build = _time.perf_counter()
    out = _qchi_and_iq(accum_I, accum_N, q_grid, chi_grid)
    if qchi_writer is not None:
        out["q_chi_frames"] = qchi_writer.dataset()
    elif build_frame_qchi:
        out["q_chi_frames"] = _stack_qchi_frames(frame_qchi)
    else:
        out["q_chi_frames"] = None
    out["iq_frames"] = _stack_iq_frames(frame_iq)

    # Detector-space xr.Dataset from the accumulated raw frames / masks.  The
    # SAXS q-maps are identical for every frame (fixed transmission geometry),
    # so store them ONCE as 2-D (row, col) rather than repeating per frame.
    if make_ds:
        ds_q2d, ds_qh2d, ds_qv2d = q2d, qh2d, qv2d
        if rotate_cw_90:
            ds_images = np.rot90(ds_images, k=-1, axes=(-2, -1))
            ds_q2d = np.rot90(ds_q2d, k=-1)
            ds_qh2d = np.rot90(ds_qh2d, k=-1)
            ds_qv2d = np.rot90(ds_qv2d, k=-1)
            ds_masks = np.rot90(ds_masks, k=-1, axes=(-2, -1))
        out["ds"] = xr.Dataset(
            {
                "intensity": (("frame", "row", "col"), ds_images),
                "q_abs": (("row", "col"), ds_q2d),
                "q_horizontal": (("row", "col"), ds_qh2d),
                "q_vertical": (("row", "col"), ds_qv2d),
                "mask": (("frame", "row", "col"), ds_masks),
            },
            coords={"frame": np.arange(n_frames, dtype=int)},
        )
    else:
        out["ds"] = None
    print(f"  [integrate_saxs] output build: "
          f"{_time.perf_counter() - _t_build:.3f}s")
    print(f"  [integrate_saxs] TOTAL: "
          f"{_time.perf_counter() - _t0_saxs:.3f}s")
    return out


# ===================================================================
# WAXS integration
# ===================================================================

def integrate_waxs(
    waxs_raw: xr.DataArray,
    mask_fn,
    n_q: int = 1000,
    n_chi: int = 360,
    cal: WAXSCalibration | None = None,
    solid_angle_correction: bool = False,
    flip_horizontal: bool = False,
    qx_shift_nm: float = 0.0,
    qy_shift_nm: float = 0.0,
    dezinger_threshold: float | None = None,
    dezinger_kernel: int = 5,
    cache_geometry: bool = True,
    pixel_splitting: int = 1,
    build_detector_ds: bool | None = None,
    build_frame_qchi: bool = True,
    frame_qchi_store: "str | Path | None" = None,
    progress: ProgressCallback | None = None,
) -> dict[str, Any]:
    """WAXS reduction via MultiPanelArcDetector per arc-angle frame.

    When *frame_qchi_store* is given, per-frame ``(q, chi)`` maps are streamed
    to that zarr path and ``out['q_chi_frames']`` is returned lazily (dask-
    backed), keeping peak memory at a single frame.  See :class:`_ZarrFrameWriter`.
    """
    _t0_waxs = _time.perf_counter()
    attrs = waxs_raw.attrs
    if cal is None:
        cal = WAXSCalibration(**_DEFAULT_CAL)

    # Lazy image source (dask array or ndarray); frames are pulled in blocks so
    # the full stack is never materialized at once.
    _img_src = waxs_raw.data
    _is_2d = waxs_raw.ndim == 2
    n_frames = 1 if _is_2d else int(waxs_raw.shape[0])

    def _load_block(s: int, e: int) -> np.ndarray:
        if _is_2d:
            return np.asarray(_img_src, dtype=float)[np.newaxis, :, :]
        return np.asarray(_img_src[s:e], dtype=float)

    arc_angles = np.asarray(
        waxs_raw.coords[waxs_raw.dims[0]].values, dtype=float
    )
    bsx_per_frame = np.asarray(
        attrs.get("smi_waxs_bsx_per_frame", [0.0] * n_frames),
        dtype=float,
    )

    # Per-frame energy (eV) — used for wavelength in q-map computation
    energy_per_frame_ev = np.asarray(
        attrs.get("smi_energy_per_frame_ev", [cal.energy_kev * 1000.0] * n_frames),
        dtype=float,
    )

    # Rotated shape from frame 0 only (no full materialization).
    img_0_rot, _ = rotate_image_and_mask(_load_block(0, 1)[0], k=cal.rotation_k)
    rot_shape = img_0_rot.shape

    def build_detector_for_angle(theta_deg: float, wavelength_nm: float | None = None):
        bc = cal.beam_center_at_angle(float(theta_deg))
        wl = wavelength_nm if wavelength_nm is not None else cal.wavelength_nm
        return MultiPanelArcDetector(
            image_shape=rot_shape,
            panel_specs=cal.make_panel_specs(),
            wavelength_nm=wl,
            pixel_size_mm=cal.pixel_size_mm,
            sample_distance_mm=cal.sample_distance_mm,
            beam_center_px=bc,
            theta_zero_deg=cal.theta_zero_deg,
            sample_offset_x_mm=cal.sample_offset_x_mm,
            sample_offset_z_mm=cal.sample_offset_z_mm,
        )

    # Pre-compute global q/chi range across all angles
    _all_q_min, _all_q_max = [], []

    # Use module-level persistent cache if requested
    _cache_key = None
    if cache_geometry:
        _cache_key = _waxs_cache_key(cal, rot_shape, flip_horizontal, qx_shift_nm, qy_shift_nm)
        _geo_cache = _WAXS_GEOMETRY_CACHE.setdefault(_cache_key, {})
    else:
        _geo_cache = {}

    # Determine per-frame wavelength (nm) from energy
    _HC_EV_NM = 1239.84198  # eV·nm
    wavelength_per_frame_nm = _HC_EV_NM / energy_per_frame_ev

    for fi, theta_val in enumerate(arc_angles):
        theta_f = float(theta_val)
        wl_nm = float(wavelength_per_frame_nm[fi])
        # Cache key includes both arc angle and wavelength
        key = (round(theta_f, 6), round(wl_nm, 8))
        if key in _geo_cache:
            # Still need q-range info even from cached entries
            qabs_px = _geo_cache[key][0]
            chi_px = _geo_cache[key][3]
            finite = np.isfinite(qabs_px) & np.isfinite(chi_px)
            if finite.any():
                _all_q_min.append(float(np.nanmin(qabs_px[finite])))
                _all_q_max.append(float(np.nanmax(qabs_px[finite])))
            continue
        det = build_detector_for_angle(theta_f, wavelength_nm=wl_nm)
        qds = det.qmap(theta_f)
        qx_px = cal.q_horizontal_sign * np.asarray(qds["qx"].values, dtype=float)
        qy_px = cal.q_vertical_sign * np.asarray(qds["qy"].values, dtype=float)
        qabs_px = np.asarray(qds["qabs"].values, dtype=float)
        sa_px = np.asarray(qds["solid_angle"].values, dtype=float)
        if flip_horizontal:
            qx_px = np.fliplr(qx_px)
            qy_px = np.fliplr(qy_px)
            qabs_px = np.fliplr(qabs_px)
            sa_px = np.fliplr(sa_px)
        qx_px += qx_shift_nm
        qy_px += qy_shift_nm
        chi_px = np.rad2deg(np.arctan2(qx_px, qy_px))
        _geo_cache[key] = (qabs_px, qx_px, qy_px, chi_px, sa_px)

        finite = np.isfinite(qabs_px) & np.isfinite(chi_px)
        if finite.any():
            _all_q_min.append(float(np.nanmin(qabs_px[finite])))
            _all_q_max.append(float(np.nanmax(qabs_px[finite])))

    if _all_q_min:
        q_edges = np.linspace(min(_all_q_min), max(_all_q_max), n_q + 1)
    else:
        q_edges = np.linspace(0, 20, n_q + 1)
    q_grid = 0.5 * (q_edges[:-1] + q_edges[1:])
    chi_edges = np.linspace(-180.0, 180.0, n_chi + 1)
    chi_grid = 0.5 * (chi_edges[:-1] + chi_edges[1:])
    print(f"  [integrate_waxs] geometry precompute ({len(arc_angles)} angles): "
          f"{_time.perf_counter() - _t0_waxs:.3f}s")

    accum_I = np.zeros((n_q, n_chi), dtype=float)
    accum_N = np.zeros((n_q, n_chi), dtype=float)
    frame_qchi: list[xr.Dataset] = []
    frame_iq: list[xr.Dataset] = []

    ds_int_frames, ds_qabs_frames = [], []
    ds_qh_frames, ds_qv_frames, ds_mask_frames = [], [], []

    # Per-geometry bin-plan cache.  WAXS geometry varies with arc angle, but
    # frames that share an (angle, wavelength) key share a q/chi map and thus a
    # pixel->bin mapping.  For a fast raster at fixed arc there is exactly one
    # key, so the plan is built once and reused for every frame — the same big
    # win as SAXS.  Keyed identically to ``_geo_cache``.
    _bin_plans: dict[tuple, _SplitBinPlan] = {}

    # Detector-space ``ds`` is large (per-frame intensity + geometry + mask) and
    # never consumed by the reduction; skip it (and its per-frame list
    # accumulation) for big scans unless explicitly requested.
    make_ds = _resolve_build_detector_ds(build_detector_ds, len(arc_angles), "waxs")

    qchi_writer = (
        _ZarrFrameWriter(frame_qchi_store, len(arc_angles), q_grid, chi_grid)
        if frame_qchi_store is not None else None
    )

    _t_loop = _time.perf_counter()
    _t_waxs_dez = 0.0
    _t_waxs_mask = 0.0
    _t_waxs_hist = 0.0
    _t_waxs_qchi = 0.0
    _t_waxs_plan = 0.0
    _t_waxs_io = 0.0

    # Double-buffered I/O: prefetch the next block while processing the current one.
    from concurrent.futures import ThreadPoolExecutor as _TPE_waxs
    _waxs_block_starts = list(range(0, n_frames, IMAGE_BLOCK_FRAMES))
    _waxs_prefetch_exec = _TPE_waxs(max_workers=1)
    _waxs_prefetch_future = _waxs_prefetch_exec.submit(
        _load_block, _waxs_block_starts[0],
        min(_waxs_block_starts[0] + IMAGE_BLOCK_FRAMES, n_frames))
    for _wbi, _bstart in enumerate(_waxs_block_starts):
        _bend = min(_bstart + IMAGE_BLOCK_FRAMES, n_frames)
        _tio = _time.perf_counter()
        _block = _waxs_prefetch_future.result()
        if _wbi + 1 < len(_waxs_block_starts):
            _next_s = _waxs_block_starts[_wbi + 1]
            _next_e = min(_next_s + IMAGE_BLOCK_FRAMES, n_frames)
            _waxs_prefetch_future = _waxs_prefetch_exec.submit(_load_block, _next_s, _next_e)
        _t_waxs_io += _time.perf_counter() - _tio
        for _j in range(_bend - _bstart):
            fi = _bstart + _j
            theta_f = float(arc_angles[fi])
            wl_nm = float(wavelength_per_frame_nm[fi])
            img_raw = _block[_j]
            bsx = float(bsx_per_frame[fi]) if fi < len(bsx_per_frame) else 0.0
            img_rot, _ = rotate_image_and_mask(img_raw, k=cal.rotation_k)

            _tm = _time.perf_counter()
            mask_rot = None
            if mask_fn is not None:
                try:
                    mask_rot = mask_fn(img_raw.shape, theta_f, bsx)
                except Exception as exc:
                    warnings.warn(
                        f"mask_fn failed for frame {fi}: {exc}", stacklevel=2
                    )
            _t_waxs_mask += _time.perf_counter() - _tm

            # Dezinger: flag hot pixels on the rotated image
            _td = _time.perf_counter()
            if dezinger_threshold is not None:
                dz_mask = dezinger(img_rot, kernel_size=dezinger_kernel,
                                   threshold=dezinger_threshold)
                if mask_rot is not None:
                    mask_rot = mask_rot & dz_mask
                else:
                    mask_rot = dz_mask
            _t_waxs_dez += _time.perf_counter() - _td

            if flip_horizontal:
                img_rot = np.fliplr(img_rot)
                if mask_rot is not None:
                    mask_rot = np.fliplr(mask_rot)

            geo_key = (round(theta_f, 6), round(wl_nm, 8))
            qabs_px, qx_px, qy_px, chi_px, sa_px = _geo_cache[geo_key]

            if solid_angle_correction:
                with np.errstate(divide="ignore", invalid="ignore"):
                    valid_sa = (
                        (mask_rot if mask_rot is not None else np.ones(img_rot.shape, bool))
                        & np.isfinite(sa_px)
                        & (sa_px > 0)
                    )
                    img_rot = np.where(valid_sa, img_rot / sa_px, np.nan)

            valid = np.isfinite(qabs_px) & np.isfinite(chi_px) & np.isfinite(img_rot)
            if mask_rot is not None:
                valid &= mask_rot

            mask_use = valid if mask_rot is None else mask_rot
            if make_ds:
                ds_int_frames.append(np.where(mask_use, img_rot, np.nan))
                ds_qabs_frames.append(qabs_px)
                ds_qh_frames.append(qx_px)
                ds_qv_frames.append(qy_px)
                ds_mask_frames.append(mask_use.astype(bool))

            _tp = _time.perf_counter()
            plan = _bin_plans.get(geo_key)
            if plan is None:
                plan = _SplitBinPlan(qabs_px, chi_px, q_edges, chi_edges,
                                     pixel_splitting=pixel_splitting)
                _bin_plans[geo_key] = plan
            _t_waxs_plan += _time.perf_counter() - _tp

            _th = _time.perf_counter()
            I_hist, N_hist = plan.integrate_frame(img_rot, valid)
            _t_waxs_hist += _time.perf_counter() - _th
            accum_I += I_hist
            accum_N += N_hist
            _tq = _time.perf_counter()
            frame_out = _qchi_and_iq(I_hist, N_hist, q_grid, chi_grid)
            if qchi_writer is not None:
                qchi_writer.write(fi, frame_out["q_chi"]["intensity"].values,
                                  frame_out["q_chi"]["counts"].values)
            elif build_frame_qchi:
                frame_qchi.append(frame_out["q_chi"])
            frame_iq.append(frame_out["iq"])
            _t_waxs_qchi += _time.perf_counter() - _tq

            if progress is not None:
                progress("waxs_integrate", fi + 1, n_frames)

    _waxs_prefetch_exec.shutdown(wait=False)
    print(f"  [integrate_waxs] per-frame loop ({len(arc_angles)} frames): "
          f"{_time.perf_counter() - _t_loop:.3f}s "
          f"(io={_t_waxs_io:.3f}s, mask={_t_waxs_mask:.3f}s, dez={_t_waxs_dez:.3f}s, "
          f"plan={_t_waxs_plan:.3f}s, hist={_t_waxs_hist:.3f}s, "
          f"qchi={_t_waxs_qchi:.3f}s)")
    _t_build = _time.perf_counter()
    out = _qchi_and_iq(accum_I, accum_N, q_grid, chi_grid)
    if qchi_writer is not None:
        out["q_chi_frames"] = qchi_writer.dataset()
    elif build_frame_qchi:
        out["q_chi_frames"] = _stack_qchi_frames(frame_qchi)
    else:
        out["q_chi_frames"] = None
    out["iq_frames"] = _stack_iq_frames(frame_iq)
    if make_ds:
        out["ds"] = xr.Dataset(
            {
                "intensity": (
                    ("frame", "row", "col"),
                    np.asarray(ds_int_frames, dtype=float),
                ),
                "q_abs": (
                    ("frame", "row", "col"),
                    np.asarray(ds_qabs_frames, dtype=float),
                ),
                "q_horizontal": (
                    ("frame", "row", "col"),
                    np.asarray(ds_qh_frames, dtype=float),
                ),
                "q_vertical": (
                    ("frame", "row", "col"),
                    np.asarray(ds_qv_frames, dtype=float),
                ),
                "mask": (
                    ("frame", "row", "col"),
                    np.asarray(ds_mask_frames, dtype=bool),
                ),
                "waxs_arc": (("frame",), np.asarray(arc_angles, dtype=float)),
            },
            coords={"frame": np.arange(len(arc_angles), dtype=int)},
        )
    else:
        out["ds"] = None
    print(f"  [integrate_waxs] output build: "
          f"{_time.perf_counter() - _t_build:.3f}s")
    print(f"  [integrate_waxs] TOTAL: "
          f"{_time.perf_counter() - _t0_waxs:.3f}s")
    return out


# ===================================================================
# Combined reduction entry point
# ===================================================================

def reduce_smi_combined(
    uid: str,
    tiled_uri: str = "https://tiled.nsls2.bnl.gov",
    catalog: str = "smi/migration",
    n_q: int = 2000,
    n_chi: int = 360,
    solid_angle_correction: bool = True,
    saxs_mask: "str | Path | dict | None" = None,
    waxs_mask: "str | Path | dict | None" = None,
    saxs_mask_path: "str | Path | None" = None,
    waxs_mask_path: "str | Path | None" = None,
    saxs_kwargs: dict[str, Any] | None = None,
    waxs_kwargs: dict[str, Any] | None = None,
    backend_options: dict[str, Any] | None = None,
    geometry: str = "transmission",
    incident_angle_deg: float = 0.0,
    saxs_beam_delta_px: Tuple[float, float] | None = None,
    waxs_beam_delta_px: Tuple[float, float] | None = None,
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
    frame_qchi_store: "str | Path | None" = "auto",
    image_cache_path: str | Path | None = "auto",
    populate_disk_cache: bool = True,
    progress: ProgressCallback | None = None,
    # ---- Optional derived-analysis stages (smi_tiled.derived) ----
    virtual_axes: "Any | None" = None,
    line_cuts: "Sequence[Any] | None" = None,
    peak_fits: "Sequence[Any] | None" = None,
) -> CombinedReductionResult:
    """
    Full SAXS + WAXS reduction pipeline.

    Parameters
    ----------
    uid : str
        Tiled run UID.
    tiled_uri, catalog : str
        Tiled connection parameters.
    n_q, n_chi : int
        Output grid dimensions.
    solid_angle_correction : bool
        Apply solid-angle correction to intensities.
    saxs_mask, waxs_mask : str, Path, dict, or None
        Mask specification.  Accepts:

        * ``None`` (default) — bundled SMI default mask
          (see :mod:`smi_tiled.defaults`).
        * ``str`` or ``Path`` — JSON file path; same schema as the
          bundled masks.
        * ``dict`` — already-parsed polygon dict (the JSON contents).
          Useful when composing or editing masks in a notebook without
          writing a temp file.  See
          :func:`make_saxs_mask_from_dict` /
          :func:`make_waxs_mask_callable_from_dict` for the schema.
    saxs_mask_path, waxs_mask_path : str or Path, optional
        Deprecated alias for ``saxs_mask`` / ``waxs_mask`` (accepts
        Path or str only).  Provided for backward compatibility with
        callers written before in-memory dict input was supported.
        When both are supplied, ``saxs_mask`` / ``waxs_mask`` wins.
    saxs_kwargs, waxs_kwargs : dict, optional
        Extra options passed to SAXS / WAXS integrators.
    backend_options : dict, optional
        Options: ``saxs_rotate_cw_90``, ``waxs_flip_horizontal``,
        ``waxs_qx_shift_nm``, ``waxs_qy_shift_nm``.
    geometry : str
        ``'transmission'`` or ``'grazing_incidence'``.
    incident_angle_deg : float
        Incident angle for GI geometry.
    saxs_beam_delta_px : (delta_row, delta_col) or None
        Additive correction to the SAXS beam center read from metadata.
        None uses the built-in defaults from SMISWAXSLoader.
    waxs_beam_delta_px : (delta_row, delta_col) or None
        Additive correction to the WAXS beam center read from metadata.
        None uses the built-in defaults from SMISWAXSLoader.
    saxs_distance_delta_mm : float or None
        Additive correction to the SAXS sample-detector distance read from
        the motor (mm). None uses the built-in default from SMISWAXSLoader.
    saxs_q_cutoff : float or None
        Explicit SAXS q cutoff in nm⁻¹. Overrides silver-behenate-based
        auto-calculation. None (default) uses silver behenate rings.
    saxs_agbh_ring_order : int
        Silver behenate ring order used for auto q cutoff (default 5).
        Ignored if ``saxs_q_cutoff`` is set.
    saxs_q_margin_fraction : float
        Fractional margin above the selected AgBh ring for q cutoff
        (default 0.08, i.e. 8%). Ignored if ``saxs_q_cutoff`` is set.
    dezinger_threshold : float or None
        Sigma threshold for median-filter hot-pixel rejection. Applied to
        both SAXS and WAXS per frame. None (default) disables dezingering.
    dezinger_kernel : int
        Kernel size for the dezinger median filter (default 5).
    cache_geometry : bool
        If True (default), cache precomputed q-maps in a module-level dict
        so that subsequent calls with the same geometry parameters skip the
        expensive pixel-position trigonometry.  Safe across scans that share
        calibration.  Call :func:`clear_geometry_cache` to free memory or
        after programmatically changing calibration parameters.
    pixel_splitting : int
        Number of sub-pixel divisions per axis for fractional pixel
        splitting during histogram binning.  1 (default) disables splitting
        (each pixel contributes to a single bin).  Values > 1 subdivide each
        pixel into an NxN grid and distribute intensity fractionally across
        bins using gradient-based interpolation of the q/chi maps.  Typical
        values are 2–4.
    build_detector_ds : bool or None
        Whether to build the per-frame *detector-space* dataset
        (``result.saxs['ds']`` / ``result.waxs['ds']`` — full image stack plus
        q-maps and masks).  This is large and not used by the reduction itself.
        ``None`` (default) is auto: built for small scans, skipped (with a
        warning) above ~50 frames to bound memory.  ``True``/``False`` force it.
    build_frame_qchi : bool
        Whether to produce the per-frame ``(q, chi)`` stacks at all (default
        True).  Set False to skip them entirely (the merged products and
        per-frame I(q) are unaffected).
    frame_qchi_store : str, Path, None, or "auto"
        Where to keep the per-frame ``(q, chi)`` stacks.  ``"auto"`` (default)
        streams them to a per-uid zarr store under the disk cache dir for large
        scans (> ~50 frames) and returns ``result.*['q_chi_frames']`` as a
        **lazy, dask-backed** dataset — so the multi-GB ``(frame, q, chi)`` array
        never lives in RAM — while small scans stay in memory.  Pass an explicit
        directory to always stream there, or ``None`` to always keep in memory.
    image_cache_path : str, Path, None, or "auto"
        Path to a pre-cached HDF5 image file (SMI Browser disk cache).
        If ``"auto"`` (default), automatically checks
        ``$SMI_BROWSER_CACHE_DIR/<uid>.h5`` (or ``$TMPDIR/smi_browser_cache/``).
        If a cache file is found, images are read from it instead of tiled;
        any missing detector fields fall back to tiled transparently.
        Pass ``None`` to disable cache lookup entirely.
    populate_disk_cache : bool
        If True (default) and the cache file does not already exist, write
        the fetched images, primary scalars, and baseline to a new HDF5
        cache file after loading from tiled.  This speeds up subsequent
        reductions of the same scan (e.g. with different parameters).
    progress : callable or None
        Optional callback ``(stage: str, current: int, total: int) -> None``
        invoked to report progress.  Stages: ``"load"``, ``"saxs_setup"``,
        ``"saxs_integrate"``, ``"waxs_setup"``, ``"waxs_integrate"``,
        ``"merge"``.  *current* is 1-based; *total* is the number of steps
        in that stage.

    Returns
    -------
    CombinedReductionResult
    """
    import time as _time
    from smi_tiled.loader import (
        TiledSMISWAXSLoader,
        _auto_cache_path,
        clear_baseline_cache,
        infer_detectors_and_steps,
        populate_cache,
        resolve_saxs_geometry,
        resolve_waxs_geometry,
    )

    saxs_kw = dict(saxs_kwargs or {})
    waxs_kw = dict(waxs_kwargs or {})
    opts = dict(backend_options or {})
    t0 = _time.perf_counter()

    # -- Build a normalizing progress wrapper --
    # Translate per-stage callbacks into an overall (stage, current, total)
    # where total = n_frames_saxs + n_frames_waxs + overhead_steps.
    # The wrapper is wired to integrate_saxs/waxs so the user's callback
    # receives a monotonically increasing `current` across the whole pipeline.
    _progress_user = progress
    _progress_offset = 0
    _progress_total = 0  # set once n_frames is known

    def _progress_wrap(stage: str, current: int, total: int) -> None:
        """Relay per-stage progress with an overall offset."""
        nonlocal _progress_offset
        if _progress_user is None:
            return
        _progress_user(stage, _progress_offset + current, _progress_total)

    def _progress_advance_stage(stage: str, steps: int) -> None:
        """Advance the offset after a stage completes, emit one callback."""
        nonlocal _progress_offset
        _progress_offset += steps
        if _progress_user is not None:
            _progress_user(stage, _progress_offset, _progress_total)

    # Resolve mask inputs: ``saxs_mask`` / ``waxs_mask`` (new, supports
    # dict) wins over ``saxs_mask_path`` / ``waxs_mask_path`` (legacy,
    # Path-only).  The mask builders downstream accept either form via
    # ``_resolve_mask_spec``, so we just pick the right one here.
    if saxs_mask is None:
        saxs_mask = saxs_mask_path
    if waxs_mask is None:
        waxs_mask = waxs_mask_path

    # Resolve image cache path
    _cache_was_missing = False
    if image_cache_path == "auto":
        image_cache_path = _auto_cache_path(uid)
        if image_cache_path is None:
            _cache_was_missing = True
    elif image_cache_path is not None:
        image_cache_path = Path(image_cache_path) if not isinstance(image_cache_path, Path) else image_cache_path
        if not image_cache_path.exists():
            _cache_was_missing = True
            image_cache_path = None

    # Load raw data — reuse a single loader (and its tiled session) for
    # everything so we don't call from_uri / authenticate twice.
    loader = TiledSMISWAXSLoader(tiled_uri=tiled_uri, catalog=catalog)
    run = loader._get_run(uid)

    # Avoid run["primary"].read() — that pulls every variable in the primary
    # stream including the multi-frame detector arrays, which can trigger
    # an HTTP 500 from the tiled backend.  infer_detectors_and_steps now
    # introspects the tiled containers directly.
    scan_info = infer_detectors_and_steps(run, None, cache_path=image_cache_path)

    saxs_raw = loader.loadSingleImage(uid, detector="saxs", image_cache_path=image_cache_path)
    waxs_raw = loader.loadSingleImage(uid, detector="waxs", image_cache_path=image_cache_path)
    has_saxs = saxs_raw is not None
    has_waxs = waxs_raw is not None
    t_load = _time.perf_counter()

    # Extract waxs_arc for SAXS dynamic masking (WAXS shadow).
    # integrate_saxs can auto-discover waxs_arc from its own DataArray's
    # coordinate axis, but for single-frame (count) scans the SAXS array
    # is squeezed to 2-D and has no waxs_arc dimension.  The WAXS
    # DataArray retains it, so prefer that; fall back to baseline.
    _saxs_waxs_arc = None
    if has_waxs and waxs_raw.ndim > 2 and waxs_raw.dims[0] == "waxs_arc":
        _saxs_waxs_arc = np.asarray(
            waxs_raw.coords["waxs_arc"].values, dtype=float
        )
    elif has_saxs and saxs_raw.ndim > 2 and saxs_raw.dims[0] == "waxs_arc":
        _saxs_waxs_arc = np.asarray(
            saxs_raw.coords["waxs_arc"].values, dtype=float
        )

    # Compute overall progress total now that we know frame counts.
    # Total = 2 (setup stages) + n_saxs_frames + n_waxs_frames + 1 (merge)
    _n_saxs_frames = (saxs_raw.shape[0] if has_saxs and saxs_raw.ndim > 2
                      else (1 if has_saxs else 0))
    _n_waxs_frames = (waxs_raw.shape[0] if has_waxs and waxs_raw.ndim > 2
                      else (1 if has_waxs else 0))
    _progress_total = 2 + _n_saxs_frames + _n_waxs_frames + 1
    _progress_offset = 1  # "load" already done
    if _progress_user is not None:
        _progress_user("load", 1, _progress_total)

    # Resolve per-frame q-chi streaming store.  ``"auto"`` (default) streams the
    # per-frame (q, chi) stacks to a per-uid zarr under the disk cache dir for
    # large scans — keeping the giant (frame, q, chi) arrays off the heap and
    # returning them lazily (dask-backed) — while small scans stay in memory.
    # ``None`` forces in-memory; an explicit path always streams there.
    def _frame_counts() -> int:
        n = 0
        if has_saxs:
            n = max(n, saxs_raw.shape[0] if saxs_raw.ndim > 2 else 1)
        if has_waxs:
            n = max(n, waxs_raw.shape[0] if waxs_raw.ndim > 2 else 1)
        return n

    _store_dir = _resolve_frame_store_dir(
        uid, frame_qchi_store, image_cache_path, _frame_counts()
    )
    saxs_qchi_store = str(_store_dir / "saxs_qchi.zarr") if _store_dir is not None else None
    waxs_qchi_store = str(_store_dir / "waxs_qchi.zarr") if _store_dir is not None else None

    # Populate disk cache for future runs if data was fetched from tiled
    if populate_disk_cache and _cache_was_missing:
        _t_debug = _time.perf_counter()
        try:
            populate_cache(uid, run, include_images=True)
        except Exception:
            pass  # cache write is best-effort
        print(f"[mask_setup] populate_cache (write): "
              f"{_time.perf_counter() - _t_debug:.3f}s")

    # -- SAXS branch --
    saxs_result: dict[str, Any] | None = None
    saxs_geo = None
    t_saxs_start = t_saxs_end = _time.perf_counter()
    if has_saxs:
        _t_debug = _time.perf_counter()
        _saxs_geo_kw: dict[str, Any] = {}
        if saxs_beam_delta_px is not None:
            _saxs_geo_kw["beam_delta_row_px"] = saxs_beam_delta_px[0]
            _saxs_geo_kw["beam_delta_col_px"] = saxs_beam_delta_px[1]
        if saxs_distance_delta_mm is not None:
            _saxs_geo_kw["distance_delta_mm"] = saxs_distance_delta_mm
        saxs_geo = resolve_saxs_geometry(run, **_saxs_geo_kw)
        print(f"[mask_setup] resolve_saxs_geometry: "
              f"{_time.perf_counter() - _t_debug:.3f}s")

        # Update saxs_raw attrs with corrected geometry
        _pixel1 = float(saxs_raw.attrs["pixel1"])
        _pixel2 = float(saxs_raw.attrs["pixel2"])
        new_attrs = dict(saxs_raw.attrs)
        new_attrs["poni1"] = saxs_geo.beam_center_row_px * _pixel1
        new_attrs["poni2"] = saxs_geo.beam_center_col_px * _pixel2
        new_attrs["dist"] = saxs_geo.dist_m
        saxs_raw.attrs.update(new_attrs)

        # SAXS mask
        _t_debug = _time.perf_counter()
        if saxs_mask is None:
            saxs_mask = saxs_kw.pop("mask_path", None)
        # When mask is a dict, skip the path-only resolver (which would
        # only check for file existence anyway); otherwise let it apply
        # bundled-default fallback for None / bare-filename inputs.
        from smi_tiled.defaults import resolve_mask_path
        if not isinstance(saxs_mask, dict):
            saxs_mask = resolve_mask_path(saxs_mask, detector="saxs")
        saxs_mask_array = None
        if saxs_mask is not None:
            saxs_mask_array = make_saxs_mask_from_spec(
                image_shape=saxs_raw.shape[-2:],
                mask_path=saxs_mask,                       # accepts dict or Path
                active_beamstop=saxs_geo.active_beamstop,
                beamstop_pos_mm=saxs_geo.beamstop_pos_mm,
                beam_center_px=(
                    saxs_geo.beam_center_row_px,
                    saxs_geo.beam_center_col_px,
                ),
            )
        print(f"[mask_setup] SAXS mask creation: "
              f"{_time.perf_counter() - _t_debug:.3f}s")

        _progress_advance_stage("saxs_setup", 1)

        # Integrate SAXS
        t_saxs_start = _time.perf_counter()
        _dyn_kw = dict(saxs_kw.get("dynamic_saxs_kwargs") or {})
        _ap = dict(_dyn_kw.pop("aperture", {}) or {})
        _ap.setdefault("agbh_ring_order", saxs_agbh_ring_order)
        _ap.setdefault("q_margin_fraction", saxs_q_margin_fraction)
        if saxs_q_cutoff is not None:
            _ap["q_cutoff"] = saxs_q_cutoff
        _dyn_kw["aperture"] = _ap

        saxs_result = integrate_saxs(
            saxs_raw=saxs_raw,
            mask=saxs_mask_array,
            n_q=n_q,
            n_chi=n_chi,
            solid_angle_correction=solid_angle_correction,
            rotate_cw_90=bool(opts.get("saxs_rotate_cw_90", False)),
            waxs_arc=_saxs_waxs_arc,
            beam_center_col_px=saxs_geo.beam_center_col_px,
            dynamic_saxs_mask=bool(saxs_kw.get("dynamic_saxs_mask", False)),
            dynamic_saxs_kwargs=_dyn_kw,
            dezinger_threshold=dezinger_threshold,
            dezinger_kernel=dezinger_kernel,
            cache_geometry=cache_geometry,
            pixel_splitting=pixel_splitting,
            build_detector_ds=build_detector_ds,
            build_frame_qchi=build_frame_qchi,
            frame_qchi_store=saxs_qchi_store,
            progress=_progress_wrap,
        )
        _progress_offset = 2 + _n_saxs_frames  # advance past saxs frames
        t_saxs_end = _time.perf_counter()

    # -- WAXS branch --
    waxs_result: dict[str, Any] | None = None
    t_waxs_start = t_waxs_end = _time.perf_counter()
    if has_waxs:
        _t_debug = _time.perf_counter()
        _waxs_geo_kw: dict[str, Any] = {}
        if waxs_beam_delta_px is not None:
            _waxs_geo_kw["beam_delta_row_px"] = waxs_beam_delta_px[0]
            _waxs_geo_kw["beam_delta_col_px"] = waxs_beam_delta_px[1]
        waxs_geo = resolve_waxs_geometry(run, **_waxs_geo_kw)
        print(f"[mask_setup] resolve_waxs_geometry: "
              f"{_time.perf_counter() - _t_debug:.3f}s")

        # WAXS mask callable
        _t_debug = _time.perf_counter()
        waxs_mask_fn = None
        if waxs_mask is None:
            waxs_mask = waxs_kw.pop("mask_path", None)
        from smi_tiled.defaults import resolve_mask_path
        if not isinstance(waxs_mask, dict):
            waxs_mask = resolve_mask_path(waxs_mask, detector="waxs")
        if waxs_mask is not None:
            waxs_bsx_pf = np.asarray(
                waxs_raw.attrs.get("smi_waxs_bsx_per_frame", []),
                dtype=float,
            )
            # waxs_bsx_ref is the bsx position where the mask polygon was
            # drawn (typically arc ≈ 0°).  If the scan started at a different
            # arc angle the first-frame bsx will be offset and we must NOT
            # use it as the reference.  Prefer an explicit value from
            # waxs_kwargs; fall back to computing the arc-0 bsx from the
            # known linear bsx-vs-arc relationship (~-4.4 mm/deg at SMI).
            _BSX_PER_ARC_DEG = -4.39  # mm/deg, SMI mechanical linkage
            if "waxs_bsx_ref" in waxs_kw:
                waxs_bsx_ref = float(waxs_kw.pop("waxs_bsx_ref"))
            elif waxs_bsx_pf.size >= 2:
                arc_pf = np.asarray(
                    waxs_raw.coords[waxs_raw.dims[0]].values, dtype=float
                )
                if arc_pf.shape[0] == waxs_bsx_pf.shape[0] and (arc_pf.max() - arc_pf.min()) > 0.5:
                    # Arc was scanned — fit slope and extrapolate to arc=0
                    slope = np.polyfit(arc_pf, waxs_bsx_pf, 1)[0]
                    waxs_bsx_ref = float(waxs_bsx_pf[0] - slope * arc_pf[0])
                else:
                    # Fixed arc with multiple frames — use known slope
                    arc_val = float(arc_pf[0])
                    waxs_bsx_ref = float(
                        waxs_bsx_pf[0] - _BSX_PER_ARC_DEG * arc_val
                    )
            elif waxs_bsx_pf.size == 1:
                # Single-frame fixed arc — use known slope
                arc_val = float(
                    waxs_raw.coords[waxs_raw.dims[0]].values[0]
                )
                waxs_bsx_ref = float(
                    waxs_bsx_pf[0] - _BSX_PER_ARC_DEG * arc_val
                )
            else:
                waxs_bsx_ref = 0.0
            waxs_mask_fn = make_waxs_mask_callable(
                waxs_mask,                                  # accepts dict or Path
                waxs_bsx_ref=waxs_bsx_ref,
                beamstop_max_abs_arc_deg=waxs_kw.pop(
                    "beamstop_max_abs_arc_deg", 15.0
                ),
            )
        print(f"[mask_setup] WAXS mask creation: "
              f"{_time.perf_counter() - _t_debug:.3f}s")

        # Build WAXS calibration
        cal_dict: dict[str, Any] = dict(_DEFAULT_CAL)
        cal_dict["beam_center_row"] = waxs_geo.beam_center_row_px
        cal_dict["beam_center_col"] = waxs_geo.beam_center_col_px
        cal_dict["energy_kev"] = waxs_geo.energy_ev / 1000.0
        # sample_distance_mm: use _DEFAULT_CAL (273 mm) rather than the
        # motor reading; the calibrated value is more accurate for SMI.
        # Users can still override via waxs_kwargs.
        if waxs_beam_col_per_arc_deg != 0:
            cal_dict["beam_col_per_arc_deg"] = waxs_beam_col_per_arc_deg
        cal_override_keys = set(WAXSCalibration.__dataclass_fields__.keys())
        for k in list(waxs_kw.keys()):
            if k in cal_override_keys:
                cal_dict[k] = waxs_kw.pop(k)
        waxs_cal = WAXSCalibration(**cal_dict)

        _progress_advance_stage("waxs_setup", 0)  # emit current position

        t_waxs_start = _time.perf_counter()
        waxs_result = integrate_waxs(
            waxs_raw=waxs_raw,
            mask_fn=waxs_mask_fn,
            n_q=n_q,
            n_chi=n_chi,
            cal=waxs_cal,
            solid_angle_correction=solid_angle_correction,
            flip_horizontal=bool(opts.get("waxs_flip_horizontal", False)),
            qx_shift_nm=float(opts.get("waxs_qx_shift_nm", 0.0)),
            qy_shift_nm=float(opts.get("waxs_qy_shift_nm", 0.0)),
            dezinger_threshold=dezinger_threshold,
            dezinger_kernel=dezinger_kernel,
            cache_geometry=cache_geometry,
            pixel_splitting=pixel_splitting,
            build_detector_ds=build_detector_ds,
            build_frame_qchi=build_frame_qchi,
            frame_qchi_store=waxs_qchi_store,
            progress=_progress_wrap,
        )
        _progress_offset = 2 + _n_saxs_frames + _n_waxs_frames
        t_waxs_end = _time.perf_counter()

    t_mask = _time.perf_counter()

    # Merge (handles None gracefully)
    t_merge_start = _time.perf_counter()
    saxs_qchi = saxs_result["q_chi"] if saxs_result else None
    waxs_qchi = waxs_result["q_chi"] if waxs_result else None
    saxs_iq = saxs_result["iq"] if saxs_result else None
    waxs_iq = waxs_result["iq"] if waxs_result else None

    merged_qchi = merge_q_chi_weighted(saxs_qchi, waxs_qchi, n_q=n_q, n_chi=n_chi)
    merged_iq = merge_iq_profiles(merged_qchi, saxs_iq, waxs_iq)
    per_frame_iq = _build_per_frame_iq(merged_iq, saxs_result, waxs_result, scan_info=scan_info)
    t_merge_end = _time.perf_counter()

    _progress_advance_stage("merge", 1)

    timing = {
        "total": t_merge_end - t0,
        "tiled_load": t_load - t0,
        "mask_setup": t_mask - t_load,
        "saxs_integrate": t_saxs_end - t_saxs_start,
        "waxs_integrate": t_waxs_end - t_waxs_start,
        "merge": t_merge_end - t_merge_start,
    }

    # Free cached baseline data for this run
    clear_baseline_cache()

    # Promote per-frame qchi from the per-detector dicts to a public
    # ``per_frame_qchi`` mapping for downstream consumers (line cuts,
    # upload schema).  The dict is keyed by detector name; absent
    # detectors are omitted.
    per_frame_qchi_map: dict[str, xr.Dataset] | None = None
    _pfq: dict[str, xr.Dataset] = {}
    if saxs_result and saxs_result.get("q_chi_frames") is not None:
        _pfq["saxs"] = saxs_result["q_chi_frames"]
    if waxs_result and waxs_result.get("q_chi_frames") is not None:
        _pfq["waxs"] = waxs_result["q_chi_frames"]
    if _pfq:
        per_frame_qchi_map = _pfq

    result = CombinedReductionResult(
        uid=uid,
        scan_info=scan_info,
        saxs=saxs_result,
        waxs=waxs_result,
        merged_qchi=merged_qchi,
        merged_iq=merged_iq,
        per_frame_iq=per_frame_iq,
        timing=timing,
        geometry=geometry,
        incident_angle_deg=incident_angle_deg,
        per_frame_qchi=per_frame_qchi_map,
    )

    # ---- Optional derived-analysis stages -------------------------------
    # Each is opt-in via a kwarg.  ``virtual_axes`` defaults to running
    # with the default config so ``fn:*`` axes from per-frame strings
    # always appear when source data is available.
    from .derived import apply_virtual_axes, apply_line_cuts, apply_peak_fits
    from .derived import VirtualAxesConfig
    va_cfg = virtual_axes if virtual_axes is not None else VirtualAxesConfig()
    apply_virtual_axes(result, va_cfg)
    if line_cuts:
        apply_line_cuts(result, list(line_cuts))
    if peak_fits:
        apply_peak_fits(result, list(peak_fits))
    return result
