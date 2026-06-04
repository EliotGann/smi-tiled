"""Line cuts (cross sections) of 2-D reduction products.

A *line cut* is a band-averaged 1-D slice through a 2-D image — either
horizontal (``kind="h"`` → average over a y-band → ``I(x)``) or vertical
(``kind="v"`` → average over an x-band → ``I(y)``).

The core math lives in :func:`compute_cross_section` (ported from
``smi_browser.figures.cuts``).  :func:`apply_line_cuts` is the
result-aware driver that walks a
:class:`~smi_tiled.integrator.CombinedReductionResult` (or GI result)
and packs all requested cuts into a ``dict[name, xr.Dataset]`` carried
on ``result.line_cuts``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Sequence

import numpy as np
import xarray as xr


# ---------------------------------------------------------------------------
# Core math
# ---------------------------------------------------------------------------

def compute_cross_section(cut: dict, x, y, image,
                          x_label: str = "", y_label: str = ""):
    """Compute a 1-D cross section through a 2-D image.

    Parameters
    ----------
    cut : dict
        ``{"kind": "h" | "v", "center": float, "width": float}``.
    x, y : 1-D arrays
        Coordinate axes (length matches columns / rows of *image*).
    image : 2-D array
        Indexed as ``image[row_y, col_x]`` (i.e. shape ``(len(y), len(x))``).
    x_label, y_label : str, optional
        Axis labels returned to the caller for plotting.

    Returns
    -------
    tuple
        ``(axis, intensity, axis_label)`` or ``None`` if inputs are
        missing.  An ``h`` cut returns ``(x, I(x), x_label)`` — band
        averaged over y.  A ``v`` cut returns ``(y, I(y), y_label)``.
        When the requested band is too narrow to capture any pixel,
        the nearest single row/column is used.
    """
    if x is None or y is None or image is None:
        return None
    c = float(cut["center"])
    w = float(cut.get("width", 0.0)) or 0.0
    half = max(w / 2.0, 0.0)
    if cut["kind"] == "h":
        mask = (y >= c - half) & (y <= c + half)
        if not np.any(mask):
            idx = int(np.argmin(np.abs(y - c)))
            section = image[idx, :].astype(float)
        else:
            section = np.nanmean(image[mask, :], axis=0)
        return x, section, x_label
    mask = (x >= c - half) & (x <= c + half)
    if not np.any(mask):
        idx = int(np.argmin(np.abs(x - c)))
        section = image[:, idx].astype(float)
    else:
        section = np.nanmean(image[:, mask], axis=1)
    return y, section, y_label


# ---------------------------------------------------------------------------
# Spec
# ---------------------------------------------------------------------------

#: Supported targets — names of the 2-D product a cut is drawn against.
#: ``merged_qchi`` is the single merged image; ``saxs_qchi`` /
#: ``waxs_qchi`` are per-frame stacks from the individual detectors;
#: ``qxy_qz`` is the GI per-frame product.
LineCutTarget = Literal["merged_qchi", "saxs_qchi", "waxs_qchi", "qxy_qz"]


@dataclass(frozen=True)
class LineCutSpec:
    """A single line cut to compute.

    Parameters
    ----------
    kind : {"h", "v"}
        Cut orientation.  ``h`` averages over a band along the y-axis
        and returns ``I(x)``; ``v`` averages over a band along the
        x-axis and returns ``I(y)``.
    center, width : float
        Band centre and full width (same units as the corresponding
        axis on the target image).
    target : str
        Which 2-D product to cut.  Defaults to ``"merged_qchi"``.
    name : str, optional
        Stable identifier (used as the dict key in
        ``result.line_cuts``).  Auto-generated from kind/center if
        omitted.
    """

    kind: Literal["h", "v"]
    center: float
    width: float
    target: str = "merged_qchi"
    name: str | None = None

    def resolved_name(self) -> str:
        if self.name:
            return self.name
        return f"{self.target}:{self.kind}_{self.center:.6g}_w{self.width:.6g}"

    def to_provenance(self) -> dict[str, Any]:
        """Stable dict form for hashing into ``reduction_hash``."""
        return {
            "kind": str(self.kind),
            "center": float(self.center),
            "width": float(self.width),
            "target": str(self.target),
            "name": self.resolved_name(),
        }


# ---------------------------------------------------------------------------
# Source resolution
# ---------------------------------------------------------------------------

def _axes_for_target(target: str) -> tuple[str, str]:
    """Return ``(x_dim, y_dim)`` of the 2-D image for *target*."""
    if target == "qxy_qz":
        return ("qxy", "qz")
    return ("q", "chi")


def _resolve_image_stack(
    result: Any, target: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, bool]:
    """Pull ``(x, y, images, per_frame)`` for *target* from *result*.

    ``images`` is shape ``(n_frames, n_y, n_x)`` when ``per_frame`` is
    True, else ``(n_y, n_x)``.  Raises ``ValueError`` if the requested
    product is unavailable on this result.
    """
    x_dim, y_dim = _axes_for_target(target)

    if target == "merged_qchi":
        ds = getattr(result, "merged_qchi", None)
        if ds is None or "intensity" not in ds:
            raise ValueError("merged_qchi unavailable on this result")
        x = np.asarray(ds[x_dim].values, dtype=float)
        y = np.asarray(ds[y_dim].values, dtype=float)
        # merged_qchi dims are (q, chi) → transpose to (y=chi, x=q).
        img = np.asarray(ds["intensity"].transpose(y_dim, x_dim).values, dtype=float)
        return x, y, img, False

    if target in ("saxs_qchi", "waxs_qchi"):
        side = target.split("_")[0]
        det = getattr(result, side, None)
        if not det or det.get("q_chi_frames") is None:
            raise ValueError(f"{target} unavailable on this result")
        ds = det["q_chi_frames"]
        x = np.asarray(ds[x_dim].values, dtype=float)
        y = np.asarray(ds[y_dim].values, dtype=float)
        # frames have dims (frame, q, chi) → (frame, chi, q).
        img = np.asarray(
            ds["intensity"].transpose("frame", y_dim, x_dim).values,
            dtype=float,
        )
        return x, y, img, True

    if target == "qxy_qz":
        # GI result.  Per-frame from frames list (or q_chi_frames Dataset).
        frames = getattr(result, "frames", None)
        qxy_grid = getattr(result, "qxy_grid", None)
        qz_grid = getattr(result, "qz_grid", None)
        if frames is None or qxy_grid is None or qz_grid is None:
            raise ValueError("qxy_qz unavailable on this result")
        x = np.asarray(qxy_grid, dtype=float)
        y = np.asarray(qz_grid, dtype=float)
        # GI frames are (qxy, qz) → (qz, qxy).
        img = np.asarray(np.stack(frames, axis=0), dtype=float).transpose(0, 2, 1)
        return x, y, img, True

    raise ValueError(f"Unknown LineCutSpec.target={target!r}")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def _cut_to_dataset(
    cut: LineCutSpec,
    x: np.ndarray,
    y: np.ndarray,
    images: np.ndarray,
    per_frame: bool,
) -> xr.Dataset:
    """Compute one cut into a single-cut xr.Dataset.

    Resulting dims are ``(frame, axis)`` for per-frame sources and
    ``(axis,)`` for the merged source.  ``axis`` is the coordinate name
    appropriate to the cut (e.g. ``"q"`` for an h-cut of qchi).
    """
    x_dim, y_dim = _axes_for_target(cut.target)
    axis_name = x_dim if cut.kind == "h" else y_dim
    axis_vals = x if cut.kind == "h" else y

    cut_dict = {"kind": cut.kind, "center": float(cut.center), "width": float(cut.width)}

    if per_frame:
        out = np.empty((images.shape[0], axis_vals.size), dtype=float)
        for fi in range(images.shape[0]):
            sec = compute_cross_section(cut_dict, x, y, images[fi])
            if sec is None:
                out[fi] = np.nan
            else:
                out[fi] = sec[1]
        ds = xr.Dataset(
            {"intensity": (("frame", axis_name), out)},
            coords={
                "frame": np.arange(images.shape[0], dtype=int),
                axis_name: axis_vals,
            },
        )
    else:
        sec = compute_cross_section(cut_dict, x, y, images)
        section = np.full(axis_vals.size, np.nan) if sec is None else sec[1]
        ds = xr.Dataset(
            {"intensity": ((axis_name,), section)},
            coords={axis_name: axis_vals},
        )

    ds.attrs.update({
        "kind": cut.kind,
        "center": float(cut.center),
        "width": float(cut.width),
        "target": cut.target,
        "axis_label": axis_name,
        "name": cut.resolved_name(),
    })
    return ds


def apply_line_cuts(
    result: Any,
    cuts: Sequence[LineCutSpec],
    *,
    per_frame: bool = True,
) -> Any:
    """Compute every cut in *cuts* and attach as ``result.line_cuts``.

    Parameters
    ----------
    result : CombinedReductionResult or GIReductionResult
        Receives the new product.  Frozen dataclasses are mutated via
        ``object.__setattr__``.
    cuts : sequence of LineCutSpec
        Cuts to compute.  Each is keyed in the output dict by
        :meth:`LineCutSpec.resolved_name`.
    per_frame : bool
        When ``True`` (the default) and the target supports per-frame
        data, the resulting Dataset carries a ``frame`` dim.  When
        ``False`` only the merged 2-D image is cut (``frame``-less
        result).

    Returns
    -------
    result
        The same *result* object, with ``result.line_cuts`` populated
        (or left ``None`` if *cuts* is empty).
    """
    if not cuts:
        return result

    out: dict[str, xr.Dataset] = {}
    for cut in cuts:
        try:
            x, y, images, has_frames = _resolve_image_stack(result, cut.target)
        except ValueError:
            continue  # silently skip cuts whose target is unavailable
        if has_frames and not per_frame:
            # Collapse to a single image by averaging frames.
            images = np.nanmean(images, axis=0)
            has_frames = False
        out[cut.resolved_name()] = _cut_to_dataset(cut, x, y, images, has_frames)

    if out:
        object.__setattr__(result, "line_cuts", out)
    return result
