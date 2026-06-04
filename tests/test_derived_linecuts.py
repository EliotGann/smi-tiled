"""Tests for `smi_tiled.derived.linecuts`."""
from __future__ import annotations

import numpy as np
import xarray as xr

from smi_tiled.derived.linecuts import (
    LineCutSpec,
    apply_line_cuts,
    compute_cross_section,
)


def test_compute_cross_section_h():
    x = np.linspace(0, 1, 5)
    y = np.linspace(0, 1, 4)
    # image[row_y, col_x] — uniform along x, increasing along y.
    img = np.tile(y[:, None], (1, x.size))
    axis, sec, label = compute_cross_section(
        {"kind": "h", "center": 0.5, "width": 0.4}, x, y, img,
        x_label="q", y_label="chi",
    )
    np.testing.assert_array_equal(axis, x)
    assert label == "q"
    # average of y values within [0.3, 0.7] (=> y=1/3 and y=2/3) = 0.5
    np.testing.assert_allclose(sec, np.full_like(x, 0.5))


def test_compute_cross_section_v_nearest_when_band_empty():
    x = np.linspace(0, 1, 5)
    y = np.linspace(0, 1, 4)
    img = np.tile(x[None, :], (y.size, 1))
    axis, sec, label = compute_cross_section(
        {"kind": "v", "center": 0.5, "width": 0.0}, x, y, img,
        x_label="q", y_label="chi",
    )
    np.testing.assert_array_equal(axis, y)
    assert label == "chi"
    # nearest column to x=0.5 is x[2]=0.5 → constant 0.5 over y.
    np.testing.assert_allclose(sec, np.full_like(y, 0.5))


def test_apply_line_cuts_merged_qchi():
    q = np.linspace(0.0, 2.0, 6)
    chi = np.linspace(-180.0, 180.0, 5)
    # Intensity = chi value (constant in q), so an h cut returns chi_center.
    intensity = np.tile(chi[None, :], (q.size, 1))   # dims (q, chi)
    counts = np.ones_like(intensity)
    merged_qchi = xr.Dataset(
        {"intensity": (("q", "chi"), intensity),
         "counts": (("q", "chi"), counts)},
        coords={"q": q, "chi": chi},
    )

    class _Result:
        pass

    r = _Result()
    r.merged_qchi = merged_qchi
    r.saxs = r.waxs = None
    r.line_cuts = None

    cuts = [
        LineCutSpec(kind="h", center=90.0, width=10.0, target="merged_qchi", name="cA"),
        LineCutSpec(kind="v", center=1.0, width=0.5, target="merged_qchi", name="cB"),
    ]
    apply_line_cuts(r, cuts)
    assert set(r.line_cuts) == {"cA", "cB"}
    da = r.line_cuts["cA"]["intensity"]
    assert da.dims == ("q",)
    # Mean over chi values in [85, 95] → only chi=90.0 → result 90.
    np.testing.assert_allclose(da.values, np.full(q.size, 90.0))

    db = r.line_cuts["cB"]["intensity"]
    assert db.dims == ("chi",)
    np.testing.assert_allclose(db.values, chi)


def test_apply_line_cuts_per_frame_saxs():
    q = np.array([0.0, 1.0, 2.0])
    chi = np.array([-90.0, 0.0, 90.0])
    # Two frames with intensity = frame_index everywhere.
    frames = np.stack([
        np.full((q.size, chi.size), 1.0),
        np.full((q.size, chi.size), 7.0),
    ])
    ds = xr.Dataset(
        {"intensity": (("frame", "q", "chi"), frames)},
        coords={"frame": [0, 1], "q": q, "chi": chi},
    )

    class _Result:
        pass

    r = _Result()
    r.merged_qchi = None
    r.saxs = {"q_chi_frames": ds}
    r.waxs = None
    r.line_cuts = None
    cuts = [LineCutSpec(kind="h", center=0.0, width=1.0,
                        target="saxs_qchi", name="band")]
    apply_line_cuts(r, cuts)
    out = r.line_cuts["band"]
    assert out["intensity"].dims == ("frame", "q")
    np.testing.assert_allclose(out["intensity"].values[0], [1.0, 1.0, 1.0])
    np.testing.assert_allclose(out["intensity"].values[1], [7.0, 7.0, 7.0])
    assert out.attrs["target"] == "saxs_qchi"


def test_apply_line_cuts_unknown_target_silently_skipped():
    class _Result:
        merged_qchi = None
        saxs = None
        waxs = None
        line_cuts = None

    r = _Result()
    apply_line_cuts(r, [LineCutSpec(kind="h", center=0.0, width=1.0,
                                    target="merged_qchi")])
    # No cuts produced, line_cuts stays None.
    assert r.line_cuts is None
