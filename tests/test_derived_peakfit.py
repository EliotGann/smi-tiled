"""Tests for `smi_tiled.derived.peakfit` (ported from smi-browser)."""
from __future__ import annotations

import threading

import numpy as np
import pytest
import xarray as xr

from smi_tiled.derived.peakfit import (
    PeakDef,
    apply_peak_fits,
    fit_peak_across_frames,
)

_GAUSS_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))


def _gauss(q, amp, mu, sigma):
    return amp * np.exp(-0.5 * ((q - mu) / sigma) ** 2)


@pytest.fixture
def synthetic():
    q = np.linspace(0.0, 10.0, 400)
    n = 30
    amps = np.linspace(1.0, 5.0, n)
    mus = np.linspace(4.0, 6.0, n)
    sigma = 0.3
    iq = np.stack([_gauss(q, a, m, sigma) for a, m in zip(amps, mus)])
    return q, iq, amps, mus, sigma


def test_recovers_gaussian_params(synthetic):
    q, iq, amps, mus, sigma = synthetic
    peak = PeakDef("p", q_min=3.0, q_max=7.0, model="gaussian", baseline="none")
    res = fit_peak_across_frames(q, iq, peak)
    assert res["success"].all()
    np.testing.assert_allclose(res["amplitude"], amps, rtol=1e-2)
    np.testing.assert_allclose(res["center"], mus, atol=2e-2)
    np.testing.assert_allclose(res["fwhm"], sigma * _GAUSS_FWHM, rtol=2e-2)
    np.testing.assert_allclose(
        res["area"], amps * sigma * np.sqrt(2 * np.pi), rtol=2e-2,
    )


def test_linear_baseline_subtraction():
    q = np.linspace(0.0, 10.0, 400)
    amp, mu, sigma = 3.0, 5.0, 0.3
    slope, intercept = 0.5, 2.0
    iq = (_gauss(q, amp, mu, sigma) + slope * q + intercept)[None, :]
    peak = PeakDef("p", q_min=3.0, q_max=7.0, model="gaussian", baseline="linear")
    res = fit_peak_across_frames(q, iq, peak)
    assert res["success"][0]
    assert res["amplitude"][0] == pytest.approx(amp, rel=2e-2)
    assert res["center"][0] == pytest.approx(mu, abs=2e-2)
    peak_nb = PeakDef("p", q_min=3.0, q_max=7.0, model="gaussian", baseline="none")
    res_nb = fit_peak_across_frames(q, iq, peak_nb)
    assert abs(res_nb["amplitude"][0] - amp) > abs(res["amplitude"][0] - amp)


def test_lorentzian_model():
    q = np.linspace(0.0, 10.0, 600)
    amp, mu, gamma = 4.0, 5.0, 0.4
    iq = (amp * gamma**2 / ((q - mu) ** 2 + gamma**2))[None, :]
    peak = PeakDef("p", q_min=3.0, q_max=7.0, model="lorentzian", baseline="none")
    res = fit_peak_across_frames(q, iq, peak)
    assert res["success"][0]
    assert res["amplitude"][0] == pytest.approx(amp, rel=2e-2)
    assert res["center"][0] == pytest.approx(mu, abs=2e-2)
    assert res["fwhm"][0] == pytest.approx(2 * gamma, rel=3e-2)


def test_nan_on_unfittable(synthetic):
    q, iq, *_ = synthetic
    iq = iq.copy()
    iq[5] = np.nan
    peak = PeakDef("p", q_min=3.0, q_max=7.0, model="gaussian", baseline="none")
    res = fit_peak_across_frames(q, iq, peak)
    assert not res["success"][5]
    assert np.isnan(res["amplitude"][5])
    assert res["success"][0]


def test_range_too_narrow_returns_all_nan(synthetic):
    q, iq, *_ = synthetic
    peak = PeakDef("p", q_min=4.99, q_max=5.0, model="gaussian", baseline="none")
    res = fit_peak_across_frames(q, iq, peak)
    assert not res["success"].any()
    assert np.isnan(res["amplitude"]).all()


def test_cancellation_returns_early(synthetic):
    q, iq, *_ = synthetic
    cancel = threading.Event()
    cancel.set()
    peak = PeakDef("p", q_min=3.0, q_max=7.0, model="gaussian", baseline="none")
    res = fit_peak_across_frames(q, iq, peak, cancel=cancel, cancel_check_every=1)
    assert not res["success"].any()


def test_progress_reports_completion(synthetic):
    q, iq, *_ = synthetic
    seen = []
    peak = PeakDef("p", q_min=3.0, q_max=7.0, model="gaussian", baseline="none")
    fit_peak_across_frames(q, iq, peak, progress=lambda d, t: seen.append((d, t)))
    assert seen[-1] == (iq.shape[0], iq.shape[0])


def test_no_peak_frame_reports_zero():
    rng = np.random.default_rng(0)
    q = np.linspace(0.0, 10.0, 400)
    peak_frame = _gauss(q, 5.0, 5.0, 0.3) + 0.01 * rng.standard_normal(q.size)
    flat_frame = 1.0 + 0.01 * rng.standard_normal(q.size)
    iq = np.stack([peak_frame, flat_frame])
    peak = PeakDef("p", q_min=4.0, q_max=6.0, model="gaussian",
                   baseline="linear", link="independent")
    res = fit_peak_across_frames(q, iq, peak)
    assert res["success"][0]
    assert not res["success"][1]
    assert res["amplitude"][1] == 0.0
    assert res["area"][1] == 0.0
    assert np.isnan(res["center"][1])
    assert np.isnan(res["fwhm"][1])


def test_width_cannot_exceed_drawn_range():
    q = np.linspace(0.0, 10.0, 400)
    iq = _gauss(q, 5.0, 5.0, 2.0)[None, :]
    peak = PeakDef("p", q_min=4.5, q_max=5.5, model="gaussian",
                   baseline="none", link="independent")
    res = fit_peak_across_frames(q, iq, peak)
    if res["success"][0]:
        assert res["fwhm"][0] <= (5.5 - 4.5) + 1e-9


def test_linked_shares_center_and_width():
    q = np.linspace(0.0, 10.0, 400)
    amps = np.array([1.0, 2.0, 3.0, 4.0])
    mu, sigma = 5.0, 0.3
    iq = np.stack([_gauss(q, a, mu, sigma) for a in amps])
    peak = PeakDef("p", q_min=3.0, q_max=7.0, model="gaussian",
                   baseline="none", link="linked")
    res = fit_peak_across_frames(q, iq, peak)
    assert res["success"].all()
    centers = res["center"][res["success"]]
    assert np.allclose(centers, centers[0])
    assert res["center"][0] == pytest.approx(mu, abs=2e-2)
    np.testing.assert_allclose(res["amplitude"], amps, rtol=5e-2)


def test_linked_falls_back_when_no_aggregate_peak():
    rng = np.random.default_rng(1)
    q = np.linspace(0.0, 10.0, 200)
    iq = 1.0 + 0.01 * rng.standard_normal((5, q.size))
    peak = PeakDef("p", q_min=4.0, q_max=6.0, model="gaussian",
                   baseline="linear", link="linked")
    res = fit_peak_across_frames(q, iq, peak)
    assert "note" in res
    assert not res["success"].any()


def test_bg_factor_widens_baseline_window():
    q = np.linspace(0.0, 10.0, 600)
    amp, mu, sigma = 3.0, 5.0, 0.3
    slope, intercept = 0.5, 2.0
    iq = (_gauss(q, amp, mu, sigma) + slope * q + intercept)[None, :]
    peak = PeakDef("p", q_min=4.5, q_max=5.5, model="gaussian",
                   baseline="linear", link="independent", bg_factor=3.0)
    res = fit_peak_across_frames(q, iq, peak)
    assert res["success"][0]
    assert res["amplitude"][0] == pytest.approx(amp, rel=5e-2)
    assert res["center"][0] == pytest.approx(mu, abs=2e-2)


def test_peakdef_key_includes_link_and_bg():
    a = PeakDef("p", 3.0, 7.0, link="linked", bg_factor=2.0)
    b = PeakDef("p", 3.0, 7.0, link="independent", bg_factor=2.0)
    c = PeakDef("p", 3.0, 7.0, link="linked", bg_factor=3.0)
    assert a.key() != b.key()
    assert a.key() != c.key()


# --- apply_peak_fits driver ------------------------------------------------

def test_apply_peak_fits_packs_dataset(synthetic):
    q, iq, *_ = synthetic
    per_frame_iq = xr.Dataset(
        {"I": (("frame", "q"), iq)},
        coords={"frame": np.arange(iq.shape[0]), "q": q},
    )

    class _Result:
        per_frame_iq = None
        peak_fits = None

    r = _Result()
    r.per_frame_iq = per_frame_iq
    peaks = [
        PeakDef("p1", q_min=3.0, q_max=7.0, baseline="none"),
        PeakDef("p2", q_min=3.5, q_max=6.5, baseline="none"),
    ]
    apply_peak_fits(r, peaks)
    assert r.peak_fits is not None
    ds = r.peak_fits
    assert ds.sizes == {"peak": 2, "frame": iq.shape[0]}
    assert set(ds.data_vars) >= {"amplitude", "center", "fwhm", "area",
                                  "success", "peak_key", "note"}
    assert list(ds["peak"].values) == ["p1", "p2"]
    assert ds["success"].values.all()
    assert isinstance(ds.attrs["peaks"], list)
    assert len(ds.attrs["peaks"]) == 2


def test_apply_peak_fits_disambiguates_duplicate_names():
    q = np.linspace(0.0, 10.0, 200)
    iq = _gauss(q, 1.0, 5.0, 0.3)[None, :]
    per_frame_iq = xr.Dataset(
        {"I": (("frame", "q"), iq)},
        coords={"frame": [0], "q": q},
    )

    class _Result:
        per_frame_iq = None
        peak_fits = None

    r = _Result()
    r.per_frame_iq = per_frame_iq
    apply_peak_fits(r, [
        PeakDef("dup", 3.0, 7.0, baseline="none"),
        PeakDef("dup", 4.0, 6.0, baseline="none"),
    ])
    assert list(r.peak_fits["peak"].values) == ["dup", "dup#2"]


def test_apply_peak_fits_empty_is_noop():
    class _Result:
        per_frame_iq = None
        peak_fits = "sentinel"

    r = _Result()
    apply_peak_fits(r, [])
    assert r.peak_fits == "sentinel"
