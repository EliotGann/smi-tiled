"""Tests for `smi_tiled.derived.virtual_axes` (ported from smi-browser)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from smi_tiled.derived.virtual_axes import (
    VirtualAxesConfig,
    apply_virtual_axes,
    derive_virtual_columns,
    parse_label_number_tokens,
)


class TestParseLabelNumberTokens:
    def test_real_target_file_name(self):
        got = parse_label_number_tokens(
            "Lucas_sample2_pos1_2450.00eV_ai0.50_wa9_bpm1.995_degC100.0")
        assert got == {
            "sample": 2.0, "pos": 1.0, "eV": 2450.0, "ai": 0.5,
            "wa": 9.0, "bpm": 1.995, "degC": 100.0,
        }

    def test_bare_number_and_textonly_excluded(self):
        assert parse_label_number_tokens("Lucas_120_run") == {}

    def test_unit_suffix_used_when_no_prefix(self):
        assert parse_label_number_tokens("2450.00eV") == {"eV": 2450.0}

    def test_negative_and_non_string(self):
        assert parse_label_number_tokens("x-3.5")["x"] == -3.5
        assert parse_label_number_tokens(None) == {}
        assert parse_label_number_tokens("") == {}

    def test_bytes_from_cache_are_decoded(self):
        assert parse_label_number_tokens(b"x_ai0.50_eV2450.00") == {
            "ai": 0.5, "eV": 2450.0}
        assert parse_label_number_tokens(np.bytes_(b"x_ai0.50")) == {"ai": 0.5}
        assert parse_label_number_tokens(np.str_("x_ai0.50")) == {"ai": 0.5}


class TestDeriveVirtualColumns:
    def _frame(self):
        return pd.DataFrame({
            "target_file_name": [
                "Lucas_pos1_2450.00eV_ai0.50_degC100.0",
                "Lucas_pos2_2460.00eV_ai4.00",
            ],
            "energy_energy": [2450.1, 2460.2],
            "ts_target_file_name": [1.0, 2.0],
            "att2_9_status": ["+- 120 uA", "+- 120 uA"],
        })

    def test_adds_prefixed_columns(self):
        out = derive_virtual_columns(self._frame())
        assert {"fn:ai", "fn:eV", "fn:pos", "fn:degC"} <= set(out.columns)
        np.testing.assert_allclose(out["fn:ai"], [0.5, 4.0])

    def test_missing_token_is_nan(self):
        out = derive_virtual_columns(self._frame())
        assert out["fn:degC"].iloc[0] == 100.0
        assert np.isnan(out["fn:degC"].iloc[1])

    def test_ignores_numeric_and_ts_and_noise(self):
        out = derive_virtual_columns(self._frame())
        assert not any(c.startswith("fn:") and "energy" in c for c in out.columns)
        assert "fn:120" not in out.columns

    def test_min_fill_drops_sparse_columns(self):
        df = pd.DataFrame({"target_file_name": [
            "run_ai0.50", "run", "run", "run",
        ]})
        out = derive_virtual_columns(df, min_fill=0.5)
        assert "fn:ai" not in out.columns

    def test_collision_across_sources_is_qualified(self):
        df = pd.DataFrame({
            "target_file_name": ["a_ai0.5", "a_ai0.6"],
            "other": ["ai9", "ai8"],
        })
        out = derive_virtual_columns(df)
        assert "fn:ai" in out.columns
        assert "fn:other:ai" in out.columns

    def test_empty_frame_returned_unchanged(self):
        empty = pd.DataFrame()
        assert derive_virtual_columns(empty).empty

    def test_bytes_object_column(self):
        df = pd.DataFrame({"target_file_name": np.array(
            [b"r_ai0.50_eV2450", b"r_ai4.00_eV2460"], dtype=object)})
        out = derive_virtual_columns(df)
        assert {"fn:ai", "fn:eV"} <= set(out.columns)
        np.testing.assert_allclose(out["fn:ai"], [0.5, 4.0])


class TestApplyVirtualAxes:
    def _result(self):
        per_frame_iq = xr.Dataset(
            {
                "I": (("frame", "q"), np.zeros((2, 3))),
                # String per-frame field already attached as object dtype.
                "target_file_name": (("frame",), np.array(
                    ["x_ai0.50_eV2450", "x_ai4.00_eV2460"], dtype=object,
                )),
            },
            coords={"frame": [0, 1], "q": [0.1, 0.2, 0.3]},
        )

        class _Result:
            pass

        r = _Result()
        r.per_frame_iq = per_frame_iq
        r.scan_info = {}
        return r

    def test_attaches_fn_axes_as_data_vars(self):
        r = self._result()
        apply_virtual_axes(r)
        assert "fn:ai" in r.per_frame_iq.data_vars
        np.testing.assert_allclose(r.per_frame_iq["fn:ai"].values, [0.5, 4.0])
        np.testing.assert_allclose(r.per_frame_iq["fn:eV"].values, [2450.0, 2460.0])

    def test_disabled_is_noop(self):
        r = self._result()
        before = set(r.per_frame_iq.data_vars)
        apply_virtual_axes(r, VirtualAxesConfig(enabled=False))
        assert set(r.per_frame_iq.data_vars) == before

    def test_primary_strings_from_scan_info(self):
        # When string fields live only in scan_info["primary_strings"]
        # (the loader path), they should still be parsed.
        per_frame_iq = xr.Dataset(
            {"I": (("frame", "q"), np.zeros((2, 3)))},
            coords={"frame": [0, 1], "q": [0.1, 0.2, 0.3]},
        )

        class _Result:
            pass

        r = _Result()
        r.per_frame_iq = per_frame_iq
        r.scan_info = {
            "primary_strings": {
                "target_file_name": np.array(
                    ["x_ai0.50", "x_ai4.00"], dtype=object,
                ),
            },
        }
        apply_virtual_axes(r)
        assert "fn:ai" in r.per_frame_iq.data_vars
