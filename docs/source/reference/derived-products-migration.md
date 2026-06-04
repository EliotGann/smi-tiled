# Derived analysis products migrated into `smi_tiled`

This page documents the analysis-product migration completed in response to
the audit at `smi-browser/docs/analysis_products_audit.md`. The goal: the
browser front-end no longer owns derived-from-reduction analysis. All such
products are produced inside `smi_tiled` and travel through the
`CombinedReductionResult` / `GIReductionResult` so they can be uploaded to
Tiled alongside the core reduction.

## What moved

Three derived-product modules previously lived in `smi-browser` and have been
ported here verbatim (with light packaging glue):

| Capability | Old home (`smi-browser`) | New home (`smi-tiled`) |
|---|---|---|
| Virtual axes from per-frame label-number strings (e.g. `fn:T_C`, `fn:V_kV`) | `smi_browser.data.scalars` | [`smi_tiled.derived.virtual_axes`](../../../src/smi_tiled/derived/virtual_axes.py) |
| Line cuts of 2D maps (q-cut at fixed χ, χ-cut at fixed q, qxy/qz cuts) | `smi_browser.figures.cuts` | [`smi_tiled.derived.linecuts`](../../../src/smi_tiled/derived/linecuts.py) |
| Per-frame peak fitting (Gaussian/Lorentzian, independent / linked / tracked) | `smi_browser.models.peakfit` | [`smi_tiled.derived.peakfit`](../../../src/smi_tiled/derived/peakfit.py) |

Public surface — all importable from `smi_tiled.derived`:

- Virtual axes: `VirtualAxesConfig`, `apply_virtual_axes`,
  `derive_virtual_columns`, `parse_label_number_tokens`, `VIRTUAL_PREFIX`
- Line cuts: `LineCutSpec`, `apply_line_cuts`, `compute_cross_section`
- Peak fits: `PeakDef`, `apply_peak_fits`, `fit_peak_across_frames`,
  `FIT_PARAMS`, `MIN_SNR`, `MIN_R2`

Behaviour was preserved. In particular, `peakfit` is a verbatim port (same
acceptance gate: SNR ≥ 3.0, R² ≥ 0.2, fitted width ≤ 97% of allowed bound),
and the `compute_cross_section` semantics for q-cut / chi-cut / qxy-cut /
qz-cut match the browser source 1:1.

## Wiring in the reduction pipeline

`reduce_smi_combined` and `reduce_smi_gi` accept three new optional kwargs:

```python
result = reduce_smi_combined(
    ...,
    virtual_axes=VirtualAxesConfig(sources=("sample_name",)),  # default config
    line_cuts=[
        LineCutSpec(kind="q_cut", center=0.0, width=5.0, target="merged_qchi"),
        LineCutSpec(kind="chi_cut", center=1.2, width=0.05),
    ],
    peak_fits=[
        PeakDef(kind="gaussian", x0=1.2, fwhm=0.05, mode="independent"),
    ],
)
```

All three default to `None` / no-op, so existing callers see no behavioural
change. When supplied, the pipeline:

1. Builds `per_frame_iq` as before.
2. **Promotes** `per_frame_qchi` onto the result as a
   `dict[str, xr.Dataset]` keyed by detector (`"saxs"`, `"waxs"`); previously
   this was discarded after merging.
3. Calls `apply_virtual_axes(result, va_cfg)` to add `fn:*` data_vars to
   `per_frame_iq` derived from string-typed columns + `scan_info["primary_strings"]`.
4. Calls `apply_line_cuts(result, line_cuts)` and stores the result as
   `result.line_cuts: dict[str, xr.Dataset]` keyed by cut name.
5. Calls `apply_peak_fits(result, peak_fits)` and stores the result as
   `result.peak_fits: xr.Dataset` with dims `(peak, frame)`.

Because `CombinedReductionResult` / `GIReductionResult` are frozen
dataclasses, the wire-up uses `object.__setattr__` to attach the optional
fields.

## Loader change required for virtual axes

`smi_tiled.loader._infer_from_cache` (HDF5 path, bulk-DataFrame path, and
fallback per-field path) now collects fields that *cannot* be cast to float
into a new dict `primary_strings: dict[str, np.ndarray]`. The reduction
exposes this as `scan_info["primary_strings"]`, which is what
`apply_virtual_axes` needs to parse `label=number` tokens from per-frame
strings (file names, sample names, …).

## Result-object schema additions

`CombinedReductionResult` has three new optional fields:

```python
per_frame_qchi: dict[str, xr.Dataset] | None = None  # {"saxs": ds, "waxs": ds}
line_cuts:      dict[str, xr.Dataset] | None = None  # {cut_name: ds}
peak_fits:      xr.Dataset | None = None             # dims: (peak, frame)
```

`GIReductionResult` gains the same `line_cuts` and `peak_fits` fields
(`per_frame_qchi` is N/A for GI because the merged map is `qxy_qz`, not
qchi-stacked).

## Upload schema additions

`smi_tiled.upload.schema`:

- `REDUCED_DATA_KEYS` adds: `"per_frame_qchi"`, `"line_cuts"`, `"peak_fits"`.
- `PROVENANCE_FIELDS` adds: `"virtual_axes_spec_hash"`,
  `"line_cuts_spec_hash"`, `"peak_fits_spec_hash"`.

`smi_tiled.upload.session`:

- New retrieval helpers: `get_per_frame_qchi(uid)`, `get_line_cuts(uid)`,
  `get_peak_fits(uid)`.
- `_provenance_dict` now includes the three derived spec hashes when the
  corresponding products are present, and folds them into `reduction_hash`.

The `upload(...)` body is still `NotImplementedError` (no behaviour change to
how data is written; only the schema and retrieval surface were extended).

## Tests added

- [`tests/test_derived_virtual_axes.py`](../../../tests/test_derived_virtual_axes.py) — 15 tests
- [`tests/test_derived_linecuts.py`](../../../tests/test_derived_linecuts.py) — 5 tests
- [`tests/test_derived_peakfit.py`](../../../tests/test_derived_peakfit.py) — 16 tests

All 36 new tests pass under `pixi run -- pytest tests/test_derived_*.py`. The
full suite is green.

---

## What `smi-browser` should now clean up

For the browser-side cleanup pass, the scope-of-deletion candidates are:

1. **`smi_browser/data/scalars.py`** — virtual-axis derivation. Browser code
   should now read `result.per_frame_iq` and use the `fn:*` data_vars that
   `smi_tiled` already produced. Or, if the browser receives a result that
   lacks them, call `smi_tiled.derived.apply_virtual_axes` instead of
   re-implementing the logic.

2. **`smi_browser/figures/cuts.py`** (`compute_cross_section` and the
   per-figure cut wrappers) — replace with
   `smi_tiled.derived.linecuts.compute_cross_section` for ad-hoc cuts, or
   read pre-computed cuts from `result.line_cuts[name]` when the reduction
   was run with `line_cuts=[…]`.

3. **`smi_browser/models/peakfit.py`** — delete in favour of
   `smi_tiled.derived.peakfit`. The acceptance gate, model functions, and
   per-frame/linked solvers are byte-equivalent. UI code should consume
   `result.peak_fits` (an `xr.Dataset` on dims `(peak, frame)` with vars
   `amplitude / center / fwhm / area / success / peak_key / note`) rather
   than re-fitting client-side.

4. **Front-end pipelines that materialise these products on demand** —
   migrate to passing `virtual_axes=`, `line_cuts=`, `peak_fits=` through to
   `reduce_smi_combined` / `reduce_smi_gi` so the products are computed
   once, attached to the reduction, and uploaded to Tiled with provenance
   hashes for reproducibility. The browser becomes a viewer of those
   products, not their producer.

5. **`primary_strings`** is now part of the reduction's `scan_info`. Any
   browser code that re-parsed string-typed primary fields out of the raw
   tiled run can drop that path and consume `scan_info["primary_strings"]`
   from the result.

The audit document remains the authoritative tracker — this page is the
contract for what `smi-tiled` now provides so that audit items 1–5 can be
closed on the browser side.
