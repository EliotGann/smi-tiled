# Reduction (`reduce_smi_combined`)

The top-level entry point is {func}`smi_tiled.reduce_smi_combined`:

```python
from smi_tiled import reduce_smi_combined

result = reduce_smi_combined(
    uid="6e61b977-…",
    tiled_uri="https://tiled.nsls2.bnl.gov",
    catalog="smi/migration",
    n_q=2000, n_chi=360,
    solid_angle_correction=True,
    geometry="transmission",
)
```

This drives the full pipeline:

1. **Load** raw SAXS + WAXS frames via {class}`~smi_tiled.TiledSMISWAXSLoader`.
2. **Resolve geometry** ({func}`~smi_tiled.resolve_saxs_geometry`,
   {func}`~smi_tiled.resolve_waxs_geometry`).
3. **Build masks** — bundled instrument geometry + beamstop + dynamic
   per-frame (see {doc}`masks`).
4. **Dezinger** hot pixels with a configurable σ threshold.
5. **Compute per-pixel q-maps** — exact 3-D pixel positions for both
   detectors (multi-panel arc model for WAXS).
6. **Bin** into `(q, χ)` via `histogram2d` with optional sub-pixel
   splitting.
7. **Merge** SAXS + WAXS into a common `(q, χ)` grid with count
   weighting; azimuthally average to merged `I(q)`.

## The result object

{class}`smi_tiled.CombinedReductionResult`:

| Field | Type | Dims / vars |
|---|---|---|
| `uid` | `str` | source run UID |
| `scan_info` | `dict` | from {func}`~smi_tiled.infer_detectors_and_steps` |
| `saxs`, `waxs` | `dict` or `None` | per-detector intermediates |
| `merged_qchi` | `xr.Dataset` or `None` | `(q, chi)` with `intensity, counts, saxs_intensity, waxs_intensity, saxs_counts, waxs_counts` |
| `merged_iq` | `xr.Dataset` or `None` | `(q,)` with `I, counts, saxs_I, waxs_I` |
| `per_frame_iq` | `xr.Dataset` or `None` | `(frame, q)` with `I, saxs_I, waxs_I` + per-frame motor scalars |
| `timing` | `dict` | per-stage wall-clock seconds |
| `geometry` | `"transmission"` or `"grazing_incidence"` | |
| `incident_angle_deg` | `float` | from sample_name or override |

For PyHyperScattering accessors that need a `DataArray`:

```python
da = result.to_dataarray("merged_iq")     # or "merged_qchi", "per_frame_iq"
```

See {meth}`~smi_tiled.CombinedReductionResult.to_dataarray`.

## Common parameter recipes

### Long-SDD (USAXS-style) scan

```python
reduce_smi_combined(
    uid="…",
    n_q=3000,              # higher q-resolution at long SDD
    pixel_splitting=3,     # smoother histograms
    cache_geometry=False,  # one-shot, don't pay the cache fill cost
)
```

### Disabling dynamic SAXS masks

The WAXS-shadow occlusion and AgBh aperture cut are on by default.
To turn one or both off:

```python
reduce_smi_combined(
    uid="…",
    saxs_kwargs={
        "dynamic_saxs_kwargs": {
            "aperture":    {"enabled": False},
            "waxs_shadow": {"enabled": False},
        },
    },
)
```

### Custom q-cutoff for SAXS

```python
reduce_smi_combined(
    uid="…",
    saxs_q_cutoff=0.6,     # nm⁻¹
)
# Or by AgBh ring order:
reduce_smi_combined(uid="…", saxs_agbh_ring_order=8)
```

### Manual beam-center override

When the calibration is off (e.g. an unfamiliar motor configuration),
bypass the resolver:

```python
reduce_smi_combined(
    uid="…",
    saxs_kwargs={
        "beam_center_row_px": 1107.0,
        "beam_center_col_px": 744.0,
    },
)
# Or add a small delta to whatever the resolver produced:
reduce_smi_combined(uid="…", saxs_beam_delta_px=(2.0, -1.5))
```

### Custom masks (in-memory dict)

See {doc}`masks` for the full schema.

```python
import json
from smi_tiled.defaults import default_saxs_mask_path

spec = json.load(open(default_saxs_mask_path()))
spec["static_regions"]["my_extra_blob"] = [[400, 400], [500, 400], [500, 500], [400, 500]]

reduce_smi_combined(uid="…", saxs_mask=spec)
```

## Grazing-incidence variant

For grazing-incidence WAXS:

```python
from smi_tiled import reduce_smi_gi

gi = reduce_smi_gi(
    uid="…",
    incident_angle_deg=0.12,   # or auto-detect from sample_name
    n_qxy=500, n_qz=500,
)
print(gi.summed_ds)            # I(qxy, qz)
```

See {class}`smi_tiled.GIReductionResult`.

## Performance: geometry caching

Computing per-pixel q-maps is the most expensive single step
(~0.1 s SAXS, several seconds WAXS).  `cache_geometry=True` (the default)
keys precomputed maps on the full calibration tuple and reuses them
across scans with identical geometry:

```python
from smi_tiled import geometry_cache_info, clear_geometry_cache

# After processing many scans:
info = geometry_cache_info()
print(f"{info['waxs_entries']} WAXS keys, ~{info['estimated_mb']:.1f} MB")

# Manually free memory (e.g. in a long-running GUI):
clear_geometry_cache()
```

The cache key includes every geometry-affecting parameter, so changing
the energy, distance, beam center, panel offsets, mask, or beam-center
deltas produces a new entry — there's no stale-data risk, only a memory
cost.

## Internals (when to dive into the integrator)

For most users `reduce_smi_combined` is enough.  If you want finer
control — separate SAXS and WAXS pipelines, custom merging, etc. — drop
down to:

- {func}`smi_tiled.integrate_saxs` — single-detector SAXS pipeline
- {func}`smi_tiled.integrate_waxs` — single-detector WAXS pipeline
- {class}`smi_tiled.MultiPanelArcDetector` — exact 3-D WAXS geometry
- {func}`smi_tiled.merge_q_chi_weighted` — weighted SAXS+WAXS merge
- {func}`smi_tiled.merge_multiple_qchi` — N-scan merge (e.g. for
  tile/mosaic averaging)
