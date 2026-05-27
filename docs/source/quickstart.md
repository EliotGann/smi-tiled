(quickstart)=
# Quick start

This page walks through reducing one SMI scan from a Tiled UID to a
1-D `I(q)` profile.

## Prerequisites

- `smi-tiled` installed with the `[tiled]` extra (see {ref}`installation`).
- Network access to the Tiled server (`https://tiled.nsls2.bnl.gov` for the
  production catalog).
- A run UID — e.g. by browsing the catalog from {meth}`smi_tiled.TiledSMISWAXSLoader.searchCatalog`,
  or from another tool like `smi-browser`.

## Minimal example

```python
from smi_tiled import reduce_smi_combined

result = reduce_smi_combined(
    uid="6e61b977-5761-489a-9a33-f3b6fdbdfc49",
    tiled_uri="https://tiled.nsls2.bnl.gov",
    catalog="smi/migration",
    n_q=2000,
    n_chi=360,
)

print(result.merged_iq)
print(result.merged_qchi)
```

The output is a {class}`smi_tiled.CombinedReductionResult` with three
xarray Datasets:

| Field | Dims | Variables |
|---|---|---|
| `merged_iq` | `(q,)` | `I, counts, saxs_I, waxs_I` |
| `merged_qchi` | `(q, chi)` | `intensity, counts, saxs_intensity, …` |
| `per_frame_iq` | `(frame, q)` | `I, saxs_I, waxs_I` (+ per-frame motor scalars) |

Plus per-detector intermediate products at `result.saxs` and `result.waxs`.

## Plotting

```python
import matplotlib.pyplot as plt
import numpy as np

iq = result.merged_iq
plt.loglog(iq["q"], iq["I"])
plt.xlabel(r"$q$ [nm$^{-1}$]")
plt.ylabel("Intensity")
plt.title(f"Scan {result.uid[:8]}…")
plt.show()
```

For 2-D plots:

```python
qchi = result.merged_qchi
plt.pcolormesh(qchi["q"], qchi["chi"], qchi["intensity"].T, norm="log")
plt.xlabel(r"$q$ [nm$^{-1}$]")
plt.ylabel(r"$\chi$ [deg]")
plt.colorbar(label="I")
```

## Using PyHyperScattering accessors (optional)

If you also have PyHyperScattering installed, extract a DataArray view
for the `rsoxs` and `fit` accessors:

```python
da = result.to_dataarray("merged_iq")    # or "merged_qchi", "per_frame_iq"
da.rsoxs.slice_q(0.5, q_width=0.05)      # PyHyper accessor
```

## Common parameters

`reduce_smi_combined` accepts many parameters; the most commonly used:

| Parameter | Default | Meaning |
|---|---|---|
| `n_q`, `n_chi` | `2000`, `360` | Output grid resolution |
| `cache_geometry` | `True` | Reuse precomputed q-maps across scans with identical geometry |
| `pixel_splitting` | `1` | Sub-pixel fractional binning (2-4 for cleanest histograms) |
| `image_cache_path` | `"auto"` | Use an HDF5 cache file if present (cuts tiled round-trips) |
| `dezinger_threshold` | `3000.0` | Hot-pixel rejection σ; `None` disables |
| `saxs_mask`, `waxs_mask` | bundled defaults | JSON path, polygon dict, or `None` |

See {doc}`user-guide/reduction` for the full list and discussion.

## Browsing the catalog

```python
from smi_tiled import TiledSMISWAXSLoader

loader = TiledSMISWAXSLoader()
df = loader.searchCatalog(
    sample="AGB",
    since="2026-05-01",
    limit=20,
)
print(df[["scan_id", "sample_name", "plan_name", "uid"]])
```

## Next steps

- {doc}`user-guide/loading` — geometry resolution, baseline/primary fallback chain
- {doc}`user-guide/reduction` — mask layers, integrator internals
- {doc}`user-guide/masks` — bundled masks, in-memory polygon dicts
- {doc}`user-guide/calibration` — re-fit geometry constants from an AGB grid scan
- {doc}`user-guide/upload` — push reduced results to a Tiled sandbox
