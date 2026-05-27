# Migration from `PyHyperScattering.SMISWAXSLoader` → `smi_tiled`

If you have existing code using `PyHyperScattering.SMISWAXSLoader` (or the SMI
parts of `PyHyperScattering.SMISWAXSIntegrator` / `smi_defaults`), this is what
changes.

## TL;DR

| Old | New |
|---|---|
| `from PyHyperScattering.SMISWAXSLoader import TiledSMISWAXSLoader` | `from smi_tiled import TiledSMISWAXSLoader` |
| `from PyHyperScattering.SMISWAXSIntegrator import reduce_smi_combined` | `from smi_tiled import reduce_smi_combined` |
| `from PyHyperScattering import smi_defaults` | `from smi_tiled import defaults` |
| `from PyHyperScattering import load; load.TiledSMISWAXSLoader(...)` | `from smi_tiled import TiledSMISWAXSLoader; TiledSMISWAXSLoader(...)` |

**Class and function names are unchanged.** Only the import paths move.

## Module renames

| Old path | New path |
|---|---|
| `PyHyperScattering.SMISWAXSLoader` | `smi_tiled.loader` |
| `PyHyperScattering.SMISWAXSIntegrator` | `smi_tiled.integrator` |
| `PyHyperScattering.smi_defaults` | `smi_tiled.defaults` |

The top-level `smi_tiled` package re-exports the common entry points
(`TiledSMISWAXSLoader`, `reduce_smi_combined`, `CombinedReductionResult`, mask
helpers, etc.) so most callers can stay on `from smi_tiled import …`.

## Data files

| Old | New |
|---|---|
| `PyHyperScattering.data.smi.masks.*` | `smi_tiled.data.masks.*` |
| `PyHyperScattering/data/smi/saxs_calibration.json` | `smi_tiled/data/saxs_calibration.json` |

If you were passing bundled-mask paths explicitly, switch to the helpers:

```python
from smi_tiled.defaults import default_saxs_mask_path, default_waxs_mask_path
saxs_path = default_saxs_mask_path()
waxs_path = default_waxs_mask_path()
```

## PyHyperScattering compatibility (optional)

`smi_tiled` doesn't depend on `pyhyperscattering`, but its DataArrays still
carry the PyHyper attrs contract (`dist, poni1, poni2, pixel1, pixel2, energy,
wavelength`). So if you have both installed, the SAXS output is consumable by
`PFGeneralIntegrator(geomethod='template_xr')`. WAXS is NOT — always use
`smi_tiled.integrator` for WAXS (the 3-panel arc isn't a flat detector).

For the xarray accessors (`da.rsoxs.*`, `da.fit.*`) that need a single
DataArray rather than the `CombinedReductionResult` Dataset, use:

```python
da = result.to_dataarray("merged_iq")     # or "merged_qchi", "per_frame_iq"
```

## smi-browser / downstream tools

If you maintain `smi-browser` or another package that imported from
`PyHyperScattering.SMISWAXSLoader`, the migration is mechanical:

```bash
# In your downstream package:
find . -name "*.py" -exec sed -i \
    -e 's/from PyHyperScattering\.SMISWAXSLoader/from smi_tiled.loader/g' \
    -e 's/from PyHyperScattering\.SMISWAXSIntegrator/from smi_tiled.integrator/g' \
    -e 's/from PyHyperScattering\.smi_defaults/from smi_tiled.defaults/g' \
    -e 's/from PyHyperScattering import smi_defaults/from smi_tiled import defaults/g' \
    {} \;
```

Add `smi-tiled[tiled]` to your dependency list. You can drop the SMI-specific
`PyHyperScattering` install (e.g. `PyHyperScattering[smi]`); keep it only if
you still use RSoXS or other PyHyper features.
