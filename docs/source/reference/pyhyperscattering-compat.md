# PyHyperScattering compatibility

`smi-tiled` was split out from
[PyHyperScattering](https://github.com/EliotGann/PyHyperScattering) in
May 2026.  It has **no runtime dependency** on PyHyperScattering, but
its output is designed to remain consumable by PyHyper's tooling for
users who have both packages installed.

## What's shared

- The `xr.DataArray` geometry-attrs contract:
  `dist, poni1, poni2, rot1, rot2, rot3, pixel1, pixel2, energy,
  wavelength`.
- The "raw frame" dim convention: `pix_y`, `pix_x` for the spatial
  axes, with scan axes (e.g. `frame`, `waxs_arc`) prepended.
- Bundled mask polygon JSON format.

## What's not shared

- **No subclass of `FileLoader`** — `smi_tiled.TiledSMISWAXSLoader` is
  Tiled-native, indexed by run UID rather than file path.
- **No pyFAI usage** — `smi_tiled.integrator` is self-contained
  (numpy + scipy + xarray + skimage).
- **WAXS geometry** — the 3-panel folded arc detector requires
  `MultiPanelArcDetector`; `PyHyperScattering.PFGeneralIntegrator`
  can't handle it (it assumes a single flat panel).

## Interop recipes

### SAXS frames through PFGeneralIntegrator

```python
import smi_tiled
import PyHyperScattering as phs

saxs = smi_tiled.TiledSMISWAXSLoader().loadSingleImage(
    uid="…", detector="saxs",
)
# SAXS DataArray has the full geometry attrs:
integrator = phs.PFGeneralIntegrator(
    geomethod='template_xr',
    template_xr=saxs,
)
result = integrator.integrateImageStack(saxs)
```

### Using PyHyper's RSoXS / Fitting accessors

Both accessors operate on `xr.DataArray`.  The smi-tiled
`reduce_smi_combined` output is `CombinedReductionResult` (a dataclass
with `xr.Dataset` members), so extract a DataArray view first:

```python
import smi_tiled

result = smi_tiled.reduce_smi_combined(uid="…")
da = result.to_dataarray("merged_iq")     # or "merged_qchi", "per_frame_iq"

# Now PyHyper accessors work:
da.rsoxs.slice_q(0.5, q_width=0.05)
da.fit.apply(my_fit_func)
```

See {meth}`~smi_tiled.CombinedReductionResult.to_dataarray`.

### Reading the wavelength attr

```{warning}
The wavelength attr is in **Ångstroms** in smi-tiled, and in **metres**
in PyHyperScattering's `SST1RSoXSLoader`.  Code that reads the attr
directly across both packages will see a 10¹⁰ scale mismatch.
```

PyHyper's `PFGeneralIntegrator` doesn't read the wavelength attr (it
derives wavelength from `energy`), so this inconsistency doesn't break
interop in practice — but be aware if you write code that reads
`da.attrs['wavelength']`.

## Migrating from PyHyperScattering.SMISWAXSLoader

See [MIGRATION.md](https://github.com/EliotGann/smi-tiled/blob/main/MIGRATION.md)
in the repo for a sed-based recipe.  Summary:

| Old | New |
|---|---|
| `from PyHyperScattering.SMISWAXSLoader import …` | `from smi_tiled.loader import …` |
| `from PyHyperScattering.SMISWAXSIntegrator import …` | `from smi_tiled.integrator import …` |
| `from PyHyperScattering.smi_defaults import …` | `from smi_tiled.defaults import …` |
| `from PyHyperScattering import smi_defaults` | `from smi_tiled import defaults` |
| `from PyHyperScattering import load; load.TiledSMISWAXSLoader(...)` | `from smi_tiled import TiledSMISWAXSLoader; TiledSMISWAXSLoader(...)` |

Class and function names are unchanged.

## Why we split

The SMI subsystem grew to ~6500 lines with zero functional dependencies
on PyHyperScattering core.  Splitting gave us:

- Independent release cadence (recalibrations ship without waiting for
  a PyHyper release).
- A natural home for the Tiled upload subsystem.
- No need to ask a PyHyper maintainer to absorb beamline-specific
  code.

PyHyperScattering still hosts `SMIRSoXSLoader.py` — the legacy
**file-based** RSoXS loader that *is* a `FileLoader` subclass.  That
belongs there.
