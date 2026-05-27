# smi-tiled

Tiled-native loader, integrator, and uploader for the NSLS-II **SMI** (beamline
12-ID) WAXS + SAXS instrument.

## What it does

- **Load** raw images from a [Tiled](https://blueskyproject.io/tiled/) catalog
  by run UID, for both detectors:
  - Pilatus 2M (SAXS, flat-panel)
  - Pilatus 900KW (WAXS, 3-panel folded arc)
- **Reduce** to `I(q)` and `I(q, χ)` via a pyFAI-independent integrator with
  exact multi-panel WAXS geometry, dynamic per-frame masks (beamstop tracking,
  WAXS-shadow occlusion on SAXS, AgBh-anchored aperture), and weighted SAXS+WAXS
  merging.
- **Upload** reduced results back into a writable Tiled sandbox so subsequent
  consumers can fetch reduced data without re-computing.

## Origin

This package was split out from
[`PyHyperScattering`](https://github.com/EliotGann/PyHyperScattering) in May 2026
once it became clear that the SMI-specific code (loader, integrator, masks,
calibration, beamstop centering) was self-contained and didn't share runtime
dependencies with the rest of PyHyperScattering. See [MIGRATION.md](MIGRATION.md)
for users moving from `PyHyperScattering.SMISWAXSLoader`.

`smi-tiled` has **no runtime dependency on `pyhyperscattering`**. The
`xr.DataArray` attrs convention (`dist, poni1, poni2, pixel1, pixel2, energy,
wavelength`) is shared so the output can be consumed by
`PyHyperScattering.PFGeneralIntegrator(geomethod='template_xr')` *if* the user
also installs PyHyperScattering — but neither package depends on the other.

## Installation

### With pixi (recommended)

```bash
git clone https://github.com/EliotGann/smi-tiled.git
cd smi-tiled
pixi install         # default env: tiled + upload + dev
pixi run test        # run the test suite
pixi shell           # drop into the env
```

### With pip

```bash
pip install smi-tiled[tiled,upload]
# or for development:
pip install -e ".[tiled,upload,dev]"
```

### Optional feature groups

| Group | What you get |
|---|---|
| `tiled` | `tiled[client]` + `bluesky-tiled-plugins` — needed for the loader |
| `upload` | `tiled[client]` + `httpx` + `zarr` — needed for the uploader |
| `dev` | `pytest`, `coverage`, `flake8` |
| `docs` | `sphinx`, `pydata-sphinx-theme` |
| `all` | everything |

## Quick start

```python
from smi_tiled import reduce_smi_combined

result = reduce_smi_combined(
    uid="6e61b977-…",
    tiled_uri="https://tiled.nsls2.bnl.gov",
    catalog="smi/migration",
    n_q=2000, n_chi=360,
)
# Result is a CombinedReductionResult; merged_iq and merged_qchi are
# xr.Datasets, per_frame_iq has the per-scan-step I(q).
print(result.merged_iq)

# For PyHyperScattering accessor compatibility:
da = result.to_dataarray("merged_iq")
da.rsoxs.slice_q(0.5)   # if PyHyperScattering is also installed
```

## Mask architecture

See `src/smi_tiled/integrator.py` module docstring for the full description.
Three layers, ANDed together:

1. **Fixed instrument geometry** (gaps, bad modules) — bundled JSON masks.
2. **Beamstop polygon** (rod or pin), motor-position-tracking.
3. **Dynamic per-frame occlusion** (WAXS shadow + AgBh aperture, SAXS only).

All public mask entry points accept either a JSON file path or an in-memory
polygon dict.

## Calibration

The geometry calibration constants
(`_SAXS_DEFAULT_DISTANCE_DELTA_MM`, `_SAXS_BEAM_*_PX_PER_MOTOR_*_MM`, …) live as
module-level defaults plus an optional override file at
`src/smi_tiled/data/saxs_calibration.json`. Re-calibrating against a fresh AGB
grid scan:

```bash
pixi run python scripts/calibrate_smi_z_scan.py <UID>
# emits /tmp/saxs_calibration.json (constants block applied at import time)
# and  /tmp/smi_beamstop_offsets.json (spline + raw table for collection code)
```

Beamstop centering offsets are emitted as both per-point measurements and a
canonical scipy B-spline (knots + coefficients + 200-pt dense grid) so the
collection code can reconstruct without `scipy` if needed.

## License

NIST public-domain (see [LICENSE](LICENSE)).
