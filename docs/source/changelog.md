# Changelog

This page tracks user-visible changes to `smi-tiled`.  Versions follow
[Semantic Versioning](https://semver.org).  Until 1.0, the API may
change in minor releases.

## Unreleased

### Added

- Initial release: package split from PyHyperScattering 0.x.
- `TiledSMISWAXSLoader` — Tiled-native loader for SMI Pilatus 2M (SAXS)
  and 900KW (WAXS).
- `reduce_smi_combined` / `reduce_smi_gi` — full reduction pipelines.
- `CombinedReductionResult` / `GIReductionResult` — result dataclasses
  with `.to_dataarray()` helpers for PyHyperScattering accessor interop.
- Three-layer masking architecture (fixed gaps, motor-tracked beamstop,
  dynamic per-frame).
- All public mask functions accept either JSON paths **or** in-memory
  polygon dicts.
- `saxs_calibration.json` JSON override mechanism applied at import time.
- `smi_beamstop_offsets.json` with spline-interpolated beamstop centering
  offsets vs `motor_z` for the collection code.
- `scripts/calibrate_smi_saxs.py` — fits BC vs (motor_x, motor_y, piezo_z).
- `scripts/calibrate_smi_z_scan.py` — fits SDD vs motor_z + beamstop offsets.
- `smi_tiled.upload` — Tiled write-back scaffolding (`UploadSession`,
  `reduction_hash`, `smi-upload` CLI).  Implementations are stubbed
  pending sandbox provisioning.
- HDF5 disk-cache fallback for tiled HTTP-500 issues.
- Chunked + parallel frame reads with retry for large multi-frame scans.
- 110 CI-safe tests; no live tiled needed for the test suite.

### Breaking changes vs PyHyperScattering.SMISWAXSLoader

- Module paths changed: see {doc}`reference/pyhyperscattering-compat`
  or `MIGRATION.md`.
- The bundled mask path moved from `PyHyperScattering/data/smi/masks/`
  to `smi_tiled/data/masks/`.
- The calibration JSON path moved from
  `PyHyperScattering/data/smi/saxs_calibration.json` to
  `smi_tiled/data/saxs_calibration.json`.
- Class and function names are **unchanged** — only import paths shift.
