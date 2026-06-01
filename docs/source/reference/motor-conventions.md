# Motor conventions and pixel mapping

This page documents the sign conventions and the relationship between
SMI's motor positions and pixel positions on the detector.

## Coordinate systems

### Detector (pixel) frame

The raw Pilatus 2M image is stored as a 2-D array with shape
`(n_rows, n_cols) = (1679, 1475)`.

- **row** axis: vertical, increasing downward (image-array convention).
- **col** axis: horizontal, increasing rightward.

In docstrings and JSON polygons we typically use `[col, row]` order
(matching SVG / canvas convention).  In numpy index expressions
(`array[row, col]`), the order is reversed.

### Q-space

Q is computed in nm⁻¹.  `q_x`, `q_y`, `q_z` are in the lab frame:

- `+x` lab points along the increasing-col direction projected onto the
  beam-orthogonal plane.
- `+y` lab points opposite the increasing-row direction (so up in
  physical space).
- `+z` lab points along the beam (toward the detector).

`q_abs = sqrt(qx² + qy² + qz²)` is what we histogram.

`chi` is defined as `atan2(qy, qx)`, in degrees, with the conventional
sign convention `chi=0` along `+x` and `chi=90°` along `+y`.

## Motor → pixel mapping (SAXS)

The EPICS beam-center PVs (`pil2M_beam_center_x_px`, `pil2M_beam_center_y_px`)
**already track the detector translation motors** — the PV value itself moves
with `motor_x`/`motor_y` at roughly the pixel pitch (~5.8/6.0 px/mm). The loader
therefore applies **no** additional motor_x/y correction by default:

```python
# Default constants in smi_tiled.loader:
_SAXS_BEAM_COL_PX_PER_MOTOR_X_MM = 0.0        # PV already tracks motor_x → no extra correction
_SAXS_BEAM_ROW_PX_PER_MOTOR_Y_MM = 0.0        # PV already tracks motor_y → no extra correction
_SAXS_BEAM_COL_PX_PER_MOTOR_Z_MM = +0.000265  # negligible — beam well-aligned to z
_SAXS_BEAM_ROW_PX_PER_MOTOR_Z_MM = -0.000350  # negligible
```

### Why the x/y slopes are 0.0

A historical AGB grid scan (`b0f165c4-…`) regressed `beam_center ≈ a + 5.82·motor_x`
/ `a + 5.99·motor_y`. That slope is the rate at which the **PV itself** moves —
not an extra term to add. Adding `(motor − ref) · slope` on top of a PV that
already encodes the motor position double-counts: for an AGB scan at
`motor_y = −18.5 mm` it shoved the beam center ~120 px off (row 877 instead of
the true ~1003, verified by AgBh ring fitting). Set these slopes non-zero **only**
for a future detector whose beam-center PV is a static, position-independent
calibration value. The per-call overrides
(`beam_col_px_per_motor_x_mm`, `beam_row_px_per_motor_y_mm`) remain available.

## SDD → motor_z mapping (SAXS)

```python
# Active constants (from b900e711-… calibration):
_SAXS_DEFAULT_DISTANCE_DELTA_MM = -13.7       # SDD = motor_z + delta
_SAXS_SDD_DELTA_MM_PER_PIEZO_Z_UM = +0.000904  # piezo_z effect
```

Reading: `motor_z_mm = 9300` → actual SDD ≈ `9286 mm`.

The `piezo_z` correction is small (~1 μm/mm per μm of piezo, i.e.
linear with the same sign as motor_z because both move the sample
relative to the detector along the beam).

## Beamstop motors

The pin beamstop is moved by `saxs_beamstop_x_pin` (column direction on
the detector) and `saxs_beamstop_y_pin` (row direction).  Sign mapping
depends on the specific stage geometry; verify physically before
automating.

Typical (current) positions: `(-227.0, +6.8) mm`.

## WAXS-specific

### `waxs_arc`

The 900KW detector pivots on an arc — `waxs_arc_deg` is the rotation
angle.  `chi = 0°` corresponds to the centre panel directly opposite
the sample.

### `waxs_bsx`

The WAXS beamstop position.  Coupled to `waxs_arc` via the SMI
mechanical linkage:
```
waxs_bsx_ref = waxs_bsx − BSX_PER_ARC_DEG × waxs_arc_deg
BSX_PER_ARC_DEG = -4.39 mm/deg
```

The beamstop is auto-disabled in the mask when `|waxs_arc| > 15°`
because the beamstop has cleared the active area.

## Energy / wavelength

Resolved from baseline `energy_energy` (eV).  The loader stores both
`attrs["energy"]` (eV) and `attrs["wavelength"]` (**Ångstroms**) on the
output DataArray.

```{note}
The wavelength attr is in **Ångstroms**, NOT metres.  This is the
convention `smi_tiled.integrator` consumes (line `attrs['wavelength']
× 1e-10` to get metres).  This differs from PyHyperScattering's
`SST1RSoXSLoader`, which stores wavelength in metres.

PyHyperScattering's `PFGeneralIntegrator` ignores the wavelength attr
(it derives wavelength from `energy`), so this difference doesn't
affect interop — but be aware if you read the attr directly.
```

## See also

- {doc}`../user-guide/calibration` — how to re-fit these constants
- {doc}`calibration-format` — the JSON override schema
- `smi_tiled.loader.resolve_saxs_geometry` (API ref)
