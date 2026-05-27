# Masking architecture

SMI masking is built up in **three layers**, applied as a logical AND
in this order.  Each layer can be inspected, swapped, or disabled
independently.

```{mermaid}
flowchart LR
  L1["Layer 1: Fixed gaps<br/>(bundled JSON, never moves)"] --> AND
  L2["Layer 2: Beamstop polygon<br/>(motor-position-tracked)"] --> AND
  L3["Layer 3: Dynamic per-frame<br/>(WAXS shadow + AgBh aperture, SAXS only)"] --> AND
  AND[["logical AND"]] --> M["Final per-frame mask"]
```

## Layer 1 — Fixed instrument geometry

Inter-module gaps and bad-pixel regions on each detector.  These are
physically wired into the hardware and never move.  Shipped as JSON
polygons:

- **SAXS** (`src/smi_tiled/data/masks/pil2M_mask_polygons.json`) — 9 polygons
  under `static_regions`: `gap_col_1`, `gap_col_2`, `gap_row_1` …
  `gap_row_7`.
- **WAXS** (`src/smi_tiled/data/masks/900KW_mask_polygons.json`) — flat
  schema with `bad_module`, `gap_upper`, `gap_lower`, `gap_left`,
  `gap_right`.

The bundled paths are exposed through
{func}`smi_tiled.defaults.default_saxs_mask_path` and
{func}`smi_tiled.defaults.default_waxs_mask_path`.

## Layer 2 — Beamstop polygon, motor-tracked

### SAXS

Two beamstop variants under `beamstops`:

- `"rod"` — the vertical bar (default for transmission)
- `"pin"` — the small disk

Which one is active is read from baseline `pil2M_active_beamstop`.
Position is computed either:

- **Modern** (`polygon_offsets_from_beam`) — offsets `[d_col, d_row]`
  relative to the resolved beam center.  The beam center already
  tracks `pil2M_motor_x/y` and `piezo_z`, so this auto-follows.
- **Legacy** (`polygon` + `reference_mm` + `pixels_per_mm`) — absolute
  pixel coordinates, shifted by
  `(motor − reference) × pixels_per_mm`.  Kept for backwards
  compatibility with older mask files.

### WAXS

A single beamstop polygon that shifts vertically by

```
bsx_shift_px = (waxs_bsx − waxs_bsx_ref) / pixel_size_mm × 1.088
```

The reference is derived from the SMI mechanical linkage:
`waxs_bsx_ref = waxs_bsx − BSX_PER_ARC_DEG × waxs_arc`
where `BSX_PER_ARC_DEG = −4.39 mm/deg`.  Auto-disabled when
`|waxs_arc| > 15°` (the beamstop has cleared the active area).

## Layer 3 — Dynamic per-frame (SAXS only)

Two extra masks computed during integration:

### WAXS-shadow mask

The WAXS detector physically blocks part of the SAXS detector's view.
The boundary column moves with `waxs_arc`.  See
`smi_tiled.integrator._make_waxs_shadow_mask`.

### AgBh aperture mask

Q-cutoff anchored to a chosen AgBh ring order (default 5) at the
current SDD.  Auto-disables when the cutoff falls beyond the detector
edge (long SDD, so the WAXS detector doesn't physically occlude
anything visible).  See `smi_tiled.integrator._make_aperture_mask`.

## Override knobs

Layers 2 and 3 are independently togglable via
{func}`smi_tiled.reduce_smi_combined`:

```python
reduce_smi_combined(
    uid="…",
    saxs_mask=my_dict_or_path_or_None,   # Layer 1+2 source
    waxs_mask=my_dict_or_path_or_None,   # Layer 1+2 source (WAXS)
    saxs_kwargs={
        "dynamic_saxs_kwargs": {
            "aperture":    {"enabled": False},   # Layer 3a off
            "waxs_shadow": {"enabled": False},   # Layer 3b off
        },
    },
    saxs_q_cutoff=0.6,           # force aperture cutoff (nm⁻¹)
    saxs_agbh_ring_order=8,      # change anchor ring
)
```

## Mask inputs: path or dict

All public mask entry points accept **either** a JSON file path **or**
an already-parsed polygon dict with the same schema.  The dict form
lets notebook users compose and edit masks in memory:

```python
import json
from smi_tiled.defaults import default_saxs_mask_path

spec = json.load(open(default_saxs_mask_path()))
spec["static_regions"]["my_extra_blob"] = [
    [100, 200], [150, 200], [150, 250], [100, 250],
]
result = reduce_smi_combined(uid="…", saxs_mask=spec)
```

Functions taking either form:

- {func}`~smi_tiled.make_saxs_mask_from_spec` /
  {func}`~smi_tiled.make_saxs_mask_from_dict`
- {func}`~smi_tiled.make_waxs_mask_callable` /
  {func}`~smi_tiled.make_waxs_mask_callable_from_dict`
- {func}`~smi_tiled.mask_for_frame`
  (via `mask_path=` accepting dict)
- {func}`~smi_tiled.reduce_smi_combined`
  (`saxs_mask=` / `waxs_mask=` kwargs)
- {func}`~smi_tiled.reduce_smi_gi` (`waxs_mask=` kwarg)

## Mask schema

See {doc}`../reference/mask-schema` for the full JSON schema, including
the legacy flat WAXS format and the modern nested format.
