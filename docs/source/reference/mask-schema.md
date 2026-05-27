# Mask JSON schema

Two schemas are accepted by every mask entry point in `smi-tiled`.

## Nested (preferred, used by the bundled SAXS mask)

```json
{
  "image_shape": [1679, 1475],
  "static_regions": {
    "gap_col_1": [[col, row], [col, row], ...],
    "gap_col_2": [...],
    "gap_row_1": [...],
    ...
  },
  "beamstops": {
    "rod": {
      "polygon_offsets_from_beam": [[d_col, d_row], ...],
      "x_motor_key": "saxs_beamstop_x_rod",
      "y_motor_key": "saxs_beamstop_y_rod"
    },
    "pin": {
      "polygon_offsets_from_beam": [[d_col, d_row], ...],
      "x_motor_key": "saxs_beamstop_x_pin",
      "y_motor_key": "saxs_beamstop_y_pin"
    }
  }
}
```

### `image_shape`

`[rows, cols]` — informational; the loader uses the actual raw frame
shape from the detector.

### `static_regions`

A dict mapping a polygon name → polygon (list of `[col, row]` pairs).
Polygons are 2-D in pixel coordinates.  Order: first coordinate is the
**column** (x), second is the **row** (y), to match the
`(x, y)`-style convention used in vector-drawing tools.

### `beamstops`

A dict keyed by beamstop variant (`"rod"`, `"pin"`).  Each entry has:

#### Modern form: `polygon_offsets_from_beam`

```json
{
  "polygon_offsets_from_beam": [[d_col, d_row], ...]
}
```

Offsets are added to the **resolved beam center** (which already
tracks `pil2M_motor_x/y` and `piezo_z`).  This is the preferred form
because it auto-follows the beam.

A list of polygons (each itself a list of offsets) is also accepted —
useful when one beamstop variant has multiple polygon regions:

```json
{
  "polygon_offsets_from_beam": [
    [[d_col, d_row], ...],
    [[d_col, d_row], ...]
  ]
}
```

#### Legacy form: absolute `polygon` + motor reference

```json
{
  "polygon": [[col, row], ...],
  "reference_mm": {"x": 0.0, "y": 0.0},
  "pixels_per_mm": {"x": 5.8140, "y": 5.8140}
}
```

The polygon's absolute pixel coords are shifted by
`(motor − reference_mm) × pixels_per_mm` for the active beamstop
motor positions.  Retained for backwards compatibility with older mask
files.

## Flat (legacy, used by the bundled WAXS mask)

```json
{
  "beamstop": [[col, row], ...],
  "bad_module": [[col, row], ...],
  "gap_upper":  [[col, row], ...],
  "gap_lower":  [[col, row], ...],
  "gap_left":   [[col, row], ...],
  "gap_right":  [[col, row], ...]
}
```

Top-level keys are polygon names.  The `"beamstop"` key (if present)
is the moving region; everything else is treated as a static gap.
{func}`~smi_tiled.make_waxs_mask_callable_from_dict` accepts both this
flat schema and the nested one above.

## Validation

The package doesn't enforce a formal JSON Schema, but you can sanity-check
a mask in Python:

```python
import json
from smi_tiled.defaults import default_saxs_mask_path, load_mask_polygons

# load_mask_polygons normalizes either schema into the nested form:
norm = load_mask_polygons(default_saxs_mask_path())
print(list(norm["static_regions"]))    # gap names
print(list(norm["beamstops"]))         # variant names
```

See {func}`~smi_tiled.defaults.load_mask_polygons`,
{func}`~smi_tiled.defaults.save_mask_polygons`.

## In-memory polygon dicts

All public mask functions accept a parsed dict directly — no JSON
file required.  See {doc}`../user-guide/masks` for examples.

## See also

- {doc}`../user-guide/masks` — three-layer mask architecture
- {func}`~smi_tiled.make_saxs_mask_from_dict`
- {func}`~smi_tiled.make_waxs_mask_callable_from_dict`
