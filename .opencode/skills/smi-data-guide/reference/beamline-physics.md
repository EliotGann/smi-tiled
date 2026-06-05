# SMI beamline physics

Code-agnostic description of the SMI beamline (NSLS-II 12-ID) instrumentation, geometry, and the physical relationships an analysis must respect. Field names below are the literal names used in the Tiled data store.

## Beamline at a glance

| Property | Value |
|---|---|
| Facility | NSLS-II (Brookhaven National Lab) |
| Sector | 12-ID |
| Name | SMI — Soft Matter Interfaces |
| Source | Undulator |
| Operating energy range | ≈ 14 – 17.5 keV (typical), set by monochromator |
| Operating modes | Transmission SAXS+WAXS, Grazing-Incidence SAXS+WAXS |
| Detectors | Pilatus 2M (SAXS, flat panel) + Pilatus 900KW (WAXS, 3-panel folded arc) |
| Sample environment | Multi-axis sample stage with piezo-fine motion, plus optional in-situ environments (heating, humidity, etc.) |

## Energy ↔ wavelength

Photons are characterized by either energy `E` (eV) or wavelength `λ` (m). The conversion is exact:

```
E · λ = h · c = 1.23984193 keV · nm
        = 1239.84193 eV · nm
        = 1.23984193 × 10⁻⁶ eV · m
```

So at 16.10 keV: `λ = 1.23984193 / 16.10 = 0.07701 nm = 0.7701 Å = 7.701 × 10⁻¹¹ m`.

Energy in SMI metadata appears under the field name `energy_energy` (eV). Some legacy / encoded sources use keV. The `keV` vs `eV` distinction is a recurring units gotcha — see `acquisition-upgrades.md`.

The momentum-transfer magnitude `q` is

```
q = (4π / λ) · sin(θ)
```

where `2θ` is the scattering angle from the direct beam to the pixel. For a flat detector at distance `D` from the sample, with the pixel at `r` mm from the beam center:

```
2θ = arctan(r / D)        (small-angle: 2θ ≈ r / D)
```

## Geometry: transmission vs grazing-incidence

Two operating modes; they differ in the sample mount and how `q` decomposes.

### Transmission

- The beam passes through a thin sample (a film, a solution capillary, a window).
- Scattering is azimuthally isotropic for an isotropic sample.
- `α_i = 0` by convention (the beam is normal to the sample mounting reference, but for analysis purposes the incident-angle term drops).
- Reduce to `I(q)` (1-D) or `I(q, χ)` (2-D, with χ the azimuthal angle on the detector).

### Grazing-incidence (GI-SAXS / GI-WAXS)

- The beam is glancing off the sample surface at a small angle `α_i` (degrees, typically 0.1°–0.5° for thin films).
- The sample is at near-grazing incidence; the diffuse scattering from the sample surface and structures is what's measured.
- `q` decomposes into:
  - `q_∥` (in-plane component, parallel to the surface)
  - `q_⊥` (out-of-plane component, perpendicular to the surface)
- `α_i` matters: it shifts where the specular peak is and where the Yoneda line sits. **Without α_i, GI data cannot be reduced correctly.**
- α_i is encoded in metadata in any of: a baseline scalar, a per-frame primary value, or `_ai<digits>` token in `sample_name`.

In an SMI run, `start.plan_name` rarely tells you which mode it was; you infer it from the sample mount (visible only in operator notes) or from whether `α_i ≠ 0` was recorded. A modern run would carry this explicitly — see `acquisition-upgrades.md`.

## Detectors

### Pilatus 2M (SAXS) — `pil2M`

A Dectris Pilatus3 2M detector positioned downstream of the sample to capture the small-angle scattering.

| Property | Value |
|---|---|
| Image field name | `pil2M_image` |
| Image shape (rows × cols) | 1679 × 1475 pixels |
| Pixel size | 0.172 mm (square) |
| Active area | ~ 289 × 254 mm |
| Layout | Single flat panel |
| Sample-detector distance (SDD) | ~1.7 – 9.3 m, set per experiment via `pil2M_motor_z` |
| Beam center | Per-run baseline scalars `pil2M_beam_center_x_px` (col) and `pil2M_beam_center_y_px` (row); typically near image center but not exactly |

The Pilatus has internal gaps between sub-modules — these appear as straight horizontal/vertical bands of dead pixels and are static across all SMI scans. Mask polygons covering them ship with any reduction package.

### Pilatus 900KW (WAXS) — `pil900KW`

A custom 3-panel "folded arc" Pilatus detector, mounted close to the sample on an arc-arm (`waxs_arc`) to capture wide-angle scattering.

| Property | Value |
|---|---|
| Image field name | `pil900KW_image` |
| Raw shape (rows × cols) | 195 × 619 (three panels concatenated) |
| Per-panel shape | 195 × 195 (with small inter-panel gaps, total width 619 across the 3) |
| Pixel size | 0.172 mm (same as SAXS) |
| Layout | 3 flat panels arranged in an arc, panel angles ≈ −7°, 0°, +7° relative to the arc axis |
| SDD (calibrated) | ≈ 0.273 m (much closer than SAXS; this is what makes it WAXS) |
| Coverage | An angular slice of `(2θ, χ)`-space; the arc swing covers a wider chi range by stepping `waxs_arc` |

The 3-panel arc is **not** a single flat panel. Any tool that assumes a flat detector (e.g. naïve pyFAI usage) will produce subtly wrong q-maps near the panel-panel seams. A multi-panel-aware reducer is required.

### WAXS arc and beamstop motion

The WAXS detector is mounted on a rotational arc:
- `waxs_arc` — the per-frame arc angle in degrees.
- `waxs_bsx` — the WAXS beamstop X position in mm.

There is a fixed mechanical linkage between them:

```
waxs_bsx ≈ waxs_bsx_0 + (-4.39 mm/deg) · waxs_arc
```

i.e. `waxs_bsx` moves **−4.39 mm per degree** of `waxs_arc`. The beamstop is on a separate motion axis but is geometrically tied to the detector arc so it stays in front of the direct beam as the arc rotates.

Reducers that mask the WAXS beamstop shadow per-frame must account for this; the shadow position depends on `waxs_arc` and `waxs_bsx` together.

## Beam stops (SAXS)

A beam stop is a small absorber that blocks the direct (un-scattered) beam from reaching the detector. Without it, the direct beam saturates a region of the image and bleeds into surrounding pixels (parasitic scatter, charge-spreading, after-pulses).

SMI's SAXS detector has **two switchable beam stops**:

### Rod beam stop

- A vertical rod, blocking a thin strip extending **upward** from the beam center.
- Used for normal transmission SAXS: the rod blocks the direct beam plus the specular reflection from below, while leaving most of the chi range usable.
- In the bundled mask polygons, the rod region is a tall, narrow rectangle aligned with the column axis.

### Pin beam stop

- A small disk near the beam center, offset by approximately **+5 pixels in the column direction** from the beam center, with radius ≈ 22 pixels (≈ 3.8 mm).
- Used when the rod's coverage is insufficient (e.g. very low q, or specific geometries where the specular reflection lands inside the rod's mask).
- Critical: the pin **does not cover the same region as the rod**. A reduction that thinks the rod is in but the pin is actually in will leave the pin's small disk un-masked, producing a spurious bump in the low-q / chi ≈ +90° region of the reduced curve.

### Which is in beam?

The state is recorded in the field `pil2M_active_beamstop`, a string equal to either `"rod"` or `"pin"`. This is the single most consequential string in SMI metadata for SAXS reduction correctness.

> **WARNING.** When `pil2M_active_beamstop` is missing from baseline / conf / primary, every standard reducer defaults to `"rod"` silently. If the actual measurement used the pin and this field is absent, the pin disk is never masked, and the reduced `I(q)` has an unphysical bump at low q. See `acquisition-upgrades.md` for the proposed fix.

The per-axis position of each beam stop is also recorded:
- `saxs_beamstop_x_rod`, `saxs_beamstop_y_rod` — rod position (mm).
- `saxs_beamstop_x_pin`, `saxs_beamstop_y_pin` — pin position (mm).
- Each has a `_user_setpoint` companion (the setpoint, vs the readback).

These positions feed into per-frame mask construction: the beam stop's polygon is shifted into pixel space using `(x, y)_mm` and the beam center.

A separate "bsx" / "bsy" pair refers to whichever beamstop is currently in beam (live position):
- `bsx` — SAXS beamstop X (mm), live.
- `bsy` — SAXS beamstop Y (mm), live.

## Sample stage motors

A typical SMI sample is mounted on a stack:

| Motor | Type | Typical range | Purpose |
|---|---|---|---|
| `stage_th` | Rotational | full | Coarse sample tilt / theta |
| `piezo_th` | Rotational, fine | small | Fine theta (alignment) |
| `piezo_x`, `piezo_y` | Translational, fine (μm-scale) | mm-range | Fine X/Y translation, used as scan axes for grid scans |
| `piezo_z` | Translational, fine (μm-scale) | mm-range | Fine Z (sample-to-detector along beam, in addition to detector translation) |

`piezo_z` is in **micrometers** in some versions of metadata and **millimeters** in others; check the units explicitly. Conversion depends on the era of the run.

## Detector motors

### SAXS

- `pil2M_motor_x`, `pil2M_motor_y` — detector translation in the plane perpendicular to the beam (mm). These shift the beam center on the detector; the EPICS PVs `pil2M_beam_center_x_px` and `pil2M_beam_center_y_px` typically already track these motions, so the beam center reported in baseline is the live, motion-corrected value.
- `pil2M_motor_z` — detector translation along the beam (mm). This sets the SAXS sample-detector distance.

### WAXS

- `pil900KW_motor_z` — WAXS detector translation along the beam (mm), fine adjustment.
- `waxs_arc`, `waxs_bsx` — see "WAXS arc and beamstop motion" above.

## Calibration

### Silver behenate (AgB) — the SMI calibrant

The community-standard SAXS calibrant is **silver behenate** (AgC₂₂H₄₃O₂), a layered crystal with a known repeat distance:

```
d (AgB) = 5.838 nm
```

By Bragg's law, the first-order ring is at:

```
q₁ = 2π / d = 2π / 5.838 nm = 1.076 nm⁻¹ = 0.1076 Å⁻¹
```

with subsequent orders at integer multiples (`q_n = n · q₁`).

A run with `sample_name` containing `AGB` or `AgB` is almost always a calibration scan. Reducing it should yield a sharp first ring at exactly `q₁`. If the observed first-ring `q` is off by more than ~0.5%, the SDD or beam center is mis-recorded — either re-calibrate (use a known SDD to back-fit the beam center, or use a known beam center to back-fit the SDD) or re-locate the calibration scan.

A typical calibration `sample_name` looks like:

```
EG_AGB_16.10keV_wa14.5_sdd2.0m
```

— operator initials, sample, energy, WAXS arc reference angle, and SDD encoded as tokens.

### Other calibrants

For specific experiments, alternative standards may be used (silicon, glassy carbon for absolute intensity). These are not in the bundled mask/calibration JSONs and require manual handling.

## Geometry summary

The minimum geometry needed to reduce a SAXS frame:

```
λ                       (from energy)
SDD                     (= |pil2M_motor_z| − reference offset)
beam_center_row_px      (= pil2M_beam_center_y_px)
beam_center_col_px      (= pil2M_beam_center_x_px)
active_beamstop         ("rod" or "pin")
beamstop_pos_mm         ({rod:{x,y}, pin:{x,y}}, mm)
```

For WAXS, additionally:

```
waxs_arc                (per-frame, deg)
waxs_bsx                (per-frame, mm — but linked to waxs_arc via -4.39 mm/deg)
WAXS panel calibration  (which 3-panel layout was used; from a calibration scan)
```

For grazing incidence, additionally:

```
α_i                     (incident angle, deg)
```

For the data-side of *where each of these lives*, see `data-organization.md`.
