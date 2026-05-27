---
sd_hide_title: true
---

# smi-tiled documentation

::::{grid} 2
:gutter: 3
:margin: 3 4 0 0

:::{grid-item}
:columns: 12

# smi-tiled

**Tiled-native loader, integrator, and uploader for the NSLS-II SMI WAXS+SAXS instrument**
(beamline 12-ID).

```{button-ref} installation
:color: primary
:expand:
:click-parent:

Get started → Installation
```

```{button-ref} quickstart
:color: secondary
:expand:
:click-parent:

Quick start
```

:::

::::

---

## What `smi-tiled` does

- **Load** raw images from a [Tiled](https://blueskyproject.io/tiled/) catalog by
  run UID, for both detectors at SMI: the Pilatus 2M (SAXS, flat-panel) and the
  Pilatus 900KW (WAXS, 3-panel folded arc).
- **Reduce** to `I(q)` and `I(q, χ)` via a pyFAI-independent integrator that
  models the multi-panel arc geometry exactly, applies dynamic per-frame masks
  (beamstop tracking, WAXS-shadow occlusion on SAXS, AgBh-anchored aperture),
  and merges SAXS+WAXS into a single count-weighted dataset.
- **Calibrate** geometry from AGB grid scans — a complete pipeline produces
  a JSON override file that ships with the package.
- **Upload** reduced results to a writable Tiled sandbox so notebooks and
  GUIs can pull them back without re-computing.

## Documentation layout

::::{grid} 2 2 3 3
:gutter: 3

:::{grid-item-card} {fas}`rocket;sd-text-info` Quick start
:link: quickstart
:link-type: ref

Reduce your first scan in under a minute.  Connect to Tiled, load by UID,
get back `I(q)` and `I(q, χ)`.
:::

:::{grid-item-card} {fas}`book;sd-text-success` User guide
:link: user-guide/index
:link-type: doc

Loading, reduction, masks, calibration, beamstop centering, upload — the
narrative-style guide to every major feature.
:::

:::{grid-item-card} {fas}`code;sd-text-warning` API reference
:link: api/index
:link-type: doc

Every public function, class, and constant.  Auto-generated from the
source code docstrings.
:::

:::{grid-item-card} {fas}`gears;sd-text-secondary` Calibration
:link: user-guide/calibration
:link-type: doc

How to re-fit the geometry constants against a fresh AGB grid scan; the
format of the calibration JSON.
:::

:::{grid-item-card} {fas}`upload;sd-text-info` Tiled upload
:link: user-guide/upload
:link-type: doc

Push reduced results back to a writable Tiled catalog.  Covers the
schema, hash-based staleness, and the `smi-upload` CLI.
:::

:::{grid-item-card} {fas}`bookmark;sd-text-muted` Reference
:link: reference/index
:link-type: doc

File formats (masks, calibration, beamstop offsets), sign and unit
conventions, and the relationship to PyHyperScattering.
:::

::::

## Indices

```{toctree}
:hidden:
:caption: Getting started

installation
quickstart
```

```{toctree}
:hidden:
:caption: User guide

user-guide/index
```

```{toctree}
:hidden:
:caption: API reference

api/index
```

```{toctree}
:hidden:
:caption: Reference

reference/index
changelog
```

* {ref}`genindex`
* {ref}`modindex`
* {ref}`search`
