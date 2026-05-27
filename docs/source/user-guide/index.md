# User guide

The user guide is organized by feature.  Each page covers one slice of
the package end-to-end — what it does, how to drive it, and the things
to watch out for.

```{toctree}
:maxdepth: 2

loading
reduction
masks
calibration
beamstop-centering
upload
```

## How to read this

If you're new, start with {doc}`loading` (how raw data gets out of
Tiled) and {doc}`reduction` (how raw frames become `I(q)` /
`I(q, χ)`).  Both pages have working examples you can paste into a
notebook.

If you're operating the beamline, the operationally important pages are
{doc}`calibration` (re-fitting geometry constants from an AGB grid scan)
and {doc}`beamstop-centering` (extracting the offset table that the
collection code feeds into EPICS).

The {doc}`masks` page explains the three-layer mask architecture; the
{doc}`upload` page covers Tiled write-back for cached reduction
results.
