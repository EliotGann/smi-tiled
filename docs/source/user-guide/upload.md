# Tiled upload

```{note}
The upload subsystem is **scaffolded** but not fully implemented at
this time.  The class shape, schema, hash, and CLI are final; the
write/read methods raise `NotImplementedError("sandbox not wired")`.
Provisioning the sandbox Tiled server is the next concrete piece of
work.  See `ROADMAP` in the README.
```

The goal of the upload subsystem is to **cache reduced data** in a
writable Tiled catalog so downstream consumers can fetch results
without re-running the reduction.

## Use case

A single reduction at SMI takes a few seconds (or longer for large
multi-frame scans).  Running 10–50 of these from a notebook or GUI is
slow enough to be annoying.  Caching the result in a Tiled sandbox
means subsequent fetches are just an HTTP read.

The same cache also lets multiple users / multiple tools share the
results of a "blessed" reduction without each re-running it.

## High-level API

```python
from smi_tiled.upload import UploadSession

session = UploadSession(
    tiled_uri="https://tiled-sandbox.nsls2.bnl.gov",
    catalog="smi-reduced",
    proposal_id="2026-1",         # optional grouping
)

# Reduce + upload in one call:
result = session.reduce_and_upload(
    uid="6e61b977-…",
    n_q=2000, n_chi=360,
    overwrite="if_stale",
)

# Or upload a pre-computed result:
from smi_tiled import reduce_smi_combined
r = reduce_smi_combined(uid="…")
session.upload(r)

# Later: fetch without re-computing
iq = session.get_merged_iq("6e61b977-…")
qchi = session.get_merged_qchi("6e61b977-…")
```

See {class}`smi_tiled.upload.UploadSession`.

## Catalog hierarchy

```
smi-reduced/
  <proposal_id>/         ← optional top-level grouping
    <uid>/               ← one node per source scan
      merged_iq          ← xr.Dataset on (q,)
      merged_qchi        ← xr.Dataset on (q, chi)
      per_frame_iq       ← xr.Dataset on (frame, q)    [if present]
```

If `proposal_id=None`, the flat layout `<catalog>/<uid>/<product>` is
used.

Each `<uid>` node carries the full reduction provenance as metadata
(see {data}`smi_tiled.upload.schema.PROVENANCE_FIELDS`), including a
**reduction hash** for staleness detection.

## Staleness via reduction_hash

The hash is computed from the canonical JSON of all reduction
parameters that affect the output (n_q, n_chi, pixel_splitting, mask
contents, dezinger threshold, calibration version, etc.).  Two calls
with semantically equivalent parameters produce identical hashes.

```python
from smi_tiled.upload import reduction_hash

h = reduction_hash({
    "uid": "…",
    "n_q": 2000, "n_chi": 360,
    "pixel_splitting": 3,
    "saxs_q_cutoff": 0.6,
    # … any reduction param
})
# → "sha256:015abd7f5cc57a2d"
```

The `overwrite` policy on {meth}`~smi_tiled.upload.UploadSession.upload`
uses this hash:

| Policy | Behavior |
|---|---|
| `"never"` | Skip if node exists; return existing URL |
| `"always"` | Overwrite unconditionally |
| `"if_stale"` (default) | Overwrite only if the stored hash differs |

## CLI

```bash
# Upload one scan
smi-upload upload --uid abc123 --tiled-uri https://… --catalog smi-reduced

# Batch (search the source catalog, reduce, upload each)
smi-upload batch --sample MyFilm --since 2026-05-01 --tiled-uri https://…

# Check whether a scan is cached and whether it's stale
smi-upload status --uid abc123 --tiled-uri https://…
```

See {mod}`smi_tiled.upload.cli`.

## Sandbox setup

The sandbox needs to be a **writable** Tiled server backed by Zarr or
HDF5.  Minimal server config:

```yaml
# tiled_config.yml (sandbox server)
trees:
  - path: /smi-reduced
    tree: tiled.adapters.zarr:ZarrGroupAdapter
    args:
      path: /data/smi-reduced-zarr
    access_control:
      access_policy: tiled.access_policies:SimpleAccessPolicy
      args:
        access_lists:
          write: ["smi-operators"]
          read: ["authenticated"]
```

Run with `tiled serve config tiled_config.yml`.

## Authentication

Set `TILED_API_KEY` in the environment, or pass `api_key=…` when
constructing the `UploadSession`.  Interactive logins are also
supported via the Tiled client's `from_uri(uri).login()` flow.

## Roadmap

The skeleton currently has working hash + schema + CLI parsing.  The
write path needs:

1. Zarr-backed write of each xr.Dataset under the `<uid>` node.
2. Provenance dict serialization as node metadata (with the
   reduction_hash and timestamp).
3. Existence + staleness check before overwrite.

The read path needs the mirror: lookup node, fetch each Dataset, return
the xr.Dataset directly (no re-typing needed since Tiled handles
xarray natively).
