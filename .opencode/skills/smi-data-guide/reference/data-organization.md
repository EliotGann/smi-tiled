# SMI data organization in Tiled

Code-agnostic reference for how an SMI run is laid out in the Tiled data store, where each analysis-critical scalar is recorded, and the algorithm to find it. Field names and stream names below are the actual identifiers in Tiled.

## The data store

| Item | Value |
|---|---|
| Tiled URL | `https://tiled.nsls2.bnl.gov` |
| Catalog | `smi/migration` (legacy + post-bluesky-tiled-plugins data unified) |
| Authentication | Public anonymous read for most data; `api_key` required for proprietary or in-progress runs |
| Run primary key | `uid` (UUID string) |
| Human-friendly key | `scan_id` (integer counter) |

A "run" is one Bluesky scan (one start–stop event). It is shaped like a `BlueskyRun` node in Tiled — internally, several streams (key-value tables) plus start/stop metadata dicts.

## Run anatomy

```
<uid> = "c849dd2c-…"   (run node)
├── metadata
│   ├── start         — the scan's start document (one dict)
│   └── stop          — the scan's stop document (one dict)
├── primary           — per-frame data table
│   ├── data          — per-frame scalars and image arrays
│   └── config        — per-detector configuration
│       ├── pil2M     — SAXS detector config
│       └── pil900KW  — WAXS detector config
└── baseline          — instrument PV snapshot (start + end of scan)
```

### `start` document

Written once when the scan starts. Always present. Common keys:

| Key | Type | Purpose |
|---|---|---|
| `scan_id` | int | Human counter. |
| `uid` | str | The same UUID that keys the run. |
| `time` | float | Unix timestamp of scan start. |
| `plan_name` | str | What kind of plan (`count`, `scan`, `grid_scan`, `fly`, …). |
| `plan_args` | dict | Verbatim arguments to the plan. |
| `sample_name` | str | Free-form operator string (see below). |
| `proposal_id` | str | Proposal/proposal number. |
| `cycle` | str | Run cycle (e.g. `"2024-2"`). |
| `user_name` | str | Operator. |
| `institution` | str | Operator's institution. |
| `detectors` | list[str] | Detectors the plan asked for. |
| `motors` | list[str] | Motors the plan was asked to vary. |
| `num_points` | int | Expected frame count. |

Anything else the operator or plan added at start time (e.g. `experiment_id`, free-form notes) is also under `start`. **There is no schema enforcement** on what gets written here; key presence varies run-to-run.

### `stop` document

Written once when the scan ends. Always present (assuming the scan actually finished — partial / aborted runs may have a `stop` written by the scanner but with `exit_status="abort"` or `"fail"`). Common keys:

| Key | Type | Purpose |
|---|---|---|
| `exit_status` | `"success"` / `"abort"` / `"fail"` | Whether the scan completed cleanly. |
| `time` | float | Unix timestamp of scan end. |
| `num_events` | dict | Frame count per stream (e.g. `{"primary": 200, "baseline": 2}`). |
| `reason` | str | Free-form abort/fail reason if applicable. |

Always check `exit_status` — a `"fail"` or `"abort"` run can still have partial data, but you should know it was incomplete.

### `primary` stream

The per-frame table. One row per frame. **Wide**: it includes the image arrays in addition to scalars.

For each frame, the columns include:

| Field family | Examples | Notes |
|---|---|---|
| Image arrays | `pil2M_image`, `pil900KW_image` | Lazy. Loading the column reads tiled chunks on demand. |
| Detector motors | `pil2M_motor_x`, `pil2M_motor_y`, `pil2M_motor_z`, `pil900KW_motor_z` | mm. |
| WAXS motion | `waxs_arc` (deg), `waxs_bsx` (mm) | Linked: `Δwaxs_bsx ≈ −4.39 · Δwaxs_arc`. |
| Sample stage | `stage_th`, `piezo_th` (deg), `piezo_x`, `piezo_y`, `piezo_z` (mm or μm) | Per-frame positions. |
| SAXS beamstop (live) | `bsx`, `bsy` (mm) | Whichever beamstop is in beam. |
| Active beamstop tag | `pil2M_active_beamstop` | String `"rod"` or `"pin"`. *Critical.* |
| Beam | `energy_energy` (eV), `bpm1`, `bpm2`, … | Beam diagnostics. |
| Counters / timing | `time` (per-frame), `seq_num` | Bluesky-internal. |
| Operator-provided | `target_file_name`, custom signals | Per-frame string fields, often parsed for virtual axes. |

> **Note on order.** Tiled does not guarantee row order matches `seq_num` order, especially for runs assembled from parallel writers. If you index into primary by row position, sort by `seq_num` first.

### `primary.config[<detector>]`

A small dict written once at scan start for each detector. Examples for `pil2M`:

| Field | Purpose |
|---|---|
| `pil2M_cam_acquire_time` | Per-frame exposure time (s). |
| `pil2M_cam_num_images` | Frames per trigger. |
| `pil2M_threshold_energy` | Detector threshold (eV) — matters for absolute intensity. |
| `pil2M_beam_center_x_px`, `pil2M_beam_center_y_px` | Beam center — sometimes mirrored from baseline. |
| `pil2M_active_beamstop` | Sometimes mirrored from baseline. |

Whether each conf field is present depends on the bluesky device definition for the era of the run.

### `baseline` stream

A snapshot of every monitored EPICS PV at scan start (and usually scan end). One row at start, one at end. **This is where ~all instrument-state scalars live.**

The baseline column set is wide (often hundreds of fields). The ones that matter for analysis:

#### Beam center (SAXS)

| Field | Type | Purpose |
|---|---|---|
| `pil2M_beam_center_x_px` | float, px | Beam center column (X). Tracks `pil2M_motor_x` automatically. |
| `pil2M_beam_center_y_px` | float, px | Beam center row (Y). Tracks `pil2M_motor_y` automatically. |
| `beam_center_x`, `beam_center_y` | float, px | Older / alternate names; sometimes present. |

#### Sample-detector distance

| Field | Type | Purpose |
|---|---|---|
| `pil2M_motor_z_user_setpoint` | float, mm | Reference SDD setpoint. |
| `pil2M_motor_z` | float, mm | Reference SDD readback (also in primary per-frame). |

The "SDD in mm" needs an interpretation: it is usually the distance from a fixed instrument reference, not necessarily from the sample. Conversion to true sample-detector distance requires a calibration offset (run an AgB calibration; back-fit).

#### Active beamstop and positions

| Field | Type | Purpose |
|---|---|---|
| `pil2M_active_beamstop` | str | `"rod"` or `"pin"`. **The most critical SAXS scalar.** |
| `saxs_beamstop_x_rod`, `saxs_beamstop_y_rod` | float, mm | Rod beamstop position. |
| `saxs_beamstop_x_pin`, `saxs_beamstop_y_pin` | float, mm | Pin beamstop position. |
| `saxs_beamstop_x_rod_user_setpoint`, etc. | float, mm | Setpoint companions of the readbacks. Read first when present (the readback can be momentarily wrong during motion). |

#### Energy

| Field | Type | Purpose |
|---|---|---|
| `energy_energy` | float, eV | Photon energy. |
| `dcm_bragg`, `dcm_e` | float | Mono diagnostics; secondary. |

#### Sample stage reference

The piezo/stage motors appear here in both `_user_setpoint` and readback forms; baseline gives the start-of-scan reference position, primary gives per-frame.

#### Sample identity

| Field | Type | Purpose |
|---|---|---|
| `target_file_name` | str | Last-saved target filename for the sample, often encoding sample tokens. Sometimes per-frame in primary instead. |

## The metadata-resolution decision tree

Almost every analysis-critical scalar (energy, SDD, beam center, active beamstop, beamstop position, incident angle, sample identity) follows the same fallback chain. To resolve any one of them:

```
1. Caller-supplied override                   STOP if present.
2. baseline scalar (PV at scan start)         STOP if present.
3. primary.config[detector] scalar            STOP if present.
4. primary first-row scalar                   STOP if present.
5. start doc field                            STOP if present.
6. parse_sample_name_geometry(start.sample_name)  STOP if a token matched.
7. bundled / hard-coded default               (silent — often wrong!)
```

This ordering has two important properties:

- **baseline is preferred over primary.** Baseline is a one-shot PV snapshot; if it recorded the value, that value is more trustworthy than reading primary's first row. Primary's first row could be a transient or a value being acquired.
- **`sample_name` parsing is last.** It exists as a fallback for the era when geometry was not always recorded as a structured PV. It silently does nothing if the sample name doesn't contain a recognizable token.

### Per-scalar canonical sources

| Scalar | Best (1st) source | Common (2nd) | Last-resort | Default if absent |
|---|---|---|---|---|
| `energy` (eV) | `baseline.energy_energy` | `primary.energy_energy` | `parse(sample_name)._keV` | `None` (some readers raise; some assume 16.1 keV) |
| `sdd` (mm) | `baseline.pil2M_motor_z_user_setpoint` | `primary.pil2M_motor_z` | `parse(sample_name)._sdd<X>m` | `None` |
| `beam_center_x_px` | `baseline.pil2M_beam_center_x_px` | `conf.pil2M_beam_center_x_px` | `baseline.beam_center_x` | center of image (rarely correct) |
| `beam_center_y_px` | `baseline.pil2M_beam_center_y_px` | `conf.pil2M_beam_center_y_px` | `baseline.beam_center_y` | center of image |
| `active_beamstop` | `baseline.pil2M_active_beamstop` | `conf.pil2M_active_beamstop` | `primary.pil2M_active_beamstop` | `"rod"` (silent — see WARNING) |
| `bs_x_rod` (mm) | `baseline.saxs_beamstop_x_rod_user_setpoint` | `baseline.saxs_beamstop_x_rod` | — | `0.0` |
| `bs_x_pin` (mm) | `baseline.saxs_beamstop_x_pin_user_setpoint` | `baseline.saxs_beamstop_x_pin` | — | `0.0` |
| `incident_angle` (deg) | `baseline.<various>` (no consistent name) | `primary.<various>` | `parse(sample_name)._ai<X>` | `0.0` (silent — assumed transmission) |
| `waxs_arc` (deg, per-frame) | `primary.waxs_arc` | — | — | required (no fallback) |
| `waxs_bsx` (mm, per-frame) | `primary.waxs_bsx` | `baseline.waxs_bsx` | — | required |

> **WARNING — silent default to `"rod"`.** If `pil2M_active_beamstop` is absent from baseline, conf, AND primary, every standard reducer falls through to default `"rod"`. The pin polygon is then never applied to a pin run, producing a spurious bump in the reduced low-q curve. This is the single most common silent-correctness bug on SMI data.
>
> **Verification:** before reduction, explicitly check `baseline["pil2M_active_beamstop"]`. If empty, check `primary` and `conf`. If still empty, the run does not record the beamstop — escalate to operator notes or refuse to reduce.

> **WARNING — silent default to `incident_angle = 0`.** Grazing-incidence runs without an explicit α_i recorded fall through to 0° (transmission). Reduction then puts the specular peak in the wrong place. Always verify α_i is recorded somewhere before treating a run as GI.

### Verification recipe (no specific package required)

In any Tiled-aware Python session:

```python
from tiled.client import from_uri
client = from_uri("https://tiled.nsls2.bnl.gov")
catalog = client["smi/migration"]

# Look up by scan_id
from tiled.queries import Key
node = catalog.search(Key("scan_id") == 1130041)
uid, run = next(iter(node.items()))

# Each tier separately
start = run.metadata.get("start", {})
stop = run.metadata.get("stop", {})

# Tier 2: baseline (handles both old and new layouts)
try:
    baseline = run["baseline"]["internal"].read()  # new layout
except Exception:
    baseline = run["baseline"].read()              # legacy layout
print("baseline active_beamstop:", baseline.get("pil2M_active_beamstop"))
print("baseline beam_center:",     baseline.get("pil2M_beam_center_x_px"),
                                    baseline.get("pil2M_beam_center_y_px"))
print("baseline energy:",          baseline.get("energy_energy"))

# Tier 3: primary.config
conf = run["primary"]["config"]["pil2M"].read()
print("conf active_beamstop:",     conf.get("pil2M_active_beamstop"))

# Tier 4: primary, first frame
primary = run["primary"]["data"]
print("primary[0] active_beamstop:", primary["pil2M_active_beamstop"][0])
print("primary[0] energy:",          primary["energy_energy"][0])

# Tier 5: start doc
print("start sample_name:",        start.get("sample_name"))
print("start plan_name:",          start.get("plan_name"))
```

(Field accessors above are illustrative — the exact API depends on the `tiled` version. Use whichever syntax your environment supports; the scalar values are what matter.)

## Image arrays

### SAXS (`pil2M_image`)

Per-frame image array on the `primary` stream. Dimensions depend on plan:

| Plan | Effective dims | Shape (per frame) |
|---|---|---|
| `count` (single frame) | `(rows, cols)` | `(1679, 1475)` |
| `scan` / `rel_scan` (N frames) | `(frame, rows, cols)` | `(N, 1679, 1475)` |
| `grid_scan` (N×M) | reshaped to `(frame, rows, cols)` with `frame = N·M` | (N·M, 1679, 1475) |
| WAXS-arc-stepped | indexed by `waxs_arc` step | `(N_arc, 1679, 1475)` |

Reading is **lazy** — Tiled returns a dask-backed view; bytes are pulled only when the column is sliced and computed.

### WAXS (`pil900KW_image`)

Per-frame image array. Raw shape is `(195, 619)` — three Pilatus panels concatenated horizontally. Some readers rotate it (`np.fliplr(np.rot90(image, k=3))`) to match the SAXS coordinate convention before integration; others keep it raw.

When stepping `waxs_arc`, you get one WAXS image per arc step: `(N_arc, 195, 619)`.

> The 3-panel arc cannot be treated as a single flat detector. Tools that assume flat-panel geometry will produce wrong q-maps near the panel-panel seams.

## Sample-name token grammar

`start.sample_name` is a free-form string. Two conventions overlap:

### Geometry tokens (legacy fallback)

Used by reducers as the bottom-of-stack metadata source:

| Token pattern | Captured | Example |
|---|---|---|
| `_<digits>keV` | photon energy in keV | `_16.10keV` → 16.10 keV |
| `_<digits>eV` | photon energy in eV (converted to keV by /1000) | `_4064.00eV` → 4.064 keV |
| `_sdd<digits>m` (or `_sdd<digits>`) | sample-detector distance in m | `_sdd2.0m` → 2.0 m |
| `_wa<digits>` | WAXS arc reference angle | `_wa14.5` → 14.5° |
| `_ai<digits>` | incident angle | `_ai0.50` → 0.50° |
| `_th<digits>` | sample theta | `_th30` → 30° |

Example calibration: `EG_AGB_16.10keV_wa14.5_sdd2.0m`.

Example GI: `Lucas_film_16.10keV_sdd2.0m_ai0.20_wa14.5`.

### Per-sample tokens (virtual axes)

A more general `<label><number>[<unit>]` pattern, used to derive per-frame "virtual axes" (often prefixed `fn:`) from string fields like `target_file_name` (which can be set per-frame to capture the current sample/condition):

```
Lucas_sample2_pos1_2450.00eV_ai0.50_wa9_bpm1.995_degC100.0
```

parses as:

```
{ "sample": 2, "pos": 1, "eV": 2450.0, "ai": 0.5,
  "wa": 9, "bpm": 1.995, "degC": 100.0 }
```

Rule of thumb for the parser:
- Token is `<optional alpha prefix><signed number><optional alpha unit>`.
- The label is the prefix when present, else the unit.
- Bare numbers with no adjacent letters are ignored.
- First occurrence of a label within the string wins.

This is how a `grid_scan` over 100 sample positions becomes a `(sample, pos)` axis pair on the reduced data — the operator encodes them in the per-frame target-file string.

> **Fragility.** Both conventions are brittle: a typo in `sample_name` (e.g. `_ai0..50`) silently changes the parsed value or drops the token entirely. The proposed fix in `acquisition-upgrades.md` is to record these as structured fields directly.

## Run discovery

Common discovery patterns (any Tiled-aware client):

| Question | Query shape |
|---|---|
| Get a run by `scan_id` | `catalog.search(Key("scan_id") == 1130041)` |
| All runs of a sample | `catalog.search(Key("sample_name") == "AgB_16.10keV")` |
| Runs in a date range | `catalog.search(TimeRange(since="2024-01-01", until="2024-12-31"))` (note: `TimeRange` may be missing on some `tiled` versions; fall back to `Key("time") >= ...`) |
| Runs by plan | `catalog.search(Key("plan_name") == "grid_scan")` |
| All AgB calibrations | `catalog.search(Regex("sample_name", r"(?i)agb"))` |

The catalog is queryable; you don't have to enumerate. But each query returns a subnode whose entries are still lazy — the runs themselves don't load until you iterate.

## Disk-cache convention

Many SMI workflows cache one HDF5 file per `uid` to avoid re-fetching from Tiled:

```
<cache_dir>/<uid>.h5
├── primary/
│   ├── pil2M_motor_x      (1-D, n_frames)
│   ├── … (every primary scalar)
│   ├── pil2M_image        (3-D, n × 1679 × 1475)
│   └── pil900KW_image     (3-D, n × 195 × 619)
├── baseline/              (table of baseline scalars)
└── attrs                  (start metadata mirrored as group attrs)
```

The cache directory is conventionally:

1. `$SMI_BROWSER_CACHE_DIR` if set.
2. `$TMPDIR/smi_browser_cache/` (or `/tmp/smi_browser_cache/`).
3. `~/.cache/smi_browser/`.

The cache is best-effort: write failures (read-only filesystem, disk full) should be caught and ignored — the data is still in Tiled.

## Where to go next

- For physical interpretation of these fields (what `waxs_arc` means, what the rod vs pin do), see `beamline-physics.md`.
- For proposed schema improvements that would eliminate most of the silent-default failure modes documented above, see `acquisition-upgrades.md`.
