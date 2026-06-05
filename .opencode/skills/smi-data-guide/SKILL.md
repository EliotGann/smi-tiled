---
name: smi-data-guide
description: Use when reading, exploring, or interpreting SMI beamline (NSLS-II 12-ID) scan data stored in Tiled, independently of any analysis package. Covers what an SMI scan physically is, what an SMI "run" looks like in the Tiled data store, where each analysis-critical scalar (energy, sample-detector distance, beam center, active beamstop, incident angle, sample identity) is recorded, and the decision algorithm for locating each one. Trigger on terms like "SMI data", "12-ID data", "what's in this scan", "what does this run contain", "sample_name parsing", "baseline vs primary", "active_beamstop", "Pilatus 2M / 900KW raw data", "AgB calibration scan", "what do I need to reduce this", or "incident angle of this frame" — when the question is about the data itself, not a particular Python package. Sibling reference `acquisition-upgrades.md` covers proposed metadata-schema improvements.
---

# SMI data guide (code-agnostic)

This skill describes SMI beamline data as it exists in the Tiled data store, independent of any reduction package. It is the conceptual companion to `smi-tiled` (which captures one specific reader/integrator). Use this skill when the question is about the data itself: what was measured, where the metadata lives, what is needed to interpret a frame, and where the schema is fragile.

## Trigger

Use this skill when:
- A user asks what is in an SMI run, what fields a scan has, or how to interpret SMI scan data without a specific tool in mind.
- A user is exploring SMI data in Tiled directly (e.g. through `tiled` CLI, a notebook, or a third-party reducer) and needs to know where energy / SDD / beam center / beamstop / incident angle live.
- A user mentions: `tiled.nsls2.bnl.gov`, `smi/migration`, `pil2M_image`, `pil900KW_image`, `pil2M_active_beamstop`, `waxs_arc`, `waxs_bsx`, `bsx`/`bsy`, `target_file_name`, `sample_name`, `scan_id`, `start doc`, `baseline stream`, `primary stream`, `silver behenate` / `AgB`, `grazing incidence` / `transmission`.
- A user is designing or critiquing the SMI metadata schema, or proposing changes to what gets recorded per scan.

Do NOT use this skill when:
- The question is specifically about the `smi-tiled` package internals (function names, code paths, masking helpers). Use the `smi-tiled` skill instead.
- The data source is not SMI (other beamlines have different conventions).

The two skills are complementary: `smi-tiled` says *how the package handles SMI data*; `smi-data-guide` says *what the data is and where it lives*.

## What SMI is, in one paragraph

SMI is the **Soft Matter Interfaces** beamline at NSLS-II (sector 12-ID), a hard-x-ray beamline (typical operating range ≈ 14–17.5 keV) instrumented for combined Small-Angle and Wide-Angle X-ray Scattering (SAXS+WAXS) on soft and hybrid materials. A scan produces simultaneous frames on two area detectors (a flat-panel SAXS detector and a 3-panel folded-arc WAXS detector), plus per-frame motor positions, plus a one-shot snapshot of every monitored EPICS PV at start of scan ("baseline"). All of this is stored as a `BlueskyRun` in Tiled.

For physical detail (detector geometry, beam stops, motors, calibration), see `reference/beamline-physics.md`.

## What an SMI scan looks like in Tiled

The data store: `https://tiled.nsls2.bnl.gov`, catalog `smi/migration`. Each scan is one **run**, identified by:
- `uid` — UUID, primary key.
- `scan_id` — integer counter, the human-friendly key.

Every run has the standard Bluesky stream layout:

| Stream | What's there |
|---|---|
| `start` (metadata) | One dict written when the scan starts: `scan_id`, `sample_name`, `plan_name`, `proposal_id`, `cycle`, `user_name`, plus user-supplied keys. |
| `stop` (metadata) | One dict written when the scan ends: `exit_status` (`"success"`/`"abort"`/`"fail"`), `time`. |
| `primary` | The per-frame table. Contains the detector image arrays (`pil2M_image`, `pil900KW_image`) and every motor / signal that was either varied (the scan axis) or read once per frame. One row per frame. |
| `primary.config[<detector>]` | Per-detector configuration scalars captured once at scan start (e.g. exposure time, gain mode). |
| `baseline` | A snapshot of every monitored PV taken at scan start (and usually scan end). One row per snapshot. This is where almost all instrument-state scalars live. |

For the full anatomy of these streams and the field-name catalog, see `reference/data-organization.md`.

## The 5 things you need to interpret any SMI frame

Independent of the reduction tool you use, every SAXS/WAXS frame requires the same 5 scalars to convert pixels into physics:

| # | Scalar | Why | Typical units |
|---|---|---|---|
| 1 | **Photon energy** (or wavelength) | Sets `q = (4π/λ) sin(θ)`. Without it, every reduced curve is wrong by the energy ratio. | eV (sometimes keV — see units gotcha) |
| 2 | **Sample-detector distance (SDD)** | Sets the pixel-to-angle map. | mm or m |
| 3 | **Beam center on the detector** (row, col in pixels) | Origin for `q` and `chi`. A 1-pixel error costs ~1% in q at the edge. | pixels |
| 4 | **Which beam stop is in beam** (`rod` or `pin`, for SAXS) | Determines which polygon must be masked. The wrong choice produces phantom data near `q_min`. | string token |
| 5 | **Incident angle** `α_i` (only for grazing-incidence runs) | Splits scattering into specular vs Yoneda components. Transmission scans have α_i = 0 by convention. | degrees |

For SMI, you also typically need the **WAXS arc angle** (`waxs_arc`, deg) per frame to interpret the WAXS panel positions, and **`waxs_bsx`** to mask the WAXS beamstop shadow.

The hard part is not computing on these scalars; it is **finding them**. The schema does not give a single canonical location for any of them.

## The metadata-resolution decision tree (must-internalize)

For each scalar above, the same 6-tier search applies. Walk top-to-bottom, take the first non-empty answer:

```
1. Caller-supplied override            (explicit kwarg / config)
2. baseline stream                     (instrument state at scan start)
3. primary.config[detector]            (one-shot detector config)
4. primary stream                      (per-frame value, take first row)
5. start doc                           (scan_id, sample_name, plan_name)
6. parse_sample_name_geometry()        (regex-extract from sample_name)
7. bundled / hard-coded default        (silent fallback)
```

The two pathologies built into this chain:

> **WARNING — silent fallback to default.** If no source records a scalar, tiers 1–6 all miss, and tier 7 (a hard-coded default) wins without warning. This is the single most common silent-correctness failure on SMI data.
>
> **Concrete example:** `pil2M_active_beamstop` is the string `"rod"` or `"pin"`. If neither baseline, conf, nor primary recorded it for a given run, every standard reducer defaults to `"rod"` — and any pin polygon that should have masked the small disk near the beam center is silently absent. The reduced curve then has a smooth, unphysical bump at low q.

> **WARNING — sample_name as fallback metadata.** SMI's older convention encodes geometry into the `sample_name` string (e.g. `EG_AGB_16.10keV_wa14.5_sdd2.0m_ai0.50`). Tier 6 parses these tokens. This works when present, but a typo in the sample name silently changes the geometry. Modern best practice: don't rely on the sample-name parse if the scalar lives in baseline.

For the full decision tree per-scalar, the regex tokens parsed out of `sample_name`, and how to verify each tier separately, see `reference/data-organization.md`.

## Run identity and scan typology

A run's `start.plan_name` tells you what kind of measurement it was. The common SMI plans:

| `plan_name` (typical) | Meaning |
|---|---|
| `count` | Single frame at fixed motor positions. |
| `scan` / `rel_scan` / `list_scan` | Step scan over one or more motors; one frame per step. |
| `grid_scan` | 2-D raster (e.g. `piezo_x` × `piezo_y`). |
| `fly` | Continuous motion with on-the-fly readout. |
| WAXS-arc plan (custom) | Steps `waxs_arc`, with `waxs_bsx` linked: see `reference/beamline-physics.md` on the −4.39 mm/deg coupling. |

For any plan, the per-frame motor positions on `primary` tell you the actual sample/detector configuration of each frame. Don't assume positions from the plan name.

## Sample identity: free string vs structured tokens

SMI's `start.sample_name` is a free-form string written by the operator. By convention, tokens of the form `<label><number>[<unit>]` carry per-sample metadata:

```
Lucas_sample2_pos1_2450.00eV_ai0.50_wa9_bpm1.995_degC100.0
```

A label-number-unit token parser (used by `smi-tiled` to derive virtual axes prefixed `fn:`) yields:

```
{ "sample": 2, "pos": 1, "eV": 2450.0, "ai": 0.5,
  "wa": 9, "bpm": 1.995, "degC": 100.0 }
```

These per-frame "virtual axes" are how SMI users typically pivot a multi-position grid scan into per-sample / per-position curves. The convention is fragile (rename the field and the axis disappears); see `reference/acquisition-upgrades.md` for the proposed structured-sample replacement.

A second convention encodes geometry into `sample_name`:
- `_<digits>keV` or `_<digits>eV` — photon energy.
- `_sdd<digits>m` — sample-detector distance in meters.
- `_wa<digits>` — WAXS arc reference angle.
- `_ai<digits>` — incident angle (grazing incidence).
- `_th<digits>` — sample theta.

This is the bottom-of-stack fallback; baseline always wins when present.

## Reference catalog

When the task needs depth on a topic, read the matching reference file under `reference/`.

| File | Topic |
|---|---|
| `reference/beamline-physics.md` | The instrument: detectors (Pilatus 2M flat / Pilatus 900KW 3-panel arc), pixel sizes, energy range, beamstops (rod/pin geometry and purpose), motor inventory, mechanical linkages, calibration standards (AgB), energy↔wavelength relation, transmission vs grazing-incidence geometry. |
| `reference/data-organization.md` | The data store: `start` / `stop` / `primary` / `primary.config[…]` / `baseline` stream layout and contents; canonical field names (`pil2M_beam_center_*`, `saxs_beamstop_x_{rod,pin}`, `pil2M_active_beamstop`, `waxs_arc`, `waxs_bsx`, …); the full per-scalar metadata-resolution decision tree; sample-name token grammar; image array dimensions and lazy-load behavior. |
| `reference/acquisition-upgrades.md` | Proposed schema improvements: explicit `active_beamstop` (no silent default); explicit per-frame `incident_angle_deg`; structured sample dictionary; unit suffixes everywhere; schema version field; mask-provenance and calibration-uid linkage. Migration path that preserves backwards compatibility with existing runs. |

## When to use which file

- **"What does this scan contain?"** → `reference/data-organization.md`.
- **"Which beam stop was in beam?" / "What's the beam center?"** → `reference/data-organization.md` (decision tree).
- **"What is `waxs_arc` and why does `waxs_bsx` move with it?"** → `reference/beamline-physics.md`.
- **"Why is the metadata so fragile?" / "What should we record going forward?"** → `reference/acquisition-upgrades.md`.
- **"How do I reduce this with `smi-tiled`?"** → use the `smi-tiled` skill (this skill does not assume you are using that package).

## House style

- This skill is **code-agnostic**. It describes the data and the physics, not any specific reader. Field names that appear (`pil2M_active_beamstop`, `target_file_name`, etc.) are the literal Tiled field names — they exist in the data, not in any package.
- When extending: add data-level facts and cross-references to other reference files; do not introduce code citations or function names. Code-specific knowledge belongs in the `smi-tiled` skill.
