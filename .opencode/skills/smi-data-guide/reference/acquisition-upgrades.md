# Acquisition metadata upgrades for SMI

A design document, not a description of the current state. This file proposes schema improvements to SMI's per-run metadata that would eliminate the silent-correctness pathologies documented in `data-organization.md`. Read `data-organization.md` first; this file presupposes you understand the current schema and its failure modes.

The intent is that this guide can guide future bsui plan / device updates and inform conversations with the beamline team about what should change at the **acquisition** layer (not the analysis layer).

## Why upgrade?

The current schema has worked for years, but it has accumulated structural problems that show up as silent-correctness failures in reduction:

1. **Silent fallback-chain defaults.** Multiple analysis-critical scalars (active beamstop, incident angle, energy in some old runs) silently fall through to a hard-coded default when no source recorded them. The reducer produces a result that looks fine but is wrong.

2. **Geometry encoded in free-form strings.** `sample_name` is parsed for `_ai<X>`, `_sdd<X>m`, `_<X>keV` tokens because there is no structured per-frame field for these. A typo silently corrupts the analysis.

3. **No schema version.** A run from 2019 and a run from 2025 have different field sets, different semantics, and no version tag. Readers have to detect the era heuristically.

4. **Mixed units.** Energy is sometimes eV, sometimes keV. Distances are sometimes mm, sometimes m, sometimes μm (piezo). Without unit suffixes on field names or unit metadata, every reader hard-codes assumptions.

5. **Sample identity is a free string.** `sample_name = "Lucas_sample2_pos1_2450eV_ai0.50_wa9"` mixes operator initials, sample ID, position, energy, incident angle, and WAXS arc into one unparseable column. Searching for "all samples named X" is regex-matching across runs.

6. **Beamstop position only at scan start.** When the rod or pin is moved mid-scan (e.g. some calibration sequences), only baseline records the start position. Per-frame movements are not captured.

7. **No mask provenance.** The mask used for a reduction is implicit: it's whatever the reducer's bundled JSON happens to be at run time. There's no record on the run pointing to which mask version was *intended* for it.

8. **No calibration-uid linkage.** A scattering run was reduced using a calibration scan from yesterday, but nothing on the run points at that calibration. If you re-reduce later, you have to remember which calibration it was paired with.

The proposals below address each, ordered roughly by impact.

## Proposal 1: kill silent defaults — explicit `active_beamstop`

**Problem.** `pil2M_active_beamstop` (`"rod"` or `"pin"`) is the single most consequential SAXS scalar. When absent from baseline / conf / primary, every reducer falls through to `"rod"`. Pin-beamstop runs without this field record produce an unphysical bump in the reduced low-q curve.

**Proposal.** Make `pil2M_active_beamstop` a **required** field in the start document, written by the bsui plan that configures the SAXS detector. If the plan cannot determine which beamstop is in beam, it should refuse to run, not write `None`.

**Schema:**

```python
start = {
    "scan_id": ...,
    "uid": ...,
    "saxs": {                   # nested namespace, see proposal 4
        "active_beamstop": "rod" | "pin",     # REQUIRED
        "beamstop_pos_mm": {
            "rod": {"x": 0.0, "y": 0.0},      # mm
            "pin": {"x": 0.0, "y": 0.0},      # mm
        },
        "beam_center_px": {
            "row": 833.0, "col": 690.0,        # px
        },
        "sdd_mm": 2730.0,                       # mm
    },
    ...
}
```

**Migration.** Old runs (no `start.saxs.active_beamstop`) keep working through the existing fallback chain. New runs hard-fail if it's missing, surfacing the error at the *acquisition* layer where the operator can fix it (not at the analysis layer hours/days later).

**Acceptance check.** A simple plan-time validator: before the run starts, assert `start.saxs.active_beamstop in {"rod", "pin"}`. If the assertion fails, the operator sees the problem immediately.

## Proposal 2: explicit per-frame `incident_angle_deg`

**Problem.** `α_i` is critical for grazing-incidence runs. Currently, it lives in:
- A per-frame primary value (in some plans).
- A baseline scalar (in some setups).
- A `_ai<X>` token in `sample_name` (in many older runs).

A run with no recorded α_i is silently treated as `α_i = 0` (transmission), which puts the specular peak in the wrong place.

**Proposal.** Add a structured field `incident_angle_deg` to the per-frame primary stream for any GI plan, and a baseline scalar `incident_angle_deg_setpoint` for the start-of-scan reference. The plan declares whether it is `"transmission"` or `"grazing_incidence"`:

```python
start = {
    "geometry": "transmission" | "grazing_incidence",     # REQUIRED
    "incident_angle_deg_nominal": 0.0,                    # if GI, the nominal α_i
    ...
}

primary = {
    "incident_angle_deg": [0.50, 0.50, 0.50, ...],        # per-frame, if GI
    ...
}

baseline = {
    "incident_angle_deg_setpoint": 0.50,                  # if GI
    ...
}
```

For transmission runs, `primary.incident_angle_deg` may be absent or a constant 0.0; the `start.geometry` field tells the reader whether to apply GI math.

**Migration.** Old runs continue to use the `sample_name` parse fallback. New runs always carry `start.geometry` explicitly.

## Proposal 3: schema version field

**Problem.** Different eras have different field sets. There is no way for a reader to know which conventions apply.

**Proposal.** Add `start.smi_schema_version` (string, e.g. `"2.0.0"`). Bump it when:
- A field's name changes.
- A field's units change.
- A new required field is added.
- The semantics of a stream changes.

```python
start = {
    "smi_schema_version": "2.0.0",
    ...
}
```

The version string follows semver: incompatible reads are major bumps. Readers can dispatch on the version explicitly.

**Migration.** Runs without `smi_schema_version` are version `0.x` (legacy). Readers default to legacy conventions when the field is absent.

## Proposal 4: namespaced, structured top-level start dict

**Problem.** The current `start` document has flat fields with implicit grouping by prefix (`saxs_*`, `pil2M_*`, etc.). Operator-supplied keys mix freely with system-supplied keys.

**Proposal.** Group start-doc fields under namespaces:

```python
start = {
    "smi_schema_version": "2.0.0",
    "scan_id": ...,
    "uid": ...,
    "time": ...,
    "geometry": "transmission",
    "plan_name": "grid_scan",
    "plan_args": { ... },

    "sample": {                     # see proposal 5
        "id": "EG_film_034",
        "name_human": "PS-PMMA blend, batch 4",
        "composition": { ... },
        "tokens": { ... },
    },

    "user": {
        "name": "...",
        "institution": "...",
        "proposal_id": "...",
        "cycle": "2024-2",
    },

    "saxs": {                       # see proposal 1
        "active_beamstop": "rod",
        ...
    },

    "waxs": {
        "arc_deg_nominal": 14.5,
        "bsx_mm_nominal": -22.0,
        "panels": [ ... ],
        ...
    },

    "calibration": {                # see proposal 8
        "saxs_uid": "abc-123-...",
        "saxs_sample": "AgB",
        "waxs_uid": "def-456-...",
    },

    "operator_notes": "free-form ...",
}
```

**Migration.** Readers with a fallback layer can look in both `start.saxs.active_beamstop` (new) and `start.pil2M_active_beamstop` (legacy). Top-level operator-supplied keys (`experiment_id`, etc.) move under a `start.user.extra` dict to keep namespaces clean.

## Proposal 5: structured sample identity

**Problem.** `sample_name = "Lucas_sample2_pos1_2450.00eV_ai0.50_wa9"` is doing the work of:
- Operator identifier.
- Sample identifier.
- Per-position label.
- Energy.
- Incident angle.
- WAXS arc setting.

A typo in any token silently breaks the analysis. Searching for "all samples of polymer X" is regex-matching across thousands of free strings.

**Proposal.** Replace `sample_name` (string) with `sample` (dict):

```python
start.sample = {
    "id": "EG_PSPMMA_034",                  # unique, machine-readable
    "name_human": "PS-PMMA blend, batch 4", # human-readable label
    "composition": {                         # structured
        "PS": {"weight_fraction": 0.5},
        "PMMA": {"weight_fraction": 0.5},
    },
    "preparation": {
        "method": "spin-coating",
        "substrate": "Si",
        "thickness_nm": 100,
    },
    "metadata_uri": "...",                   # optional pointer to ELN
    "tokens": {                              # legacy-compat virtual axes
        "sample": 2,
        "pos": 1,
        "eV": 2450.0,
        "ai": 0.50,
        "wa": 9,
    },
    "name_legacy": "Lucas_sample2_pos1_2450.00eV_ai0.50_wa9",  # original string
}
```

**Migration.** A migration tool reads `sample_name` from old runs and populates the new `sample` dict (with `tokens` from the existing parser, and `name_legacy` as the original). Old readers can still read `sample_name`.

**Per-frame variation.** When the operator wants per-frame sample identification (a grid scan over 100 positions), the per-frame string field stays in primary (`target_file_name`), but the parser becomes:

```python
primary.frame_metadata = [
    {"sample_id": "...", "position_idx": 1, "tokens": {...}},
    {"sample_id": "...", "position_idx": 2, "tokens": {...}},
    ...
]
```

— a structured per-frame dict, not a parsed-from-string dict.

## Proposal 6: units in the schema, not in conventions

**Problem.** `energy_energy` is eV. Some older runs have it in keV. Distances are mm in baseline, m in some `sample_name` tokens, and μm for piezo motors. Every reader hard-codes a unit assumption per field.

**Proposal.** Two complementary changes:

(a) **Field-name suffixes** — when a field's units cannot be standardized:

```
energy_ev                # was: energy_energy
sdd_mm                   # was: pil2M_motor_z
incident_angle_deg
piezo_z_um               # explicit μm
beam_center_row_px       # explicit pixel
beamstop_x_mm
```

(b) **Optional unit metadata** on every numeric field:

```python
primary.metadata["energy_ev"] = {
    "units": "eV",
    "description": "Photon energy",
}
```

(stored as Tiled column-level metadata, not as separate columns).

**Migration.** Old field names map through an alias table maintained by the new readers. Schema versions ≥ 2.0.0 use the new names; versions < 2.0.0 use the alias map.

## Proposal 7: per-frame beamstop and beam-center

**Problem.** When the rod or pin is moved during a scan (rare but it happens — some calibration sequences use a `bs_scan` plan), only baseline records the *start* position. Per-frame movement is silently dropped.

Likewise, if the detector is translated during a scan (e.g. to extend q-coverage), the beam center moves with it. The EPICS PV is supposed to track this, but if the PV is broken or out-of-date, the baseline value is wrong for late frames.

**Proposal.** Promote beamstop position and beam center to per-frame fields when they are known to vary:

```python
primary = {
    "saxs_beamstop_x_mm": [22.0, 22.0, 23.0, ...],   # per-frame, if scanned
    "saxs_beamstop_y_mm": [...],
    "saxs_beam_center_row_px": [...],                # per-frame, if detector moves
    "saxs_beam_center_col_px": [...],
    "saxs_active_beamstop": ["rod", "rod", "pin", ...],  # per-frame, allows mid-scan switching
    ...
}
```

For scans where these don't change, primary can omit the field; the baseline scalar is then authoritative.

**Migration.** Existing readers' resolution chain (baseline → primary first row) already handles per-frame fields when present; this proposal makes per-frame the *recommended* layer, not just allowed.

## Proposal 8: explicit calibration linkage

**Problem.** A scattering run was reduced using AgB calibration scan `1130039` from earlier in the session. That information lives only in the operator's notebook. If the run is re-reduced six months later by a different person, the calibration linkage is lost.

**Proposal.** Store calibration uids on the run:

```python
start.calibration = {
    "saxs_uid": "abc-...-...",       # the AgB run
    "saxs_sample": "AgB",            # the calibrant material
    "saxs_d_nm": 5.838,              # the known d-spacing
    "saxs_q1_inv_nm": 1.076,         # the predicted first ring
    "waxs_uid": "def-...-...",       # the WAXS calibration run
    "waxs_sample": "Si_NIST_640d",
}
```

These are **pointers**, not derived values. The reducer follows the pointer to fetch the calibration parameters fresh.

**Acceptance check.** A reducer that needs calibration walks: caller-override → `start.calibration.saxs_uid` → "no calibration" (refuse, don't silently use a default). The third tier is the failure mode that should not exist; this proposal makes it visible.

## Proposal 9: mask provenance

**Problem.** The mask used in a reduction is whatever the reducer's bundled JSON happens to be. There is no record on the run of which mask was *intended*.

This shows up most painfully when a mask file is updated (e.g. a new pin polygon is committed). All previous reductions are now silently inconsistent with the new mask, with no way to tell which was used for any given output.

**Proposal.** Store mask provenance on the run alongside calibration:

```python
start.mask = {
    "saxs_mask_id": "smi_pil2M_v3",       # named version
    "saxs_mask_hash": "sha256:...",        # content-addressed
    "saxs_mask_uri": "tiled://...",        # optional, where to fetch
    "waxs_mask_id": "smi_pil900KW_v2",
    ...
}
```

**Migration.** Old runs default to `"smi_pil2M_legacy"` (the as-shipped bundled mask at the time of acquisition). New runs always carry an explicit ID and content hash.

## Proposal 10: per-detector configuration completeness

**Problem.** `primary.config[pil2M]` is sparse: some runs record exposure time, some don't. Threshold energy (which affects absolute intensity) is rarely recorded.

**Proposal.** Define a required per-detector config schema:

```python
primary.config["pil2M"] = {
    "acquire_time_s": 0.5,                 # REQUIRED
    "num_images": 1,                       # REQUIRED (frames per trigger)
    "threshold_energy_ev": 8500,           # REQUIRED — affects absolute I
    "active_beamstop": "rod",              # mirrored from start.saxs (DRY)
    "image_shape_px": [1679, 1475],
    "pixel_size_m": [0.000172, 0.000172],
    ...
}

primary.config["pil900KW"] = {
    "acquire_time_s": 0.5,
    "num_images": 1,
    "threshold_energy_ev": 8500,
    "panel_count": 3,
    "panel_offsets_deg": [-7.0, 0.0, 7.0],
    "panel_shape_px": [195, 195],
    ...
}
```

**Migration.** Validators run at scan-start to assert all required fields are present.

## Implementation roadmap

A concrete order of operations for rolling out the above:

1. **Schema version** (Proposal 3) — zero-cost, instantly disambiguates eras.
2. **Calibration linkage** (Proposal 8) and **mask provenance** (Proposal 9) — record-only fields, no reduction-side change required.
3. **Explicit `start.geometry`** (Proposal 2, partial) — stops the silent transmission default for GI runs.
4. **Explicit `start.saxs.active_beamstop`** (Proposal 1) — kills the most consequential silent default.
5. **Structured sample dict** (Proposal 5) — biggest UX improvement; coexists with `sample_name` initially.
6. **Per-frame beamstop / beam center** (Proposal 7) — limited to plans that actually need it; a strict superset of current behavior.
7. **Unit suffixes everywhere** (Proposal 6) — needs a deprecation cycle on field names.
8. **Per-detector config completeness** (Proposal 10) — needs bsui device-definition updates.

## What this preserves

These proposals do **not** change:

- The Tiled URL or catalog name (`smi/migration` stays).
- The Bluesky `BlueskyRun` shape (start / stop / primary / baseline streams).
- Image array layouts (`pil2M_image` is still 1679×1475, `pil900KW_image` is still 195×619).
- The use of `uid` and `scan_id` as identifiers.
- The legacy `sample_name` string, which continues to be written for backwards compatibility.

A reducer that supports schema version 2.0.0 must also support legacy. The migration is one-way: old runs stay old, new runs follow the new schema.

## Why now

The cost of these proposals is moderate — most are device-definition and plan-time validator changes — and the cost of *not* doing them compounds with every new run. Each new pin run without `pil2M_active_beamstop` recorded is a future analyst, who didn't observe the acquisition, silently reducing it as a rod run.

The schema upgrades shift the locus of correctness from "the analyst remembers / has good operator notes" to "the run, on its own, contains everything needed to interpret it." That is the property that makes data archives useful at decade scale.
