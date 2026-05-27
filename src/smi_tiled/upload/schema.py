"""Catalog hierarchy and metadata schema for SMI reduced-data Tiled sandbox.

Layout
------
    smi-reduced/
      <proposal_id>/         ← optional top-level grouping
        <uid>/               ← one node per source scan
          merged_iq          ← xr.Dataset on (q,)
          merged_qchi        ← xr.Dataset on (q, chi)
          per_frame_iq       ← xr.Dataset on (frame, q)

Each ``<uid>`` node carries the full reduction provenance as metadata
(see :data:`PROVENANCE_FIELDS`) plus a ``reduction_hash`` for
stale-cache detection.
"""
from __future__ import annotations

#: Names of the xr.Dataset products that may be stored under a ``<uid>``
#: node.  ``per_frame_iq`` is optional (single-frame scans omit it).
REDUCED_DATA_KEYS: tuple[str, ...] = ("merged_iq", "merged_qchi", "per_frame_iq")

#: Fields copied from the source run + reduction params into per-node metadata.
#: Used both for inspection and for ``reduction_hash`` computation.
PROVENANCE_FIELDS: tuple[str, ...] = (
    # source run identity
    "uid", "scan_id", "sample_name", "plan_name",
    # source geometry
    "energy_kev", "sdd_mm",
    "detector_kinds",          # ('saxs',) / ('waxs',) / ('saxs', 'waxs')
    "geometry", "incident_angle_deg",
    # reduction parameters
    "n_q", "n_chi",
    "pixel_splitting", "dezinger_threshold", "dezinger_kernel",
    "saxs_q_cutoff", "saxs_agbh_ring_order", "saxs_q_margin_fraction",
    "waxs_beam_col_per_arc_deg", "solid_angle_correction",
    "saxs_beam_delta_px", "waxs_beam_delta_px", "saxs_distance_delta_mm",
    "saxs_mask_hash", "waxs_mask_hash",
    # package identity
    "smi_tiled_version",
    "reduction_hash",
    "upload_timestamp",
)


def catalog_path_for(uid: str, proposal_id: str | None = None) -> tuple[str, ...]:
    """Return the slash-segment catalog path for a given UID.

    Two-level layout if a proposal_id is supplied; otherwise flat.
    """
    if proposal_id:
        return (str(proposal_id), str(uid))
    return (str(uid),)
