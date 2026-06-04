"""Optional, composable derived-analysis stages for a reduction result.

These helpers attach further analysis products onto an existing
:class:`~smi_tiled.integrator.CombinedReductionResult` (or
:class:`~smi_tiled.integrator.GIReductionResult`) without re-running the
heavy integration pipeline:

* :mod:`smi_tiled.derived.virtual_axes` — parse numeric tokens out of
  per-frame string fields (``target_file_name`` etc.) and attach them as
  ``fn:*`` data variables on ``per_frame_iq``.
* :mod:`smi_tiled.derived.linecuts` — frame-by-frame cross sections of
  the merged or per-frame qchi maps.
* :mod:`smi_tiled.derived.peakfit` — per-frame Gaussian/Lorentzian peak
  fits across ``per_frame_iq``.

Each helper has the form ``apply_*(result, spec, ...) -> result``: the
returned result is the same object decorated with the new product (the
underlying frozen dataclass is mutated in-place via ``object.__setattr__``
to keep the result identity stable).

The same stages are invoked automatically from
:func:`smi_tiled.integrator.reduce_smi_combined` when the matching
``virtual_axes`` / ``line_cuts`` / ``peak_fits`` kwargs are supplied.
"""
from __future__ import annotations

from .virtual_axes import (
    VIRTUAL_PREFIX,
    VirtualAxesConfig,
    apply_virtual_axes,
    derive_virtual_columns,
    parse_label_number_tokens,
)
from .linecuts import (
    LineCutSpec,
    apply_line_cuts,
    compute_cross_section,
)
from .peakfit import (
    FIT_PARAMS,
    MIN_R2,
    MIN_SNR,
    PeakDef,
    apply_peak_fits,
    fit_peak_across_frames,
)

__all__ = [
    "VIRTUAL_PREFIX",
    "VirtualAxesConfig",
    "apply_virtual_axes",
    "derive_virtual_columns",
    "parse_label_number_tokens",
    "LineCutSpec",
    "apply_line_cuts",
    "compute_cross_section",
    "FIT_PARAMS",
    "MIN_R2",
    "MIN_SNR",
    "PeakDef",
    "apply_peak_fits",
    "fit_peak_across_frames",
]
