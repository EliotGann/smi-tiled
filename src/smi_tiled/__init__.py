"""smi-tiled — Tiled-native SMI (NSLS-II 12-ID) WAXS+SAXS loader/integrator/uploader.

Quick start
-----------
>>> from smi_tiled import reduce_smi_combined, TiledSMISWAXSLoader
>>> result = reduce_smi_combined(uid="6e61b977-…")
>>> da = result.to_dataarray("merged_iq")   # for PyHyperScattering accessors

See :mod:`smi_tiled.loader` (raw image loading), :mod:`smi_tiled.integrator`
(reduction pipeline), :mod:`smi_tiled.defaults` (bundled masks + calibration),
and :mod:`smi_tiled.upload` (Tiled write-back, optional).
"""
from __future__ import annotations

# Top-level convenience re-exports.  Detailed APIs live under the submodules.
from .loader import (
    TiledSMISWAXSLoader,
    SAXSGeometry,
    WAXSGeometry,
    WAXSPanelGeometry,
    resolve_saxs_geometry,
    resolve_waxs_geometry,
    parse_sample_name_geometry,
    infer_detectors_and_steps,
    load_saxs_raw,
    load_waxs_raw,
    populate_cache,
    clear_baseline_cache,
    SAXS_IMAGE_FIELD,
    WAXS_IMAGE_FIELD,
    WAXS_ARC_FIELD,
    WAXS_BSX_FIELD,
)

from .integrator import (
    reduce_smi_combined,
    reduce_smi_gi,
    CombinedReductionResult,
    GIReductionResult,
    integrate_saxs,
    integrate_waxs,
    MultiPanelArcDetector,
    WAXSCalibration,
    PanelSpec,
    ProgressCallback,
    make_saxs_mask_from_spec,
    make_saxs_mask_from_dict,
    make_waxs_mask_callable,
    make_waxs_mask_callable_from_dict,
    polygons_to_mask,
    mask_for_frame,
    clear_geometry_cache,
    geometry_cache_info,
    merge_q_chi_weighted,
    merge_iq_profiles,
    merge_multiple_qchi,
    merge_reduction_results,
)

from . import defaults

try:
    from ._version import __version__  # populated by setuptools-scm
except ImportError:  # source checkout without an editable install
    __version__ = "0.0.0+unknown"

__all__ = [
    "__version__",
    # Loader
    "TiledSMISWAXSLoader",
    "SAXSGeometry", "WAXSGeometry", "WAXSPanelGeometry",
    "resolve_saxs_geometry", "resolve_waxs_geometry",
    "parse_sample_name_geometry", "infer_detectors_and_steps",
    "load_saxs_raw", "load_waxs_raw",
    "populate_cache", "clear_baseline_cache",
    "SAXS_IMAGE_FIELD", "WAXS_IMAGE_FIELD", "WAXS_ARC_FIELD", "WAXS_BSX_FIELD",
    # Integrator
    "reduce_smi_combined", "reduce_smi_gi",
    "CombinedReductionResult", "GIReductionResult",
    "integrate_saxs", "integrate_waxs",
    "MultiPanelArcDetector", "WAXSCalibration", "PanelSpec",
    "ProgressCallback",
    # Masks
    "make_saxs_mask_from_spec", "make_saxs_mask_from_dict",
    "make_waxs_mask_callable", "make_waxs_mask_callable_from_dict",
    "polygons_to_mask", "mask_for_frame",
    # Caches
    "clear_geometry_cache", "geometry_cache_info",
    # Merge helpers
    "merge_q_chi_weighted", "merge_iq_profiles",
    "merge_multiple_qchi", "merge_reduction_results",
    # Submodules
    "defaults",
]
