"""Tiled write-back for reduced SMI data.

The upload subsystem is optional — install ``smi-tiled[upload]`` to pull
in ``tiled[client]``, ``httpx``, and ``zarr``.

Typical use:

>>> from smi_tiled.upload import UploadSession
>>> session = UploadSession(
...     tiled_uri="https://tiled-sandbox.nsls2.bnl.gov",
...     catalog="smi-reduced",
... )
>>> session.upload(result)            # result: CombinedReductionResult
>>> iq = session.get_merged_iq(uid)   # later retrieval, no re-compute
"""
from .session import UploadSession
from .hash import reduction_hash
from .schema import (
    catalog_path_for,
    REDUCED_DATA_KEYS,
    PROVENANCE_FIELDS,
)

__all__ = [
    "UploadSession",
    "reduction_hash",
    "catalog_path_for",
    "REDUCED_DATA_KEYS",
    "PROVENANCE_FIELDS",
]
