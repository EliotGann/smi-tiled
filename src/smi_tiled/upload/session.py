"""``UploadSession`` — write reduced SMI data to a writable Tiled catalog.

This is the central object users interact with.  It encapsulates:

  - the writable Tiled client (authenticated as needed),
  - the catalog hierarchy convention (see :mod:`smi_tiled.upload.schema`),
  - upload / retrieval / staleness logic.

.. note::
   This module is a **skeleton** — the public API and docstrings are
   final, but the upload/get methods raise ``NotImplementedError`` until
   the sandbox Tiled server is provisioned.  See ``ROADMAP`` in the
   package README.
"""
from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass, field
from typing import Any, Literal

import xarray as xr

from .hash import reduction_hash
from .schema import (
    PROVENANCE_FIELDS,
    REDUCED_DATA_KEYS,
    catalog_path_for,
)


OverwritePolicy = Literal["never", "always", "if_stale"]


@dataclass
class UploadResult:
    """Outcome of a single :meth:`UploadSession.upload` call."""

    uid: str
    catalog_url: str
    new_node: bool          # True if the node was created on this call
    reduction_hash: str
    skipped: bool = False   # True if overwrite='never' and node existed
    written_keys: tuple[str, ...] = ()


@dataclass
class UploadSession:
    """Writable Tiled session for SMI reduced data.

    Parameters
    ----------
    tiled_uri : str
        URI of a *writable* Tiled server (typically a sandbox).
    catalog : str
        Slash-separated catalog path within the server (e.g.
        ``"smi-reduced"``).
    api_key : str, optional
        Tiled API key.  Falls back to ``$TILED_API_KEY`` env var.
    proposal_id : str, optional
        If supplied, scans are stored under
        ``<catalog>/<proposal_id>/<uid>``.  Otherwise the flat layout
        ``<catalog>/<uid>`` is used.
    """

    tiled_uri: str
    catalog: str
    api_key: str | None = None
    proposal_id: str | None = None

    _client: Any = field(default=None, init=False, repr=False)
    _catalog_node: Any = field(default=None, init=False, repr=False)

    # ----------------------------------------------------------------
    # Tiled client lifecycle
    # ----------------------------------------------------------------

    def _connect(self) -> Any:
        if self._client is not None:
            return self._client
        # Lazy import — only needed when actually talking to a server.
        from tiled.client import from_uri  # type: ignore[import-not-found]
        import os
        kwargs: dict[str, Any] = {}
        api_key = self.api_key or os.environ.get("TILED_API_KEY")
        if api_key:
            kwargs["api_key"] = api_key
        self._client = from_uri(self.tiled_uri, **kwargs)
        return self._client

    def _get_catalog(self) -> Any:
        if self._catalog_node is not None:
            return self._catalog_node
        node = self._connect()
        for segment in self.catalog.split("/"):
            if segment:
                node = node[segment]
        self._catalog_node = node
        return node

    # ----------------------------------------------------------------
    # Write paths
    # ----------------------------------------------------------------

    def upload(
        self,
        result: Any,                          # CombinedReductionResult
        *,
        overwrite: OverwritePolicy = "if_stale",
        extra_params: dict[str, Any] | None = None,
    ) -> UploadResult:
        """Write a :class:`CombinedReductionResult` to the sandbox.

        Parameters
        ----------
        result : CombinedReductionResult
            From :func:`smi_tiled.integrator.reduce_smi_combined`.
        overwrite : {'never', 'always', 'if_stale'}
            How to handle an existing node for the same UID.
        extra_params : dict, optional
            Additional reduction parameters to feed into the
            ``reduction_hash`` (e.g. mask file paths, env metadata).
        """
        raise NotImplementedError(
            "UploadSession.upload — sandbox Tiled server not yet wired up. "
            "See smi-tiled-upload PLAN.md."
        )

    def reduce_and_upload(
        self,
        uid: str,
        *,
        overwrite: OverwritePolicy = "if_stale",
        **reduce_kwargs: Any,
    ) -> UploadResult:
        """Convenience: reduce a UID and upload in one call.

        Equivalent to::

            from smi_tiled import reduce_smi_combined
            result = reduce_smi_combined(uid=uid, **reduce_kwargs)
            session.upload(result, overwrite=overwrite)
        """
        from smi_tiled import reduce_smi_combined
        result = reduce_smi_combined(uid=uid, **reduce_kwargs)
        return self.upload(result, overwrite=overwrite,
                           extra_params=reduce_kwargs)

    # ----------------------------------------------------------------
    # Read paths
    # ----------------------------------------------------------------

    def get_merged_iq(self, uid: str) -> xr.Dataset:
        """Fetch the cached merged I(q) for *uid*."""
        return self._get_product(uid, "merged_iq")

    def get_merged_qchi(self, uid: str) -> xr.Dataset:
        """Fetch the cached merged I(q, chi) for *uid*."""
        return self._get_product(uid, "merged_qchi")

    def get_per_frame_iq(self, uid: str) -> xr.Dataset:
        """Fetch the cached per-frame I(q) for *uid*."""
        return self._get_product(uid, "per_frame_iq")

    def _get_product(self, uid: str, key: str) -> xr.Dataset:
        if key not in REDUCED_DATA_KEYS:
            raise ValueError(
                f"key must be one of {REDUCED_DATA_KEYS}, got {key!r}"
            )
        raise NotImplementedError(
            f"UploadSession._get_product({uid!r}, {key!r}) — sandbox not wired."
        )

    # ----------------------------------------------------------------
    # Inspection
    # ----------------------------------------------------------------

    def status(self, uid: str) -> dict[str, Any]:
        """Return cache status for *uid* (exists / hash / timestamp / size)."""
        raise NotImplementedError("UploadSession.status — sandbox not wired.")

    # ----------------------------------------------------------------
    # Helpers
    # ----------------------------------------------------------------

    def _node_path(self, uid: str) -> tuple[str, ...]:
        return catalog_path_for(uid, proposal_id=self.proposal_id)

    @staticmethod
    def _provenance_dict(
        result: Any,
        reduction_params: dict[str, Any] | None,
    ) -> dict[str, Any]:
        """Build the metadata payload attached to the <uid> node."""
        scan_info = getattr(result, "scan_info", None) or {}
        rd: dict[str, Any] = {
            "uid": result.uid,
            "scan_id": scan_info.get("scan_id"),
            "sample_name": scan_info.get("sample_name"),
            "plan_name": scan_info.get("plan_name"),
            "geometry": getattr(result, "geometry", None),
            "incident_angle_deg": getattr(result, "incident_angle_deg", None),
            "upload_timestamp": _dt.datetime.now(
                _dt.timezone.utc
            ).isoformat(timespec="seconds"),
        }
        if reduction_params:
            rd.update({k: v for k, v in reduction_params.items()
                       if k in PROVENANCE_FIELDS})
        rd["reduction_hash"] = reduction_hash(reduction_params or {})
        try:
            from smi_tiled import __version__ as _v
            rd["smi_tiled_version"] = _v
        except Exception:
            pass
        return rd
