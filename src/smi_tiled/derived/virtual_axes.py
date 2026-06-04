"""Virtual per-frame axes parsed out of structured string fields.

This is the smi-tiled home of what used to live in
``smi_browser.data.scalars`` — the regex-based parser that turns
structured per-frame strings (such as
``"Lucas_sample2_pos1_2450.00eV_ai0.50_degC100.0"``) into numeric
``fn:*`` columns that can drive any axis selector.

The classic usage is :func:`derive_virtual_columns` on a per-frame
``pandas.DataFrame``.  The reduction pipeline calls the higher-level
:func:`apply_virtual_axes` on a
:class:`~smi_tiled.integrator.CombinedReductionResult` so the resulting
``fn:*`` axes ride along on ``result.per_frame_iq`` as 1-D data variables
on the ``frame`` dim.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr

#: Matches a ``label`` + ``number`` (+ ``unit``) token, where at least one of
#: the alphabetic groups touches the number, e.g. ``ai0.50`` (prefix ``ai``),
#: ``2450.00eV`` (unit ``eV``), ``degC100.0`` (prefix ``degC``).  Bare numbers
#: with no adjacent letters (``120``) do not yield a label and are ignored.
_TOKEN_RE = re.compile(
    r"([A-Za-z][A-Za-z%°µ/]*)?(-?\d+(?:\.\d+)?)([A-Za-z%°µ/]+)?"
)

#: Default prefix marking columns derived from a string field.
VIRTUAL_PREFIX = "fn:"


@dataclass(frozen=True)
class VirtualAxesConfig:
    """Knobs for :func:`apply_virtual_axes` / :func:`derive_virtual_columns`.

    Parameters
    ----------
    sources : tuple of str, optional
        Names of string-typed per-frame fields to parse.  ``None``
        (default) means auto-discover every non-numeric, non-``ts_``
        per-frame data var.
    min_fill : float
        A derived column is kept only when its non-NaN fraction is at
        least this value.  Filters out noise from free-text status
        fields.  Default 0.5.
    prefix : str
        Prepended to every derived column.  Default ``"fn:"``.
    enabled : bool
        Master switch.  When ``False`` the pipeline skips parsing
        entirely.
    """

    sources: tuple[str, ...] | None = None
    min_fill: float = 0.5
    prefix: str = VIRTUAL_PREFIX
    enabled: bool = True

    def to_provenance(self) -> dict[str, Any]:
        """Stable dict form for hashing into ``reduction_hash``."""
        return {
            "sources": list(self.sources) if self.sources else None,
            "min_fill": float(self.min_fill),
            "prefix": str(self.prefix),
            "enabled": bool(self.enabled),
        }


def parse_label_number_tokens(s: Any) -> dict[str, float]:
    """Extract ``{label: number}`` pairs from a structured string.

    Tokens are ``label`` + ``number`` (+ ``unit``) runs; the column label
    is the alphabetic *prefix* when present, else the trailing *unit*.
    A bare number with no adjacent letters yields nothing.

    Examples
    --------
    ``"Lucas_sample2_pos1_2450.00eV_ai0.50_wa9_bpm1.995_degC100.0"`` →
    ``{"sample": 2, "pos": 1, "eV": 2450.0, "ai": 0.5, "wa": 9,
       "bpm": 1.995, "degC": 100.0}``.
    """
    out: dict[str, float] = {}
    # Strings cached in HDF5 come back as ``bytes`` (h5py variable-length
    # strings); tiled returns ``str``.  Accept both, plus numpy string scalars.
    if isinstance(s, (bytes, bytearray, np.bytes_)):
        try:
            s = bytes(s).decode("utf-8", "replace")
        except Exception:
            return out
    elif not isinstance(s, str):
        if isinstance(s, np.str_):
            s = str(s)
        else:
            return out
    if not s:
        return out
    for prefix, number, unit in _TOKEN_RE.findall(s):
        label = prefix or unit
        if not label:
            continue  # bare number — no usable axis label
        try:
            value = float(number)
        except ValueError:
            continue
        # First occurrence of a label within the string wins.
        out.setdefault(label, value)
    return out


def derive_virtual_columns(
    df: pd.DataFrame,
    *,
    prefix: str = VIRTUAL_PREFIX,
    sources: list[str] | tuple[str, ...] | None = None,
    min_fill: float = 0.5,
) -> pd.DataFrame:
    """Append numeric columns parsed from structured string columns.

    For each *source* string column, every cell is parsed with
    :func:`parse_label_number_tokens`; the union of labels becomes one
    float column each (``NaN`` where a frame lacks that token).  Columns
    whose non-NaN fraction is below ``min_fill`` are dropped.

    Parameters
    ----------
    sources : sequence of column names, optional
        If ``None`` (default), every object/string per-frame column is
        scanned (numeric and ``ts_`` timestamp columns are skipped).
    prefix : str
        Prepended to every derived column (default ``"fn:"``).  On a
        label collision across sources, the later one is qualified as
        ``prefix + source + ":" + label``.

    Returns the original frame unchanged if nothing is derived.
    """
    if df is None or df.empty:
        return df

    if sources is None:
        sources = [
            c for c in df.columns
            if not str(c).startswith("ts_")
            and not pd.api.types.is_numeric_dtype(df[c])
        ]

    n = len(df)
    new_cols: dict[str, np.ndarray] = {}
    for src in sources:
        if src not in df.columns:
            continue
        parsed = [parse_label_number_tokens(v) for v in df[src].tolist()]
        labels: list[str] = []
        for d in parsed:
            for lab in d:
                if lab not in labels:
                    labels.append(lab)
        for lab in labels:
            col = np.array([d.get(lab, np.nan) for d in parsed], dtype=float)
            if np.isfinite(col).sum() < min_fill * n:
                continue
            name = f"{prefix}{lab}"
            if name in df.columns or name in new_cols:
                name = f"{prefix}{src}:{lab}"  # disambiguate across sources
            if name in df.columns or name in new_cols:
                continue
            new_cols[name] = col

    if not new_cols:
        return df
    out = df.copy()
    for name, col in new_cols.items():
        out[name] = col
    return out


def _per_frame_string_dataframe(
    per_frame_iq: xr.Dataset,
    scan_info: dict[str, Any] | None,
    sources: tuple[str, ...] | None,
) -> pd.DataFrame:
    """Build a per-frame DataFrame of string-typed candidate columns.

    Sources are merged from two places, in this order (first wins):

    1. ``per_frame_iq`` data variables on the ``frame`` dim whose dtype
       is object / string / bytes.
    2. ``scan_info["primary_strings"]`` — a ``{name: 1-D array}`` dict
       written by :func:`smi_tiled.loader.infer_detectors_and_steps`
       for per-frame string columns it skipped while building the
       numeric ``step_candidates``.
    """
    n_frames = int(per_frame_iq.sizes.get("frame", 0))
    if n_frames == 0:
        return pd.DataFrame()

    columns: dict[str, np.ndarray] = {}

    def _is_stringy(arr: np.ndarray) -> bool:
        if arr.dtype == object:
            return True
        kind = getattr(arr.dtype, "kind", "")
        return kind in ("U", "S", "O")

    # 1) data vars on per_frame_iq
    for name, var in per_frame_iq.data_vars.items():
        if var.dims != ("frame",):
            continue
        values = np.asarray(var.values)
        if values.shape[0] != n_frames:
            continue
        if not _is_stringy(values):
            continue
        columns[str(name)] = values

    # 2) primary_strings dict from scan_info
    if scan_info is not None:
        primary_strings = scan_info.get("primary_strings") or {}
        for name, arr in primary_strings.items():
            if name in columns:
                continue
            values = np.asarray(arr)
            if values.shape[0] != n_frames:
                continue
            if not _is_stringy(values):
                continue
            columns[str(name)] = values

    if sources is not None:
        columns = {k: v for k, v in columns.items() if k in sources}

    if not columns:
        return pd.DataFrame()
    return pd.DataFrame(columns)


def apply_virtual_axes(
    result: Any,
    config: VirtualAxesConfig | None = None,
) -> Any:
    """Attach ``fn:*`` per-frame data variables onto ``result.per_frame_iq``.

    Parses every string-typed per-frame field reachable from *result*
    (``per_frame_iq`` data vars + ``scan_info["primary_strings"]``)
    according to *config*, then rebuilds ``result.per_frame_iq`` with
    the resulting numeric columns attached as 1-D data variables on the
    ``frame`` dim.

    Returns *result* (mutated in place via ``object.__setattr__`` for
    frozen dataclasses).  A no-op when ``config.enabled`` is ``False``
    or ``per_frame_iq`` is missing/empty.
    """
    config = config or VirtualAxesConfig()
    if not config.enabled:
        return result
    per_frame_iq = getattr(result, "per_frame_iq", None)
    if per_frame_iq is None:
        return result

    df = _per_frame_string_dataframe(
        per_frame_iq, getattr(result, "scan_info", None), config.sources,
    )
    if df.empty:
        return result

    derived = derive_virtual_columns(
        df,
        prefix=config.prefix,
        sources=list(df.columns),
        min_fill=config.min_fill,
    )
    new_cols = [c for c in derived.columns if c not in df.columns]
    if not new_cols:
        return result

    updates: dict[str, tuple[tuple[str, ...], np.ndarray]] = {}
    for name in new_cols:
        if name in per_frame_iq.data_vars:
            continue  # never clobber existing
        updates[name] = (("frame",), np.asarray(derived[name].values, dtype=float))

    if not updates:
        return result

    new_per_frame_iq = per_frame_iq.assign({
        n: xr.DataArray(v, dims=d) for n, (d, v) in updates.items()
    })
    object.__setattr__(result, "per_frame_iq", new_per_frame_iq)
    return result
