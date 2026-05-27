"""Reduction-parameter hash for stale-cache detection."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


def _stable_repr(value: Any) -> Any:
    """Normalize a value into something json.dumps-stable and hashable.

    - Path → str(path)
    - dict → dict with sorted keys
    - tuple/list → list of stable-repr'd entries
    - numpy scalars → Python scalars
    - other → str() fallback
    """
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _stable_repr(v) for k, v in sorted(value.items())}
    if isinstance(value, (list, tuple)):
        return [_stable_repr(v) for v in value]
    # numpy scalars
    item = getattr(value, "item", None)
    if callable(item):
        try:
            return item()
        except Exception:
            pass
    return str(value)


def reduction_hash(params: dict[str, Any], *, algo: str = "sha256",
                   length: int = 16) -> str:
    """Stable short hex digest of a reduction-parameter dict.

    Two calls with semantically equivalent params (any key order, any
    Path-vs-str representation, etc.) produce identical hashes.

    Parameters
    ----------
    params : dict
        Reduction parameters.  Values that are not JSON-native are
        coerced (Path → str, numpy scalars → Python scalars, …).
    algo : str
        Any algorithm name accepted by :mod:`hashlib`.  Default sha256.
    length : int
        How many leading hex characters to return.  Default 16
        (= 64 bits, enough for cache-key purposes).

    Returns
    -------
    str
        ``<algo>:<hex>`` so the algorithm choice is recorded with the
        hash itself.
    """
    normalized = _stable_repr(params)
    blob = json.dumps(normalized, separators=(",", ":"), sort_keys=True)
    digest = hashlib.new(algo, blob.encode("utf-8")).hexdigest()[:length]
    return f"{algo}:{digest}"
