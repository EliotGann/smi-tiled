"""CLI entry point for ``smi-upload``.

Usage:
    smi-upload upload  --uid <UID> [--n-q 2000 --tiled-uri …]
    smi-upload batch   --sample "MyFilm" --since 2026-05-01 [...]
    smi-upload status  --uid <UID>
    smi-upload get-iq  --uid <UID> --out merged_iq.nc
"""
from __future__ import annotations

import argparse
import sys
from typing import Sequence


def _add_common_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--tiled-uri", required=True,
                   help="Writable Tiled server URI.")
    p.add_argument("--catalog", default="smi-reduced",
                   help="Catalog path within the server.")
    p.add_argument("--proposal-id", default=None,
                   help="Optional proposal-id sub-group.")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="smi-upload",
                                     description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_upload = sub.add_parser("upload", help="Upload one scan.")
    _add_common_args(p_upload)
    p_upload.add_argument("--uid", required=True)
    p_upload.add_argument("--n-q", type=int, default=2000)
    p_upload.add_argument("--n-chi", type=int, default=360)
    p_upload.add_argument("--overwrite",
                          choices=["never", "always", "if_stale"],
                          default="if_stale")

    p_batch = sub.add_parser("batch", help="Search and upload many scans.")
    _add_common_args(p_batch)
    p_batch.add_argument("--sample", default=None)
    p_batch.add_argument("--since", default=None)
    p_batch.add_argument("--until", default=None)
    p_batch.add_argument("--limit", type=int, default=None)
    p_batch.add_argument("--overwrite",
                         choices=["never", "always", "if_stale"],
                         default="if_stale")

    p_status = sub.add_parser("status", help="Show cache status for a UID.")
    _add_common_args(p_status)
    p_status.add_argument("--uid", required=True)

    args = parser.parse_args(argv)

    # Lazy import so `smi-upload --help` works without tiled installed.
    from .session import UploadSession
    session = UploadSession(
        tiled_uri=args.tiled_uri,
        catalog=args.catalog,
        proposal_id=args.proposal_id,
    )

    if args.cmd == "upload":
        result = session.reduce_and_upload(
            uid=args.uid,
            overwrite=args.overwrite,
            n_q=args.n_q,
            n_chi=args.n_chi,
        )
        print(f"OK  uid={result.uid}  url={result.catalog_url}  "
              f"new={result.new_node}  hash={result.reduction_hash}")
        return 0

    if args.cmd == "batch":
        print("smi-upload batch — not yet implemented (sandbox not wired)")
        return 1

    if args.cmd == "status":
        print(session.status(args.uid))
        return 0

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
