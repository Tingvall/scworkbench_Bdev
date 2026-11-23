#!/usr/bin/env python3
"""Create a recoverable active-development clone of the frozen HVP bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


EXPECTED_SOURCE_MANIFEST = "1f9e399600a9ab74b6037301af8a6b521bf4669b8cfd4b76f80f2d07ba437b98"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-bundle", type=Path, required=True)
    parser.add_argument("--clone", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()

    observed = sha256(args.source_bundle / "manifest.json")
    if observed != EXPECTED_SOURCE_MANIFEST:
        raise RuntimeError(f"Frozen source manifest changed: {observed}")
    if args.clone.exists():
        raise RuntimeError(f"Refusing to overwrite existing clone: {args.clone}")
    args.clone.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([
        "cp", "--recursive", "--reflink=auto", "--",
        str(args.source_bundle), str(args.clone),
    ], check=True)
    if sha256(args.clone / "manifest.json") != EXPECTED_SOURCE_MANIFEST:
        raise RuntimeError("Clone manifest is not initially identical")
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt = {
        "status": "complete",
        "source_bundle": str(args.source_bundle),
        "source_manifest_sha256": EXPECTED_SOURCE_MANIFEST,
        "clone": str(args.clone),
        "clone_initial_manifest_sha256": EXPECTED_SOURCE_MANIFEST,
        "copy_mode": "cp_recursive_reflink_auto",
        "source_mutation": False,
    }
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")


if __name__ == "__main__":
    main()
