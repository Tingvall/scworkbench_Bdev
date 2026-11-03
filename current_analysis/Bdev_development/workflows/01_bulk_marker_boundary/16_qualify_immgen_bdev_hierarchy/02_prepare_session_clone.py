#!/usr/bin/env python3
"""Create the copy-on-write bundle required for unpromoted scWorkbench sessions."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import yaml


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).resolve().read_text())
    inputs = yaml.safe_load(Path(config["inputs_manifest"]).read_text())
    parent = inputs["parent_execution"]
    source = Path(parent["execution_clone"])
    destination = Path(config["storage"]["session_bundle"])
    if destination.exists():
        raise RuntimeError(f"Session bundle already exists: {destination}")
    if sha256_file(source / "manifest.json") != parent["manifest_sha256"]:
        raise RuntimeError("Parent bundle manifest identity changed")
    source_manifest = json.loads((source / "manifest.json").read_text())
    if source_manifest["current_state"]["fingerprint"] != parent["bundle_fingerprint"]:
        raise RuntimeError("Parent bundle fingerprint changed")
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["cp", "--archive", "--reflink=auto", "--", str(source), str(destination)], check=True)
    if sha256_file(destination / "manifest.json") != parent["manifest_sha256"]:
        raise RuntimeError("Session clone manifest differs from parent")
    receipt = Path(config["storage"]["results_root"]) / "06_manifests" / "02_session_clone_receipt.json"
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps({
        "schema_version": "bdev.hierarchical_session_clone_receipt.v1",
        "status": "copied_for_unpromoted_session_execution",
        "parent": str(source),
        "session_bundle": str(destination),
        "fingerprint_before_execution": parent["bundle_fingerprint"],
        "manifest_sha256_before_execution": parent["manifest_sha256"],
        "copy_command": ["cp", "--archive", "--reflink=auto"],
        "promotion_planned": False,
    }, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
