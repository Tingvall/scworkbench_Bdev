#!/usr/bin/env python3
"""Create an isolated reflink-capable bundle clone for the bounded Bdev test."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.resolve().read_text())
    source = Path(config["inputs"]["source_bundle"]["path"])
    destination = Path(config["storage"]["intermediate_root"]) / "05_objects" / "01_session_bundle"
    if destination.exists():
        raise RuntimeError(f"Refusing to overwrite bundle clone: {destination}")
    if sha256(source / "manifest.json") != config["inputs"]["source_bundle"]["manifest_sha256"]:
        raise RuntimeError("Source bundle manifest checksum changed")
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["cp", "--archive", "--reflink=auto", "--", str(source), str(destination)], check=True)
    cloned_manifest = destination / "manifest.json"
    state = json.loads(cloned_manifest.read_text())
    if state["current_state"]["fingerprint"] != config["inputs"]["source_bundle"]["fingerprint"]:
        raise RuntimeError("Cloned bundle fingerprint changed")
    receipt = {
        "schema_version": "bdev.stable_immgen_broad_bundle_clone.v1",
        "status": "complete",
        "source": str(source),
        "destination": str(destination),
        "manifest_sha256": sha256(cloned_manifest),
        "fingerprint": state["current_state"]["fingerprint"],
        "copy_command": ["cp", "--archive", "--reflink=auto"],
    }
    output = Path(config["storage"]["results_root"]) / "06_manifests" / "02_bundle_clone_receipt.json"
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
