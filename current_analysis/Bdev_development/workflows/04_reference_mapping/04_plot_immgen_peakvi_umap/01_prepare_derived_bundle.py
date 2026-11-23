#!/usr/bin/env python3
"""Create a copy-on-write bundle for the derived ImmGen PeakVI UMAP."""

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


def merge(base: dict, override: dict) -> dict:
    result = dict(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(path: Path) -> dict:
    current = yaml.safe_load(path.read_text())
    base_name = current.pop("base_config", None)
    if base_name is None:
        return current
    return merge(load_config((path.parent / base_name).resolve()), current)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = load_config(config_path)

    results = Path(config["storage"]["results_root"])
    for category in ("01_tables", "02_figures", "03_reports", "04_logs", "05_objects", "06_manifests"):
        (results / category).mkdir(parents=True, exist_ok=True)
    Path(config["storage"]["sessions"]).mkdir(parents=True, exist_ok=True)

    source = Path(config["inputs"]["mapping_bundle"]["path"])
    source_manifest = source / "manifest.json"
    expected = config["inputs"]["mapping_bundle"]["manifest_sha256"]
    if sha256(source_manifest) != expected:
        raise RuntimeError("Completed ImmGen mapping-bundle manifest checksum changed")
    mapping_table = Path(config["inputs"]["mapping_table"]["path"])
    if sha256(mapping_table) != config["inputs"]["mapping_table"]["sha256"]:
        raise RuntimeError("ImmGen native-L4 mapping table checksum changed")
    l4_payload = Path(config["inputs"]["authenticated_l4"]["path"])
    if sha256(l4_payload) != config["inputs"]["authenticated_l4"]["sha256"]:
        raise RuntimeError("Authenticated fitted-L4 payload checksum changed")

    destination = Path(config["storage"]["derived_bundle"])
    if destination.exists():
        raise RuntimeError(f"Refusing to overwrite derived bundle: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["cp", "--recursive", "--reflink=auto", "--", str(source), str(destination)],
        check=True,
    )
    if sha256(destination / "manifest.json") != expected:
        raise RuntimeError("Derived-bundle manifest differs before UMAP execution")

    receipt = {
        "schema_version": "bdev.immgen_peakvi_umap_preparation.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "source_bundle": str(source),
        "source_manifest_sha256": expected,
        "derived_bundle": str(destination),
        "derived_manifest_sha256_before": sha256(destination / "manifest.json"),
        "mapping_table_sha256": sha256(mapping_table),
        "authenticated_l4_sha256": sha256(l4_payload),
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "source_mutated": False,
        "paper_analysis": False,
    }
    output = results / "06_manifests" / "01_preparation_receipt.json"
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
