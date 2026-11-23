#!/usr/bin/env python3
"""Prepare the isolated bundle and frozen Parse38 GMT for 03c."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
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
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())

    source = Path(config["inputs"]["source_bundle"]["path"])
    source_manifest = source / "manifest.json"
    if sha256(source_manifest) != config["inputs"]["source_bundle"]["manifest_sha256"]:
        raise RuntimeError("Archived source-bundle manifest checksum changed")

    gene_sets = Path(config["inputs"]["gene_sets"]["path"])
    if sha256(gene_sets) != config["inputs"]["gene_sets"]["sha256"]:
        raise RuntimeError("Parse38 GMT checksum changed")
    lines = [line for line in gene_sets.read_text().splitlines() if line]
    memberships = sum(max(0, len(line.split("\t")) - 2) for line in lines)
    if len(lines) != int(config["inputs"]["gene_sets"]["expected_sets"]):
        raise RuntimeError("Parse38 GMT set count changed")
    if memberships != int(config["inputs"]["gene_sets"]["expected_memberships"]):
        raise RuntimeError("Parse38 GMT membership count changed")

    results = Path(config["storage"]["results_root"])
    for name in ("01_tables", "02_figures", "03_reports", "04_logs", "05_objects", "06_manifests"):
        (results / name).mkdir(parents=True, exist_ok=True)
    frozen_gmt = results / "01_tables" / "01_parse38_cytokine_gene_sets.gmt"
    if frozen_gmt.exists():
        if sha256(frozen_gmt) != sha256(gene_sets):
            raise RuntimeError("Existing frozen GMT differs from declared input")
    else:
        shutil.copyfile(gene_sets, frozen_gmt)

    destination = Path(config["storage"]["execution_clone"])
    if destination.exists():
        raise RuntimeError(f"Refusing to overwrite execution clone: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    Path(config["storage"]["sessions"]).mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["cp", "--recursive", "--reflink=auto", "--", str(source), str(destination)],
        check=True,
    )
    if sha256(destination / "manifest.json") != sha256(source_manifest):
        raise RuntimeError("Execution-clone manifest differs before migration")

    receipt = {
        "schema_version": "bdev.parse38_depth_adjusted_preparation.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "source_bundle": str(source),
        "source_manifest_sha256": sha256(source_manifest),
        "execution_clone": str(destination),
        "frozen_gmt": str(frozen_gmt),
        "frozen_gmt_sha256": sha256(frozen_gmt),
        "gene_sets": len(lines),
        "memberships": memberships,
        "archive_mutated": False,
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
    }
    output = results / "06_manifests" / "01_preparation_receipt.json"
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
