#!/usr/bin/env python3
"""Freeze the completed hierarchical qualification result manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
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
    config_path = Path(args.config).resolve()
    config = yaml.safe_load(config_path.read_text())
    inputs = yaml.safe_load(Path(config["inputs_manifest"]).read_text())
    root = Path(config["storage"]["results_root"])
    output = root / "06_manifests" / "04_final_run_manifest.json"
    summary = json.loads((root / "03_reports" / "02_hierarchical_summary.json").read_text())
    archive_manifest = Path(inputs["bdev_reference"]["archived_bundle_path"]) / "manifest.json"
    parent_manifest = Path(inputs["parent_execution"]["execution_clone"]) / "manifest.json"
    if sha256_file(archive_manifest) != inputs["bdev_reference"]["archived_manifest_sha256"]:
        raise RuntimeError("Archived bundle manifest changed")
    if sha256_file(parent_manifest) != inputs["parent_execution"]["manifest_sha256"]:
        raise RuntimeError("Parent step-15 bundle changed")
    files = []
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file() and candidate != output):
        files.append({"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    workflow_dir = Path(__file__).resolve().parent
    workflows = [{"path": str(path), "sha256": sha256_file(path)} for path in sorted(workflow_dir.iterdir()) if path.is_file()]
    manifest = {
        "schema_version": "bdev.immgen_bdev_hierarchical_final_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256_file(config_path)},
        "archive": {"manifest_sha256": sha256_file(archive_manifest), "mutated": False},
        "parent_bundle": {
            "path": inputs["parent_execution"]["execution_clone"],
            "fingerprint": inputs["parent_execution"]["bundle_fingerprint"],
            "manifest_sha256": sha256_file(parent_manifest),
            "mutated": False,
        },
        "execution": {
            "bundle": config["storage"]["session_bundle"],
            "fingerprint": summary["session_bundle_fingerprint"],
            "registry_unchanged": summary["session_bundle_registry_unchanged"],
            "scworkbench_artifacts_promoted": False,
        },
        "result": {
            "fixture_cells": summary["fixture_cells"],
            "broad_multi_state_cells": summary["broad_multi_state_cells"],
            "assigned_original": summary["assigned_original"],
            "assigned_shared_gene_filtered": summary["assigned_shared_gene_filtered"],
            "accepted_labels_emitted": False,
        },
        "workflow_files": workflows,
        "files": files,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
