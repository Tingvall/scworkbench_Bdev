#!/usr/bin/env python3
"""Hash the completed bounded-labeling evidence and freeze its final manifest."""

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
    files = []
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        if path == output:
            continue
        files.append(
            {
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    archive_manifest = Path(inputs["bdev_reference"]["bundle_path"]) / "manifest.json"
    summary = json.loads((root / "03_reports" / "02_labeling_summary.json").read_text())
    workflow_dir = Path(__file__).resolve().parent
    workflow_files = [
        {
            "path": str(path),
            "sha256": sha256_file(path),
        }
        for path in sorted(workflow_dir.iterdir())
        if path.is_file() and path.name != "04_final_run_manifest.json"
    ]
    manifest = {
        "schema_version": "bdev.immgen_bdev_labeling_final_manifest.v1",
        "status": "complete_flat_labeling_not_accepted",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256_file(config_path)},
        "archive": {
            "path": str(Path(inputs["bdev_reference"]["bundle_path"])),
            "manifest_sha256_expected": inputs["bdev_reference"]["manifest_sha256"],
            "manifest_sha256_observed_after_run": sha256_file(archive_manifest),
            "mutated": False,
        },
        "execution_clone": {
            "path": config["storage"]["execution_clone"],
            "final_bundle_fingerprint": summary["bundle_fingerprint_after_promotion"],
        },
        "result": {
            "collection_sets": summary["collection_sets"],
            "scored_sets": summary["scored_sets"],
            "fixture_cells": summary["fixture_cells"],
            "assigned_cells": summary["assigned_cells"],
            "unassigned_cells": summary["unassigned_cells"],
            "exact_pt_preserved": summary["exact_pt_preserved"],
            "accepted_labels_emitted": False,
        },
        "workflow_files": workflow_files,
        "files": files,
    }
    if manifest["archive"]["manifest_sha256_observed_after_run"] != manifest[
        "archive"
    ]["manifest_sha256_expected"]:
        raise RuntimeError("Archived bundle manifest changed during the run")
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
