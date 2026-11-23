#!/usr/bin/env python3
"""Inventory the completed Bdev-only full-reference qualification."""

from __future__ import annotations

import argparse
import hashlib
import json
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
    root = Path(config["storage"]["results_root"])
    output = root / "06_manifests" / "04_run_file_inventory.json"
    required = {
        "01_tables/01_peak_scores.tsv.gz",
        "01_tables/02_hybrid20k_peak_manifest.tsv",
        "01_tables/05_candidate_metrics.tsv",
        "03_reports/01_feature_selection_summary.json",
        "03_reports/02_candidate_qualification_summary.json",
        "06_manifests/01_runtime_install.json",
        "06_manifests/02_feature_selection_receipt.json",
    }
    qualification_path = root / "03_reports" / "02_candidate_qualification_summary.json"
    qualification = json.loads(qualification_path.read_text())
    if qualification["status"] == "passed":
        required |= {
            "03_reports/03_frozen_reference_manifest.json",
            "06_manifests/03_freeze_execution_receipt.json",
        }
    files = [
        {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(root.rglob("*")) if path.is_file() and path != output
    ]
    missing = sorted(required - {entry["path"] for entry in files})
    if missing:
        raise RuntimeError(f"Required outputs are missing: {missing}")
    manifest = {
        "schema_version": "bdev.full_peakvi_run_inventory.v1",
        "status": "complete" if qualification["status"] == "passed" else "complete_gate_failed",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "files": files,
        "file_count": len(files),
        "external_query_use": False,
        "real_immgen_use": False,
        "paper_analysis": False,
        "implementation_unlocked": False,
        "sif_built": False,
        "candidate_qualification_status": qualification["status"],
        "selected_candidate": qualification.get("selected_candidate"),
        "frozen_reference_created": qualification["status"] == "passed",
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "complete", "files": len(files), "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
