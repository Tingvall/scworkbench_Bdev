#!/usr/bin/env python3
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
    config_path = parser.parse_args().config.resolve()
    config = yaml.safe_load(config_path.read_text())
    sif = Path(config["build"]["output_sif"])
    apptainer = config["build"]["apptainer"]
    inspected = subprocess.check_output(
        [apptainer, "inspect", "--json", str(sif)], text=True)
    labels = json.loads(inspected)["data"]["attributes"]["labels"]
    expected = {
        "CoreCommit": "4aca71d1c75b4e1532911c7ce68d945f2a319ca8",
        "Profile": "bdev-development-20261118",
        "Stack": "bdev-qualified-eight-pack-candidate",
        "BaseSIFSHA256": config["base_sif"]["sha256"],
    }
    observed = {key: labels.get(key) for key in expected}
    mismatches = {key: {"expected": expected[key], "observed": observed[key]}
                  for key in expected if observed[key] != expected[key]}
    if not mismatches:
        raise RuntimeError("candidate labels unexpectedly match; rejection is not justified")
    root = Path(config["build"]["results_root"])
    report = {
        "schema_version": "bdev.candidate_sif_rejection.v1",
        "status": "rejected",
        "run_identity": config["run_identity"],
        "reason": "inherited_project_labels_not_replaced",
        "sif_path": str(sif),
        "sif_bytes": sif.stat().st_size,
        "sif_sha256": sha256(sif),
        "label_mismatches": mismatches,
        "scientific_content_qualification_run": False,
        "implementation_unlocked": False,
    }
    report_path = root / "03_reports" / "01_stale_inherited_labels_rejection.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    receipt = {
        "schema_version": "bdev.candidate_sif_rejection_receipt.v1",
        "status": "complete_rejected",
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "rejection_report_sha256": sha256(report_path),
        "sif_sha256": report["sif_sha256"],
    }
    receipt_path = root / "06_manifests" / "02_rejected_candidate_receipt.json"
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

