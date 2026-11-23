#!/usr/bin/env python3
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
    config_path = parser.parse_args().config.resolve()
    config = yaml.safe_load(config_path.read_text())
    root = Path(config["storage"]["results_root"])
    report_path = root / "03_reports" / "01_gated_continuous_transfer_qualification.json"
    runtime_path = root / "06_manifests" / "01_runtime_and_clone_receipt.json"
    execution_path = root / "06_manifests" / "02_execution_receipt.json"
    required_tables = [
        "01_reference_fold_stratum_counts.tsv.gz",
        "02_reference_fold_sample_domain_counts.tsv",
        "03_balanced_fold_heldout_predictions.tsv.gz",
        "04_calibration_metrics.tsv",
        "05_support_thresholds.tsv",
        "06_calibration_coverage.tsv",
        "07_gated_fixture_mapping.tsv.gz",
        "08_continuous_acceptance.tsv",
    ]
    for path in [report_path, runtime_path, execution_path] + [
            root / "01_tables" / name for name in required_tables]:
        if not path.exists():
            raise RuntimeError(f"missing qualification artifact: {path}")
    report = json.loads(report_path.read_text())
    execution = json.loads(execution_path.read_text())
    if not report.get("gate_passed") or report.get("status") != "passed":
        raise RuntimeError("gated continuous-transfer qualification did not pass")
    if report.get("latent_dimensions") != 20:
        raise RuntimeError("qualification did not use 20 latent dimensions")
    if report.get("fold_count") != 5:
        raise RuntimeError("qualification did not use five folds")
    if report.get("raw_estimate_rows") != report.get("query_cells"):
        raise RuntimeError("raw estimate accounting is incomplete")
    if not execution.get("source_manifest_unchanged"):
        raise RuntimeError("qualification parent manifest changed")

    inventory_path = root / "06_manifests" / "03_run_inventory.json"
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path != inventory_path:
            files.append({
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    inventory = {
        "schema_version": "bdev.gated_continuous_transfer_inventory.v1",
        "status": "complete_passed",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "qualification_report_sha256": sha256(report_path),
        "runtime_receipt_sha256": sha256(runtime_path),
        "execution_receipt_sha256": sha256(execution_path),
        "implementation_unlocked": False,
        "full_query_mapping_performed": False,
        "sif_built": False,
        "files": files,
    }
    inventory_path.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
