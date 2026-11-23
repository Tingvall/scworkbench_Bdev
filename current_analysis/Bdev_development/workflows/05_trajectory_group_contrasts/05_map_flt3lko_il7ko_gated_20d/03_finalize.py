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
    report_path = root / "03_reports" / "01_full_gated_20d_mapping.json"
    runtime_path = root / "06_manifests" / "01_runtime_and_clone_receipt.json"
    execution_path = root / "06_manifests" / "02_execution_receipt.json"
    required = [report_path, runtime_path, execution_path,
        root / "01_tables" / "01_full_query_slot_map.tsv.gz",
        root / "01_tables" / "02_full_gated_20d_mapping.tsv.gz",
        root / "01_tables" / "03_full_gated_20d_group_counts.tsv",
        root / "01_tables" / "04_qualified_support_thresholds.tsv"]
    for path in required:
        if not path.exists():
            raise RuntimeError(f"missing full-mapping artifact: {path}")
    report = json.loads(report_path.read_text())
    execution = json.loads(execution_path.read_text())
    if report.get("status") != "passed" or not report.get("gate_passed"):
        raise RuntimeError("full gated 20D mapping did not pass")
    if report.get("query_cells") != config["scope"]["expected_query_cells"]:
        raise RuntimeError("full query accounting is incomplete")
    if report.get("raw_estimate_rows") != report.get("query_cells"):
        raise RuntimeError("raw estimate accounting is incomplete")
    if report.get("latent_dimensions") != 20 or report.get("model_refit"):
        raise RuntimeError("model contract was violated")
    if not execution.get("qualified_source_manifest_unchanged"):
        raise RuntimeError("qualified source bundle changed")

    inventory_path = root / "06_manifests" / "03_run_inventory.json"
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path != inventory_path:
            files.append({"path": str(path.relative_to(root)),
                "bytes": path.stat().st_size, "sha256": sha256(path)})
    inventory = {
        "schema_version": "bdev.flt3lko_il7ko_full_gated_20d_inventory.v1",
        "status": "complete_passed", "run_identity": config["run_identity"],
        "config_path": str(config_path), "config_sha256": sha256(config_path),
        "report_sha256": sha256(report_path),
        "runtime_receipt_sha256": sha256(runtime_path),
        "execution_receipt_sha256": sha256(execution_path),
        "implementation_unlocked": False, "genotype_contrast_performed": False,
        "sif_built": False, "files": files}
    inventory_path.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
