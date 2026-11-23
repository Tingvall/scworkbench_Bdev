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
    report_path = root / "03_reports" / "01_immgen_bulk_gated_20d_mapping.json"
    prepare_path = root / "06_manifests" / "01_prepare_receipt.json"
    execution_path = root / "06_manifests" / "02_execution_receipt.json"
    required = [
        report_path,
        prepare_path,
        execution_path,
        root / "01_tables" / "01_bulk_query_slot_map.tsv",
        root / "01_tables" / "02_immgen_bulk_gated_20d_mapping.tsv",
        root / "01_tables" / "03_immgen_bulk_gated_20d_neighbors.tsv.gz",
        root / "01_tables" / "04_immgen_bulk_fixed_umap_placements.tsv",
        root / "01_tables" / "05_qualified_support_thresholds.tsv",
        root / "02_figures" / "01_immgen_bulk_gated_20d_fixed_atac_umap.png",
        root / "02_figures" / "01_immgen_bulk_gated_20d_fixed_atac_umap.pdf",
    ]
    for path in required:
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"missing or empty output: {path}")

    report = json.loads(report_path.read_text())
    execution = json.loads(execution_path.read_text())
    if report.get("status") != "passed" or not report.get("gate_passed"):
        raise RuntimeError("bulk mapping report did not pass")
    if report.get("query_libraries") != 28 or report.get("raw_estimate_rows") != 28:
        raise RuntimeError("bulk query accounting is incomplete")
    if report.get("fixed_umap_placements") != 28:
        raise RuntimeError("fixed-UMAP placement accounting is incomplete")
    if report.get("latent_dimensions") != 20 or report.get("selected_features") != 20000:
        raise RuntimeError("selected model dimensions changed")
    if report.get("model_refit") or report.get("umap_refit"):
        raise RuntimeError("a prohibited refit was reported")
    if not execution.get("qualified_source_manifest_unchanged"):
        raise RuntimeError("qualified source manifest changed")
    if not execution.get("bulk_source_manifest_unchanged"):
        raise RuntimeError("bulk source manifest changed")
    if not execution.get("fixed_umap_source_manifest_unchanged"):
        raise RuntimeError("fixed-UMAP source manifest changed")
    if not execution.get("model_payload_unchanged"):
        raise RuntimeError("model payload changed")

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
        "schema_version": "bdev.immgen_bulk_gated_20d_inventory.v1",
        "status": "complete_passed",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "report_sha256": sha256(report_path),
        "prepare_receipt_sha256": sha256(prepare_path),
        "execution_receipt_sha256": sha256(execution_path),
        "implementation_unlocked": False,
        "paper_analysis": False,
        "sif_built": False,
        "files": files,
    }
    inventory_path.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
