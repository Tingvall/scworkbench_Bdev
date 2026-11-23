#!/usr/bin/env python3
"""Validate and inventory the read-only frozen-model visualization run."""

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

    summary_path = root / "03_reports" / "01_fixed_umap_comparison_summary.json"
    receipt_path = root / "06_manifests" / "01_execution_receipt.json"
    summary = json.loads(summary_path.read_text())
    receipt = json.loads(receipt_path.read_text())
    if summary["status"] != "complete_read_only_visualization":
        raise RuntimeError("visualization summary is not complete")
    if not summary["source_inputs_unchanged"] or summary["model_refit"]:
        raise RuntimeError("read-only visualization contract was not satisfied")
    if summary["trajectory_output_status"] != "diagnostic_unaccepted":
        raise RuntimeError("trajectory output is not marked diagnostic/unaccepted")
    if receipt["config_sha256"] != sha256(config_path):
        raise RuntimeError("execution receipt config checksum mismatch")
    workflow_path = Path(receipt["workflow_path"])
    if receipt["workflow_sha256"] != sha256(workflow_path):
        raise RuntimeError("execution receipt workflow checksum mismatch")

    required = [
        "02_figures/01_four_model_fixed_atac_umap_overview.png",
        "02_figures/01_four_model_fixed_atac_umap_overview.pdf",
        "02_figures/02_hybrid20k_fixed_atac_umap_detail.png",
        "02_figures/03_cluster20k_fixed_atac_umap_detail.png",
        "02_figures/04_overall20k_fixed_atac_umap_detail.png",
        "02_figures/05_full191433_fixed_atac_umap_detail.png",
        "02_figures/06_four_model_diagnostic_native_l4.png",
        "02_figures/06_four_model_diagnostic_native_l4.pdf",
        "01_tables/01_fixed_umap_candidate_placements.tsv",
        "01_tables/03_authenticated_native_l4_curve.tsv",
    ]
    missing = [relative for relative in required if not (root / relative).is_file()]
    if missing:
        raise RuntimeError(f"missing required visualization outputs: {missing}")

    inventory_path = root / "06_manifests" / "02_run_file_inventory.json"
    records = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p != inventory_path):
        records.append({
            "path": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    inventory = {
        "schema_version": "bdev.frozen_peakvi_fixed_umap_inventory.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "files": records,
    }
    inventory_path.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
