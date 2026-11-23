#!/usr/bin/env python3
"""Inventory the derived ImmGen PeakVI UMAP result files."""

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
    root = Path(config["storage"]["results_root"])
    output = root / "06_manifests" / "03_run_file_inventory.json"
    files = [
        {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(root.rglob("*")) if path.is_file() and path != output
    ]
    required = {
        "01_tables/01_bdev_immgen_peakvi_umap_coordinates.tsv",
        "02_figures/01_immgen_population_peakvi_umap.pdf",
        "02_figures/01_immgen_population_peakvi_umap.png",
        "02_figures/02_native_l4_peakvi_umap.pdf",
        "02_figures/02_native_l4_peakvi_umap.png",
        "03_reports/01_immgen_peakvi_umap_summary.json",
        "04_logs/02_fit_transform_and_plot.log",
        "06_manifests/01_preparation_receipt.json",
        "06_manifests/02_execution_receipt.json",
    }
    observed = {entry["path"] for entry in files}
    missing = sorted(required - observed)
    if missing:
        raise RuntimeError(f"Required outputs are missing: {missing}")
    source_manifest = Path(config["inputs"]["mapping_bundle"]["path"]) / "manifest.json"
    if sha256(source_manifest) != config["inputs"]["mapping_bundle"]["manifest_sha256"]:
        raise RuntimeError("Source mapping bundle changed during derived UMAP run")
    design = config["umap"]["design"]
    inventory = {
        "schema_version": "bdev.immgen_peakvi_umap_inventory.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "file_count": len(files),
        "files": files,
        "design": design,
        "reference_only_fit": design == "reference_only_umap_fit_then_frozen_model_query_transform",
        "frozen_model_query_transform": design == "reference_only_umap_fit_then_frozen_model_query_transform",
        "authenticated_union_fit": design == "authenticated_union_fit",
        "source_mutated": False,
        "paper_analysis": False,
        "implementation_unlocked": False,
        "sif_built": False,
    }
    output.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "complete", "files": len(files), "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
