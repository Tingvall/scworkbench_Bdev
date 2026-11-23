#!/usr/bin/env python3
"""Create the immutable result-file inventory for ImmGen trajectory mapping."""

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
    output = root / "06_manifests" / "06_run_file_inventory.json"
    files = [
        {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(root.rglob("*")) if path.is_file() and path != output
    ]
    required = {
        "01_tables/01_immgen_library_run_manifest.tsv",
        "01_tables/02_bdev_exact_peak_axis.bed",
        "01_tables/03_immgen_preprocessing_qc.tsv",
        "01_tables/04_immgen_library_l4_mapping.tsv",
        "01_tables/05_immgen_population_l4_summary.tsv",
        "02_figures/01_immgen_library_native_l4.pdf",
        "02_figures/01_immgen_library_native_l4.png",
        "02_figures/02_immgen_reference_support.pdf",
        "02_figures/02_immgen_reference_support.png",
        "03_reports/01_immgen_atac_trajectory_mapping_summary.json",
        "06_manifests/01_runtime_install.json",
        "06_manifests/02_ena_resolution.json",
        "06_manifests/03_peak_axis_export.json",
        "06_manifests/04_preprocessing_receipt.json",
        "06_manifests/05_execution_receipt.json",
    }
    observed = {entry["path"] for entry in files}
    missing = sorted(required - observed)
    if missing:
        raise RuntimeError(f"Required outputs are missing: {missing}")
    manifest = {
        "schema_version": "bdev.immgen_atac_trajectory_run_inventory.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "file_count": len(files),
        "files": files,
        "real_immgen_projected": True,
        "categorical_cluster_transfer": False,
        "paper_analysis": False,
        "implementation_unlocked": False,
        "sif_built": False,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "complete", "files": len(files), "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
