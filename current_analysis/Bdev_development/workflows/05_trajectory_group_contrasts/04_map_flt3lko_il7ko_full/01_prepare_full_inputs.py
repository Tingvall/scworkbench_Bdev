#!/usr/bin/env python3
"""Freeze all QC-selected KO cells and the exact frozen HVP-depth20k axis."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import yaml


EXPECTED = {
    "FLT3LKO_BM_WT": ("Flt3l", "WT"),
    "FLT3LKO_BM_KO": ("Flt3l", "KO"),
    "IL7KO_BM_WT": ("Il7", "WT"),
    "IL7KO_BM_KO": ("Il7", "KO"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_hash(path: Path, expected: str) -> None:
    observed = sha256(path)
    if observed != expected:
        raise RuntimeError(f"Checksum mismatch for {path}: {observed}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    cfg = yaml.safe_load(config_path.read_text())
    root = Path(cfg["storage"]["results_root"])
    if root.exists():
        raise RuntimeError(f"Refusing to overwrite result root: {root}")
    for category in ("01_tables", "02_figures", "03_reports", "04_logs",
                     "05_objects", "06_manifests"):
        (root / category).mkdir(parents=True, exist_ok=False)

    inputs = cfg["inputs"]
    for key in ("qc_metadata", "selected_axis", "reciprocal_liftover_map"):
        require_hash(Path(inputs[key]["path"]), inputs[key]["sha256"])
    require_hash(Path(cfg["runtime"]["base_sif"]),
                 cfg["runtime"]["base_sif_sha256"])
    require_hash(Path(inputs["frozen_model_bundle"]["path"]) / "manifest.json",
                 inputs["frozen_model_bundle"]["manifest_sha256"])
    require_hash(Path(inputs["fixed_bdev_bundle"]["path"]) / "manifest.json",
                 inputs["fixed_bdev_bundle"]["manifest_sha256"])

    cells = []
    with gzip.open(inputs["qc_metadata"]["path"], "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["qc_use_downstream"].upper() != "TRUE":
                continue
            library = row["library_id"]
            if library not in EXPECTED:
                raise RuntimeError(f"Unexpected selected library: {library}")
            if (row["knockout_target"], row["genotype"]) != EXPECTED[library]:
                raise RuntimeError(f"Identity mismatch for {row['cell_id']}")
            cells.append(row)
    cells.sort(key=lambda row: (list(EXPECTED).index(row["library_id"]),
                                row["cell_id"]))
    if len(cells) != cfg["scope"]["expected_query_cells"]:
        raise RuntimeError(f"Expected 23,252 selected cells, found {len(cells)}")
    if len({row["cell_id"] for row in cells}) != len(cells):
        raise RuntimeError("QC-selected cell IDs are not unique")
    counts = Counter(row["library_id"] for row in cells)

    with gzip.open(inputs["selected_axis"]["path"], "rt", newline="") as handle:
        axis = [row for row in csv.DictReader(handle, delimiter="\t")
                if row["candidate"] == inputs["selected_axis"]["candidate"]]
    if len(axis) != inputs["selected_axis"]["expected_features"]:
        raise RuntimeError(f"Expected 20,000 selected features, found {len(axis)}")
    if [int(row["feature_order"]) for row in axis] != list(range(1, 20001)):
        raise RuntimeError("Frozen feature order is not exact")

    cell_columns = [
        "cell_id", "library_id", "knockout_target", "genotype", "tissue",
        "rna_umi", "detected_genes", "mitochondrial_percent",
        "atac_fragments", "tss_enrichment", "nucleosome_signal",
    ]
    cells_path = root / "01_tables" / "01_all_qc_selected_cells.tsv.gz"
    with gzip.open(cells_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cell_columns, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(cells)
    axis_path = root / "01_tables" / "02_hvp_depth20k_axis.tsv.gz"
    with gzip.open(axis_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle,
                                fieldnames=["feature_order", "feature_id", "selection_rank"],
                                delimiter="\t", extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(axis)

    report = {
        "schema_version": "bdev.flt3lko_il7ko_full_input_audit.v1",
        "status": "passed",
        "run_identity": cfg["run_identity"],
        "query_cells": len(cells),
        "library_counts": dict(counts),
        "features": len(axis),
        "model_refit": False,
        "umap_refit": False,
        "paper_analysis": False,
        "sif_built": False,
    }
    (root / "03_reports" / "01_full_input_audit.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    receipt = {
        **report,
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "outputs": {
            "cells": {"path": str(cells_path), "sha256": sha256(cells_path)},
            "axis": {"path": str(axis_path), "sha256": sha256(axis_path)},
        },
    }
    (root / "06_manifests" / "01_input_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
