#!/usr/bin/env python3
"""Inventory the completed 03c up-signature GSEA UMAP visualization run."""

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
    output = root / "06_manifests" / "02_run_file_inventory.json"
    files = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p != output):
        files.append({
            "path": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    figures = root / "02_figures"
    pngs = sorted(figures.glob("*.png"))
    pdfs = sorted(figures.glob("*.pdf"))
    expected = int(config["acceptance"]["expected_up_signatures"])
    if len(pngs) != expected or len(pdfs) != expected:
        raise RuntimeError(
            f"Expected {expected} PNG/PDF pairs; observed {len(pngs)}/{len(pdfs)}"
        )
    required = {
        "01_tables/01_up_geneset_gsea_umap_index.tsv",
        "01_tables/02_up_geneset_gsea_cell_assignments.tsv",
        "03_reports/01_up_geneset_gsea_umap_summary.json",
        "04_logs/01_up_geneset_gsea_umaps.log",
        "06_manifests/01_plotting_receipt.json",
    }
    observed = {entry["path"] for entry in files}
    missing = sorted(required - observed)
    if missing:
        raise RuntimeError(f"Required visualization outputs are missing: {missing}")
    manifest = {
        "schema_version": "bdev.parse38_up_gsea_umaps_inventory.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "parent_run_identity": config["scope"]["parent_run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "files": files,
        "file_count": len(files),
        "figure_pairs": expected,
        "parent_run_mutated": False,
        "paper_analysis": False,
        "implementation_unlocked": False,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": "complete",
        "files": len(files),
        "figure_pairs": expected,
        "output": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
