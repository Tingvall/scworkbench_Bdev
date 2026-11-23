#!/usr/bin/env python3
"""Write the immutable file inventory for the completed 03c run."""

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
    output = root / "06_manifests" / "03_run_file_inventory.json"
    files = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p != output):
        files.append(
            {
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    required = {
        "03_reports/01_depth_adjusted_qualification_summary.json",
        "01_tables/10_depth_adjusted_scores.tsv",
        "01_tables/11_depth_adjustment_diagnostics.tsv",
        "01_tables/29_method_comparison.tsv",
        "02_figures/01_depth_adjusted_dynamics_heatmap.pdf",
        "02_figures/02_depth_adjusted_signal_interval_map.pdf",
        "02_figures/03_sct_preranked_gsea_heatmap.pdf",
        "02_figures/04_flt3l_up_depth_adjusted_four_view_summary.pdf",
        "02_figures/05_il7_up_depth_adjusted_four_view_summary.pdf",
    }
    observed = {entry["path"] for entry in files}
    missing = sorted(required - observed)
    if missing:
        raise RuntimeError(f"Required 03c outputs are missing: {missing}")
    manifest = {
        "schema_version": "bdev.parse38_depth_adjusted_run_inventory.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "files": files,
        "file_count": len(files),
        "paper_analysis": False,
        "implementation_unlocked": False,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "complete", "files": len(files), "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
