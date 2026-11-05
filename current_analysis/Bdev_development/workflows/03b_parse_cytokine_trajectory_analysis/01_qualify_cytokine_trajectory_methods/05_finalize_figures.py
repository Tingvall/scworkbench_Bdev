#!/usr/bin/env python3
"""Write a checksummed addendum manifest for the step 03b method figures."""

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
    final_run_manifest = root / "06_manifests" / "03_final_run_manifest.json"
    output = root / "06_manifests" / "04_figure_manifest.json"
    if output.exists():
        raise RuntimeError(f"Refusing to overwrite figure manifest: {output}")

    expected = [
        root / "02_figures" / f"{stem}.{suffix}"
        for stem in (
            "01_dynamics_ucell_heatmap",
            "02_signal_interval_map",
            "03_preranked_gsea_heatmap",
        )
        for suffix in ("png", "pdf")
    ]
    missing = [str(path) for path in expected if not path.is_file()]
    if missing:
        raise RuntimeError(f"Missing expected figures: {missing}")

    workflow_dir = Path(__file__).resolve().parent
    source_tables = [
        root / "01_tables" / name
        for name in (
            "09_dynamics_signal_table.tsv",
            "10_dynamics_dynamics_table.tsv",
            "14_intervals_global_test.tsv",
            "16_intervals_signal_intervals.tsv",
            "19_gsea_global_enrichment.tsv",
            "20_gsea_window_enrichment.tsv",
        )
    ]
    manifest = {
        "schema_version": "bdev.parse38_trajectory_figure_manifest.v1",
        "status": "complete_development_figures",
        "run_identity": config["run_identity"],
        "parent_run_manifest": {
            "path": str(final_run_manifest),
            "sha256": sha256(final_run_manifest),
        },
        "plot_script": {
            "path": str(workflow_dir / "04_plot_cytokine_trajectory_methods.R"),
            "sha256": sha256(workflow_dir / "04_plot_cytokine_trajectory_methods.R"),
        },
        "source_tables": [
            {"path": str(path.relative_to(root)), "sha256": sha256(path)}
            for path in source_tables
        ],
        "figures": [
            {
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in expected
        ],
        "source_bundle_mutated": False,
        "outputs_promoted": False,
        "implementation_unlocked": False,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"figure_manifest_sha256={sha256(output)}")


if __name__ == "__main__":
    main()
