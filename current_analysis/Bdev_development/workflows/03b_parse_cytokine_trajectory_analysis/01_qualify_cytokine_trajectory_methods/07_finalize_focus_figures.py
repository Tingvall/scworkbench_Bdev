#!/usr/bin/env python3
"""Write a checksummed addendum for the FLT3L-up and IL-7-up figures."""

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
    config = yaml.safe_load(args.config.resolve().read_text())
    root = Path(config["storage"]["results_root"])
    output = root / "06_manifests" / "05_focus_figure_manifest.json"
    if output.exists():
        raise RuntimeError(f"Refusing to overwrite focus figure manifest: {output}")

    stems = (
        "04_flt3l_up_four_method_summary",
        "05_il7_up_four_method_summary",
    )
    figures = [
        root / "02_figures" / f"{stem}.{suffix}"
        for stem in stems
        for suffix in ("png", "pdf")
    ]
    missing = [str(path) for path in figures if not path.is_file()]
    if missing:
        raise RuntimeError(f"Missing focused figures: {missing}")

    workflow_dir = Path(__file__).resolve().parent
    score_payload = (
        Path(config["storage"]["sessions"])
        / "parse38_cytokine_trajectory_qualification"
        / "stage/operation_pack_v2/transition_gwdcodtyyl9dkv1x"
        / "outputs/signature_scores/payload.rds"
    )
    source_tables = [
        root / "01_tables" / name
        for name in (
            "09_dynamics_signal_table.tsv",
            "10_dynamics_dynamics_table.tsv",
            "13_intervals_signal_curve.tsv",
            "14_intervals_global_test.tsv",
            "16_intervals_signal_intervals.tsv",
            "19_gsea_global_enrichment.tsv",
            "20_gsea_window_enrichment.tsv",
            "21_gsea_enrichment_intervals.tsv",
        )
    ]
    parent = root / "06_manifests" / "04_figure_manifest.json"
    script = workflow_dir / "06_plot_flt3l_il7_four_methods.R"
    bundle_manifest = Path(config["storage"]["execution_clone"]) / "manifest.json"
    manifest = {
        "schema_version": "bdev.parse38_trajectory_focus_figure_manifest.v1",
        "status": "complete_development_figures",
        "run_identity": config["run_identity"],
        "targets": ["FLT3L up", "IL-7 up"],
        "views": ["ATAC-UMAP UCell", "dynamics", "signal intervals", "preranked GSEA"],
        "parent_figure_manifest": {"path": str(parent), "sha256": sha256(parent)},
        "plot_script": {"path": str(script), "sha256": sha256(script)},
        "atac_umap_source": {
            "space_id": config["inputs"]["fit_space"],
            "bundle_manifest": str(bundle_manifest),
            "bundle_manifest_sha256": sha256(bundle_manifest),
            "read_only": True,
        },
        "ucell_score_payload": {"path": str(score_payload), "sha256": sha256(score_payload)},
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
            for path in figures
        ],
        "source_bundle_mutated": False,
        "outputs_promoted": False,
        "implementation_unlocked": False,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"focus_figure_manifest_sha256={sha256(output)}")


if __name__ == "__main__":
    main()
