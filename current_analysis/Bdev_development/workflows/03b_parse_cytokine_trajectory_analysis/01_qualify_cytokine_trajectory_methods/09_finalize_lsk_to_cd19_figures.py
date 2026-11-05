#!/usr/bin/env python3
"""Write a checksummed addendum for the LSK-to-CD19 focused figures."""

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
    output = root / "06_manifests" / "06_lsk_to_cd19_figure_manifest.json"
    if output.exists():
        raise RuntimeError(f"Refusing to overwrite direction figure manifest: {output}")

    stems = (
        "06_flt3l_up_lsk_to_cd19_four_method_summary",
        "07_il7_up_lsk_to_cd19_four_method_summary",
    )
    figures = [
        root / "02_figures" / f"{stem}.{suffix}"
        for stem in stems
        for suffix in ("png", "pdf")
    ]
    missing = [str(path) for path in figures if not path.is_file()]
    if missing:
        raise RuntimeError(f"Missing LSK-to-CD19 figures: {missing}")

    workflow_dir = Path(__file__).resolve().parent
    session_stage = (
        Path(config["storage"]["sessions"])
        / "parse38_cytokine_trajectory_qualification"
        / "stage/operation_pack_v2"
    )
    score_payload = (
        session_stage
        / "transition_gwdcodtyyl9dkv1x/outputs/signature_scores/payload.rds"
    )
    model_payload = (
        session_stage
        / "transition_7gdx4c4n3cn8iuti/outputs/model/model_tree/model.rds"
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
    parent = root / "06_manifests" / "05_focus_figure_manifest.json"
    script = workflow_dir / "08_plot_lsk_to_cd19_focus_figures.R"
    bundle_manifest = Path(config["storage"]["execution_clone"]) / "manifest.json"
    manifest = {
        "schema_version": "bdev.parse38_trajectory_direction_figure_manifest.v1",
        "status": "complete_development_figures",
        "run_identity": config["run_identity"],
        "targets": ["FLT3L up", "IL-7 up"],
        "views": [
            "ATAC-UMAP UCell with fitted L4 arrow",
            "dynamics",
            "signal intervals",
            "preranked GSEA",
        ],
        "display_semantics": {
            "analysis_trajectory": "authenticated fitted Slingshot L4 in native orientation",
            "native_cluster_path": ["5", "1", "2", "16", "13", "12", "7", "3"],
            "display_transform": "display_fraction = 1 - native_L4_fraction",
            "display_zero": "LSK",
            "display_one": "CD19",
            "directional_annotations": {
                "spearman_rho": "display_rho = -native_rho",
                "curve_class": "turns_on/off and increasing/decreasing are reversed",
                "global_monotonic_gsea_nes": "display_NES = -native_NES",
                "local_window_gsea_nes": "unchanged; only window positions are reversed",
            },
            "underlying_results_modified": False,
        },
        "parent_focus_figure_manifest": {"path": str(parent), "sha256": sha256(parent)},
        "plot_script": {"path": str(script), "sha256": sha256(script)},
        "atac_umap_source": {
            "space_id": config["inputs"]["fit_space"],
            "bundle_manifest": str(bundle_manifest),
            "bundle_manifest_sha256": sha256(bundle_manifest),
            "read_only": True,
        },
        "fitted_l4_model": {"path": str(model_payload), "sha256": sha256(model_payload)},
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
    print(f"direction_figure_manifest_sha256={sha256(output)}")


if __name__ == "__main__":
    main()
