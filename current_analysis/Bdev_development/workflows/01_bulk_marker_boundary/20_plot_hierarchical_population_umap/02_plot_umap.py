#!/usr/bin/env python3
"""Plot hierarchical population policies on the Bdev WNN UMAP."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/bdev_hierarchical_umap_mpl")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    require(config["schema_version"] == "bdev.hierarchical_population_umap.v1", "Unexpected config schema")
    inputs = {name: Path(spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        require(sha256(path) == config["inputs"][name]["sha256"], f"Checksum mismatch: {name}")

    root = Path(config["storage"]["results_root"])
    tables = root / "01_tables"
    figures = root / "02_figures"
    reports = root / "03_reports"
    manifests = root / "06_manifests"
    for directory in (figures, reports, manifests):
        directory.mkdir(parents=True, exist_ok=True)
    coordinate_path = tables / "01_bdev_umap_atac_coordinates.tsv.gz"
    receipt_path = manifests / "01_coordinate_export_receipt.json"
    require(coordinate_path.is_file() and receipt_path.is_file(), "Coordinate export is incomplete")

    coordinates = pd.read_csv(coordinate_path, sep="\t")
    decisions = pd.read_csv(inputs["decisions"], sep="\t")
    require(len(coordinates) == 1746 and coordinates["cell_id"].is_unique, "Unexpected coordinate inventory")
    require(len(decisions) == 1746 * 3, "Unexpected decision inventory")
    merged = decisions.merge(coordinates, on="cell_id", how="left", validate="many_to_one")
    require(merged[["UMAP_1", "UMAP_2"]].notna().all().all(), "Missing UMAP coordinates")

    status_colors = {"unassigned": "#CFCFCF", "broad_only": "#80B1D3", "fine_label": "#08519C"}
    status_names = {"unassigned": "Fully unassigned", "broad_only": "Broad state only", "fine_label": "Fine population call"}
    policies = config["plot"]["comparison_policies"]
    policy_titles = {
        "no_margin": "No margin",
        "fixed_margin_0.1": "Fixed margin = 0.10",
        "population_calibrated": "Population-specific margins",
    }
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4), sharex=True, sharey=True)
    for axis, policy in zip(axes, policies):
        data = merged[merged["policy"] == policy].copy()
        for status in ("unassigned", "broad_only", "fine_label"):
            part = data[data["decision"] == status]
            axis.scatter(part["UMAP_1"], part["UMAP_2"], s=7, c=status_colors[status],
                         alpha=0.82, linewidths=0, rasterized=True)
        counts = data["decision"].value_counts()
        axis.set_title(
            f"{policy_titles[policy]}\n"
            f"fine {counts.get('fine_label', 0):,} | broad only {counts.get('broad_only', 0):,} | unassigned {counts.get('unassigned', 0):,}",
            fontsize=11,
        )
        axis.set_xlabel("ATAC UMAP 1")
        axis.set_aspect("equal", adjustable="box")
        axis.spines[["top", "right", "left", "bottom"]].set_visible(False)
        axis.set_xticks([]); axis.set_yticks([])
    axes[0].set_ylabel("ATAC UMAP 2")
    legend = [Line2D([0], [0], marker="o", color="none", markerfacecolor=status_colors[key],
                     markeredgecolor="none", markersize=7, label=status_names[key])
              for key in ("fine_label", "broad_only", "unassigned")]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Bdev hierarchical ImmGen policy comparison", fontsize=15, y=0.99)
    fig.tight_layout(rect=(0, 0.055, 1, 0.96))
    status_png = figures / "01_bdev_hierarchical_policy_status_umap.png"
    status_pdf = figures / "01_bdev_hierarchical_policy_status_umap.pdf"
    fig.savefig(status_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(status_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    lineage_colors = {
        "B": "#1F78B4", "Stem&Prog": "#6A3D9A", "abT": "#FF7F00",
        "gdT": "#E31A1C", "ILC": "#33A02C", "myeloid": "#8C510A",
        "monocyte": "#E7298A", "DC": "#00A6A6", "stroma": "#7A9A01",
    }
    primary = merged[merged["policy"] == config["plot"]["primary_policy"]].copy()
    fig, axis = plt.subplots(figsize=(14.2, 8.8))
    fig.subplots_adjust(left=0.05, right=0.68, bottom=0.08, top=0.90)
    unassigned = primary[primary["decision"] == "unassigned"]
    axis.scatter(unassigned["UMAP_1"], unassigned["UMAP_2"], s=8, c="#D0D0D0",
                 alpha=0.65, linewidths=0, rasterized=True)
    broad_only = primary[primary["decision"] == "broad_only"]
    for lineage, part in broad_only.groupby("broad_top_state"):
        axis.scatter(part["UMAP_1"], part["UMAP_2"], s=8,
                     c=lineage_colors.get(lineage, "#999999"), alpha=0.18,
                     linewidths=0, rasterized=True)
    fine = primary[primary["decision"] == "fine_label"]
    for lineage, part in fine.groupby("winner_lineage"):
        axis.scatter(part["UMAP_1"], part["UMAP_2"], s=12,
                     c=lineage_colors.get(lineage, "#555555"), alpha=0.92,
                     linewidths=0, rasterized=True)
    minimum_label_cells = int(config["plot"]["minimum_cells_for_direct_population_label"])
    label_rows = []
    for label, part in fine.groupby("fine_label"):
        if len(part) >= minimum_label_cells:
            label_rows.append((label, len(part), float(part["UMAP_1"].median()), float(part["UMAP_2"].median())))
    numbered_labels = sorted(label_rows, key=lambda item: (str(fine.loc[fine["fine_label"] == item[0], "winner_lineage"].iloc[0]), item[0]))
    for number, (label, count, x, y) in enumerate(numbered_labels, 1):
        axis.text(x, y, str(number), fontsize=6.4, fontweight="bold", ha="center", va="center",
                  color="#111111", zorder=10,
                  bbox={"boxstyle": "circle,pad=0.16", "facecolor": "white", "edgecolor": "#333333", "alpha": 0.88, "linewidth": 0.45})
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=color,
                      markeredgecolor="none", markersize=7, label=lineage)
               for lineage, color in lineage_colors.items()]
    handles.extend([
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#777777", alpha=0.25,
               markeredgecolor="none", markersize=7, label="Broad-state only (pale)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#D0D0D0",
               markeredgecolor="none", markersize=7, label="Fully unassigned"),
    ])
    axis.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
                fontsize=8, title="Lineage / status", title_fontsize=9)
    key_lines = [f"{number:>2}  {label} ({count})" for number, (label, count, _, _) in enumerate(numbered_labels, 1)]
    split = (len(key_lines) + 1) // 2
    fig.text(0.715, 0.53, "Numbered population labels", fontsize=9, fontweight="bold", ha="left")
    fig.text(0.715, 0.50, "\n".join(key_lines[:split]), fontsize=7.6, family="monospace", ha="left", va="top")
    fig.text(0.855, 0.50, "\n".join(key_lines[split:]), fontsize=7.6, family="monospace", ha="left", va="top")
    axis.set_title(
        "Bdev population-specific calibrated calls\n"
        "385 provisional fine calls; numbered populations have at least 5 cells",
        fontsize=14,
    )
    axis.set_xlabel("ATAC UMAP 1"); axis.set_ylabel("ATAC UMAP 2")
    axis.set_xticks([]); axis.set_yticks([])
    axis.set_aspect("equal", adjustable="datalim")
    axis.spines[["top", "right", "left", "bottom"]].set_visible(False)
    label_png = figures / "02_bdev_population_calibrated_labels_umap.png"
    label_pdf = figures / "02_bdev_population_calibrated_labels_umap.pdf"
    fig.savefig(label_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(label_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    wide = decisions.pivot(index="cell_id", columns="policy", values="decision")
    transition = pd.crosstab(wide["fixed_margin_0.1"], wide["population_calibrated"])
    transition_path = tables / "02_fixed_vs_calibrated_status_transitions.tsv"
    transition.rename_axis("fixed_margin_0.1").to_csv(transition_path, sep="\t")
    thresholds = pd.read_csv(inputs["population_thresholds"], sep="\t")
    available = thresholds[thresholds["threshold_available"] == True]
    summary = {
        "schema_version": "bdev.hierarchical_population_umap_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "space_id": config["bundle"]["space_id"],
        "fixture_cells": len(coordinates),
        "fixed_fine_calls": int((wide["fixed_margin_0.1"] == "fine_label").sum()),
        "calibrated_fine_calls": int((wide["population_calibrated"] == "fine_label").sum()),
        "calibrated_only_fine_calls": int(((wide["population_calibrated"] == "fine_label") & (wide["fixed_margin_0.1"] != "fine_label")).sum()),
        "fixed_only_fine_calls": int(((wide["fixed_margin_0.1"] == "fine_label") & (wide["population_calibrated"] != "fine_label")).sum()),
        "both_fine_calls": int(((wide["fixed_margin_0.1"] == "fine_label") & (wide["population_calibrated"] == "fine_label")).sum()),
        "calibrated_thresholds": len(available),
        "thresholds_below_0.1": int((available["calibrated_minimum_margin"] < 0.1 - 1e-12).sum()),
        "thresholds_above_0.1": int((available["calibrated_minimum_margin"] > 0.1 + 1e-12).sum()),
        "labels_accepted": False,
        "bundle_mutated": False,
    }
    summary_path = reports / "01_umap_policy_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    report_path = reports / "02_umap_policy_explanation.md"
    report_path.write_text(
        "# Bdev hierarchical policy UMAP\n\n"
        f"The fixed 0.10 policy makes {summary['fixed_fine_calls']} fine calls. The population-calibrated policy makes {summary['calibrated_fine_calls']}: "
        f"{summary['both_fine_calls']} are shared, {summary['calibrated_only_fine_calls']} are added, and {summary['fixed_only_fine_calls']} fixed-margin calls are removed. "
        f"Of {summary['calibrated_thresholds']} available population thresholds, {summary['thresholds_below_0.1']} are below 0.10 and {summary['thresholds_above_0.1']} are above 0.10.\n\n"
        "All colors and labels are diagnostic; Bdev has no independent fine-label truth in this comparison.\n"
    )

    workflow_dir = Path(__file__).resolve().parent
    generated = [coordinate_path, receipt_path, status_png, status_pdf, label_png, label_pdf, transition_path, summary_path, report_path]
    manifest = {
        "schema_version": "bdev.hierarchical_population_umap_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in inputs.items()},
        "workflow_files": [
            {"path": str(path), "sha256": sha256(path)}
            for path in sorted(workflow_dir.iterdir()) if path.is_file()
        ],
        "outputs": [
            {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in generated
        ],
        "bundle_fingerprint": config["bundle"]["fingerprint"],
        "space_id": config["bundle"]["space_id"],
        "bundle_mutated": False,
        "labels_accepted": False,
    }
    manifest_path = manifests / "02_final_run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"manifest_sha256={sha256(manifest_path)}")


if __name__ == "__main__":
    main()
