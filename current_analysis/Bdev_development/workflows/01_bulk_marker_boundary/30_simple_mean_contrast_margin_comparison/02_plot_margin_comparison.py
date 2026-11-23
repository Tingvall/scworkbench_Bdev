#!/usr/bin/env python3
"""Plot pure argmax, margin 0.01, and margin 0.05 on the Bdev ATAC UMAP."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/bdev_simple_mean_margin_comparison_mpl")
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


def style_axis(axis: plt.Axes) -> None:
    axis.set_xlabel("ATAC UMAP 1")
    axis.set_ylabel("ATAC UMAP 2")
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_aspect("equal", adjustable="box")
    axis.spines[["top", "right", "left", "bottom"]].set_visible(False)


def scatter_labels(axis, data, field, colors, unassigned, size=8):
    background = data[data[field] == unassigned]
    axis.scatter(background["UMAP_1"], background["UMAP_2"], s=size,
                 c=colors[unassigned], alpha=0.48, linewidths=0,
                 rasterized=True, zorder=1)
    counts = data.loc[data[field] != unassigned, field].value_counts()
    for label in counts.index:
        part = data[data[field] == label]
        axis.scatter(part["UMAP_1"], part["UMAP_2"], s=size,
                     c=[colors[label]], alpha=0.88, linewidths=0,
                     rasterized=True, zorder=2)


def annotate_numbers(axis, data, field, numbers, colors, unassigned, minimum):
    for label, part in data[data[field] != unassigned].groupby(field):
        if len(part) < minimum:
            continue
        axis.text(float(part["UMAP_1"].median()), float(part["UMAP_2"].median()),
                  str(numbers[label]), ha="center", va="center", fontsize=5.7,
                  fontweight="bold", color="white", zorder=10,
                  bbox={"boxstyle": "circle,pad=0.15", "facecolor": colors[label],
                        "edgecolor": "white", "linewidth": 0.5, "alpha": 0.95})


def save_pair(fig, directory: Path, stem: str) -> list[Path]:
    png = directory / f"{stem}.png"
    pdf = directory / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return [png, pdf]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    root = Path(config["storage"]["results_root"])
    tables = root / "01_tables"
    figures = root / "02_figures"
    reports = root / "03_reports"
    manifests = root / "06_manifests"
    table_path = tables / "01_cell_margin_comparison.tsv.gz"
    require(table_path.is_file(), "Run 01_export_margin_decisions.R first")
    data = pd.read_csv(table_path, sep="\t", keep_default_na=False)
    summary_table = pd.read_csv(tables / "02_threshold_summary.tsv", sep="\t")
    rules = config["comparison"]["rules"]
    unassigned = config["comparison"]["unassigned_label"]
    minimum = int(config["comparison"]["minimum_cells_for_direct_label"])
    require(len(data) == int(config["acceptance_boundary"]["expected_cells"]),
            "Unexpected cell count")

    pure_counts = data.loc[data["label_pure_argmax"] != unassigned,
                           "label_pure_argmax"].value_counts()
    ordered_labels = list(pure_counts.index)
    cmap = plt.get_cmap("turbo")
    positions = np.linspace(0.02, 0.96, len(ordered_labels))
    colors = {label: cmap(position) for label, position in zip(ordered_labels, positions)}
    colors[unassigned] = "#D0D0D0"
    numbers = {label: number for number, label in enumerate(ordered_labels, 1)}

    generated: list[Path] = []
    fig, axes = plt.subplots(1, 3, figsize=(22.5, 7.1), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.035, right=0.77, bottom=0.08, top=0.88, wspace=0.06)
    for axis, rule in zip(axes, rules):
        field = f"label_{rule['rule_id']}"
        row = summary_table[summary_table["rule_id"] == rule["rule_id"]].iloc[0]
        scatter_labels(axis, data, field, colors, unassigned)
        annotate_numbers(axis, data, field, numbers, colors, unassigned, minimum)
        style_axis(axis)
        axis.set_title(
            f"{rule['display_label']}\n{int(row.assigned):,} assigned; "
            f"{int(row.assigned_labels)} labels", fontsize=11,
        )
    fig.suptitle("Raw UCell winner-margin comparison on the Bdev ATAC UMAP",
                 fontsize=15, y=0.97)
    handles = [Line2D(
        [0], [0], marker="o", color="none", markerfacecolor=colors[label],
        markeredgecolor="none", markersize=5.5,
        label=f"{numbers[label]:>2}  {label} ({int(pure_counts[label])})",
    ) for label in ordered_labels]
    handles.append(Line2D([0], [0], marker="o", color="none",
                          markerfacecolor=colors[unassigned], markeredgecolor="none",
                          markersize=5.5, label="unassigned"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.785, 0.89),
               ncol=2, frameon=False, fontsize=7.2,
               title="Pure-winner label (cells)", title_fontsize=9,
               handletextpad=0.3, columnspacing=0.7)
    generated.extend(save_pair(fig, figures, "01_pure_winner_margin001_margin005_atac_umap"))

    data["confidence_stratum"] = np.select(
        [data["top_tie"].astype(str).str.lower().eq("true"),
         data["score_margin"] < 0.01,
         data["score_margin"] < 0.05],
        ["top-score tie", "pure winner only (<0.01)", "passes 0.01 only (<0.05)"],
        default="passes 0.05",
    )
    stratum_order = ["top-score tie", "pure winner only (<0.01)",
                     "passes 0.01 only (<0.05)", "passes 0.05"]
    stratum_colors = {
        "top-score tie": "#C7C7C7",
        "pure winner only (<0.01)": "#FDD49E",
        "passes 0.01 only (<0.05)": "#F16913",
        "passes 0.05": "#2166AC",
    }
    stratum_counts = data["confidence_stratum"].value_counts()
    fig, axis = plt.subplots(figsize=(9.3, 8.0))
    for zorder, stratum in enumerate(stratum_order, 1):
        part = data[data["confidence_stratum"] == stratum]
        axis.scatter(part["UMAP_1"], part["UMAP_2"], s=12,
                     c=stratum_colors[stratum], alpha=0.82,
                     linewidths=0, rasterized=True, zorder=zorder)
    style_axis(axis)
    axis.set_title("Nested raw-UCell winner-margin confidence strata", fontsize=13)
    axis.legend(handles=[Line2D(
        [0], [0], marker="o", color="none", markerfacecolor=stratum_colors[s],
        markeredgecolor="none", markersize=7,
        label=f"{s} ({int(stratum_counts.get(s, 0)):,})",
    ) for s in reversed(stratum_order)], loc="upper left", frameon=False,
        fontsize=9)
    fig.tight_layout()
    generated.extend(save_pair(fig, figures, "02_margin_confidence_strata_atac_umap"))

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.6))
    x = np.arange(len(summary_table))
    labels = list(summary_table["display_label"])
    axes[0].bar(x, summary_table["assigned"], color="#2166AC")
    axes[0].bar(x, summary_table["unassigned"], bottom=summary_table["assigned"],
                color="#D0D0D0")
    axes[0].set_xticks(x, labels, rotation=15, ha="right")
    axes[0].set_ylabel("Cells")
    axes[0].set_title("Assignment coverage")
    axes[0].legend(handles=[Line2D([0], [0], color="#2166AC", lw=7,
                                   label="assigned"),
                            Line2D([0], [0], color="#D0D0D0", lw=7,
                                   label="unassigned")], frameon=False)
    axes[1].bar(x, summary_table["assigned_labels"], color="#4DAF4A")
    axes[1].set_xticks(x, labels, rotation=15, ha="right")
    axes[1].set_ylabel("Distinct winning labels")
    axes[1].set_title("Label representation")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Effect of the winner-margin rule", fontsize=14)
    fig.tight_layout()
    generated.extend(save_pair(fig, figures, "03_margin_assignment_summary"))

    report = {
        "schema_version": "bdev.simple_mean_contrast_margin_comparison_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "cells": int(len(data)),
        "exact_top_score_ties": int(data["top_tie"].astype(str).str.lower().eq("true").sum()),
        "rules": summary_table.to_dict(orient="records"),
        "confidence_strata": {
            key: int(stratum_counts.get(key, 0)) for key in stratum_order
        },
        "ties_unassigned": True,
        "labels_accepted": False,
        "artifacts_promoted": False,
        "bundle_mutated": False,
    }
    report_path = reports / "01_margin_comparison_summary.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    generated.append(report_path)

    workflow_dir = Path(__file__).resolve().parent
    manifest_path = manifests / "02_final_run_manifest.json"
    output_files = sorted(
        path for path in root.rglob("*") if path.is_file() and path != manifest_path
    )
    manifest = {
        "schema_version": "bdev.simple_mean_contrast_margin_comparison_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "inputs": {
            name: {"path": spec["path"], "sha256": sha256(Path(spec["path"]))}
            for name, spec in config["inputs"].items()
        },
        "workflow_files": [
            {"path": str(path), "sha256": sha256(path)}
            for path in sorted(workflow_dir.iterdir()) if path.is_file()
        ],
        "outputs": [
            {"path": str(path.relative_to(root)), "bytes": path.stat().st_size,
             "sha256": sha256(path)} for path in output_files
        ],
        "file_count": len(output_files),
        "bundle_mutated": False,
        "labels_accepted": False,
        "artifacts_promoted": False,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    print(f"manifest_sha256={sha256(manifest_path)}")


if __name__ == "__main__":
    main()

