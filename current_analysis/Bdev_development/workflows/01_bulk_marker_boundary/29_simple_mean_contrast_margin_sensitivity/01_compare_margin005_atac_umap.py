#!/usr/bin/env python3
"""Compare raw-UCell winner margins 0.01 and 0.05 on the Bdev ATAC UMAP."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/bdev_simple_mean_margin005_mpl")
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


def scatter_labels(
    axis: plt.Axes,
    data: pd.DataFrame,
    label_field: str,
    colors: dict[str, object],
    unassigned: str,
    size: float = 11,
) -> None:
    background = data[data[label_field] == unassigned]
    axis.scatter(
        background["UMAP_1"], background["UMAP_2"], s=size,
        c=colors[unassigned], alpha=0.48, linewidths=0, rasterized=True,
        zorder=1,
    )
    counts = data.loc[data[label_field] != unassigned, label_field].value_counts()
    for label in counts.index:
        part = data[data[label_field] == label]
        axis.scatter(
            part["UMAP_1"], part["UMAP_2"], s=size, c=[colors[label]],
            alpha=0.88, linewidths=0, rasterized=True, zorder=2,
        )


def annotate_numbers(
    axis: plt.Axes,
    data: pd.DataFrame,
    label_field: str,
    label_numbers: dict[str, int],
    colors: dict[str, object],
    unassigned: str,
    minimum_cells: int,
) -> None:
    for label, part in data[data[label_field] != unassigned].groupby(label_field):
        if len(part) < minimum_cells:
            continue
        axis.text(
            float(part["UMAP_1"].median()), float(part["UMAP_2"].median()),
            str(label_numbers[label]), ha="center", va="center", fontsize=6.2,
            fontweight="bold", color="white", zorder=10,
            bbox={"boxstyle": "circle,pad=0.17", "facecolor": colors[label],
                  "edgecolor": "white", "linewidth": 0.6, "alpha": 0.95},
        )


def save_pair(fig: plt.Figure, figures: Path, stem: str) -> list[Path]:
    png = figures / f"{stem}.png"
    pdf = figures / f"{stem}.pdf"
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
    require(
        config["schema_version"] == "bdev.simple_mean_contrast_margin_sensitivity.v1",
        "Unexpected config schema",
    )

    inputs = {name: Path(spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        require(sha256(path) == config["inputs"][name]["sha256"],
                f"Checksum mismatch: {name}")

    root = Path(config["storage"]["results_root"])
    tables = root / "01_tables"
    figures = root / "02_figures"
    reports = root / "03_reports"
    manifests = root / "06_manifests"
    for directory in (tables, figures, reports, manifests):
        directory.mkdir(parents=True, exist_ok=True)

    coordinates = pd.read_csv(inputs["coordinates"], sep="\t")
    decisions = pd.read_csv(
        inputs["decisions_margin001"], sep="\t", keep_default_na=False
    )
    expected_cells = int(config["acceptance_boundary"]["expected_cells"])
    require(len(coordinates) == expected_cells and coordinates["cell_id"].is_unique,
            "Unexpected ATAC coordinate inventory")
    require(len(decisions) == expected_cells and decisions["cell_id"].is_unique,
            "Unexpected decision inventory")
    require(set(coordinates["cell_id"]) == set(decisions["cell_id"]),
            "Coordinate and decision cell inventories differ")
    require(config["inputs"]["coordinates"]["space_id"] == "umap.atac",
            "Coordinates are not declared as umap.atac")

    data = decisions.merge(coordinates, on="cell_id", how="left", validate="one_to_one")
    require(data[["UMAP_1", "UMAP_2"]].notna().all().all(),
            "Missing ATAC UMAP coordinates")
    unassigned = str(config["comparison"]["unassigned_label"])
    baseline_margin = float(config["comparison"]["baseline_minimum_raw_margin"])
    diagnostic_margin = float(config["comparison"]["diagnostic_minimum_raw_margin"])
    require(diagnostic_margin > baseline_margin, "Diagnostic margin is not stricter")
    data["label_margin001"] = data["label"].replace("", unassigned).fillna(unassigned)
    baseline_assigned = data["label_margin001"] != unassigned
    data["label_margin005"] = np.where(
        baseline_assigned & (data["score_margin"] >= diagnostic_margin),
        data["label_margin001"], unassigned,
    )
    diagnostic_assigned = data["label_margin005"] != unassigned
    data["margin_change_status"] = np.select(
        [diagnostic_assigned, baseline_assigned & ~diagnostic_assigned],
        ["retained assignment", "newly unassigned at 0.05"],
        default="already unassigned at 0.01",
    )
    require(int(baseline_assigned.sum()) ==
            int(config["acceptance_boundary"]["baseline_expected_assigned"]),
            "Baseline assigned count changed")
    require(int(diagnostic_assigned.sum()) ==
            int(config["acceptance_boundary"]["diagnostic_expected_assigned"]),
            "Diagnostic assigned count changed")

    old_counts = data.loc[baseline_assigned, "label_margin001"].value_counts()
    new_counts = data.loc[diagnostic_assigned, "label_margin005"].value_counts()
    ordered_labels = list(old_counts.index)
    cmap = plt.get_cmap("turbo")
    positions = np.linspace(0.02, 0.96, len(ordered_labels))
    colors = {label: cmap(position) for label, position in zip(ordered_labels, positions)}
    colors[unassigned] = "#D0D0D0"
    label_numbers = {label: number for number, label in enumerate(ordered_labels, 1)}
    minimum_cells = int(config["comparison"]["minimum_cells_for_direct_label"])

    label_comparison = pd.DataFrame({
        "label": ordered_labels,
        "assigned_margin001": [int(old_counts.get(label, 0)) for label in ordered_labels],
        "assigned_margin005": [int(new_counts.get(label, 0)) for label in ordered_labels],
    })
    label_comparison["newly_unassigned"] = (
        label_comparison["assigned_margin001"] - label_comparison["assigned_margin005"]
    )
    label_comparison["retained_fraction"] = (
        label_comparison["assigned_margin005"] /
        label_comparison["assigned_margin001"]
    )
    label_path = tables / "01_label_count_comparison.tsv"
    label_comparison.to_csv(label_path, sep="\t", index=False)

    sample_field = str(config["comparison"]["sample_field"])
    sample_rows = []
    for sample_id, part in data.groupby(sample_field, sort=True):
        old_n = int((part["label_margin001"] != unassigned).sum())
        new_n = int((part["label_margin005"] != unassigned).sum())
        sample_rows.append({
            sample_field: sample_id,
            "cells": len(part),
            "assigned_margin001": old_n,
            "assigned_margin005": new_n,
            "newly_unassigned": old_n - new_n,
            "retained_fraction_of_margin001": new_n / old_n if old_n else 0,
        })
    sample_comparison = pd.DataFrame(sample_rows)
    sample_path = tables / "02_sample_count_comparison.tsv"
    sample_comparison.to_csv(sample_path, sep="\t", index=False)
    cell_path = tables / "03_cell_margin_comparison.tsv.gz"
    data[[
        "cell_id", sample_field, "label_margin001", "label_margin005",
        "top_score", "runner_up_score", "score_margin", "margin_change_status",
        "UMAP_1", "UMAP_2", "pt",
    ]].to_csv(cell_path, sep="\t", index=False, compression="gzip")

    generated: list[Path] = [label_path, sample_path, cell_path]

    fig, axis = plt.subplots(figsize=(15.5, 9.2))
    fig.subplots_adjust(left=0.05, right=0.70, bottom=0.07, top=0.91)
    scatter_labels(axis, data, "label_margin005", colors, unassigned)
    annotate_numbers(
        axis, data, "label_margin005", label_numbers, colors, unassigned,
        minimum_cells,
    )
    style_axis(axis)
    assigned005 = int(diagnostic_assigned.sum())
    axis.set_title(
        "Simple ImmGen mean-contrast labels at raw UCell margin ≥ 0.05\n"
        f"{assigned005:,} assigned; {len(data) - assigned005:,} unassigned",
        fontsize=14,
    )
    handles = [
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor=colors[label],
            markeredgecolor="none", markersize=6,
            label=f"{label_numbers[label]:>2}  {label} ({int(new_counts[label])})",
        )
        for label in ordered_labels if label in new_counts.index
    ]
    handles.append(Line2D(
        [0], [0], marker="o", color="none", markerfacecolor=colors[unassigned],
        markeredgecolor="none", markersize=6,
        label=f"    {unassigned} ({len(data) - assigned005:,})",
    ))
    fig.legend(
        handles=handles, loc="upper left", bbox_to_anchor=(0.705, 0.89),
        ncol=1, frameon=False, fontsize=8.0, title="Label (cells)",
        title_fontsize=9, handletextpad=0.35,
    )
    generated.extend(save_pair(fig, figures, "01_margin005_labels_atac_umap"))

    fig, axes = plt.subplots(1, 2, figsize=(15.2, 7.0), sharex=True, sharey=True)
    for axis, field, title in (
        (axes[0], "label_margin001", "Margin ≥ 0.01: 1,416 assigned"),
        (axes[1], "label_margin005", "Margin ≥ 0.05: 443 assigned"),
    ):
        scatter_labels(axis, data, field, colors, unassigned, size=9)
        annotate_numbers(axis, data, field, label_numbers, colors, unassigned,
                         minimum_cells)
        style_axis(axis)
        axis.set_title(title, fontsize=12)
    fig.suptitle(
        "Effect of the raw UCell winner–runner-up margin on Bdev labels",
        fontsize=15, y=0.98,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    generated.extend(save_pair(fig, figures, "02_margin001_vs_margin005_atac_umap"))

    status_colors = {
        "already unassigned at 0.01": "#D0D0D0",
        "newly unassigned at 0.05": "#E68613",
        "retained assignment": "#2166AC",
    }
    status_order = list(status_colors)
    status_counts = data["margin_change_status"].value_counts()
    fig, axis = plt.subplots(figsize=(9.2, 8.0))
    for zorder, status in enumerate(status_order, 1):
        part = data[data["margin_change_status"] == status]
        axis.scatter(
            part["UMAP_1"], part["UMAP_2"], s=12,
            c=status_colors[status], alpha=0.82 if zorder > 1 else 0.5,
            linewidths=0, rasterized=True, zorder=zorder,
        )
    style_axis(axis)
    axis.set_title(
        "Assignments removed by increasing the raw UCell margin to 0.05",
        fontsize=13,
    )
    status_handles = [Line2D(
        [0], [0], marker="o", color="none", markerfacecolor=status_colors[status],
        markeredgecolor="none", markersize=7,
        label=f"{status} ({int(status_counts.get(status, 0)):,})",
    ) for status in reversed(status_order)]
    axis.legend(handles=status_handles, loc="upper left", frameon=False,
                fontsize=9)
    fig.tight_layout()
    generated.extend(save_pair(fig, figures, "03_margin005_assignment_changes_atac_umap"))

    newly_unassigned = int((baseline_assigned & ~diagnostic_assigned).sum())
    summary = {
        "schema_version": "bdev.simple_mean_contrast_margin_sensitivity_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "cells": len(data),
        "baseline_margin": baseline_margin,
        "diagnostic_margin": diagnostic_margin,
        "baseline_assigned": int(baseline_assigned.sum()),
        "diagnostic_assigned": assigned005,
        "newly_unassigned": newly_unassigned,
        "diagnostic_unassigned": int((~diagnostic_assigned).sum()),
        "baseline_assigned_labels": int(data.loc[baseline_assigned, "label_margin001"].nunique()),
        "diagnostic_assigned_labels": int(data.loc[diagnostic_assigned, "label_margin005"].nunique()),
        "retained_fraction_of_baseline_assignments": assigned005 / int(baseline_assigned.sum()),
        "labels_accepted": False,
        "artifacts_promoted": False,
        "bundle_mutated": False,
    }
    summary_path = reports / "01_margin_sensitivity_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    generated.append(summary_path)

    workflow_dir = Path(__file__).resolve().parent
    manifest = {
        "schema_version": "bdev.simple_mean_contrast_margin_sensitivity_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "inputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in inputs.items()
        },
        "workflow_files": [
            {"path": str(path), "sha256": sha256(path)}
            for path in sorted(workflow_dir.iterdir()) if path.is_file()
        ],
        "outputs": [
            {"path": str(path.relative_to(root)), "bytes": path.stat().st_size,
             "sha256": sha256(path)} for path in generated
        ],
        "bundle_mutated": False,
        "labels_accepted": False,
        "artifacts_promoted": False,
    }
    manifest_path = manifests / "01_final_run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"manifest_sha256={sha256(manifest_path)}")


if __name__ == "__main__":
    main()

