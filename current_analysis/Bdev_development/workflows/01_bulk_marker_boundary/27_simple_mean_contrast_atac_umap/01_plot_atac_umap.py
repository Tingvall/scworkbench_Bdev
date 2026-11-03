#!/usr/bin/env python3
"""Plot simple mean-contrast diagnostic labels on the existing Bdev ATAC UMAP."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", "/tmp/bdev_simple_mean_contrast_atac_umap_mpl")
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


def scatter_labels(axis: plt.Axes, data: pd.DataFrame, colors: dict[str, object],
                   unassigned: str, size: float, alpha: float) -> None:
    background = data[data["label"] == unassigned]
    axis.scatter(background["UMAP_1"], background["UMAP_2"], s=size, c=colors[unassigned],
                 alpha=0.55, linewidths=0, rasterized=True, zorder=1)
    counts = data.loc[data["label"] != unassigned, "label"].value_counts()
    for label in counts.index:
        part = data[data["label"] == label]
        axis.scatter(part["UMAP_1"], part["UMAP_2"], s=size, c=[colors[label]],
                     alpha=alpha, linewidths=0, rasterized=True, zorder=2)


def annotate_numbers(axis: plt.Axes, data: pd.DataFrame, label_numbers: dict[str, int],
                     colors: dict[str, object], minimum_cells: int, unassigned: str) -> None:
    for label, part in data[data["label"] != unassigned].groupby("label"):
        if len(part) < minimum_cells:
            continue
        axis.text(float(part["UMAP_1"].median()), float(part["UMAP_2"].median()),
                  str(label_numbers[label]), ha="center", va="center", fontsize=6.2,
                  fontweight="bold", color="white", zorder=10,
                  bbox={"boxstyle": "circle,pad=0.17", "facecolor": colors[label],
                        "edgecolor": "white", "linewidth": 0.6, "alpha": 0.95})


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").lower()


def highlight_label(axis: plt.Axes, data: pd.DataFrame, label: str, color: object,
                    size: float, alpha: float) -> None:
    background = data[data["label"] != label]
    target = data[data["label"] == label]
    axis.scatter(background["UMAP_1"], background["UMAP_2"], s=size, c="#D7D7D7",
                 alpha=0.48, linewidths=0, rasterized=True, zorder=1)
    axis.scatter(target["UMAP_1"], target["UMAP_2"], s=size + 3, c=[color],
                 alpha=alpha, linewidths=0, rasterized=True, zorder=2)
    axis.set_title(f"{label} (n={len(target):,})", fontsize=10)
    style_axis(axis)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    require(config["schema_version"] == "bdev.simple_mean_contrast_atac_umap.v1",
            "Unexpected config schema")

    inputs = {name: Path(spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        require(sha256(path) == config["inputs"][name]["sha256"], f"Checksum mismatch: {name}")

    root = Path(config["storage"]["results_root"])
    tables = root / "01_tables"
    figures = root / "02_figures"
    reports = root / "03_reports"
    manifests = root / "06_manifests"
    for directory in (tables, figures, reports, manifests):
        directory.mkdir(parents=True, exist_ok=True)

    coordinates = pd.read_csv(inputs["coordinates"], sep="\t")
    decisions = pd.read_csv(inputs["decisions"], sep="\t", keep_default_na=False)
    require(len(coordinates) == 1746 and coordinates["cell_id"].is_unique,
            "Unexpected ATAC coordinate inventory")
    require(len(decisions) == 1746 and decisions["cell_id"].is_unique,
            "Unexpected decision inventory")
    require(set(coordinates["cell_id"]) == set(decisions["cell_id"]),
            "Coordinate and decision cell inventories differ")
    require(config["inputs"]["coordinates"]["space_id"] == "umap.atac",
            "Coordinates are not declared as umap.atac")

    data = decisions.merge(coordinates, on="cell_id", how="left", validate="one_to_one")
    require(data[["UMAP_1", "UMAP_2"]].notna().all().all(), "Missing ATAC UMAP coordinates")
    unassigned = str(config["plot"]["unassigned_label"])
    data["label"] = data["label"].replace("", unassigned).fillna(unassigned)
    sample_field = str(config["plot"]["sample_field"])
    data[sample_field] = data[sample_field].astype(str)

    assigned_counts = data.loc[data["label"] != unassigned, "label"].value_counts()
    ordered_labels = list(assigned_counts.index)
    cmap = plt.get_cmap("turbo")
    positions = np.linspace(0.02, 0.96, len(ordered_labels))
    colors = {label: cmap(position) for label, position in zip(ordered_labels, positions)}
    colors[unassigned] = "#D0D0D0"
    label_numbers = {label: number for number, label in enumerate(ordered_labels, 1)}

    label_counts = (data.groupby([sample_field, "label"], observed=True).size()
                    .rename("cells").reset_index())
    label_counts["fraction_sample"] = label_counts["cells"] / label_counts.groupby(sample_field)["cells"].transform("sum")
    counts_path = tables / "01_label_counts_by_sample.tsv"
    label_counts.to_csv(counts_path, sep="\t", index=False)
    joined_path = tables / "02_atac_umap_labels.tsv.gz"
    data[["cell_id", sample_field, "label", "UMAP_1", "UMAP_2", "top_score",
          "runner_up_score", "score_margin", "abstain_reason", "pt"]].to_csv(
              joined_path, sep="\t", index=False, compression="gzip")

    minimum_cells = int(config["plot"]["minimum_cells_for_direct_label"])
    alpha = float(config["plot"]["point_alpha"])
    fig, axis = plt.subplots(figsize=(15.5, 9.2))
    fig.subplots_adjust(left=0.05, right=0.70, bottom=0.07, top=0.91)
    scatter_labels(axis, data, colors, unassigned,
                   float(config["plot"]["point_size_overall"]), alpha)
    annotate_numbers(axis, data, label_numbers, colors, minimum_cells, unassigned)
    style_axis(axis)
    axis.set_title(
        "Simple ImmGen mean-contrast labels on the Bdev ATAC UMAP\n"
        f"{len(data) - int((data['label'] == unassigned).sum()):,} assigned; "
        f"{int((data['label'] == unassigned).sum()):,} unassigned",
        fontsize=14,
    )
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[label],
                      markeredgecolor="none", markersize=6,
                      label=f"{label_numbers[label]:>2}  {label} ({assigned_counts[label]})")
               for label in ordered_labels]
    handles.append(Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[unassigned],
                          markeredgecolor="none", markersize=6,
                          label=f"    {unassigned} ({int((data['label'] == unassigned).sum())})"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.705, 0.89), ncol=2,
               frameon=False, fontsize=7.4, title="Label (total cells)", title_fontsize=9,
               handletextpad=0.35, columnspacing=0.8)
    overall_png = figures / "01_simple_mean_contrast_labels_atac_umap.png"
    overall_pdf = figures / "01_simple_mean_contrast_labels_atac_umap.pdf"
    fig.savefig(overall_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(overall_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    sample_ids = sorted(data[sample_field].unique())
    require(len(sample_ids) == 4, "Expected four sample IDs")
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 10.4), sharex=True, sharey=True)
    for axis, sample_id in zip(axes.flat, sample_ids):
        part = data[data[sample_field] == sample_id]
        scatter_labels(axis, part, colors, unassigned,
                       float(config["plot"]["point_size_facet"]), alpha)
        annotate_numbers(axis, part, label_numbers, colors, minimum_cells, unassigned)
        assigned = int((part["label"] != unassigned).sum())
        axis.set_title(f"sampleId {sample_id}: {assigned}/{len(part)} assigned", fontsize=11)
        style_axis(axis)
    fig.suptitle("Simple ImmGen mean-contrast labels by sample on the ATAC UMAP", fontsize=14, y=0.98)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    facet_png = figures / "02_simple_mean_contrast_labels_by_sample_atac_umap.png"
    facet_pdf = figures / "02_simple_mean_contrast_labels_by_sample_atac_umap.pdf"
    fig.savefig(facet_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(facet_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    highlight_threshold = int(config["plot"]["highlight_label_minimum_cells_exclusive"])
    highlight_labels = [label for label in ordered_labels if int(assigned_counts[label]) > highlight_threshold]
    highlight_dir = figures / "labels_over20"
    highlight_dir.mkdir(parents=True, exist_ok=True)
    highlight_outputs = []
    ncols = 3
    nrows = int(np.ceil(len(highlight_labels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12.6, 4.0 * nrows), sharex=True, sharey=True)
    for axis, label in zip(np.asarray(axes).flat, highlight_labels):
        highlight_label(axis, data, label, colors[label],
                        float(config["plot"]["point_size_facet"]), alpha)
    for axis in np.asarray(axes).flat[len(highlight_labels):]:
        axis.set_axis_off()
    fig.suptitle(
        f"Labels with more than {highlight_threshold} cells on the Bdev ATAC UMAP\n"
        "Target label colored; all other cells grey",
        fontsize=14, y=0.995,
    )
    fig.tight_layout(rect=(0, 0.01, 1, 0.965))
    highlights_png = figures / "03_each_label_over20_atac_umap.png"
    highlights_pdf = figures / "03_each_label_over20_atac_umap.pdf"
    fig.savefig(highlights_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(highlights_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    highlight_outputs.extend([highlights_png, highlights_pdf])

    for number, label in enumerate(highlight_labels, 1):
        fig, axis = plt.subplots(figsize=(7.2, 7.2))
        highlight_label(axis, data, label, colors[label],
                        float(config["plot"]["point_size_overall"]), alpha)
        fig.tight_layout()
        stem = f"{number:02d}_{safe_name(label)}_atac_umap"
        png_path = highlight_dir / f"{stem}.png"
        pdf_path = highlight_dir / f"{stem}.pdf"
        fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
        fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        highlight_outputs.extend([png_path, pdf_path])

    summary = {
        "schema_version": "bdev.simple_mean_contrast_atac_umap_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "space_id": "umap.atac",
        "cells": int(len(data)),
        "assigned": int((data["label"] != unassigned).sum()),
        "unassigned": int((data["label"] == unassigned).sum()),
        "assigned_labels": int(len(ordered_labels)),
        "highlight_threshold_exclusive": highlight_threshold,
        "highlighted_labels": [
            {"label": label, "cells": int(assigned_counts[label])} for label in highlight_labels
        ],
        "sample_ids": sample_ids,
        "labels_accepted": False,
        "artifacts_promoted": False,
        "bundle_mutated": False,
    }
    summary_path = reports / "01_atac_umap_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    generated = [counts_path, joined_path, overall_png, overall_pdf, facet_png, facet_pdf,
                 highlights_png, highlights_pdf, *highlight_outputs[2:], summary_path]
    workflow_dir = Path(__file__).resolve().parent
    manifest = {
        "schema_version": "bdev.simple_mean_contrast_atac_umap_manifest.v1",
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
        "space_id": "umap.atac",
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
