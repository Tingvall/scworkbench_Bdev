#!/usr/bin/env python3
"""Summarize bootstrap label stability and plot it on the ATAC UMAP."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/bdev_mean_contrast_bootstrap_mpl")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def style_umap(axis: plt.Axes) -> None:
    axis.set_xlabel("ATAC UMAP 1"); axis.set_ylabel("ATAC UMAP 2")
    axis.set_xticks([]); axis.set_yticks([])
    axis.set_aspect("equal", adjustable="box")
    axis.spines[["top", "right", "left", "bottom"]].set_visible(False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    cfg = yaml.safe_load(config_path.read_text())
    for name in ("baseline_decisions", "atac_coordinates"):
        path = Path(cfg["inputs"][name]["path"])
        if sha256(path) != cfg["inputs"][name]["sha256"]:
            raise RuntimeError(f"Input checksum changed: {name}")
    if cfg["inputs"]["atac_coordinates"]["space_id"] != "umap.atac":
        raise RuntimeError("Coordinates are not declared as umap.atac")

    root = Path(cfg["storage"]["results_root"])
    tables, figures, reports, manifests = [root / name for name in ("01_tables", "02_figures", "03_reports", "06_manifests")]
    for directory in (tables, figures, reports, manifests):
        directory.mkdir(parents=True, exist_ok=True)
    baseline = pd.read_csv(cfg["inputs"]["baseline_decisions"]["path"], sep="\t", keep_default_na=False)
    boot = pd.read_csv(tables / "05_bootstrap_bdev_decisions.tsv.gz", sep="\t", keep_default_na=False)
    coordinates = pd.read_csv(cfg["inputs"]["atac_coordinates"]["path"], sep="\t")
    rounds = int(cfg["bootstrap"]["rounds"])
    unknown = str(cfg["scoring"]["unknown_label"])
    if len(boot) != len(baseline) * rounds or boot["bootstrap_round"].nunique() != rounds:
        raise RuntimeError("Bootstrap decision inventory is incomplete")

    base = baseline[["cell_id", "label", "sampleId", "pt"]].rename(columns={"label": "baseline_label"})
    joined = boot.merge(base[["cell_id", "baseline_label"]], on="cell_id", validate="many_to_one")
    joined["same_as_baseline"] = joined["label"] == joined["baseline_label"]
    joined["assigned"] = joined["label"] != unknown
    cell_rows = []
    for cell_id, group in joined.groupby("cell_id", sort=False):
        counts = group["label"].value_counts()
        modes = sorted(counts[counts == counts.max()].index)
        base_row = base.loc[base["cell_id"] == cell_id].iloc[0]
        cell_rows.append({
            "cell_id": cell_id,
            "sampleId": str(base_row["sampleId"]),
            "pt": base_row["pt"],
            "baseline_label": base_row["baseline_label"],
            "bootstrap_mode_label": modes[0],
            "same_as_baseline_fraction": group["same_as_baseline"].mean(),
            "assigned_fraction": group["assigned"].mean(),
            "mode_fraction": counts.max() / rounds,
            "unique_bootstrap_labels": counts.size,
        })
    cells = pd.DataFrame(cell_rows)
    threshold = float(cfg["stability"]["stable_cell_agreement_threshold"])
    cells["stable_at_threshold"] = cells["same_as_baseline_fraction"] >= threshold
    cells.to_csv(tables / "06_per_cell_label_stability.tsv.gz", sep="\t", index=False, compression="gzip")

    label_rows = []
    for label, group in cells.groupby("baseline_label"):
        values = group["same_as_baseline_fraction"]
        label_rows.append({
            "baseline_label": label, "cells": len(group),
            "mean_agreement": values.mean(), "median_agreement": values.median(),
            "p10_agreement": values.quantile(0.1), "p25_agreement": values.quantile(0.25),
            "p75_agreement": values.quantile(0.75),
            "stable_cells": int(group["stable_at_threshold"].sum()),
            "stable_fraction": group["stable_at_threshold"].mean(),
            "mean_assigned_fraction": group["assigned_fraction"].mean(),
        })
    labels = pd.DataFrame(label_rows).sort_values(["cells", "baseline_label"], ascending=[False, True])
    labels.to_csv(tables / "07_baseline_label_stability.tsv", sep="\t", index=False)
    transitions = (joined.groupby(["baseline_label", "label"]).size().rename("cell_rounds").reset_index())
    transitions["fraction_within_baseline_label"] = transitions["cell_rounds"] / transitions.groupby("baseline_label")["cell_rounds"].transform("sum")
    transitions.to_csv(tables / "08_label_transition_frequencies.tsv", sep="\t", index=False)

    population_rounds = pd.read_csv(tables / "03_population_geneset_stability_by_round.tsv", sep="\t")
    population_summary = (population_rounds.groupby(["population", "population_libraries", "stability_evidence"])
                          .agg(mean_overlap_top20=("overlap_baseline_top20", "mean"),
                               minimum_overlap_top20=("overlap_baseline_top20", "min"),
                               mean_jaccard_top20=("jaccard_baseline_top20", "mean"),
                               minimum_jaccard_top20=("jaccard_baseline_top20", "min"))
                          .reset_index().sort_values(["stability_evidence", "mean_jaccard_top20", "population"], ascending=[True, True, True]))
    population_summary.to_csv(tables / "09_population_geneset_stability_summary.tsv", sep="\t", index=False)

    plot = cells.merge(coordinates, on="cell_id", validate="one_to_one")
    fig, axis = plt.subplots(figsize=(8.3, 8.0))
    points = axis.scatter(plot["UMAP_1"], plot["UMAP_2"], c=plot["same_as_baseline_fraction"],
                          cmap="viridis", vmin=0, vmax=1, s=12, alpha=0.9, linewidths=0, rasterized=True)
    style_umap(axis)
    axis.set_title(f"Bootstrap agreement with baseline labels on ATAC UMAP\n{rounds} ImmGen bootstrap rounds", fontsize=13)
    colorbar = fig.colorbar(points, ax=axis, fraction=0.035, pad=0.02)
    colorbar.set_label("Fraction retaining baseline label")
    fig.tight_layout()
    stability_png = figures / "01_per_cell_bootstrap_agreement_atac_umap.png"
    stability_pdf = figures / "01_per_cell_bootstrap_agreement_atac_umap.pdf"
    fig.savefig(stability_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stability_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.5), sharex=True, sharey=True)
    for axis, sample_id in zip(axes.flat, sorted(plot["sampleId"].astype(str).unique())):
        part = plot[plot["sampleId"].astype(str) == sample_id]
        points = axis.scatter(part["UMAP_1"], part["UMAP_2"], c=part["same_as_baseline_fraction"],
                              cmap="viridis", vmin=0, vmax=1, s=11, alpha=0.9, linewidths=0, rasterized=True)
        style_umap(axis)
        axis.set_title(f"sampleId {sample_id}; median={part['same_as_baseline_fraction'].median():.2f}")
    fig.suptitle("Bootstrap label agreement by sample on the ATAC UMAP", fontsize=14)
    colorbar = fig.colorbar(points, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
    colorbar.set_label("Fraction retaining baseline label")
    facet_png = figures / "02_per_cell_bootstrap_agreement_by_sample_atac_umap.png"
    facet_pdf = figures / "02_per_cell_bootstrap_agreement_by_sample_atac_umap.pdf"
    fig.savefig(facet_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(facet_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    min_cells = int(cfg["stability"]["label_plot_minimum_cells_exclusive"])
    shown = labels[(labels["cells"] > min_cells) & (labels["baseline_label"] != unknown)].copy().sort_values("median_agreement")
    fig, axis = plt.subplots(figsize=(9.0, max(4.5, len(shown) * 0.42)))
    y = np.arange(len(shown))
    axis.barh(y, shown["median_agreement"], color="#3182BD", alpha=0.9)
    axis.hlines(y, shown["p25_agreement"], shown["p75_agreement"], color="#08306B", linewidth=3)
    axis.axvline(threshold, color="#CB181D", linestyle="--", linewidth=1.2, label=f"stability threshold {threshold:.2f}")
    axis.set_yticks(y, [f"{label} (n={cells:,})" for label, cells in zip(shown["baseline_label"], shown["cells"])])
    axis.set_xlim(0, 1.01); axis.set_xlabel("Per-cell agreement with baseline (median; IQR line)")
    axis.set_title(f"Label stability across {rounds} ImmGen bootstrap rounds")
    axis.spines[["top", "right"]].set_visible(False); axis.legend(frameon=False)
    fig.tight_layout()
    label_png = figures / "03_label_bootstrap_stability.png"
    label_pdf = figures / "03_label_bootstrap_stability.pdf"
    fig.savefig(label_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(label_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    assigned_baseline = cells[cells["baseline_label"] != unknown]
    summary = {
        "schema_version": "bdev.simple_mean_contrast_bootstrap_stability_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": cfg["run_identity"],
        "bootstrap_rounds": rounds,
        "cells": len(cells),
        "overall_cell_round_agreement": joined["same_as_baseline"].mean(),
        "median_per_cell_agreement": cells["same_as_baseline_fraction"].median(),
        "stable_cell_threshold": threshold,
        "stable_cells": int(cells["stable_at_threshold"].sum()),
        "stable_cell_fraction": cells["stable_at_threshold"].mean(),
        "baseline_assigned_cells": len(assigned_baseline),
        "baseline_assigned_stable_fraction": assigned_baseline["stable_at_threshold"].mean(),
        "bootstrap_assessable_populations": int((population_summary["stability_evidence"] == "bootstrap_assessable").sum()),
        "singleton_populations_structurally_unassessable": int((population_summary["stability_evidence"] != "bootstrap_assessable").sum()),
        "bdev_used_for_gene_selection": False,
        "default_changed_by_run": False,
        "artifacts_promoted": False,
        "accepted_labels_emitted": False,
    }
    summary_path = reports / "02_bootstrap_stability_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    top_labels = labels[(labels["cells"] > min_cells)].sort_values("cells", ascending=False)
    report_lines = [
        "# Simple mean-contrast bootstrap stability", "",
        f"Across {rounds} seeded within-population ImmGen bootstrap rounds, the overall cell-round agreement with the baseline was {summary['overall_cell_round_agreement']:.3f}.",
        f"{summary['stable_cells']}/{summary['cells']} cells ({summary['stable_cell_fraction']:.1%}) retained the baseline label in at least {threshold:.0%} of rounds.",
        f"Among baseline-assigned cells, the stable fraction was {summary['baseline_assigned_stable_fraction']:.1%}.", "",
        "Eight singleton ImmGen populations are structurally unassessable by within-population resampling and must not be interpreted as validated by invariant bootstrap sets.", "",
        "## Labels with more than 20 baseline cells", "",
        "| Label | Cells | Median agreement | Stable fraction |", "|---|---:|---:|---:|",
    ]
    for row in top_labels.itertuples():
        report_lines.append(f"| {row.baseline_label} | {row.cells} | {row.median_agreement:.3f} | {row.stable_fraction:.1%} |")
    report_lines += ["", "This is diagnostic evidence only. Bdev was used to measure transfer stability, not to select genes; the default and source bundle were not changed."]
    report_path = reports / "03_bootstrap_stability_report.md"
    report_path.write_text("\n".join(report_lines) + "\n")

    generated = sorted([path for directory in (tables, figures, reports, manifests) for path in directory.rglob("*") if path.is_file() and path.name != "02_final_run_manifest.json"])
    workflow_dir = Path(__file__).resolve().parent
    manifest = {
        "schema_version": "bdev.simple_mean_contrast_bootstrap_stability_manifest.v1",
        "status": "complete_diagnostic_only", "run_identity": cfg["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "inputs": {name: {"path": spec["path"], "sha256": spec["sha256"]}
                   for name, spec in cfg["inputs"].items() if "sha256" in spec},
        "workflow_files": [{"path": str(path), "sha256": sha256(path)} for path in sorted(workflow_dir.iterdir()) if path.is_file()],
        "outputs": [{"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in generated],
        "bundle_fingerprint_unchanged": True, "default_changed_by_run": False,
        "artifacts_promoted": False, "accepted_labels_emitted": False,
    }
    manifest_path = manifests / "02_final_run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"manifest_sha256={sha256(manifest_path)}")


if __name__ == "__main__":
    main()
