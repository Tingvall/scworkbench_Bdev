#!/usr/bin/env python3
from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.stats import spearmanr


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_config(path: Path) -> dict:
    current = yaml.safe_load(path.read_text())
    base = current.pop("base_config", None)
    if base is None:
        return current
    merged = read_config((path.parent / base).resolve())
    def update(a, b):
        for k, v in b.items():
            if isinstance(v, dict) and isinstance(a.get(k), dict):
                update(a[k], v)
            else:
                a[k] = v
    update(merged, current)
    return merged


def entropy(weights: np.ndarray) -> float:
    x = weights[np.isfinite(weights) & (weights > 0)]
    if not len(x):
        return float("nan")
    x = x / x.sum()
    return float(-(x * np.log(x)).sum())


def centered(values: pd.Series, groups: pd.Series) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce")
    return x - x.groupby(groups).transform("median")


def safe_spearman(x, y):
    x = pd.to_numeric(pd.Series(x), errors="coerce")
    y = pd.to_numeric(pd.Series(y), errors="coerce")
    keep = x.notna() & y.notna()
    if keep.sum() < 6 or x[keep].nunique() < 2 or y[keep].nunique() < 2:
        return float("nan"), float("nan"), int(keep.sum())
    r, p = spearmanr(x[keep], y[keep])
    return float(r), float(p), int(keep.sum())


config_path = Path(__import__("sys").argv[2]).resolve()
assert __import__("sys").argv[1] == "--config"
cfg = read_config(config_path)
root = Path(cfg["storage"]["results_root"])
tables, figures, reports, logs, manifests = [root / x for x in
    ["01_tables", "02_figures", "03_reports", "04_logs", "06_manifests"]]
for p in [tables, figures, reports, logs, manifests]:
    p.mkdir(parents=True, exist_ok=True)

coords = pd.read_csv(tables / "01_position_driver_latent_coordinates.tsv", sep="\t")
latent_cols = [x for x in coords.columns if x != "cell_id"]
coords = coords.set_index("cell_id")
labels = pd.read_csv(tables / "02_provisional_label_diagnostic_cells.tsv.gz", sep="\t")
ref = pd.read_csv(tables / "03_reference_cell_technical_metadata.tsv.gz", sep="\t")
ledger = pd.read_csv(cfg["inputs"]["neighbor_driver_ledger"]["path"], sep="\t")
mapping_qc = pd.read_csv(cfg["inputs"]["immgen_qc"]["path"], sep="\t")
matrix_qc = pd.read_csv(tables / "06_immgen_matrix_qc.tsv", sep="\t")
unknown = cfg["diagnostic"]["unknown_label"]

# Stability-weighted, sample-balanced provisional RNA-label centroids.
label_centroids = []
for label, x in labels.loc[labels["label"] != unknown].groupby("label", sort=True):
    x = x.loc[x["cell_id"].isin(coords.index)].copy()
    if len(x) < int(cfg["diagnostic"]["minimum_label_cells_for_centroid"]):
        continue
    w = pd.to_numeric(x["same_as_baseline_fraction"], errors="coerce").fillna(0.25).clip(lower=0.05).to_numpy()
    samples = x["orig.ident"].astype(str)
    sample_n = samples.value_counts()
    w = w / samples.map(sample_n).to_numpy()
    w = w / w.sum()
    c = np.average(coords.loc[x["cell_id"], latent_cols].to_numpy(), axis=0, weights=w)
    label_centroids.append({"label": label, "n_cells": len(x), "n_samples": samples.nunique(),
                            **{d: c[i] for i, d in enumerate(latent_cols)}})
centroids = pd.DataFrame(label_centroids).set_index("label")
centroids.reset_index().to_csv(tables / "07_provisional_rna_label_latent_centroids.tsv", sep="\t", index=False)

# Give every mapping-reference cell a diagnostic nearest provisional label.
ref_xyz = coords.loc[ref["cell_id"], latent_cols].to_numpy()
cent_xyz = centroids[latent_cols].to_numpy()
dist = np.sqrt(((ref_xyz[:, None, :] - cent_xyz[None, :, :]) ** 2).sum(axis=2))
order = np.argsort(dist, axis=1)
ref["nearest_provisional_label"] = centroids.index.to_numpy()[order[:, 0]]
ref["nearest_label_distance"] = dist[np.arange(len(ref)), order[:, 0]]
ref["nearest_label_margin"] = dist[np.arange(len(ref)), order[:, 1]] - dist[np.arange(len(ref)), order[:, 0]]

reuse = ledger.groupby("cell_id").agg(reuse_libraries=("query_cell_id", "nunique"),
    reuse_assignments=("query_cell_id", "size"), total_neighbor_weight=("normalized_weight", "sum")).reset_index()
ref = ref.merge(reuse, on="cell_id", how="left")
ref[["reuse_libraries", "reuse_assignments", "total_neighbor_weight"]] = ref[["reuse_libraries", "reuse_assignments", "total_neighbor_weight"]].fillna(0)
ref.to_csv(tables / "08_reference_hub_label_and_technical_audit.tsv.gz", sep="\t", index=False, compression="gzip")

# Library position and exact latent-dimension pull toward observed neighbors.
library_rows, dim_rows = [], []
for gsm, x in ledger.groupby("query_cell_id", sort=False):
    pop = str(x["population"].iloc[0])
    if gsm not in coords.index:
        continue
    q = coords.loc[gsm, latent_cols].to_numpy(float)
    w = x["normalized_weight"].to_numpy(float); w = w / w.sum()
    o = np.average(coords.loc[x["cell_id"], latent_cols].to_numpy(float), axis=0, weights=w)
    provisional = ref.set_index("cell_id").loc[x["cell_id"], "nearest_provisional_label"].to_numpy()
    expected_weight = float(w[provisional == pop].sum())
    cluster6_weight = float(w[x["reference_cluster"].astype(str).to_numpy() == str(cfg["diagnostic"]["cluster_of_interest"])].sum())
    sample_weights = x.assign(_w=w).groupby("reference_sample")["_w"].sum().sort_values(ascending=False)
    has_expected = pop in centroids.index
    e = centroids.loc[pop, latent_cols].to_numpy(float) if has_expected else np.repeat(np.nan, len(latent_cols))
    pull = (q - e) ** 2 - (q - o) ** 2 if has_expected else np.repeat(np.nan, len(latent_cols))
    expected_distance = float(np.linalg.norm(q-e)) if has_expected else np.nan
    library_rows.append(dict(gsm=gsm, population=pop, support_call=x["support_call"].iloc[0], has_sufficient_expected_label_centroid=has_expected,
        expected_label_cells=int(centroids.loc[pop, "n_cells"]) if has_expected else 0, distance_to_expected_centroid=expected_distance,
        distance_to_observed_neighbor_centroid=float(np.linalg.norm(q-o)), observed_over_expected_distance_ratio=float(np.linalg.norm(q-o)/expected_distance) if has_expected else np.nan,
        expected_provisional_label_neighbor_weight=expected_weight, cluster6_neighbor_weight=cluster6_weight,
        dominant_reference_sample=str(sample_weights.index[0]), dominant_reference_sample_weight=float(sample_weights.iloc[0]),
        reference_sample_entropy=entropy(sample_weights.to_numpy()), median_neighbor_distance=float(np.median(x["distance"])),
        mapped_native_l4=float(x["query_accepted_native_l4"].iloc[0]), query_umap_1=float(x["query_UMAP_1"].iloc[0]),
        query_umap_2=float(x["query_UMAP_2"].iloc[0])))
    for j, d in enumerate(latent_cols):
        if not has_expected:
            continue
        dim_rows.append(dict(gsm=gsm, population=pop, latent_dimension=d, query_value=q[j], expected_centroid=e[j],
            observed_neighbor_centroid=o[j], squared_distance_advantage_observed_vs_expected=pull[j],
            absolute_fraction_of_pull=abs(pull[j]) / np.abs(pull).sum() if np.abs(pull).sum() else 0))
libraries = pd.DataFrame(library_rows)
dims = pd.DataFrame(dim_rows)
libraries.to_csv(tables / "09_library_label_position_diagnostic.tsv", sep="\t", index=False)
dims.to_csv(tables / "10_library_latent_dimension_drivers.tsv", sep="\t", index=False)

# Reference technical bias and cluster-6 contrasts.
technical = [x for x in cfg["diagnostic"]["technical_fields"] if x in ref.columns and x not in {"orig.ident", "Phase"}]
tech_rows = []
for metric in technical:
    r, p, n = safe_spearman(ref[metric], ref["reuse_libraries"])
    c6 = pd.to_numeric(ref.loc[ref[cfg["diagnostic"]["cluster_field"]].astype(str) == str(cfg["diagnostic"]["cluster_of_interest"]), metric], errors="coerce")
    other = pd.to_numeric(ref.loc[ref[cfg["diagnostic"]["cluster_field"]].astype(str) != str(cfg["diagnostic"]["cluster_of_interest"]), metric], errors="coerce")
    scale = np.nanmedian(np.abs(other - np.nanmedian(other)))
    tech_rows.append(dict(level="reference", predictor=metric, outcome="reuse_libraries", n=n, spearman=r, p_value=p,
        cluster6_median=float(np.nanmedian(c6)), other_median=float(np.nanmedian(other)),
        cluster6_median_difference_in_other_mad=float((np.nanmedian(c6)-np.nanmedian(other))/scale) if scale else np.nan))
    r_centered, p_centered, n_centered = safe_spearman(
        centered(ref[metric], ref[cfg["diagnostic"]["cluster_field"]].astype(str)),
        centered(ref["reuse_libraries"], ref[cfg["diagnostic"]["cluster_field"]].astype(str)))
    tech_rows.append(dict(level="reference_cluster_centered", predictor=metric, outcome="reuse_libraries", n=n_centered,
        spearman=r_centered, p_value=p_centered, cluster6_median=np.nan, other_median=np.nan,
        cluster6_median_difference_in_other_mad=np.nan))

# Sample over-representation among neighbor weights.
ref_sample_prop = ref["orig.ident"].astype(str).value_counts(normalize=True)
neighbor_sample = ledger.groupby("reference_sample")["normalized_weight"].sum() / ledger["query_cell_id"].nunique()
neighbor_sample.index = neighbor_sample.index.astype(str)
sample_audit = pd.DataFrame({"reference_fraction": ref_sample_prop, "mean_neighbor_weight": neighbor_sample}).fillna(0)
sample_audit["enrichment_ratio"] = sample_audit["mean_neighbor_weight"] / sample_audit["reference_fraction"]
sample_audit.rename_axis("sample").reset_index().to_csv(tables / "11_reference_sample_enrichment.tsv", sep="\t", index=False)

# Query technical associations after removing the population median.
q = libraries.merge(mapping_qc, on=["gsm", "population"], how="left").merge(matrix_qc, on="gsm", how="left")
if "mapped_deduplicated_blacklist_filtered_fragments" in q:
    q["peak_overlap_fraction"] = q["total_peak_overlaps"] / q["mapped_deduplicated_blacklist_filtered_fragments"]
query_predictors = [x for x in ["mapped_deduplicated_blacklist_filtered_fragments", "total_peak_overlaps", "peak_overlap_fraction", "matrix_total", "nonzero_peaks", "nonzero_fraction"] if x in q]
query_outcomes = ["distance_to_expected_centroid", "expected_provisional_label_neighbor_weight", "cluster6_neighbor_weight", "median_neighbor_distance", "mapped_native_l4", "query_umap_1", "query_umap_2"]
replicated = q["population"].map(q["population"].value_counts()) >= 2
for predictor in query_predictors:
    for outcome in query_outcomes:
        qr = q.loc[replicated]
        r, p, n = safe_spearman(centered(qr[predictor], qr["population"]), centered(qr[outcome], qr["population"]))
        tech_rows.append(dict(level="immgen_population_centered", predictor=predictor, outcome=outcome, n=n, spearman=r, p_value=p,
            cluster6_median=np.nan, other_median=np.nan, cluster6_median_difference_in_other_mad=np.nan))
q.to_csv(tables / "12_library_position_and_technical_qc.tsv", sep="\t", index=False)

# Labelled Bdev cells: within-label QC association with distance from their own label centroid.
label_coord = labels.loc[labels["label"].isin(centroids.index) & labels["cell_id"].isin(coords.index)].copy()
label_coord["distance_to_own_label_centroid"] = [np.linalg.norm(coords.loc[c, latent_cols].to_numpy(float) - centroids.loc[l, latent_cols].to_numpy(float))
    for c, l in zip(label_coord["cell_id"], label_coord["label"])]
for predictor in technical:
    if predictor not in label_coord:
        continue
    r, p, n = safe_spearman(centered(label_coord[predictor], label_coord["label"]), centered(label_coord["distance_to_own_label_centroid"], label_coord["label"]))
    tech_rows.append(dict(level="rna_label_cells_label_centered", predictor=predictor, outcome="distance_to_own_label_centroid", n=n,
        spearman=r, p_value=p, cluster6_median=np.nan, other_median=np.nan, cluster6_median_difference_in_other_mad=np.nan))
technical_audit = pd.DataFrame(tech_rows)
technical_audit.to_csv(tables / "13_technical_bias_associations.tsv", sep="\t", index=False)

# Link individual latent dimensions to technical covariates after biological-group centering.
latent_tech_rows = []
ref_latent = ref.merge(coords[latent_cols], left_on="cell_id", right_index=True, how="left")
for d in latent_cols:
    for predictor in technical:
        r, p, n = safe_spearman(
            centered(ref_latent[predictor], ref_latent[cfg["diagnostic"]["cluster_field"]].astype(str)),
            centered(ref_latent[d], ref_latent[cfg["diagnostic"]["cluster_field"]].astype(str)))
        latent_tech_rows.append(dict(level="reference_cluster_centered", latent_dimension=d, predictor=predictor, n=n, spearman=r, p_value=p))
q_latent = q.merge(coords[latent_cols], left_on="gsm", right_index=True, how="left")
q_latent_replicated = q_latent.loc[q_latent["population"].map(q_latent["population"].value_counts()) >= 2]
for d in latent_cols:
    for predictor in query_predictors:
        r, p, n = safe_spearman(centered(q_latent_replicated[predictor], q_latent_replicated["population"]),
                                centered(q_latent_replicated[d], q_latent_replicated["population"]))
        latent_tech_rows.append(dict(level="immgen_population_centered", latent_dimension=d, predictor=predictor, n=n, spearman=r, p_value=p))
label_latent = label_coord.merge(coords[latent_cols], left_on="cell_id", right_index=True, how="left")
for d in latent_cols:
    for predictor in technical:
        r, p, n = safe_spearman(centered(label_latent[predictor], label_latent["label"]), centered(label_latent[d], label_latent["label"]))
        latent_tech_rows.append(dict(level="rna_label_cells_label_centered", latent_dimension=d, predictor=predictor, n=n, spearman=r, p_value=p))
latent_technical = pd.DataFrame(latent_tech_rows)
latent_technical.to_csv(tables / "16_latent_dimension_technical_associations.tsv", sep="\t", index=False)

# Recurring candidate peaks and flags for assay-wide/technical behavior.
peaks = pd.read_csv(tables / "04_peak_driver_candidates.tsv.gz", sep="\t")
peak_summary = peaks.groupby("feature_id").agg(libraries=("gsm", "nunique"), populations=("population", "nunique"),
    mean_alignment_score=("alignment_score", "mean"), median_rank=("rank", "median"), chromosome=("chromosome", "first"),
    width=("width", "first"), reference_prevalence=("reference_prevalence", "mean"), immgen_bulk_prevalence=("immgen_bulk_prevalence", "mean")).reset_index()
peak_summary["ubiquitous_bulk_flag"] = peak_summary["immgen_bulk_prevalence"] >= 0.9
peak_summary["bulk_reference_prevalence_gap"] = peak_summary["immgen_bulk_prevalence"] - peak_summary["reference_prevalence"]
peak_summary["noncanonical_contig_flag"] = ~peak_summary["chromosome"].astype(str).str.match(r"^chr(?:[1-9]|1[0-9]|X|Y)$")
peak_summary = peak_summary.sort_values(["libraries", "populations", "mean_alignment_score"], ascending=[False, False, False])
peak_summary.to_csv(tables / "14_recurring_peak_driver_candidates.tsv.gz", sep="\t", index=False, compression="gzip")
axis_audit = pd.read_csv(tables / "15_peak_axis_prevalence_audit.tsv.gz", sep="\t")
axis_r, axis_p, axis_n = safe_spearman(axis_audit["reference_prevalence"], axis_audit["immgen_bulk_prevalence"])
axis_zero_common = int(((axis_audit["immgen_bulk_prevalence"] == 0) & (axis_audit["reference_prevalence"] >= 0.5)).sum())
axis_one_rare = int(((axis_audit["immgen_bulk_prevalence"] == 1) & (axis_audit["reference_prevalence"] <= 0.1)).sum())
recurring = peak_summary.loc[peak_summary["libraries"] >= 5]
recurring_all_or_none = int(recurring["immgen_bulk_prevalence"].isin([0, 1]).sum())
dim_summary = dims.groupby("latent_dimension").agg(mean_advantage=("squared_distance_advantage_observed_vs_expected", "mean"),
    median_absolute_fraction=("absolute_fraction_of_pull", "median")).sort_values("mean_advantage", ascending=False)
dominant_dim = str(dim_summary.index[0])
dominant_dim_tech = latent_technical.loc[(latent_technical["latent_dimension"] == dominant_dim) &
    (latent_technical["level"] == "immgen_population_centered") & latent_technical["spearman"].notna()].copy()
dominant_dim_tech = dominant_dim_tech.reindex(dominant_dim_tech["spearman"].abs().sort_values(ascending=False).index)

# Figures.
heat = dims.pivot(index="gsm", columns="latent_dimension", values="squared_distance_advantage_observed_vs_expected")
heat = heat.loc[[x for x in libraries.sort_values("population")["gsm"] if x in heat.index]]
fig, ax = plt.subplots(figsize=(10, 9)); im = ax.imshow(heat, aspect="auto", cmap="coolwarm", vmin=-np.nanmax(abs(heat.to_numpy())), vmax=np.nanmax(abs(heat.to_numpy())))
ax.set_xticks(range(len(heat.columns)), heat.columns, rotation=45, ha="right"); ax.set_yticks(range(len(heat.index)), heat.index, fontsize=7)
ax.set_title("Latent dimensions pulling ImmGen libraries toward observed neighbors\npositive = closer to observed neighbors than provisional RNA-label centroid")
fig.colorbar(im, ax=ax, label="Squared-distance advantage"); fig.tight_layout(); fig.savefig(figures / "01_latent_dimension_pull.png", dpi=220); plt.close(fig)

position_plot = libraries.loc[libraries["has_sufficient_expected_label_centroid"]].copy()
fig, ax = plt.subplots(figsize=(8, 7)); sc = ax.scatter(position_plot["distance_to_expected_centroid"], position_plot["distance_to_observed_neighbor_centroid"],
    c=position_plot["cluster6_neighbor_weight"], cmap="viridis", s=45)
limit = max(ax.get_xlim()[1], ax.get_ylim()[1]); ax.plot([0, limit], [0, limit], "--", color="grey")
for _, row in position_plot.iterrows(): ax.annotate(row["population"], (row["distance_to_expected_centroid"], row["distance_to_observed_neighbor_centroid"]), fontsize=6, alpha=.8)
ax.set(xlabel="Distance to corresponding provisional RNA-label centroid", ylabel="Distance to observed neighbor centroid", title="Expected-label versus observed-neighborhood position")
fig.colorbar(sc, ax=ax, label="Cluster 6 neighbor weight"); fig.tight_layout(); fig.savefig(figures / "02_expected_vs_observed_position.png", dpi=220); plt.close(fig)

plot_tech = technical_audit.loc[technical_audit["spearman"].notna()].copy()
plot_tech["key"] = plot_tech["level"] + ": " + plot_tech["predictor"] + " -> " + plot_tech["outcome"]
plot_tech = plot_tech.reindex(plot_tech["spearman"].abs().sort_values(ascending=False).index).head(25).sort_values("spearman")
fig, ax = plt.subplots(figsize=(10, 8)); colors = np.where(plot_tech["spearman"] >= 0, "#B2182B", "#2166AC")
ax.barh(plot_tech["key"], plot_tech["spearman"], color=colors); ax.axvline(0, color="black", lw=.7); ax.set_xlabel("Spearman correlation")
ax.set_title("Largest technical associations (diagnostic; population/label centered where indicated)"); fig.tight_layout(); fig.savefig(figures / "03_technical_bias_associations.png", dpi=220); plt.close(fig)

summary = {
    "schema_version": "bdev.immgen_position_driver_summary.v1", "status": "complete_diagnostic", "run_identity": cfg["run_identity"],
    "rna_label_centroids": int(len(centroids)), "libraries_evaluated": int(len(libraries)), "libraries_with_expected_label_centroid": int(libraries["has_sufficient_expected_label_centroid"].sum()), "reference_cells": int(len(ref)),
    "median_expected_label_neighbor_weight": float(libraries["expected_provisional_label_neighbor_weight"].median()),
    "median_cluster6_neighbor_weight": float(libraries["cluster6_neighbor_weight"].median()),
    "median_observed_over_expected_distance_ratio": float(libraries["observed_over_expected_distance_ratio"].median()),
    "reference_sample_max_enrichment": float(sample_audit["enrichment_ratio"].max()),
    "reference_sample_max_enrichment_sample": str(sample_audit["enrichment_ratio"].idxmax()),
    "recurring_peak_candidates_5plus_libraries": int((peak_summary["libraries"] >= 5).sum()),
    "recurring_peak_candidates_5plus_all_or_none_bulk": recurring_all_or_none,
    "peak_axis_reference_vs_bulk_prevalence_spearman": axis_r,
    "peak_axis_absent_all_bulk_common_reference": axis_zero_common,
    "peak_axis_present_all_bulk_rare_reference": axis_one_rare,
    "dominant_pull_latent_dimension": dominant_dim,
    "dominant_pull_dimension_median_absolute_fraction": float(dim_summary.loc[dominant_dim, "median_absolute_fraction"]),
    "dominant_pull_dimension_strongest_query_technical_predictor": str(dominant_dim_tech.iloc[0]["predictor"]) if len(dominant_dim_tech) else None,
    "dominant_pull_dimension_strongest_query_technical_spearman": float(dominant_dim_tech.iloc[0]["spearman"]) if len(dominant_dim_tech) else None,
    "labels_accepted": False, "peak_candidates_causal": False, "scientific_status": "diagnostic_only_mapping_not_accepted",
    "source_manifests_unchanged": True, "paper_analysis": False, "implementation_unlocked": False, "sif_built": False,
}
(reports / "02_position_driver_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
(logs / "02_analyze_attribution.log").write_text(json.dumps(summary, indent=2) + "\n")
(manifests / "02_analysis_receipt.json").write_text(json.dumps({**summary, "config_path": str(config_path), "config_sha256": sha256(config_path)}, indent=2) + "\n")
print(json.dumps(summary, indent=2))
