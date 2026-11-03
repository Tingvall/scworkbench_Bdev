#!/usr/bin/env python3
"""Generate simple ImmGen-only population sets and compare them post hoc."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


DIRECT_LABEL_MAP = {
    "LTHSC": "LTHSC.34-.BM", "STHSC": "STHSC.150-.BM",
    "MPP3": "MPP3.48+.BM", "MPP4": "MPP4.135+.BM",
    "CLP": "proB.CLP.BM", "FrA": "proB.FrA.BM",
    "FrBC": "proB.FrBC.BM", "FrE": "B.FrE.BM", "B1b": "B1b.PC",
    "DN1": "preT.DN1.Th", "DN2a": "preT.DN2a.Th",
    "DN2b": "preT.DN2b.Th", "DN3": "preT.DN3.Th",
    "DN4": "T.DN4.Th", "ISP": "T.ISP.Th", "DP": "T.DP.Th",
    "Treg": "Treg.4.25hi.Sp", "NKT": "NKT.Sp", "ILC2": "ILC2.SI",
    "pDC": "DC.pDC.Sp", "D_CD4": "DC.4+.Sp", "D_CD8": "DC.8+.Sp",
    "Mo_Ly6Cpos": "Mo.6C+II-.Bl", "Mo_Ly6Cneg": "Mo.6C-II-.Bl",
    "Neutrophil": "GN.BM", "BEC": "BEC.SLN", "LEC": "LEC.SLN",
    "FRC": "FRC.SLN", "IAP": "IAP.SLN",
}


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
    inputs = config["inputs"]
    for name in ("selected_counts", "metadata", "supplied_sets_for_posthoc_comparison"):
        path = Path(inputs[name]["path"])
        if sha256(path) != inputs[name]["sha256"]:
            raise RuntimeError(f"Input checksum changed: {name}")

    counts = pd.read_csv(inputs["selected_counts"]["path"], sep="\t", index_col=0)
    metadata = pd.read_csv(inputs["metadata"]["path"], sep="\t", dtype=str).fillna("")
    if list(counts.columns) != list(metadata["library_id"]):
        raise RuntimeError("Count columns and metadata are misaligned")
    if counts.index.duplicated().any():
        raise RuntimeError("Gene symbols are duplicated")
    populations = list(dict.fromkeys(metadata["population"]))
    if len(populations) != int(config["discovery"]["expected_populations"]):
        raise RuntimeError("Population inventory changed")

    library_sizes = counts.sum(axis=0).replace(0, np.nan)
    log_cpm = np.log2(counts.divide(library_sizes, axis=1) * 1_000_000 + 0.5)
    profiles = log_cpm.T.assign(population=metadata["population"].to_numpy()).groupby("population", sort=False).median().T
    floor = float(config["discovery"]["minimum_target_log2_cpm"])
    set_size = int(config["discovery"]["genes_per_population"])

    ranking_rows: list[pd.DataFrame] = []
    selected: dict[str, list[str]] = {}
    for population in populations:
        target = profiles[population]
        competitor = profiles.drop(columns=population).mean(axis=1)
        ranking = pd.DataFrame({
            "population": population,
            "gene_symbol": profiles.index,
            "target_log2_cpm": target.to_numpy(),
            "competitor_mean_log2_cpm": competitor.to_numpy(),
        })
        ranking["specificity_contrast"] = ranking["target_log2_cpm"] - ranking["competitor_mean_log2_cpm"]
        ranking = ranking[(ranking["target_log2_cpm"] >= floor) & (ranking["specificity_contrast"] > 0)].sort_values(
            ["specificity_contrast", "target_log2_cpm", "gene_symbol"],
            ascending=[False, False, True], kind="mergesort").reset_index(drop=True)
        ranking["rank"] = np.arange(1, len(ranking) + 1)
        if len(ranking) < set_size:
            raise RuntimeError(f"Too few eligible genes for {population}")
        selected[population] = list(ranking.head(set_size)["gene_symbol"])
        ranking_rows.append(ranking.head(100).copy())

    root = Path(config["storage"]["results_root"])
    tables, reports = root / "01_tables", root / "03_reports"
    tables.mkdir(parents=True, exist_ok=True); reports.mkdir(parents=True, exist_ok=True)
    rankings = pd.concat(ranking_rows, ignore_index=True)
    rankings.to_csv(tables / "01_top100_rankings.tsv.gz", sep="\t", index=False, compression="gzip")
    membership = rankings[rankings["rank"] <= set_size].copy()
    membership.to_csv(tables / "02_selected_membership.tsv", sep="\t", index=False)
    with (tables / "03_simple_mean_contrast_sets.gmt").open("w", encoding="utf-8") as handle:
        for population in populations:
            handle.write("\t".join([population, "ImmGen population-balanced mean-contrast top20", *selected[population]]) + "\n")

    supplied = pd.read_csv(inputs["supplied_sets_for_posthoc_comparison"]["path"], sep="\t", dtype=str)
    comparison_rows = []
    for supplied_label, population in DIRECT_LABEL_MAP.items():
        ranking = rankings[rankings["population"] == population].set_index("gene_symbol")["rank"]
        ranks = ranking.reindex(list(supplied[supplied_label]))
        comparison_rows.append({
            "supplied_label": supplied_label, "immgen_population": population,
            "supplied_genes": len(ranks), "overlap_top20": int((ranks <= 20).sum()),
            "overlap_top50": int((ranks <= 50).sum()), "overlap_top100": int((ranks <= 100).sum()),
            "median_rank_within_top100_or_missing": float(ranks.fillna(101).median()),
        })
    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(tables / "04_posthoc_supplied_set_similarity.tsv", sep="\t", index=False)
    summary = {
        "schema_version": "bdev.simple_mean_contrast_discovery_summary.v1",
        "status": "complete_immgen_only_discovery",
        "run_identity": config["run_identity"],
        "populations": len(populations), "genes_per_population": set_size,
        "selection_used_supplied_sets": False, "selection_used_bdev": False,
        "posthoc_direct_label_matches": len(comparison),
        "posthoc_overlap_top20": int(comparison["overlap_top20"].sum()),
        "posthoc_overlap_top50": int(comparison["overlap_top50"].sum()),
        "posthoc_overlap_top100": int(comparison["overlap_top100"].sum()),
        "posthoc_supplied_genes": int(comparison["supplied_genes"].sum()),
    }
    (reports / "01_discovery_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

