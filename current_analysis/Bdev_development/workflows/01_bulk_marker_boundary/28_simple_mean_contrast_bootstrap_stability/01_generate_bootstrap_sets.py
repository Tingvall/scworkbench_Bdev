#!/usr/bin/env python3
"""Generate deterministic within-population ImmGen bootstrap gene sets."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gmt(path: Path) -> dict[str, list[str]]:
    result = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        result[fields[0]] = fields[2:]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    cfg = yaml.safe_load(config_path.read_text())
    if cfg["schema_version"] != "bdev.simple_mean_contrast_bootstrap_stability.v1":
        raise RuntimeError("Unexpected config schema")
    for name in ("selected_counts", "metadata", "baseline_gmt"):
        path = Path(cfg["inputs"][name]["path"])
        if sha256(path) != cfg["inputs"][name]["sha256"]:
            raise RuntimeError(f"Input checksum changed: {name}")

    counts = pd.read_csv(cfg["inputs"]["selected_counts"]["path"], sep="\t", index_col=0)
    metadata = pd.read_csv(cfg["inputs"]["metadata"]["path"], sep="\t", dtype=str).fillna("")
    if list(counts.columns) != list(metadata["library_id"]):
        raise RuntimeError("Counts and metadata are misaligned")
    populations = list(dict.fromkeys(metadata["population"]))
    if len(populations) != int(cfg["discovery"]["expected_populations"]):
        raise RuntimeError("Population inventory changed")
    baseline = read_gmt(Path(cfg["inputs"]["baseline_gmt"]["path"]))
    if list(baseline) != populations or any(len(genes) != 20 for genes in baseline.values()):
        raise RuntimeError("Baseline set inventory changed")

    library_sizes = counts.sum(axis=0).replace(0, np.nan)
    log_cpm = np.log2(counts.divide(library_sizes, axis=1) * 1_000_000 + 0.5)
    pop_indices = {
        population: np.flatnonzero(metadata["population"].to_numpy() == population)
        for population in populations
    }
    rounds = int(cfg["bootstrap"]["rounds"])
    seed = int(cfg["bootstrap"]["random_seed"])
    rng = np.random.default_rng(seed)
    floor = float(cfg["discovery"]["minimum_target_log2_cpm"])
    set_size = int(cfg["discovery"]["genes_per_population"])

    membership_rows = []
    plan_rows = []
    stability_rows = []
    frequency: dict[tuple[str, str], int] = {}
    for bootstrap_round in range(1, rounds + 1):
        sampled_profiles = {}
        for population in populations:
            source = pop_indices[population]
            sampled = rng.choice(source, size=len(source), replace=True)
            sampled_profiles[population] = log_cpm.iloc[:, sampled].median(axis=1)
            plan_rows.append({
                "bootstrap_round": bootstrap_round,
                "population": population,
                "population_libraries": len(source),
                "source_library_ids": ";".join(metadata.iloc[source]["library_id"]),
                "sampled_library_ids": ";".join(metadata.iloc[sampled]["library_id"]),
                "distinct_sampled_libraries": len(np.unique(sampled)),
                "stability_evidence": "structurally_unassessable_singleton" if len(source) == 1 else "bootstrap_assessable",
            })
        profiles = pd.DataFrame(sampled_profiles, index=log_cpm.index)
        for population in populations:
            target = profiles[population]
            competitor = profiles.drop(columns=population).mean(axis=1)
            ranking = pd.DataFrame({
                "gene_symbol": profiles.index,
                "target_log2_cpm": target.to_numpy(),
                "competitor_mean_log2_cpm": competitor.to_numpy(),
            })
            ranking["specificity_contrast"] = ranking["target_log2_cpm"] - ranking["competitor_mean_log2_cpm"]
            ranking = ranking[(ranking["target_log2_cpm"] >= floor) & (ranking["specificity_contrast"] > 0)].sort_values(
                ["specificity_contrast", "target_log2_cpm", "gene_symbol"],
                ascending=[False, False, True], kind="mergesort").head(set_size).reset_index(drop=True)
            if len(ranking) != set_size:
                raise RuntimeError(f"Too few genes for {population} in round {bootstrap_round}")
            selected = list(ranking["gene_symbol"])
            overlap = len(set(selected).intersection(baseline[population]))
            union = len(set(selected).union(baseline[population]))
            n_libraries = len(pop_indices[population])
            stability_rows.append({
                "bootstrap_round": bootstrap_round,
                "population": population,
                "population_libraries": n_libraries,
                "stability_evidence": "structurally_unassessable_singleton" if n_libraries == 1 else "bootstrap_assessable",
                "overlap_baseline_top20": overlap,
                "jaccard_baseline_top20": overlap / union,
            })
            for rank, row in ranking.iterrows():
                gene = str(row["gene_symbol"])
                membership_rows.append({
                    "bootstrap_round": bootstrap_round,
                    "population": population,
                    "rank": rank + 1,
                    "gene_symbol": gene,
                    "target_log2_cpm": row["target_log2_cpm"],
                    "competitor_mean_log2_cpm": row["competitor_mean_log2_cpm"],
                    "specificity_contrast": row["specificity_contrast"],
                })
                frequency[(population, gene)] = frequency.get((population, gene), 0) + 1

    root = Path(cfg["storage"]["results_root"])
    tables = root / "01_tables"
    reports = root / "03_reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    membership = pd.DataFrame(membership_rows)
    membership.to_csv(tables / "01_bootstrap_set_membership.tsv.gz", sep="\t", index=False, compression="gzip")
    plan = pd.DataFrame(plan_rows)
    plan.to_csv(tables / "02_bootstrap_plan.tsv.gz", sep="\t", index=False, compression="gzip")
    stability = pd.DataFrame(stability_rows)
    stability.to_csv(tables / "03_population_geneset_stability_by_round.tsv", sep="\t", index=False)
    frequencies = pd.DataFrame([
        {"population": population, "gene_symbol": gene, "selected_rounds": count,
         "selection_frequency": count / rounds, "baseline_member": gene in baseline[population]}
        for (population, gene), count in frequency.items()
    ]).sort_values(["population", "selection_frequency", "gene_symbol"], ascending=[True, False, True])
    frequencies.to_csv(tables / "04_gene_selection_frequency.tsv.gz", sep="\t", index=False, compression="gzip")

    summary = {
        "schema_version": "bdev.simple_mean_contrast_bootstrap_discovery_summary.v1",
        "status": "complete_ImmGen_only_bootstrap_discovery",
        "run_identity": cfg["run_identity"],
        "rounds": rounds,
        "random_seed": seed,
        "populations": len(populations),
        "singleton_populations": int(sum(len(indices) == 1 for indices in pop_indices.values())),
        "bootstrap_assessable_populations": int(sum(len(indices) > 1 for indices in pop_indices.values())),
        "bdev_used_for_gene_selection": False,
        "default_changed_by_run": False,
    }
    (reports / "01_bootstrap_discovery_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
