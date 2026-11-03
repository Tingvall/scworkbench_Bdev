#!/usr/bin/env python3
"""Discover stable broad programs from ImmGen and evaluate them without Bdev tuning."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter, defaultdict
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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def empirical_percentile(training: np.ndarray, value: float) -> float:
    return float(
        (np.sum(training < value) + 0.5 * np.sum(training == value) + 0.5)
        / (len(training) + 1)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    require(config["schema_version"] == "bdev.stable_immgen_broad_genesets.v1", "Unexpected schema")

    inputs = config["inputs"]
    for name in ("selected_counts", "metadata", "fixture_cells"):
        path = Path(inputs[name]["path"])
        require(path.is_file(), f"Missing input {name}: {path}")
        require(sha256(path) == inputs[name]["sha256"], f"Checksum mismatch: {name}")
    bundle_manifest = Path(inputs["source_bundle"]["path"]) / "manifest.json"
    require(sha256(bundle_manifest) == inputs["source_bundle"]["manifest_sha256"], "Bundle manifest changed")

    root = Path(config["storage"]["results_root"])
    require(not root.exists(), f"Refusing to overwrite run: {root}")
    tables, reports, manifests = (root / name for name in ("01_tables", "03_reports", "06_manifests"))
    for directory in (tables, reports, manifests):
        directory.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(inputs["selected_counts"]["path"], sep="\t", index_col=0)
    metadata = pd.read_csv(inputs["metadata"]["path"], sep="\t", dtype=str).fillna("")
    require(list(counts.columns) == list(metadata["library_id"]), "Counts and metadata are misaligned")
    population_metadata = metadata.drop_duplicates("population").set_index("population")
    require(not population_metadata.index.has_duplicates, "Population metadata are not unique")
    lineages = list(config["discovery"]["lineages"])
    require(set(population_metadata["table_s1_lineage"]) == set(lineages), "Lineage inventory changed")

    max_rank = int(config["discovery"]["ucell_max_rank"])
    ranks = counts.rank(axis=0, method="average", ascending=False).clip(upper=max_rank)
    contributions = (max_rank - ranks) / (max_rank - 1)
    population_contributions = (
        contributions.T.assign(population=metadata["population"].to_numpy())
        .groupby("population", sort=False).median().T
    )
    target_quantile = float(config["discovery"]["target_consistency_quantile"])
    competitor_quantile = float(config["discovery"]["competitor_typical_quantile"])
    maximum_genes = int(config["discovery"]["maximum_genes_per_lineage"])

    def select_markers(training_populations: list[str], target: str) -> pd.DataFrame:
        training_lineage = population_metadata.loc[training_populations, "table_s1_lineage"]
        target_populations = list(training_lineage[training_lineage == target].index)
        require(bool(target_populations), f"No training populations for {target}")
        target_value = population_contributions[target_populations].quantile(target_quantile, axis=1)
        competitor_values = []
        for competitor in lineages:
            if competitor == target:
                continue
            competitor_populations = list(training_lineage[training_lineage == competitor].index)
            competitor_values.append(
                population_contributions[competitor_populations]
                .quantile(competitor_quantile, axis=1).rename(competitor)
            )
        competitor_envelope = pd.concat(competitor_values, axis=1).max(axis=1)
        ranking = pd.DataFrame({
            "gene_symbol": counts.index,
            "target_consistency": target_value.to_numpy(),
            "competitor_envelope": competitor_envelope.to_numpy(),
        })
        ranking["specificity_margin"] = ranking["target_consistency"] - ranking["competitor_envelope"]
        ranking = ranking[
            (ranking["target_consistency"] > 0)
            & (ranking["specificity_margin"] > 0)
        ].sort_values(
            ["specificity_margin", "target_consistency", "gene_symbol"],
            ascending=[False, False, True], kind="mergesort"
        ).head(maximum_genes).reset_index(drop=True)
        ranking["selection_rank"] = np.arange(1, len(ranking) + 1)
        require(len(ranking) >= 5, f"Fewer than five positive stable candidates for {target}")
        return ranking

    prediction_rows: list[dict[str, object]] = []
    selection_records: dict[str, list[tuple[str, pd.DataFrame]]] = defaultdict(list)
    all_populations = list(population_metadata.index)
    for heldout_population in all_populations:
        training_populations = [value for value in all_populations if value != heldout_population]
        fold_sets: dict[str, list[str]] = {}
        for lineage in lineages:
            ranking = select_markers(training_populations, lineage)
            fold_sets[lineage] = list(ranking["gene_symbol"])
            if population_metadata.loc[heldout_population, "table_s1_lineage"] == lineage:
                selection_records[lineage].append((heldout_population, ranking))

        training_ids = list(metadata.loc[metadata["population"] != heldout_population, "library_id"])
        heldout_ids = list(metadata.loc[metadata["population"] == heldout_population, "library_id"])
        training_scores = pd.DataFrame({
            lineage: contributions.loc[genes, training_ids].mean(axis=0)
            for lineage, genes in fold_sets.items()
        })
        heldout_scores = pd.DataFrame({
            lineage: contributions.loc[genes, heldout_ids].mean(axis=0)
            for lineage, genes in fold_sets.items()
        })
        true_lineage = population_metadata.loc[heldout_population, "table_s1_lineage"]
        for library_id in heldout_ids:
            calibrated = {
                lineage: empirical_percentile(
                    training_scores[lineage].to_numpy(), float(heldout_scores.loc[library_id, lineage])
                )
                for lineage in lineages
            }
            ordered = sorted(lineages, key=lambda value: (-calibrated[value], value))
            winner, runner = ordered[:2]
            prediction_rows.append({
                "heldout_library": library_id,
                "heldout_population": heldout_population,
                "true_lineage": true_lineage,
                "winner_lineage": winner,
                "winner_correct": winner == true_lineage,
                "top_score": calibrated[winner],
                "runner_up_score": calibrated[runner],
                "score_margin": calibrated[winner] - calibrated[runner],
                "selected_genes_winner": len(fold_sets[winner]),
                "selected_genes_true_lineage": len(fold_sets[true_lineage]),
            })
    predictions = pd.DataFrame(prediction_rows)
    predictions.to_csv(tables / "03_immgen_leave_population_out_predictions.tsv", sep="\t", index=False)

    membership_rows: list[dict[str, object]] = []
    final_sets: dict[str, list[str]] = {}
    stability_rows: list[dict[str, object]] = []
    for lineage in lineages:
        records = selection_records[lineage]
        frequencies: Counter[str] = Counter()
        ranks_by_gene: dict[str, list[int]] = defaultdict(list)
        fold_metrics: dict[str, list[tuple[float, float, float]]] = defaultdict(list)
        for heldout_population, ranking in records:
            for row in ranking.itertuples(index=False):
                frequencies[row.gene_symbol] += 1
                ranks_by_gene[row.gene_symbol].append(int(row.selection_rank))
                fold_metrics[row.gene_symbol].append(
                    (float(row.target_consistency), float(row.competitor_envelope), float(row.specificity_margin))
                )
        ordered = sorted(
            frequencies,
            key=lambda gene: (
                -frequencies[gene],
                float(np.mean(ranks_by_gene[gene])),
                gene,
            ),
        )[:maximum_genes]
        final_sets[lineage] = ordered
        folds = len(records)
        for final_rank, gene in enumerate(ordered, 1):
            values = np.asarray(fold_metrics[gene])
            membership_rows.append({
                "lineage": lineage,
                "gene_symbol": gene,
                "final_rank": final_rank,
                "target_population_folds": folds,
                "folds_selected": frequencies[gene],
                "selection_frequency": frequencies[gene] / folds,
                "mean_fold_rank": float(np.mean(ranks_by_gene[gene])),
                "mean_target_consistency": float(values[:, 0].mean()),
                "mean_competitor_envelope": float(values[:, 1].mean()),
                "mean_specificity_margin": float(values[:, 2].mean()),
            })
        stability_rows.append({
            "lineage": lineage,
            "ImmGen_populations": folds,
            "final_genes": len(ordered),
            "genes_selected_in_at_least_80pct_folds": sum(
                frequencies[gene] / folds >= 0.8 for gene in ordered
            ),
            "minimum_final_selection_frequency": min(frequencies[gene] / folds for gene in ordered),
        })
    membership = pd.DataFrame(membership_rows)
    with gzip.open(tables / "01_immgen_stable_broad_membership.tsv.gz", "wt", encoding="utf-8") as handle:
        membership.to_csv(handle, sep="\t", index=False)
    with (tables / "02_immgen_stable_broad_sets.gmt").open("w", encoding="utf-8") as handle:
        for lineage in lineages:
            set_id = "IMMGEN_STABLE_BROAD_" + "".join(c if c.isalnum() else "_" for c in lineage).upper()
            handle.write("\t".join([set_id, f"ImmGen-only stable broad lineage: {lineage}", *final_sets[lineage]]) + "\n")

    score_minimum = float(config["internal_evaluation"]["score_minimum"])
    fixed_margin = float(config["internal_evaluation"]["fixed_margin"])
    predictions["fixed_gate_accepted"] = (
        (predictions["top_score"] >= score_minimum)
        & (predictions["score_margin"] >= fixed_margin)
    )
    lineage_summary = predictions.groupby("true_lineage", sort=False).agg(
        heldout_libraries=("heldout_library", "size"),
        winner_correct=("winner_correct", "sum"),
        winner_accuracy=("winner_correct", "mean"),
        fixed_gate_accepted=("fixed_gate_accepted", "sum"),
    ).reset_index()
    correct_accepted = predictions[predictions["fixed_gate_accepted"]].groupby("true_lineage")["winner_correct"].sum()
    lineage_summary["fixed_gate_correct"] = lineage_summary["true_lineage"].map(correct_accepted).fillna(0).astype(int)
    lineage_summary["fixed_gate_precision"] = np.where(
        lineage_summary["fixed_gate_accepted"] > 0,
        lineage_summary["fixed_gate_correct"] / lineage_summary["fixed_gate_accepted"],
        np.nan,
    )
    lineage_summary.to_csv(tables / "04_immgen_internal_lineage_summary.tsv", sep="\t", index=False)
    pd.DataFrame(stability_rows).to_csv(tables / "05_final_set_stability_summary.tsv", sep="\t", index=False)

    accepted = predictions[predictions["fixed_gate_accepted"]]
    summary = {
        "schema_version": "bdev.stable_immgen_broad_internal_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "data_source": "ImmGen_only",
        "libraries": len(predictions),
        "populations": len(all_populations),
        "lineages": len(lineages),
        "nested_leave_population_out_winner_correct": int(predictions["winner_correct"].sum()),
        "nested_leave_population_out_winner_accuracy": float(predictions["winner_correct"].mean()),
        "fixed_gate_score_minimum": score_minimum,
        "fixed_gate_margin": fixed_margin,
        "fixed_gate_accepted": len(accepted),
        "fixed_gate_correct": int(accepted["winner_correct"].sum()),
        "fixed_gate_precision": float(accepted["winner_correct"].mean()) if len(accepted) else None,
        "final_set_sizes": {lineage: len(final_sets[lineage]) for lineage in lineages},
        "Bdev_used_for_selection": False,
        "accepted_gene_sets_emitted": False,
    }
    (reports / "01_immgen_internal_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    manifest = {
        "schema_version": "bdev.stable_immgen_broad_input_manifest.v1",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "inputs": {
            name: {"path": spec["path"], "sha256": sha256(Path(spec["path"]))}
            for name, spec in inputs.items() if name != "source_bundle"
        },
        "source_bundle": {
            "path": inputs["source_bundle"]["path"],
            "manifest_sha256": sha256(bundle_manifest),
            "fingerprint": inputs["source_bundle"]["fingerprint"],
        },
    }
    (manifests / "01_input_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
