#!/usr/bin/env python3
"""Build ImmGen-only lineage-anchored narrow modules and test target holdouts."""

from __future__ import annotations

import argparse
import gzip
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
    require(config["schema_version"] == "bdev.lineage_anchored_narrow_genesets.v1", "Unexpected schema")
    inputs = config["inputs"]
    for name in ("selected_counts", "metadata", "fold_plan", "fixture_cells"):
        path = Path(inputs[name]["path"])
        require(path.is_file(), f"Missing input {name}: {path}")
        require(sha256(path) == inputs[name]["sha256"], f"Checksum mismatch: {name}")
    bundle_manifest = Path(inputs["source_bundle"]["path"]) / "manifest.json"
    require(sha256(bundle_manifest) == inputs["source_bundle"]["manifest_sha256"], "Bundle manifest changed")

    root = Path(config["storage"]["results_root"])
    require(not root.exists(), f"Refusing to overwrite run: {root}")
    tables, reports, manifests = (root / value for value in ("01_tables", "03_reports", "06_manifests"))
    for directory in (tables, reports, manifests):
        directory.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(inputs["selected_counts"]["path"], sep="\t", index_col=0)
    metadata = pd.read_csv(inputs["metadata"]["path"], sep="\t", dtype=str).fillna("")
    plans = pd.read_csv(inputs["fold_plan"]["path"], sep="\t", keep_default_na=False)
    require(list(counts.columns) == list(metadata["library_id"]), "Counts and metadata are misaligned")
    population_metadata = metadata.drop_duplicates("population").set_index("population")
    populations = list(population_metadata.index)
    lineage = population_metadata["table_s1_lineage"].to_dict()
    lineages = sorted(set(lineage.values()))
    require(len(populations) == int(config["discovery"]["populations"]), "Population inventory changed")
    require(len(lineages) == int(config["discovery"]["lineages"]), "Lineage inventory changed")

    max_rank = int(config["discovery"]["ucell_max_rank"])
    ranks = counts.rank(axis=0, method="average", ascending=False).clip(upper=max_rank)
    contributions = (max_rank - ranks) / (max_rank - 1)
    population_contributions = (
        contributions.T.assign(population=metadata["population"].to_numpy())
        .groupby("population", sort=False).median().T
    )
    target_quantile = float(config["discovery"]["target_consistency_quantile"])
    anchor_size = int(config["discovery"]["anchor_genes_per_population"])
    local_size = int(config["discovery"]["local_genes_per_population"])

    def derive_modules(target: str, heldout_library: str | None = None) -> dict[str, pd.DataFrame]:
        target_ids = list(metadata.loc[metadata["population"] == target, "library_id"])
        training_ids = [value for value in target_ids if value != heldout_library]
        require(bool(training_ids), f"No target training library for {target}")
        target_value = contributions[training_ids].quantile(target_quantile, axis=1)
        target_lineage = lineage[target]
        other_lineage_values = []
        for other_lineage in lineages:
            if other_lineage == target_lineage:
                continue
            other_populations = [value for value in populations if lineage[value] == other_lineage]
            other_lineage_values.append(
                population_contributions[other_populations].median(axis=1).rename(other_lineage)
            )
        anchor_competitor = pd.concat(other_lineage_values, axis=1).max(axis=1)
        same_lineage = [
            value for value in populations if lineage[value] == target_lineage and value != target
        ]
        require(bool(same_lineage), f"No local competitors for {target}")
        local_competitor = population_contributions[same_lineage].max(axis=1)

        def rank_module(competitor: pd.Series, size: int, excluded: set[str]) -> pd.DataFrame:
            ranking = pd.DataFrame({
                "gene_symbol": counts.index,
                "target_consistency": target_value.to_numpy(),
                "competitor_envelope": competitor.to_numpy(),
            })
            ranking["specificity_margin"] = ranking["target_consistency"] - ranking["competitor_envelope"]
            ranking = ranking[
                (ranking["target_consistency"] > 0)
                & (ranking["specificity_margin"] > 0)
                & (~ranking["gene_symbol"].isin(excluded))
            ].sort_values(
                ["specificity_margin", "target_consistency", "gene_symbol"],
                ascending=[False, False, True], kind="mergesort"
            ).head(size).reset_index(drop=True)
            ranking["module_rank"] = np.arange(1, len(ranking) + 1)
            require(len(ranking) == size, f"Undersized module for {target}")
            return ranking

        anchor = rank_module(anchor_competitor, anchor_size, set())
        local = rank_module(local_competitor, local_size, set(anchor["gene_symbol"]))
        return {"anchor": anchor, "local": local}

    full_modules = {target: derive_modules(target) for target in populations}
    mapping_rows: list[dict[str, object]] = []
    membership_rows: list[dict[str, object]] = []
    set_genes: dict[str, list[str]] = {}
    for population_index, target in enumerate(populations, 1):
        for module in ("anchor", "local"):
            set_id = f"IMMGEN_NARROW_{module.upper()}_{population_index:02d}"
            ranking = full_modules[target][module]
            set_genes[set_id] = list(ranking["gene_symbol"])
            mapping_rows.append({
                "geneset_id": set_id,
                "target": target,
                "lineage": lineage[target],
                "module": module,
                "genes": len(ranking),
            })
            for row in ranking.itertuples(index=False):
                membership_rows.append({
                    "geneset_id": set_id,
                    "target": target,
                    "lineage": lineage[target],
                    "module": module,
                    "gene_symbol": row.gene_symbol,
                    "module_rank": int(row.module_rank),
                    "target_consistency": float(row.target_consistency),
                    "competitor_envelope": float(row.competitor_envelope),
                    "specificity_margin": float(row.specificity_margin),
                })
    mapping = pd.DataFrame(mapping_rows)
    membership = pd.DataFrame(membership_rows)
    mapping.to_csv(tables / "01_module_map.tsv", sep="\t", index=False)
    with gzip.open(tables / "02_module_membership.tsv.gz", "wt", encoding="utf-8") as handle:
        membership.to_csv(handle, sep="\t", index=False)
    with (tables / "03_score_collection.gmt").open("w", encoding="utf-8") as handle:
        for row in mapping.itertuples(index=False):
            handle.write("\t".join([
                row.geneset_id,
                f"ImmGen-only {row.module} module for {row.target}",
                *set_genes[row.geneset_id],
            ]) + "\n")

    raw_anchor = pd.DataFrame({
        target: contributions.loc[full_modules[target]["anchor"]["gene_symbol"]].mean(axis=0)
        for target in populations
    })
    raw_local = pd.DataFrame({
        target: contributions.loc[full_modules[target]["local"]["gene_symbol"]].mean(axis=0)
        for target in populations
    })
    candidate_rows: list[dict[str, object]] = []
    holdout_index = 0
    for plan in plans.itertuples(index=False):
        heldout_libraries = [value for value in plan.heldout_units.split(";") if value]
        for heldout_library in heldout_libraries:
            holdout_index += 1
            target = plan.target
            fold_modules = derive_modules(target, heldout_library)
            fold_anchor = contributions.loc[fold_modules["anchor"]["gene_symbol"]].mean(axis=0)
            fold_local = contributions.loc[fold_modules["local"]["gene_symbol"]].mean(axis=0)
            training_ids = [value for value in metadata["library_id"] if value != heldout_library]
            anchor_scores: dict[str, float] = {}
            local_scores: dict[str, float] = {}
            for candidate in populations:
                candidate_anchor = fold_anchor if candidate == target else raw_anchor[candidate]
                candidate_local = fold_local if candidate == target else raw_local[candidate]
                anchor_scores[candidate] = empirical_percentile(
                    candidate_anchor.loc[training_ids].to_numpy(), float(candidate_anchor.loc[heldout_library])
                )
                local_scores[candidate] = empirical_percentile(
                    candidate_local.loc[training_ids].to_numpy(), float(candidate_local.loc[heldout_library])
                )
            for candidate in populations:
                conflict = max(
                    anchor_scores[other] for other in populations
                    if lineage[other] != lineage[candidate]
                )
                candidate_rows.append({
                    "holdout_index": holdout_index,
                    "heldout_library": heldout_library,
                    "true_population": target,
                    "true_lineage": lineage[target],
                    "candidate_population": candidate,
                    "candidate_lineage": lineage[candidate],
                    "anchor_score": anchor_scores[candidate],
                    "local_score": local_scores[candidate],
                    "conflict_score": conflict,
                    "anchor_contrast": anchor_scores[candidate] - conflict,
                })
    candidates = pd.DataFrame(candidate_rows)
    require(candidates["holdout_index"].nunique() == int(config["internal_evaluation"]["heldout_libraries_expected"]), "Holdout inventory changed")
    with gzip.open(tables / "04_immgen_holdout_candidate_scores.tsv.gz", "wt", encoding="utf-8") as handle:
        candidates.to_csv(handle, sep="\t", index=False)

    ablations = [
        ("anchor1_local1_conflict0", 1.0, 1.0, 0.0),
        ("anchor1_local1_conflict0.5", 1.0, 1.0, 0.5),
        ("anchor1_local1_conflict1", 1.0, 1.0, 1.0),
        ("anchor2_local1_conflict1", 2.0, 1.0, 1.0),
    ]
    ablation_rows = []
    decisions_by_variant: dict[str, pd.DataFrame] = {}
    for variant, anchor_weight, local_weight, conflict_weight in ablations:
        scored = candidates.copy()
        scored["composite_score"] = (
            anchor_weight * scored["anchor_score"]
            + local_weight * scored["local_score"]
            - conflict_weight * scored["conflict_score"]
        )
        decision_rows = []
        for _, rows in scored.groupby("holdout_index", sort=True):
            rows = rows.sort_values(["composite_score", "candidate_population"], ascending=[False, True], kind="mergesort")
            winner, runner = rows.iloc[0], rows.iloc[1]
            decision_rows.append({
                "holdout_index": int(winner["holdout_index"]),
                "heldout_library": winner["heldout_library"],
                "true_population": winner["true_population"],
                "true_lineage": winner["true_lineage"],
                "winner_population": winner["candidate_population"],
                "winner_lineage": winner["candidate_lineage"],
                "winner_exact": winner["candidate_population"] == winner["true_population"],
                "winner_lineage_correct": winner["candidate_lineage"] == winner["true_lineage"],
                "anchor_score": float(winner["anchor_score"]),
                "local_score": float(winner["local_score"]),
                "conflict_score": float(winner["conflict_score"]),
                "anchor_contrast": float(winner["anchor_contrast"]),
                "top_composite_score": float(winner["composite_score"]),
                "runner_up_composite_score": float(runner["composite_score"]),
                "composite_margin": float(winner["composite_score"] - runner["composite_score"]),
            })
        decisions = pd.DataFrame(decision_rows)
        decisions_by_variant[variant] = decisions
        ablation_rows.append({
            "variant": variant,
            "anchor_weight": anchor_weight,
            "local_weight": local_weight,
            "conflict_weight": conflict_weight,
            "heldout_libraries": len(decisions),
            "exact_correct": int(decisions["winner_exact"].sum()),
            "exact_accuracy": float(decisions["winner_exact"].mean()),
            "lineage_correct": int(decisions["winner_lineage_correct"].sum()),
            "lineage_accuracy": float(decisions["winner_lineage_correct"].mean()),
        })
    ablation = pd.DataFrame(ablation_rows)
    ablation.to_csv(tables / "05_immgen_weight_ablation.tsv", sep="\t", index=False)

    selected_variant = "anchor2_local1_conflict1"
    decisions = decisions_by_variant[selected_variant]
    minimum_anchor = float(config["acceptance"]["minimum_anchor_score"])
    minimum_contrast = float(config["acceptance"]["minimum_anchor_contrast"])
    minimum_margin = float(config["acceptance"]["minimum_composite_margin"])
    decisions["accepted"] = (
        (decisions["anchor_score"] >= minimum_anchor)
        & (decisions["anchor_contrast"] >= minimum_contrast)
        & (decisions["composite_margin"] >= minimum_margin)
    )
    decisions.to_csv(tables / "06_immgen_selected_holdout_decisions.tsv", sep="\t", index=False)
    accepted = decisions[decisions["accepted"]]
    lineage_summary = decisions.groupby("true_lineage", sort=False).agg(
        heldout_libraries=("heldout_library", "size"),
        unthresholded_exact=("winner_exact", "sum"),
        unthresholded_lineage_correct=("winner_lineage_correct", "sum"),
        accepted=("accepted", "sum"),
    ).reset_index()
    accepted_exact = accepted.groupby("true_lineage")["winner_exact"].sum()
    accepted_lineage = accepted.groupby("true_lineage")["winner_lineage_correct"].sum()
    lineage_summary["accepted_exact"] = lineage_summary["true_lineage"].map(accepted_exact).fillna(0).astype(int)
    lineage_summary["accepted_lineage_correct"] = lineage_summary["true_lineage"].map(accepted_lineage).fillna(0).astype(int)
    lineage_summary.to_csv(tables / "07_immgen_lineage_summary.tsv", sep="\t", index=False)

    summary = {
        "schema_version": "bdev.lineage_anchored_narrow_internal_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "data_source": "ImmGen_only",
        "populations_competing": len(populations),
        "broad_gate_used": False,
        "heldout_libraries": len(decisions),
        "unthresholded_exact_correct": int(decisions["winner_exact"].sum()),
        "unthresholded_exact_accuracy": float(decisions["winner_exact"].mean()),
        "unthresholded_lineage_correct": int(decisions["winner_lineage_correct"].sum()),
        "unthresholded_lineage_accuracy": float(decisions["winner_lineage_correct"].mean()),
        "accepted": len(accepted),
        "accepted_exact_correct": int(accepted["winner_exact"].sum()),
        "accepted_exact_precision": float(accepted["winner_exact"].mean()),
        "accepted_lineage_correct": int(accepted["winner_lineage_correct"].sum()),
        "accepted_lineage_precision": float(accepted["winner_lineage_correct"].mean()),
        "weight_selection_is_independent_validation": False,
        "competitor_program_discovery_is_fully_nested": False,
        "Bdev_used_for_selection": False,
        "accepted_labels_emitted": False,
    }
    (reports / "01_immgen_internal_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    input_manifest = {
        "schema_version": "bdev.lineage_anchored_narrow_input_manifest.v1",
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
    (manifests / "01_input_manifest.json").write_text(json.dumps(input_manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
