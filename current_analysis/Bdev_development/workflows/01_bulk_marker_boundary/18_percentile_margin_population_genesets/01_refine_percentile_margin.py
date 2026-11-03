#!/usr/bin/env python3
"""Refine ImmGen programs against the unchanged flat percentile-margin reducer."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str], compressed: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if compressed:
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
                    writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(rows)
    else:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


def slug(value: str) -> str:
    encoded = value.replace("+", "_POS_").replace("-", "_NEG_").replace("&", "_AND_")
    return re.sub(r"(^_+|_+$)", "", re.sub(r"[^A-Za-z0-9]+", "_", encoded)).upper()


def fmt(value: float) -> str:
    return f"{value:.12g}" if np.isfinite(value) else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    require(
        config["schema_version"] == "bdev.immgen_percentile_margin_population_genesets.v1",
        "Unexpected config schema",
    )

    inputs = {name: Path(spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        require(sha256(path) == config["inputs"][name]["sha256"], f"Checksum mismatch for {name}")

    run_identity = config["run_identity"]
    intermediate_root = Path(config["outputs"]["intermediate_root"]) / run_identity
    results_root = Path(config["outputs"]["results_root"]) / run_identity
    trace_path = intermediate_root / config["outputs"]["trace"]
    result_paths = {name: results_root / rel for name, rel in config["outputs"]["tables"].items()}
    report_paths = {name: results_root / rel for name, rel in config["outputs"]["reports"].items()}
    manifest_path = results_root / config["outputs"]["manifest"]
    require(not results_root.exists(), f"Refusing to overwrite existing run: {results_root}")
    require(not intermediate_root.exists(), f"Refusing to overwrite existing run: {intermediate_root}")

    counts = pd.read_csv(inputs["selected_counts"], sep="\t", index_col=0)
    metadata = pd.read_csv(inputs["metadata"], sep="\t", dtype=str).fillna("")
    require(list(counts.columns) == list(metadata["library_id"]), "Count columns and metadata are misaligned")
    all_ids = list(metadata["library_id"])
    population_order = list(dict.fromkeys(metadata["population"]))
    require(len(population_order) == 64, f"Expected 64 populations, found {len(population_order)}")
    max_rank = int(config["method"]["scoring"]["max_rank"])
    capped_ranks = counts.rank(axis=0, method="average", ascending=False).clip(upper=max_rank)

    def score_genes(genes: list[str]) -> pd.Series:
        require(bool(genes), "Cannot score an empty program")
        n_genes = len(genes)
        baseline = n_genes * (n_genes + 1) / 2
        denominator = n_genes * max_rank - baseline
        require(denominator > 0, f"Invalid UCell denominator for {n_genes} genes")
        rank_sum = capped_ranks.loc[genes, all_ids].sum(axis=0)
        return (1 - (rank_sum - baseline) / denominator).clip(lower=0, upper=1)

    def calibrate(scores: pd.Series) -> pd.Series:
        return scores.rank(method="average", pct=True)

    plan_rows = read_tsv(inputs["joint_fold_plan"])
    ranking_rows = read_tsv(inputs["joint_fold_rankings"])
    plan_by_target: dict[str, list[dict[str, str]]] = defaultdict(list)
    rankings_by_target_fold: dict[tuple[str, int], list[tuple[int, str]]] = defaultdict(list)
    for row in plan_rows:
        plan_by_target[row["target"]].append(row)
    for row in ranking_rows:
        rankings_by_target_fold[(row["target"], int(row["fold"]))].append(
            (int(row["fold_rank"]), row["gene_symbol"])
        )
    require(set(plan_by_target) == set(population_order), "Fold plan does not cover all populations")

    candidate_sets: dict[str, list[str]] = {}
    candidate_sources: dict[str, str] = {}
    candidate_ranks: dict[str, dict[str, int]] = {}
    fold_presence: dict[str, dict[str, int]] = {}
    for target in population_order:
        folds = sorted({int(row["fold"]) for row in plan_by_target[target]})
        fold_lists = [
            [gene for _, gene in sorted(rankings_by_target_fold[(target, fold)])]
            for fold in folds
        ]
        if not all(fold_lists):
            candidate, source = [], "one_or_more_folds_without_joint_candidates"
        elif len(fold_lists) == 1:
            candidate, source = fold_lists[0], "single_full_data_joint_top100"
        else:
            common = set(fold_lists[0]).intersection(*fold_lists[1:])
            if common:
                candidate = sorted(
                    common,
                    key=lambda gene: (
                        np.mean([genes.index(gene) + 1 for genes in fold_lists]),
                        max(genes.index(gene) + 1 for genes in fold_lists),
                        gene,
                    ),
                )
                source = "joint_crossfold_intersection"
            else:
                union = set().union(*fold_lists)
                candidate = sorted(
                    union,
                    key=lambda gene: (
                        -sum(gene in genes for genes in fold_lists),
                        np.mean([genes.index(gene) + 1 for genes in fold_lists if gene in genes]),
                        gene,
                    ),
                )[:100]
                source = "joint_empty_intersection_ranked_union_fallback"
        candidate_sets[target] = candidate
        candidate_sources[target] = source
        candidate_ranks[target] = {gene: rank for rank, gene in enumerate(candidate, 1)}
        fold_presence[target] = {gene: sum(gene in genes for genes in fold_lists) for gene in candidate}

    require(all(candidate_sets.values()), "At least one population has no candidate program")
    final_sets = {target: list(candidate_sets[target]) for target in population_order}
    calibrated = pd.DataFrame(
        {target: calibrate(score_genes(final_sets[target])) for target in population_order},
        index=all_ids,
    )
    metadata_indexed = metadata.set_index("library_id")
    target_ids = {
        target: list(metadata.loc[metadata["population"] == target, "library_id"])
        for target in population_order
    }
    floor = int(config["method"]["refinement"]["standard_floor"])
    trace_rows: list[dict] = []
    refinement_rows: list[dict] = []

    def margin_metrics(target: str, target_scores: pd.Series, competitor_scores: pd.DataFrame, ids: list[str]) -> dict:
        strongest = competitor_scores.loc[ids].max(axis=1)
        margins = target_scores.loc[ids] - strongest
        winners = competitor_scores.loc[ids].idxmax(axis=1)
        return {
            "mean_margin": float(margins.mean()),
            "minimum_margin": float(margins.min()),
            "mean_target_score": float(target_scores.loc[ids].mean()),
            "strongest_competitor": str(winners.value_counts().sort_index().idxmax()),
        }

    def greedy_refine(
        target: str,
        genes: list[str],
        training_ids: list[str],
        competitor_scores: pd.DataFrame,
        context: str,
        coordinate_pass: int,
    ) -> tuple[list[str], pd.Series, dict, int]:
        current = list(genes)
        current_scores = calibrate(score_genes(current))
        rank_sum = capped_ranks.loc[current, all_ids].sum(axis=0).to_numpy(dtype=float)
        current_metrics = margin_metrics(target, current_scores, competitor_scores, training_ids)
        initial = dict(current_metrics)
        initial_genes = len(current)
        accepted_deletions = 0
        iteration = 0
        local_floor = min(floor, len(current))
        id_positions = [all_ids.index(library_id) for library_id in training_ids]
        while len(current) > local_floor:
            genes_after = len(current) - 1
            baseline = genes_after * (genes_after + 1) / 2
            denominator = genes_after * max_rank - baseline
            deletion_sums = rank_sum[None, :] - capped_ranks.loc[current, all_ids].to_numpy(dtype=float)
            deletion_raw = np.clip(1 - (deletion_sums - baseline) / denominator, 0, 1)
            deletion_frame = pd.DataFrame(deletion_raw).rank(axis=1, method="average", pct=True).to_numpy()
            competitor_max = competitor_scores.loc[training_ids].max(axis=1).to_numpy()
            margins = deletion_frame[:, id_positions] - competitor_max[None, :]
            mean_margins = margins.mean(axis=1)
            minimum_margins = margins.min(axis=1)
            mean_scores = deletion_frame[:, id_positions].mean(axis=1)
            best = max(
                range(len(current)),
                key=lambda index: (
                    mean_margins[index],
                    minimum_margins[index],
                    mean_scores[index],
                    candidate_ranks.get(target, {}).get(current[index], 0),
                    current[index],
                ),
            )
            iteration += 1
            improves = mean_margins[best] > current_metrics["mean_margin"] + 1e-12
            trace_rows.append({
                "context": context,
                "target": target,
                "coordinate_pass": coordinate_pass,
                "iteration": iteration,
                "genes_before": len(current),
                "tested_deletion": current[best],
                "mean_margin_before": fmt(current_metrics["mean_margin"]),
                "mean_margin_after": fmt(float(mean_margins[best])),
                "minimum_margin_after": fmt(float(minimum_margins[best])),
                "mean_target_score_after": fmt(float(mean_scores[best])),
                "accepted_deletion": str(improves).lower(),
            })
            if not improves:
                break
            selected_gene = current[best]
            current_scores = pd.Series(deletion_frame[best], index=all_ids)
            current_metrics = margin_metrics(target, current_scores, competitor_scores, training_ids)
            rank_sum = rank_sum - capped_ranks.loc[selected_gene, all_ids].to_numpy(dtype=float)
            current.remove(selected_gene)
            accepted_deletions += 1
        summary = {
            "context": context,
            "target": target,
            "coordinate_pass": coordinate_pass,
            "training_libraries": len(training_ids),
            "initial_genes": initial_genes,
            "final_genes": len(current),
            "accepted_deletions": accepted_deletions,
            "initial_mean_margin": fmt(initial["mean_margin"]),
            "final_mean_margin": fmt(current_metrics["mean_margin"]),
            "final_minimum_margin": fmt(current_metrics["minimum_margin"]),
            "final_mean_target_score": fmt(current_metrics["mean_target_score"]),
            "strongest_competitor": current_metrics["strongest_competitor"],
        }
        return current, current_scores, summary, accepted_deletions

    # Repeated coordinate passes make every program compete on its own target
    # libraries against the currently refined other 63 programs.
    coordinate_pass = 0
    while True:
        coordinate_pass += 1
        pass_deletions = 0
        for target_i, target in enumerate(population_order, 1):
            competitors = calibrated.drop(columns=target)
            refined, scores, summary, deletions = greedy_refine(
                target,
                final_sets[target],
                target_ids[target],
                competitors,
                "complete_atlas",
                coordinate_pass,
            )
            final_sets[target] = refined
            calibrated[target] = scores
            refinement_rows.append(summary)
            pass_deletions += deletions
            print(
                f"coordinate pass {coordinate_pass}, target {target_i}/64 {target}: "
                f"{summary['initial_genes']} -> {summary['final_genes']}",
                flush=True,
            )
        if pass_deletions == 0:
            break
        require(coordinate_pass < 100, "Coordinate refinement did not converge")

    minimum_score = float(config["flat_default_qualification"]["minimum_calibrated_score"])
    minimum_margin = float(config["flat_default_qualification"]["minimum_calibrated_margin"])
    tie_tolerance = float(config["flat_default_qualification"]["tie_tolerance"])

    def reduce_values(values: pd.Series) -> dict:
        ordered = values.sort_values(ascending=False, kind="mergesort")
        winner = str(ordered.index[0])
        top_score = float(ordered.iloc[0])
        runner_score = float(ordered.iloc[1])
        margin = top_score - runner_score
        if top_score < minimum_score:
            reason = "min_score"
        elif margin <= tie_tolerance:
            reason = "tie"
        elif margin < minimum_margin:
            reason = "min_margin"
        else:
            reason = "none"
        return {
            "winner": winner,
            "label": winner if reason == "none" else "unassigned",
            "top_score": top_score,
            "runner_up_score": runner_score,
            "score_margin": margin,
            "abstain_reason": reason,
        }

    decision_rows: list[dict] = []
    for library_id in all_ids:
        result = reduce_values(calibrated.loc[library_id])
        expected = metadata_indexed.loc[library_id, "population"]
        expected_score = float(calibrated.loc[library_id, expected])
        decision_rows.append({
            "library_id": library_id,
            "expected_population": expected,
            "winner_population": result["winner"],
            "label": result["label"],
            "top_score": fmt(result["top_score"]),
            "runner_up_score": fmt(result["runner_up_score"]),
            "score_margin": fmt(result["score_margin"]),
            "abstain_reason": result["abstain_reason"],
            "target_score": fmt(expected_score),
            "target_rank": 1 + int((calibrated.loc[library_id] > expected_score).sum()),
            "winner_correct": str(result["winner"] == expected).lower(),
            "assigned_correct": str(result["label"] == expected).lower(),
        })

    # Membership holdout: only the target membership is rebuilt and refined;
    # the complete-data refined competitors remain fixed.
    holdout_rows: list[dict] = []
    for target in population_order:
        for plan in sorted(plan_by_target[target], key=lambda row: int(row["fold"])):
            heldout = [value for value in plan["heldout_units"].split(";") if value]
            if not heldout:
                continue
            fold = int(plan["fold"])
            fold_genes = [gene for _, gene in sorted(rankings_by_target_fold[(target, fold)])]
            training = [library_id for library_id in target_ids[target] if library_id not in heldout]
            if not fold_genes or not training:
                for library_id in heldout:
                    holdout_rows.append({
                        "target": target,
                        "fold": fold,
                        "heldout_library": library_id,
                        "evidence_tier": plan["evidence_tier"],
                        "fold_genes": len(fold_genes),
                        "winner_population": "",
                        "label": "unassigned",
                        "top_score": "",
                        "runner_up_score": "",
                        "score_margin": "",
                        "target_score": "",
                        "target_rank": 0,
                        "abstain_reason": "no_fold_training_program",
                        "winner_correct": "false",
                        "assigned_correct": "false",
                    })
                continue
            competitors = calibrated.drop(columns=target)
            fold_refined, fold_scores, summary, _ = greedy_refine(
                target,
                fold_genes,
                training,
                competitors,
                f"target_holdout_fold_{fold}",
                0,
            )
            refinement_rows.append(summary)
            for library_id in heldout:
                values = calibrated.loc[library_id].copy()
                values.loc[target] = float(fold_scores.loc[library_id])
                result = reduce_values(values)
                target_score = float(values.loc[target])
                holdout_rows.append({
                    "target": target,
                    "fold": fold,
                    "heldout_library": library_id,
                    "evidence_tier": plan["evidence_tier"],
                    "fold_genes": len(fold_refined),
                    "winner_population": result["winner"],
                    "label": result["label"],
                    "top_score": fmt(result["top_score"]),
                    "runner_up_score": fmt(result["runner_up_score"]),
                    "score_margin": fmt(result["score_margin"]),
                    "target_score": fmt(target_score),
                    "target_rank": 1 + int((values > target_score).sum()),
                    "abstain_reason": result["abstain_reason"],
                    "winner_correct": str(result["winner"] == target).lower(),
                    "assigned_correct": str(result["label"] == target).lower(),
                })

    membership_rows: list[dict] = []
    summary_rows: list[dict] = []
    qualification_rows: list[dict] = []
    eligible_targets: list[str] = []
    minimum_holdout = float(
        config["flat_default_qualification"]["minimum_target_membership_holdout_assigned_correct_fraction"]
    )
    minimum_genes = int(config["flat_default_qualification"]["minimum_genes"])
    for target in population_order:
        plan = plan_by_target[target][0]
        genes = final_sets[target]
        target_program_scores = calibrated.loc[target_ids[target], target]
        population_means = {
            population: float(calibrated.loc[target_ids[population], target].mean())
            for population in population_order
        }
        complete_mean_rank = 1 + sum(
            score >= population_means[target]
            for population, score in population_means.items()
            if population != target
        )
        validation = [row for row in holdout_rows if row["target"] == target]
        holdout_winner = (
            sum(row["winner_correct"] == "true" for row in validation) / len(validation)
            if validation else float("nan")
        )
        holdout_assigned = (
            sum(row["assigned_correct"] == "true" for row in validation) / len(validation)
            if validation else float("nan")
        )
        reasons = []
        if int(plan["independent_units"]) < 2:
            reasons.append("singleton")
        if len(genes) < minimum_genes:
            reasons.append("fewer_than_five_genes")
        if complete_mean_rank != 1:
            reasons.append("complete_mean_rank_not_first")
        if not validation:
            reasons.append("no_target_membership_holdout")
        elif holdout_assigned < minimum_holdout:
            reasons.append("holdout_assigned_correct_fraction")
        eligible = not reasons
        if eligible:
            eligible_targets.append(target)
        status = "flat_label_candidate" if eligible else "score_only_not_flat_qualified"
        competitor_means = {population: score for population, score in population_means.items() if population != target}
        strongest_population = max(competitor_means, key=lambda population: (competitor_means[population], population))
        summary_rows.append({
            "geneset_id": "IMMGEN_POP_" + slug(target),
            "target": target,
            "source_role": plan["source_role"],
            "lineage": plan["lineage"],
            "independent_units": plan["independent_units"],
            "evidence_tier": plan["evidence_tier"],
            "candidate_source": candidate_sources[target],
            "candidate_genes": len(candidate_sets[target]),
            "final_genes": len(genes),
            "status": status,
            "complete_mean_target_rank": complete_mean_rank,
            "strongest_other_population": strongest_population,
            "target_population_mean_percentile": fmt(population_means[target]),
            "strongest_other_population_mean_percentile": fmt(competitor_means[strongest_population]),
            "target_library_mean_classifier_margin": fmt(
                float((target_program_scores - calibrated.drop(columns=target).loc[target_ids[target]].max(axis=1)).mean())
            ),
            "holdout_assigned_correct_fraction": fmt(holdout_assigned),
            "flat_label_eligible": str(eligible).lower(),
        })
        qualification_rows.append({
            "geneset_id": "IMMGEN_POP_" + slug(target),
            "target": target,
            "lineage": plan["lineage"],
            "independent_units": plan["independent_units"],
            "evidence_tier": plan["evidence_tier"],
            "final_genes": len(genes),
            "complete_mean_target_rank": complete_mean_rank,
            "holdout_libraries": len(validation),
            "holdout_winner_correct_fraction": fmt(holdout_winner),
            "holdout_assigned_correct_fraction": fmt(holdout_assigned),
            "flat_label_eligible": str(eligible).lower(),
            "ineligibility_reasons": ";".join(reasons),
        })
        for final_rank, gene in enumerate(genes, 1):
            membership_rows.append({
                "geneset_id": "IMMGEN_POP_" + slug(target),
                "target": target,
                "source_role": plan["source_role"],
                "lineage": plan["lineage"],
                "independent_units": plan["independent_units"],
                "evidence_tier": plan["evidence_tier"],
                "candidate_source": candidate_sources[target],
                "gene_symbol": gene,
                "candidate_rank": candidate_ranks[target][gene],
                "fold_presence": fold_presence[target][gene],
                "final_rank": final_rank,
                "candidate_genes": len(candidate_sets[target]),
                "final_genes": len(genes),
            })

    membership_fields = [
        "geneset_id", "target", "source_role", "lineage", "independent_units", "evidence_tier",
        "candidate_source", "gene_symbol", "candidate_rank", "fold_presence", "final_rank",
        "candidate_genes", "final_genes",
    ]
    write_tsv(result_paths["membership"], membership_rows, membership_fields, compressed=True)
    write_tsv(result_paths["target_summary"], summary_rows, list(summary_rows[0]))
    write_tsv(result_paths["full_decisions"], decision_rows, list(decision_rows[0]))
    write_tsv(result_paths["holdout_validation"], holdout_rows, list(holdout_rows[0]))
    write_tsv(result_paths["qualification"], qualification_rows, list(qualification_rows[0]))
    write_tsv(result_paths["refinement_summary"], refinement_rows, list(refinement_rows[0]))
    trace_fields = [
        "context", "target", "coordinate_pass", "iteration", "genes_before", "tested_deletion",
        "mean_margin_before", "mean_margin_after", "minimum_margin_after", "mean_target_score_after",
        "accepted_deletion",
    ]
    write_tsv(trace_path, trace_rows, trace_fields, compressed=True)
    for key, targets in (("score_gmt", population_order), ("label_gmt", eligible_targets)):
        result_paths[key].parent.mkdir(parents=True, exist_ok=True)
        with result_paths[key].open("w") as handle:
            for target in targets:
                status = next(row["status"] for row in summary_rows if row["target"] == target)
                description = f"diagnostic_only|status={status}"
                handle.write("\t".join(["IMMGEN_POP_" + slug(target), description, *final_sets[target]]) + "\n")

    complete_winner = sum(row["winner_correct"] == "true" for row in decision_rows) / len(decision_rows)
    complete_assigned = sum(row["assigned_correct"] == "true" for row in decision_rows) / len(decision_rows)
    holdout_winner = sum(row["winner_correct"] == "true" for row in holdout_rows) / len(holdout_rows)
    holdout_assigned = sum(row["assigned_correct"] == "true" for row in holdout_rows) / len(holdout_rows)
    abstentions = pd.Series([row["abstain_reason"] for row in decision_rows]).value_counts().to_dict()
    run_summary = {
        "schema_version": "bdev.immgen_percentile_margin_population_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": run_identity,
        "coordinate_passes": coordinate_pass,
        "score_programs": len(population_order),
        "flat_label_candidate_programs": len(eligible_targets),
        "complete_atlas_libraries": len(decision_rows),
        "complete_atlas_winner_accuracy": complete_winner,
        "complete_atlas_assigned_correct_fraction": complete_assigned,
        "complete_atlas_abstention_reasons": abstentions,
        "target_membership_holdout_libraries": len(holdout_rows),
        "target_membership_holdout_winner_accuracy": holdout_winner,
        "target_membership_holdout_assigned_correct_fraction": holdout_assigned,
        "single_cell_expression_read": False,
        "scworkbench_invoked": False,
        "feature_sets_registered": False,
        "artifacts_promoted": False,
        "accepted_sets_emitted": 0,
    }
    report_paths["summary"].parent.mkdir(parents=True, exist_ok=True)
    report_paths["summary"].write_text(json.dumps(run_summary, indent=2, sort_keys=True) + "\n")
    report_paths["report"].write_text(
        "# Percentile-margin ImmGen population refinement\n\n"
        f"- Score programs: {len(population_order)}/64\n"
        f"- Flat-label candidate programs: {len(eligible_targets)}/64\n"
        f"- Complete-atlas winner accuracy: {100 * complete_winner:.1f}% (optimistic)\n"
        f"- Complete-atlas assigned-correct fraction: {100 * complete_assigned:.1f}%\n"
        f"- Target-membership holdout winner accuracy: {100 * holdout_winner:.1f}%\n"
        f"- Target-membership holdout assigned-correct fraction: {100 * holdout_assigned:.1f}%\n\n"
        "The default score and margin gates were not changed. All outputs remain diagnostic and unpromoted.\n"
    )

    generated = {"trace": trace_path, **result_paths, **report_paths}
    script_path = Path(__file__).resolve()
    manifest = {
        "schema_version": "bdev.immgen_percentile_margin_population_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": run_identity,
        "inputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in {"config": config_path, **inputs}.items()
        },
        "scripts": {script_path.name: {"path": str(script_path), "sha256": sha256(script_path)}},
        "outputs": {
            name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}
            for name, path in generated.items()
        },
        "single_cell_expression_read": False,
        "scworkbench_invoked": False,
        "feature_sets_registered": False,
        "artifacts_promoted": False,
        "accepted_sets_emitted": 0,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(run_summary, indent=2, sort_keys=True))
    print(f"manifest_sha256={sha256(manifest_path)}")


if __name__ == "__main__":
    main()
