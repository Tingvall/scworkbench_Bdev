#!/usr/bin/env python3
"""Refine joint-specific ImmGen programs and qualify flat-label candidates."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
from collections import Counter, defaultdict
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


def resolve(base: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (base / path).resolve()


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
                    writer.writeheader(); writer.writerows(rows)
    else:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)


def slug(value: str) -> str:
    encoded = value.replace("+", "_POS_").replace("-", "_NEG_").replace("&", "_AND_")
    return re.sub(r"(^_+|_+$)", "", re.sub(r"[^A-Za-z0-9]+", "_", encoded)).upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    require(config["schema_version"] == "bdev.immgen_global_local_population_genesets.v1", "Unexpected config schema")
    inputs = {name: resolve(config_path.parent, spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        require(sha256(path) == config["inputs"][name]["sha256"], f"Checksum mismatch for {name}")

    intermediate_root = Path(config["outputs"]["intermediate_root"])
    results_root = Path(config["outputs"]["results_root"])
    intermediate = {name: intermediate_root / rel for name, rel in config["outputs"]["intermediate_tables"].items()}
    results = {name: results_root / rel for name, rel in config["outputs"]["result_tables"].items()}
    reports = {name: results_root / rel for name, rel in config["outputs"]["reports"].items()}
    for name in ("fold_plan", "fold_rankings", "model_audit"):
        require(intermediate[name].is_file(), f"Missing model output {name}: {intermediate[name]}")

    counts = pd.read_csv(inputs["selected_counts"], sep="\t", index_col=0)
    metadata = pd.read_csv(inputs["metadata"], sep="\t", dtype=str).fillna("")
    require(list(counts.columns) == list(metadata["library_id"]), "Count columns and metadata are misaligned")
    population_order = list(dict.fromkeys(metadata["population"]))
    require(len(population_order) == int(config["target_scope"]["expected_populations"]), "Unexpected population count")
    population_info = metadata.drop_duplicates("population").set_index("population")
    discovery_role = config["target_scope"]["discovery_role"]
    heldout_role = config["target_scope"]["held_out_role"]
    max_rank = int(config["method"]["scoring"]["max_rank"])
    capped_ranks = counts.rank(axis=0, method="average", ascending=False).clip(upper=max_rank)

    def evaluation_metadata(target: str) -> pd.DataFrame:
        selected = metadata[metadata["role"] == discovery_role].copy()
        if population_info.loc[target, "role"] == heldout_role:
            selected = pd.concat([selected, metadata[metadata["population"] == target].copy()], ignore_index=True)
        return selected

    def score_from_rank_sum(rank_sum: pd.Series, n_genes: int) -> pd.Series:
        baseline = n_genes * (n_genes + 1) / 2
        denominator = n_genes * max_rank - baseline
        require(denominator > 0, f"Invalid score denominator for {n_genes} genes")
        return (1 - (rank_sum - baseline) / denominator).clip(lower=0, upper=1)

    def score_genes(genes: list[str], ids: list[str]) -> pd.Series:
        require(bool(genes), "Cannot score an empty program")
        return score_from_rank_sum(capped_ranks.loc[genes, ids].sum(axis=0), len(genes))

    def population_metrics(target: str, library_scores: pd.Series, md: pd.DataFrame) -> dict:
        indexed = md.set_index("library_id").loc[list(library_scores.index)]
        populations = list(dict.fromkeys(indexed["population"]))
        means = {
            population: float(library_scores.loc[indexed.index[indexed["population"] == population]].mean())
            for population in populations
        }
        require(target in means, f"Target absent from score inventory: {target}")
        competitors = {population: value for population, value in means.items() if population != target}
        strongest = max(competitors, key=lambda population: (competitors[population], population))
        sd = float(np.std(list(means.values()), ddof=1))
        gap = (means[target] - competitors[strongest]) / sd if sd > 0 else float("-inf")
        target_ids = indexed.index[indexed["population"] == target]
        return {
            "gap": gap,
            "strongest_competitor": strongest,
            "target_rank": 1 + sum(value >= means[target] for value in competitors.values()),
            "target_library_win_fraction": float(np.mean(library_scores.loc[target_ids] > competitors[strongest])),
        }

    plan_rows = read_tsv(intermediate["fold_plan"])
    ranking_rows = read_tsv(intermediate["fold_rankings"])
    plan_by_target: dict[str, list[dict[str, str]]] = defaultdict(list)
    rankings_by_target_fold: dict[tuple[str, int], list[tuple[int, str]]] = defaultdict(list)
    for row in plan_rows:
        plan_by_target[row["target"]].append(row)
    for row in ranking_rows:
        rankings_by_target_fold[(row["target"], int(row["fold"]))].append((int(row["fold_rank"]), row["gene_symbol"]))
    require(set(plan_by_target) == set(population_order), "Fold plan does not cover all populations")

    candidate_sets: dict[str, list[str]] = {}
    candidate_source: dict[str, str] = {}
    candidate_rank: dict[str, dict[str, int]] = {}
    fold_presence: dict[str, dict[str, int]] = {}
    for target in population_order:
        folds = sorted({int(row["fold"]) for row in plan_by_target[target]})
        fold_lists = [[gene for _, gene in sorted(rankings_by_target_fold[(target, fold)])] for fold in folds]
        if not all(fold_lists):
            candidate, source = [], "one_or_more_folds_without_joint_candidates"
        elif len(fold_lists) == 1:
            candidate, source = fold_lists[0], "single_full_data_joint_top100"
        else:
            intersection = set(fold_lists[0]).intersection(*fold_lists[1:])
            if intersection:
                candidate = sorted(intersection, key=lambda gene: (
                    np.mean([genes.index(gene) + 1 for genes in fold_lists]),
                    max(genes.index(gene) + 1 for genes in fold_lists), gene))
                source = "joint_crossfold_intersection"
            else:
                union = set().union(*fold_lists)
                candidate = sorted(union, key=lambda gene: (
                    -sum(gene in genes for genes in fold_lists),
                    np.mean([genes.index(gene) + 1 for genes in fold_lists if gene in genes]), gene
                ))[: int(config["method"]["empty_intersection_fallback"]["maximum_genes"])]
                source = "joint_empty_intersection_ranked_union_fallback"
        candidate_sets[target] = candidate
        candidate_source[target] = source
        candidate_rank[target] = {gene: rank for rank, gene in enumerate(candidate, 1)}
        fold_presence[target] = {gene: sum(gene in genes for genes in fold_lists) for gene in candidate}

    final_sets: dict[str, list[str]] = {}
    trace_rows: list[dict] = []
    standard_floor = int(config["method"]["refinement"]["standard_floor"])
    for target_i, target in enumerate(population_order, 1):
        current = list(candidate_sets[target])
        if not current:
            final_sets[target] = []
            print(f"refined {target_i} of {len(population_order)} {target}: no joint candidates", flush=True)
            continue
        md = evaluation_metadata(target)
        ids = list(md["library_id"])
        rank_sum = capped_ranks.loc[current, ids].sum(axis=0)
        current_metrics = population_metrics(target, score_from_rank_sum(rank_sum, len(current)), md)
        floor = standard_floor if len(current) >= standard_floor else len(current)
        iteration = 0
        while len(current) > floor:
            genes_after = len(current) - 1
            baseline = genes_after * (genes_after + 1) / 2
            denominator = genes_after * max_rank - baseline
            deletion_sums = rank_sum.to_numpy()[None, :] - capped_ranks.loc[current, ids].to_numpy(dtype=float)
            deletion_scores = np.clip(1 - (deletion_sums - baseline) / denominator, 0, 1)
            aligned = md.set_index("library_id").loc[ids]
            populations = list(dict.fromkeys(aligned["population"]))
            population_scores = np.column_stack([
                deletion_scores[:, np.flatnonzero(aligned["population"].to_numpy() == population)].mean(axis=1)
                for population in populations
            ])
            target_column = populations.index(target)
            competitor_columns = [index for index, population in enumerate(populations) if population != target]
            strongest_columns = np.array([
                max(competitor_columns, key=lambda index: (population_scores[row, index], populations[index]))
                for row in range(len(current))
            ])
            row_index = np.arange(len(current))
            sd = population_scores.std(axis=1, ddof=1)
            gaps = np.divide(
                population_scores[:, target_column] - population_scores[row_index, strongest_columns], sd,
                out=np.full(len(current), -np.inf), where=sd > 0)
            best_index = max(range(len(current)), key=lambda index: (gaps[index], -candidate_rank[target][current[index]]))
            improves = np.isfinite(gaps[best_index]) and gaps[best_index] > current_metrics["gap"]
            iteration += 1
            for index, gene in enumerate(current):
                trace_rows.append({
                    "target": target, "iteration": iteration, "genes_before": len(current),
                    "tested_deletion": gene, "candidate_rank": candidate_rank[target][gene],
                    "gap_before": f"{current_metrics['gap']:.12g}", "gap_after": f"{gaps[index]:.12g}",
                    "improvement": f"{gaps[index]-current_metrics['gap']:.12g}",
                    "strongest_competitor_after": populations[strongest_columns[index]],
                    "selected_deletion": str(improves and index == best_index).lower(),
                })
            if not improves:
                break
            selected_gene = current[best_index]
            rank_sum = rank_sum - capped_ranks.loc[selected_gene, ids]
            current_metrics = {"gap": float(gaps[best_index]), "strongest_competitor": populations[strongest_columns[best_index]]}
            current.remove(selected_gene)
        final_sets[target] = current
        print(f"refined {target_i} of {len(population_order)} {target}: {len(candidate_sets[target])} -> {len(current)}", flush=True)

    all_ids = list(metadata["library_id"])
    nonempty_targets = [target for target in population_order if final_sets[target]]
    require(len(nonempty_targets) > 1, "Fewer than two nonempty programs")
    raw_scores = pd.DataFrame({target: score_genes(final_sets[target], all_ids) for target in nonempty_targets}, index=all_ids)
    calibrated = raw_scores.rank(axis=0, method="average", pct=True)
    minimum_score = float(config["flat_default_qualification"]["minimum_calibrated_score"])
    minimum_margin = float(config["flat_default_qualification"]["minimum_calibrated_margin"])
    tie_tolerance = 1e-8

    def reduce_values(values: pd.Series) -> dict:
        ordered = values.sort_values(ascending=False, kind="mergesort")
        top_label = str(ordered.index[0]); top_score = float(ordered.iloc[0])
        runner_score = float(ordered.iloc[1]); margin = top_score - runner_score
        if top_score < minimum_score:
            reason = "min_score"
        elif margin <= tie_tolerance:
            reason = "tie"
        elif margin < minimum_margin:
            reason = "min_margin"
        else:
            reason = "none"
        return {
            "winner": top_label, "label": top_label if reason == "none" else "unassigned",
            "top_score": top_score, "runner_up_score": runner_score,
            "score_margin": margin, "abstain_reason": reason,
            "winner_rank": ordered,
        }

    decision_rows: list[dict] = []
    for library_id in all_ids:
        result = reduce_values(calibrated.loc[library_id])
        expected = metadata.loc[metadata["library_id"] == library_id, "population"].iloc[0]
        target_score = float(calibrated.loc[library_id, expected]) if expected in calibrated.columns else float("nan")
        target_rank = 1 + int((calibrated.loc[library_id] > target_score).sum()) if np.isfinite(target_score) else 0
        decision_rows.append({
            "library_id": library_id, "expected_population": expected,
            "winner_population": result["winner"], "label": result["label"],
            "top_score": f"{result['top_score']:.12g}", "runner_up_score": f"{result['runner_up_score']:.12g}",
            "score_margin": f"{result['score_margin']:.12g}", "abstain_reason": result["abstain_reason"],
            "target_score": f"{target_score:.12g}" if np.isfinite(target_score) else "",
            "target_rank": target_rank,
            "winner_correct": str(result["winner"] == expected).lower(),
            "assigned_correct": str(result["label"] == expected).lower(),
        })

    full_summary_rows: list[dict] = []
    for target in population_order:
        rows = [row for row in decision_rows if row["expected_population"] == target]
        full_summary_rows.append({
            "population": target, "libraries": len(rows),
            "winner_correct": sum(row["winner_correct"] == "true" for row in rows),
            "assigned_correct": sum(row["assigned_correct"] == "true" for row in rows),
            "unassigned": sum(row["label"] == "unassigned" for row in rows),
            "winner_accuracy": f"{sum(row['winner_correct']=='true' for row in rows)/len(rows):.8f}",
            "assigned_correct_fraction": f"{sum(row['assigned_correct']=='true' for row in rows)/len(rows):.8f}",
        })

    # Target-membership holdout: replace the target column with the fold-specific
    # training membership, while complete-data competitor columns remain fixed.
    holdout_rows: list[dict] = []
    for target in population_order:
        for plan in sorted(plan_by_target[target], key=lambda row: int(row["fold"])):
            heldout = [value for value in plan["heldout_units"].split(";") if value]
            if not heldout:
                continue
            fold_genes = [gene for _, gene in sorted(rankings_by_target_fold[(target, int(plan["fold"]))])]
            if not fold_genes:
                for library_id in heldout:
                    holdout_rows.append({
                        "target": target, "fold": plan["fold"], "heldout_library": library_id,
                        "evidence_tier": plan["evidence_tier"], "fold_genes": 0,
                        "winner_population": "", "label": "unassigned", "top_score": "",
                        "runner_up_score": "", "score_margin": "", "target_score": "",
                        "target_rank": 0, "abstain_reason": "no_fold_genes",
                        "winner_correct": "false", "assigned_correct": "false",
                    })
                continue
            fold_raw = score_genes(fold_genes, all_ids)
            fold_calibrated = fold_raw.rank(method="average", pct=True)
            for library_id in heldout:
                values = calibrated.loc[library_id].copy()
                values.loc[target] = float(fold_calibrated.loc[library_id])
                result = reduce_values(values)
                target_score = float(values.loc[target])
                target_rank = 1 + int((values > target_score).sum())
                holdout_rows.append({
                    "target": target, "fold": plan["fold"], "heldout_library": library_id,
                    "evidence_tier": plan["evidence_tier"], "fold_genes": len(fold_genes),
                    "winner_population": result["winner"], "label": result["label"],
                    "top_score": f"{result['top_score']:.12g}", "runner_up_score": f"{result['runner_up_score']:.12g}",
                    "score_margin": f"{result['score_margin']:.12g}", "target_score": f"{target_score:.12g}",
                    "target_rank": target_rank, "abstain_reason": result["abstain_reason"],
                    "winner_correct": str(result["winner"] == target).lower(),
                    "assigned_correct": str(result["label"] == target).lower(),
                })

    membership_rows: list[dict] = []
    target_metrics: dict[str, dict] = {}
    for target in population_order:
        genes = final_sets[target]
        plan = plan_by_target[target][0]
        if genes:
            md = evaluation_metadata(target); ids = list(md["library_id"])
            target_metrics[target] = population_metrics(target, score_genes(genes, ids), md)
        for final_rank, gene in enumerate(genes, 1):
            membership_rows.append({
                "geneset_id": "IMMGEN_POP_" + slug(target), "target": target,
                "source_role": plan["source_role"], "lineage": plan["lineage"],
                "independent_units": plan["independent_units"], "evidence_tier": plan["evidence_tier"],
                "candidate_source": candidate_source[target], "gene_symbol": gene,
                "candidate_rank": candidate_rank[target][gene], "fold_presence": fold_presence[target][gene],
                "final_rank": final_rank, "candidate_genes": len(candidate_sets[target]), "final_genes": len(genes),
            })

    minimum_holdout = float(config["flat_default_qualification"]["minimum_target_membership_holdout_assigned_correct_fraction"])
    qualification_rows: list[dict] = []
    summary_rows: list[dict] = []
    label_eligible_targets: list[str] = []
    for target in population_order:
        plan = plan_by_target[target][0]
        validation = [row for row in holdout_rows if row["target"] == target]
        holdout_fraction = (sum(row["assigned_correct"] == "true" for row in validation) / len(validation)) if validation else float("nan")
        genes = final_sets[target]
        metrics = target_metrics.get(target)
        reasons = []
        if int(plan["independent_units"]) < 2:
            reasons.append("singleton")
        if len(genes) < int(config["method"]["small_panel_minimum_genes"]):
            reasons.append("fewer_than_five_genes")
        if metrics is None or metrics["target_rank"] != 1:
            reasons.append("complete_mean_rank_not_first")
        if not validation:
            reasons.append("no_target_membership_holdout")
        elif holdout_fraction < minimum_holdout:
            reasons.append("holdout_assigned_correct_fraction")
        eligible = not reasons
        if eligible:
            label_eligible_targets.append(target)
        qualification_rows.append({
            "geneset_id": "IMMGEN_POP_" + slug(target), "target": target,
            "lineage": plan["lineage"], "independent_units": plan["independent_units"],
            "evidence_tier": plan["evidence_tier"], "final_genes": len(genes),
            "complete_mean_target_rank": metrics["target_rank"] if metrics else 0,
            "holdout_libraries": len(validation),
            "holdout_winner_correct_fraction": f"{sum(row['winner_correct']=='true' for row in validation)/len(validation):.8f}" if validation else "",
            "holdout_assigned_correct_fraction": f"{holdout_fraction:.8f}" if validation else "",
            "flat_label_eligible": str(eligible).lower(),
            "ineligibility_reasons": ";".join(reasons),
        })
        status = "flat_label_candidate" if eligible else ("score_only_no_program" if not genes else "score_only_not_flat_qualified")
        summary_rows.append({
            "geneset_id": "IMMGEN_POP_" + slug(target), "target": target,
            "source_role": plan["source_role"], "lineage": plan["lineage"],
            "independent_units": plan["independent_units"], "evidence_tier": plan["evidence_tier"],
            "candidate_source": candidate_source[target], "candidate_genes": len(candidate_sets[target]),
            "final_genes": len(genes), "status": status,
            "ucell_gap": f"{metrics['gap']:.12g}" if metrics else "",
            "strongest_competitor": metrics["strongest_competitor"] if metrics else "",
            "target_rank": metrics["target_rank"] if metrics else 0,
            "target_library_win_fraction": f"{metrics['target_library_win_fraction']:.8f}" if metrics else "",
            "holdout_assigned_correct_fraction": f"{holdout_fraction:.8f}" if validation else "",
            "flat_label_eligible": str(eligible).lower(),
        })

    membership_fields = [
        "geneset_id", "target", "source_role", "lineage", "independent_units", "evidence_tier",
        "candidate_source", "gene_symbol", "candidate_rank", "fold_presence", "final_rank", "candidate_genes", "final_genes",
    ]
    write_tsv(results["membership"], membership_rows, membership_fields, compressed=True)
    write_tsv(results["target_summary"], summary_rows, list(summary_rows[0]))
    trace_fields = [
        "target", "iteration", "genes_before", "tested_deletion", "candidate_rank", "gap_before",
        "gap_after", "improvement", "strongest_competitor_after", "selected_deletion",
    ]
    write_tsv(intermediate["refinement_trace"], trace_rows, trace_fields, compressed=True)
    for key, targets in (("score_gmt", nonempty_targets), ("label_gmt", label_eligible_targets)):
        results[key].parent.mkdir(parents=True, exist_ok=True)
        with results[key].open("w") as handle:
            for target in targets:
                row = next(item for item in summary_rows if item["target"] == target)
                description = f"diagnostic_only|tier={row['evidence_tier']}|status={row['status']}"
                handle.write("\t".join([row["geneset_id"], description, *final_sets[target]]) + "\n")
    write_tsv(results["full_confusion"], decision_rows, list(decision_rows[0]))
    write_tsv(results["full_confusion_summary"], full_summary_rows, list(full_summary_rows[0]))
    write_tsv(results["holdout_validation"], holdout_rows, list(holdout_rows[0]))
    write_tsv(results["qualification_summary"], qualification_rows, list(qualification_rows[0]))

    full_winner_accuracy = sum(row["winner_correct"] == "true" for row in decision_rows) / len(decision_rows)
    full_assigned_accuracy = sum(row["assigned_correct"] == "true" for row in decision_rows) / len(decision_rows)
    run_summary = {
        "schema_version": "bdev.immgen_global_local_population_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "targets": len(population_order),
        "nonempty_score_programs": len(nonempty_targets),
        "flat_label_candidate_programs": len(label_eligible_targets),
        "complete_atlas_libraries": len(decision_rows),
        "complete_atlas_winner_accuracy": full_winner_accuracy,
        "complete_atlas_assigned_correct_fraction": full_assigned_accuracy,
        "complete_atlas_unassigned": sum(row["label"] == "unassigned" for row in decision_rows),
        "target_membership_holdout_libraries": len(holdout_rows),
        "target_membership_holdout_assigned_correct_fraction": sum(row["assigned_correct"] == "true" for row in holdout_rows) / len(holdout_rows),
        "single_cell_expression_read": False,
        "scworkbench_invoked": False,
        "feature_sets_registered": False,
        "artifacts_promoted": False,
        "accepted_sets_emitted": 0,
    }
    reports["final_summary"].parent.mkdir(parents=True, exist_ok=True)
    reports["final_summary"].write_text(json.dumps(run_summary, indent=2, sort_keys=True) + "\n")
    report_lines = [
        "# Global-local ImmGen population discovery",
        "",
        f"- Nonempty score programs: {len(nonempty_targets)}/64",
        f"- Flat-label candidate programs: {len(label_eligible_targets)}/64",
        f"- Complete-atlas winner accuracy: {100*full_winner_accuracy:.1f}% (optimistic)",
        f"- Complete-atlas assigned-correct fraction: {100*full_assigned_accuracy:.1f}%",
        f"- Target-membership holdout assigned-correct fraction: {100*run_summary['target_membership_holdout_assigned_correct_fraction']:.1f}%",
        "",
        "All products remain diagnostic. The score collection retains every nonempty program; the flat-label collection contains only programs passing the frozen qualification rules.",
    ]
    reports["report"].write_text("\n".join(report_lines) + "\n")

    workflow_dir = Path(__file__).resolve().parent
    generated = {**intermediate, **results, **reports}
    scripts = [workflow_dir / "01_prepare_joint_rankings.py", Path(__file__).resolve()]
    manifest = {
        "schema_version": "bdev.immgen_global_local_population_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in {"config": config_path, **inputs}.items()},
        "scripts": {path.name: {"path": str(path), "sha256": sha256(path)} for path in scripts},
        "outputs": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in generated.items()},
        "single_cell_expression_read": False,
        "scworkbench_invoked": False,
        "feature_sets_registered": False,
        "artifacts_promoted": False,
        "accepted_sets_emitted": 0,
    }
    manifest_path = results_root / config["outputs"]["manifest"]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(run_summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
