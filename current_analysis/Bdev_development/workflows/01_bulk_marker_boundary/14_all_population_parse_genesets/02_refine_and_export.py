#!/usr/bin/env python3
"""Assemble, refine, and export diagnostic ImmGen population genesets."""

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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    require(config["schema_version"] == "bdev.immgen_all_population_parse.v1", "Unexpected config schema")
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
    max_rank = int(config["method"]["scoring"]["max_rank"])
    capped_ranks = counts.rank(axis=0, method="average", ascending=False).clip(upper=max_rank)
    discovery_role = config["target_scope"]["discovery_role"]
    heldout_role = config["target_scope"]["held_out_role"]

    population_info = metadata.drop_duplicates("population").set_index("population")

    def evaluation_metadata(target: str) -> pd.DataFrame:
        info = population_info.loc[target]
        local = metadata[(metadata["table_s1_lineage"] == info["table_s1_lineage"]) & (metadata["role"] == discovery_role)].copy()
        if info["role"] == heldout_role:
            local = pd.concat([local, metadata[metadata["population"] == target].copy()], ignore_index=True)
        return local

    def score_from_rank_sum(rank_sum: pd.Series, n_genes: int) -> pd.Series:
        baseline = n_genes * (n_genes + 1) / 2
        denominator = n_genes * max_rank - baseline
        require(denominator > 0, f"Invalid score denominator for {n_genes} genes")
        return (1 - (rank_sum - baseline) / denominator).clip(lower=0, upper=1)

    def metrics(target: str, library_scores: pd.Series, md: pd.DataFrame) -> dict:
        indexed = md.set_index("library_id").loc[list(library_scores.index)]
        populations = list(dict.fromkeys(indexed["population"]))
        population_scores = {
            population: float(library_scores.loc[indexed.index[indexed["population"] == population]].mean())
            for population in populations
        }
        require(target in population_scores, f"Target absent from score inventory: {target}")
        competitors = {population: score for population, score in population_scores.items() if population != target}
        strongest = max(competitors, key=lambda population: (competitors[population], population))
        sd = float(np.std(list(population_scores.values()), ddof=1))
        gap = (population_scores[target] - competitors[strongest]) / sd if sd > 0 else float("-inf")
        target_ids = indexed.index[indexed["population"] == target]
        return {
            "gap": gap,
            "strongest_competitor": strongest,
            "target_rank": 1 + sum(value >= population_scores[target] for value in competitors.values()),
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
    require(set(plan_by_target) == set(population_order), "Fold plan does not cover every population")

    candidate_sets: dict[str, list[str]] = {}
    candidate_source: dict[str, str] = {}
    candidate_rank: dict[str, dict[str, int]] = {}
    fold_presence: dict[str, dict[str, int]] = {}
    for target in population_order:
        folds = sorted({int(row["fold"]) for row in plan_by_target[target]})
        fold_lists = [[gene for _, gene in sorted(rankings_by_target_fold[(target, fold)])] for fold in folds]
        require(all(fold_lists), f"At least one fold produced no genes for {target}")
        if len(fold_lists) == 1:
            candidate = fold_lists[0]
            source = "single_full_data_top100"
        else:
            intersection = set(fold_lists[0]).intersection(*fold_lists[1:])
            if intersection:
                candidate = sorted(
                    intersection,
                    key=lambda gene: (np.mean([genes.index(gene) + 1 for genes in fold_lists]), max(genes.index(gene) + 1 for genes in fold_lists), gene),
                )
                source = "crossfold_intersection"
            else:
                union = set().union(*fold_lists)
                candidate = sorted(
                    union,
                    key=lambda gene: (-sum(gene in genes for genes in fold_lists), np.mean([genes.index(gene) + 1 for genes in fold_lists if gene in genes]), gene),
                )[: int(config["method"]["empty_intersection_fallback"]["maximum_genes"])]
                source = "empty_intersection_ranked_union_fallback"
        require(candidate, f"No candidate genes for {target}")
        candidate_sets[target] = candidate
        candidate_source[target] = source
        candidate_rank[target] = {gene: rank for rank, gene in enumerate(candidate, 1)}
        fold_presence[target] = {gene: sum(gene in genes for genes in fold_lists) for gene in candidate}

    final_sets: dict[str, list[str]] = {}
    trace_rows: list[dict] = []
    standard_floor = int(config["method"]["refinement"]["standard_floor"])
    for target_i, target in enumerate(population_order, 1):
        current = list(candidate_sets[target])
        md = evaluation_metadata(target)
        ids = list(md["library_id"])
        rank_sum = capped_ranks.loc[current, ids].sum(axis=0)
        current_metrics = metrics(target, score_from_rank_sum(rank_sum, len(current)), md)
        floor = standard_floor if len(current) >= standard_floor else len(current)
        iteration = 0
        while len(current) > floor:
            genes_after = len(current) - 1
            baseline = genes_after * (genes_after + 1) / 2
            denominator = genes_after * max_rank - baseline
            require(denominator > 0, f"Invalid deletion-score denominator for {target}")
            deletion_sums = rank_sum.to_numpy()[None, :] - capped_ranks.loc[current, ids].to_numpy(dtype=float)
            deletion_scores = np.clip(1 - (deletion_sums - baseline) / denominator, 0, 1)
            aligned = md.set_index("library_id").loc[ids]
            populations = list(dict.fromkeys(aligned["population"]))
            population_scores = np.column_stack(
                [deletion_scores[:, np.flatnonzero(aligned["population"].to_numpy() == population)].mean(axis=1) for population in populations]
            )
            target_column = populations.index(target)
            competitor_columns = [index for index, population in enumerate(populations) if population != target]
            strongest_columns = np.array(
                [max(competitor_columns, key=lambda index: (population_scores[row, index], populations[index])) for row in range(len(current))]
            )
            row_index = np.arange(len(current))
            standard_deviation = population_scores.std(axis=1, ddof=1)
            gaps = np.divide(
                population_scores[:, target_column] - population_scores[row_index, strongest_columns],
                standard_deviation,
                out=np.full(len(current), -np.inf),
                where=standard_deviation > 0,
            )
            best_index = max(range(len(current)), key=lambda index: (gaps[index], -candidate_rank[target][current[index]]))
            improves = np.isfinite(gaps[best_index]) and gaps[best_index] > current_metrics["gap"]
            iteration += 1
            for index, gene in enumerate(current):
                gap = gaps[index]
                trace_rows.append(
                    {
                        "target": target,
                        "iteration": iteration,
                        "genes_before": len(current),
                        "tested_deletion": gene,
                        "candidate_rank": candidate_rank[target][gene],
                        "gap_before": f"{current_metrics['gap']:.12g}",
                        "gap_after": f"{gap:.12g}",
                        "improvement": f"{gap-current_metrics['gap']:.12g}",
                        "strongest_competitor_after": populations[strongest_columns[index]],
                        "selected_deletion": str(improves and index == best_index).lower(),
                    }
                )
            if not improves:
                break
            selected_gene = current[best_index]
            rank_sum = rank_sum - capped_ranks.loc[selected_gene, ids]
            current_metrics = {
                "gap": float(gaps[best_index]),
                "strongest_competitor": populations[strongest_columns[best_index]],
            }
            current.remove(selected_gene)
        final_sets[target] = current
        print(f"refined {target_i} of {len(population_order)} {target}: {len(candidate_sets[target])} -> {len(current)}", flush=True)

    membership_rows: list[dict] = []
    summary_rows: list[dict] = []
    for target in population_order:
        plan = plan_by_target[target][0]
        genes = final_sets[target]
        md = evaluation_metadata(target)
        ids = list(md["library_id"])
        rank_sum = capped_ranks.loc[genes, ids].sum(axis=0)
        result = metrics(target, score_from_rank_sum(rank_sum, len(genes)), md)
        small_panel = len(genes) < int(config["method"]["standard_minimum_genes"])
        fallback = candidate_source[target] == "empty_intersection_ranked_union_fallback"
        if fallback:
            status = "provisional_ranked_union_fallback"
        elif plan["evidence_tier"] == "provisional_singleton":
            status = "provisional_singleton"
        elif small_panel:
            status = "small_crossfold_intersection"
        else:
            status = "diagnostic_crossfold_set"
        if plan["source_role"] == heldout_role:
            status = "heldout_derived_" + status
        for final_rank, gene in enumerate(genes, 1):
            membership_rows.append(
                {
                    "geneset_id": "IMMGEN_POP_" + slug(target),
                    "target": target,
                    "source_role": plan["source_role"],
                    "lineage": plan["lineage"],
                    "independent_units": plan["independent_units"],
                    "evidence_tier": plan["evidence_tier"],
                    "candidate_source": candidate_source[target],
                    "gene_symbol": gene,
                    "candidate_rank": candidate_rank[target][gene],
                    "fold_presence": fold_presence[target][gene],
                    "final_rank": final_rank,
                    "candidate_genes": len(candidate_sets[target]),
                    "final_genes": len(genes),
                    "status": status,
                }
            )
        summary_rows.append(
            {
                "geneset_id": "IMMGEN_POP_" + slug(target),
                "target": target,
                "source_role": plan["source_role"],
                "lineage": plan["lineage"],
                "independent_units": plan["independent_units"],
                "evidence_tier": plan["evidence_tier"],
                "candidate_source": candidate_source[target],
                "candidate_genes": len(candidate_sets[target]),
                "final_genes": len(genes),
                "status": status,
                "ucell_gap": f"{result['gap']:.12g}",
                "strongest_competitor": result["strongest_competitor"],
                "target_rank": result["target_rank"],
                "target_library_win_fraction": f"{result['target_library_win_fraction']:.8f}",
                "full_data_mean_pass": str(result["target_rank"] == 1).lower(),
            }
        )

    write_tsv(results["membership"], membership_rows, list(membership_rows[0]), compressed=True)
    write_tsv(results["target_summary"], summary_rows, list(summary_rows[0]))
    write_tsv(intermediate["refinement_trace"], trace_rows, list(trace_rows[0]), compressed=True)
    results["candidate_gmt"].parent.mkdir(parents=True, exist_ok=True)
    with results["candidate_gmt"].open("w") as handle:
        for row in summary_rows:
            target = row["target"]
            description = f"diagnostic_only|role={row['source_role']}|tier={row['evidence_tier']}|status={row['status']}"
            handle.write("\t".join([row["geneset_id"], description, *final_sets[target]]) + "\n")

    status_counts = Counter(row["status"] for row in summary_rows)
    tier_counts = Counter(row["evidence_tier"] for row in summary_rows)
    summary = {
        "schema_version": "bdev.immgen_all_population_parse_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "targets": len(summary_rows),
        "targets_with_set": sum(bool(final_sets[target]) for target in population_order),
        "discovery_targets": sum(row["source_role"] == discovery_role for row in summary_rows),
        "held_out_derived_targets": sum(row["source_role"] == heldout_role for row in summary_rows),
        "full_data_mean_passes": sum(row["full_data_mean_pass"] == "true" for row in summary_rows),
        "small_panels": sum(int(row["final_genes"]) < int(config["method"]["standard_minimum_genes"]) for row in summary_rows),
        "ranked_union_fallbacks": sum(row["candidate_source"] == "empty_intersection_ranked_union_fallback" for row in summary_rows),
        "evidence_tier_counts": dict(sorted(tier_counts.items())),
        "target_status_counts": dict(sorted(status_counts.items())),
        "single_cell_expression_read": False,
        "accepted_sets_emitted": 0,
        "candidate_gmt_created": True,
        "importable_artifact_created": False,
        "scworkbench_invoked": False,
    }
    reports["final_summary"].parent.mkdir(parents=True, exist_ok=True)
    reports["final_summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    workflow_dir = Path(__file__).resolve().parent
    generated = {**intermediate, **results, **reports}
    scripts = [workflow_dir / "01_fit_population_rankings.R", Path(__file__).resolve()]
    manifest = {
        "schema_version": "bdev.immgen_all_population_parse_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in {"config": config_path, **inputs}.items()},
        "scripts": {path.name: {"path": str(path), "sha256": sha256(path)} for path in scripts},
        "outputs": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in generated.items()},
        "discovery_validation_boundary_preserved": True,
        "single_cell_expression_read": False,
        "accepted_sets_emitted": 0,
        "candidate_gmt_created": True,
        "importable_artifact_created": False,
        "scworkbench_invoked": False,
    }
    manifest_path = results_root / config["outputs"]["manifest"]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
