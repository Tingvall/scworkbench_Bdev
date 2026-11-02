#!/usr/bin/env python3
"""Refine low-replication Parse candidates and assemble an adaptive collection."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
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


def as_bool(value: object) -> bool:
    return str(value).lower() == "true"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    require(config["schema_version"] == "bdev.immgen_replication_adaptive_parse_transfer.v1", "Unexpected config schema")
    config_dir = args.config.parent.resolve()
    inputs = {name: resolve(config_dir, spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        observed = sha256(path)
        require(observed == config["inputs"][name]["sha256"], f"Checksum mismatch for {name}: {observed}")

    output_root = Path(config["outputs"]["root"])
    tables = {name: output_root / rel for name, rel in config["outputs"]["tables"].items()}
    reports = {name: output_root / rel for name, rel in config["outputs"]["reports"].items()}
    for name in ("adaptive_plan", "low_replication_rankings", "low_replication_candidates", "model_audit"):
        require(tables[name].is_file(), f"Missing model output {name}: {tables[name]}")

    counts = pd.read_csv(inputs["selected_counts"], sep="\t", index_col=0)
    metadata = pd.read_csv(inputs["metadata"], sep="\t", dtype=str).fillna("")
    require(list(counts.columns) == list(metadata["library_id"]), "Count columns and metadata are misaligned")
    max_rank = int(config["method"]["scoring"]["max_rank"])
    capped_ranks = counts.rank(axis=0, method="average", ascending=False).clip(upper=max_rank)
    discovery_meta = metadata[metadata["role"] == "discovery"].copy()
    heldout_meta = metadata[metadata["role"] == "held_out_validation"].copy()
    approved_targets: dict[str, list[str]] = defaultdict(list)
    for row in read_tsv(inputs["contrast_plan"]):
        if row["target"] not in approved_targets[row["family"]]:
            approved_targets[row["family"]].append(row["target"])

    group_column = {"superlineage": "superlineage", "broad_lineage": "table_s1_lineage", "b_branch": "b_branch", "b_local_stage": "population"}

    def family_metadata(family: str) -> pd.DataFrame:
        local = discovery_meta.copy()
        if family in {"b_branch", "b_local_stage"}:
            local = local[local["table_s1_lineage"] == "B"]
        if family == "b_local_stage":
            local = local[local["population"].isin(approved_targets[family])]
        return local

    def score_from_rank_sum(rank_sum: pd.Series, n_genes: int) -> pd.Series:
        baseline = n_genes * (n_genes + 1) / 2
        denominator = n_genes * max_rank - baseline
        require(denominator > 0, f"Invalid score denominator for {n_genes} genes")
        return (1 - (rank_sum - baseline) / denominator).clip(lower=0, upper=1)

    def gene_set_score(genes: list[str], library_ids: list[str]) -> pd.Series:
        require(genes, "Cannot score an empty set")
        require(all(gene in capped_ranks.index for gene in genes), "Membership gene absent from count matrix")
        return score_from_rank_sum(capped_ranks.loc[genes, library_ids].sum(axis=0), len(genes))

    def unit_scores_and_groups(family: str, library_scores: pd.Series, md: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
        local = md.set_index("library_id").loc[list(library_scores.index)]
        if family == "b_local_stage":
            return library_scores, local["population"]
        units = list(dict.fromkeys(local["population"]))
        scores = pd.Series({unit: float(library_scores.loc[local.index[local["population"] == unit]].mean()) for unit in units})
        unit_md = local.drop_duplicates("population").set_index("population").loc[units]
        return scores, unit_md[group_column[family]]

    def metrics(family: str, target: str, library_scores: pd.Series, md: pd.DataFrame) -> dict:
        unit_scores, groups = unit_scores_and_groups(family, library_scores, md)
        group_values = {group: float(np.median(unit_scores.loc[groups.index[groups == group]])) for group in dict.fromkeys(groups)}
        require(target in group_values, f"Target absent from scores: {family}/{target}")
        competitors = {group: value for group, value in group_values.items() if group != target}
        strongest = max(competitors, key=lambda group: (competitors[group], group))
        sd = float(np.std(list(group_values.values()), ddof=1))
        gap = (group_values[target] - competitors[strongest]) / sd if sd > 0 else float("-inf")
        target_scores = unit_scores.loc[groups.index[groups == target]]
        return {"gap": gap, "strongest_competitor": strongest, "target_rank": 1 + sum(value >= group_values[target] for value in competitors.values()), "target_unit_win_fraction": float(np.mean(target_scores > competitors[strongest]))}

    plan_rows = read_tsv(tables["adaptive_plan"])
    target_order: list[tuple[str, str]] = []
    tier: dict[tuple[str, str], str] = {}
    units: dict[tuple[str, str], int] = {}
    for row in plan_rows:
        key = (row["family"], row["target"])
        if key not in tier:
            target_order.append(key)
            tier[key] = row["evidence_tier"]
            units[key] = int(row["independent_units"])
    require(len(target_order) == 27, f"Expected 27 targets, found {len(target_order)}")
    require(Counter(tier.values()) == Counter({"high_threefold": 14, "moderate_twofold": 10, "provisional_singleton": 3}), "Unexpected evidence-tier inventory")

    exact_sets: dict[tuple[str, str], list[tuple[int, str]]] = defaultdict(list)
    exact_candidate_sizes: dict[tuple[str, str], int] = {}
    for row in read_tsv(inputs["exact_transfer_final_membership"]):
        key = (row["family"], row["target"])
        exact_sets[key].append((int(row["final_rank"]), row["gene_symbol"]))
        exact_candidate_sizes[key] = int(row["candidate_genes"])
    high_sets = {key: [gene for _, gene in sorted(rows)] for key, rows in exact_sets.items()}
    require(len(high_sets) == 14, "Expected 14 completed high-tier sets")

    candidate_sets: dict[tuple[str, str], list[str]] = defaultdict(list)
    candidate_rank: dict[tuple[str, str], dict[str, int]] = defaultdict(dict)
    construction_pass: dict[tuple[str, str], bool] = {}
    for row in read_tsv(tables["low_replication_candidates"]):
        key = (row["family"], row["target"])
        candidate_sets[key].append(row["gene_symbol"])
        candidate_rank[key][row["gene_symbol"]] = int(row["candidate_rank"])
        construction_pass[key] = as_bool(row["construction_pass"])
    require(len(candidate_sets) == 13, "Expected candidates for 13 low-replication targets")

    adaptive_sets: dict[tuple[str, str], list[str]] = dict(high_sets)
    trace_rows: list[dict] = []
    floor = int(config["method"]["refinement"]["floor"])
    for key in target_order:
        family, target = key
        if tier[key] == "high_threefold":
            continue
        if not construction_pass.get(key, False):
            adaptive_sets[key] = []
            continue
        current = list(candidate_sets[key])
        md = family_metadata(family)
        ids = list(md["library_id"])
        rank_sum = capped_ranks.loc[current, ids].sum(axis=0)
        current_scores = score_from_rank_sum(rank_sum, len(current))
        current_metrics = metrics(family, target, current_scores, md)
        iteration = 0
        while len(current) > floor:
            choices = []
            for gene in current:
                candidate_sum = rank_sum - capped_ranks.loc[gene, ids]
                candidate_scores = score_from_rank_sum(candidate_sum, len(current) - 1)
                candidate_metrics = metrics(family, target, candidate_scores, md)
                choices.append((candidate_metrics["gap"], -candidate_rank[key][gene], gene, candidate_metrics, candidate_sum))
            best = max(choices, key=lambda item: (item[0], item[1]))
            improves = np.isfinite(best[0]) and best[0] > current_metrics["gap"]
            iteration += 1
            for gap, _, gene, tested, _ in choices:
                trace_rows.append({"family": family, "target": target, "evidence_tier": tier[key], "iteration": iteration, "genes_before": len(current), "tested_deletion": gene, "candidate_rank": candidate_rank[key][gene], "gap_before": f"{current_metrics['gap']:.12g}", "gap_after": f"{gap:.12g}", "improvement": f"{gap-current_metrics['gap']:.12g}", "strongest_competitor_after": tested["strongest_competitor"], "selected_deletion": str(improves and gene == best[2]).lower()})
            if not improves:
                break
            current.remove(best[2])
            current_metrics = best[3]
            rank_sum = best[4]
        adaptive_sets[key] = current

    membership_rows: list[dict] = []
    target_rows: list[dict] = []
    for key in target_order:
        family, target = key
        genes = adaptive_sets.get(key, [])
        candidate_genes = exact_candidate_sizes[key] if tier[key] == "high_threefold" else len(candidate_sets.get(key, []))
        if genes:
            md = family_metadata(family)
            result = metrics(family, target, gene_set_score(genes, list(md["library_id"])), md)
            status = "reused_exact_threefold" if tier[key] == "high_threefold" else ("constructed_twofold" if tier[key] == "moderate_twofold" else "constructed_singleton_provisional")
        else:
            result = None
            status = "candidate_below_20"
        for final_rank, gene in enumerate(genes, 1):
            membership_rows.append({"family": family, "target": target, "independent_units": units[key], "evidence_tier": tier[key], "gene_symbol": gene, "final_rank": final_rank, "candidate_genes": candidate_genes, "final_genes": len(genes), "status": status})
        target_rows.append({"family": family, "target": target, "independent_units": units[key], "evidence_tier": tier[key], "candidate_genes": candidate_genes, "final_genes": len(genes), "status": status, "ucell_gap": f"{result['gap']:.12g}" if result else "", "strongest_competitor": result["strongest_competitor"] if result else "", "target_rank": result["target_rank"] if result else "", "target_unit_win_fraction": f"{result['target_unit_win_fraction']:.8f}" if result else "", "full_data_mean_pass": str(result["target_rank"] == 1).lower() if result else "false", "reproducibility_claim": "threefold_stability" if tier[key] == "high_threefold" else ("twofold_stability_weaker" if tier[key] == "moderate_twofold" else "none_singleton_provisional")})

    write_tsv(tables["adaptive_final_membership"], membership_rows, list(membership_rows[0]), compressed=True)
    write_tsv(tables["refinement_trace"], trace_rows, list(trace_rows[0]), compressed=True)
    write_tsv(tables["adaptive_target_summary"], target_rows, list(target_rows[0]))

    prior_comparison = read_tsv(inputs["prior_target_comparison"])
    comparison_fields = list(prior_comparison[0]) + ["evidence_tier"]
    comparison_rows = [{**row, "evidence_tier": ""} for row in prior_comparison]
    for row in target_rows:
        comparison_rows.append({"method": "replication_adaptive_parse", "family": row["family"], "target": row["target"], "independent_target_units": row["independent_units"], "genes": row["final_genes"], "status": row["status"], "ucell_gap": row["ucell_gap"], "strongest_competitor": row["strongest_competitor"], "target_rank": row["target_rank"], "target_unit_win_fraction": row["target_unit_win_fraction"], "full_data_mean_pass": row["full_data_mean_pass"], "selection_uses_full_discovery_scores": "true", "evidence_tier": row["evidence_tier"]})
    write_tsv(tables["method_target_comparison"], comparison_rows, comparison_fields)

    prior_heldout = read_tsv(inputs["prior_heldout_comparison"])
    heldout_rows = list(prior_heldout)
    broad_sets = {target: genes for (family, target), genes in adaptive_sets.items() if family == "broad_lineage" and genes}
    expected = heldout_meta.drop_duplicates("population").set_index("population")["table_s1_lineage"].to_dict()
    for population in dict.fromkeys(heldout_meta["population"]):
        ids = list(heldout_meta.loc[heldout_meta["population"] == population, "library_id"])
        scores = {target: float(gene_set_score(genes, ids).mean()) for target, genes in broad_sets.items()}
        predicted = max(scores, key=lambda target: (scores[target], target))
        heldout_rows.append({"method": "replication_adaptive_parse", "population": population, "expected_lineage": expected[population], "predicted_lineage": predicted, "correct": str(predicted == expected[population]).lower(), "winning_score": f"{scores[predicted]:.12g}", "broad_panel_sets": len(broad_sets), "scores_json": json.dumps(scores, sort_keys=True, separators=(",", ":"))})
    write_tsv(tables["heldout_broad_comparison"], heldout_rows, list(heldout_rows[0]))

    prior_methods = read_tsv(inputs["prior_method_summary"])
    method_fields = list(prior_methods[0]) + ["high_threefold_sets", "moderate_twofold_sets", "provisional_singleton_sets"]
    method_rows = [{**row, "high_threefold_sets": "", "moderate_twofold_sets": "", "provisional_singleton_sets": ""} for row in prior_methods]
    adaptive_heldout = [row for row in heldout_rows if row["method"] == "replication_adaptive_parse"]
    method_rows.append({"method": "replication_adaptive_parse", "targets": 27, "targets_with_set": sum(bool(adaptive_sets.get(key)) for key in target_order), "full_data_mean_passes": sum(row["full_data_mean_pass"] == "true" for row in target_rows), "heldout_broad_correct": sum(row["correct"] == "true" for row in adaptive_heldout), "heldout_broad_total": len(adaptive_heldout), "broad_panel_sets": len(broad_sets), "single_cell_expression_used": "false", "high_threefold_sets": sum(bool(adaptive_sets.get(key)) and tier[key] == "high_threefold" for key in target_order), "moderate_twofold_sets": sum(bool(adaptive_sets.get(key)) and tier[key] == "moderate_twofold" for key in target_order), "provisional_singleton_sets": sum(bool(adaptive_sets.get(key)) and tier[key] == "provisional_singleton" for key in target_order)})
    write_tsv(tables["method_summary"], method_rows, method_fields)

    adaptive_method = method_rows[-1]
    summary = {"schema_version": "bdev.immgen_replication_adaptive_parse_summary.v1", "status": "complete_diagnostic_only", "run_identity": config["run_identity"], "adaptive_method": adaptive_method, "target_status_counts": dict(sorted(Counter(row["status"] for row in target_rows).items())), "evidence_tier_counts": dict(sorted(Counter(row["evidence_tier"] for row in target_rows).items())), "singleton_targets": [row["target"] for row in target_rows if row["evidence_tier"] == "provisional_singleton"], "single_cell_expression_read": False, "accepted_sets_emitted": 0, "importable_artifact_created": False, "scworkbench_invoked": False}
    reports["adaptive_summary"].parent.mkdir(parents=True, exist_ok=True)
    reports["adaptive_summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    generated = {**tables, **reports}
    script_dir = Path(__file__).resolve().parent
    scripts = [script_dir / "01_fit_low_replication_rankings.R", Path(__file__).resolve()]
    manifest = {"schema_version": "bdev.immgen_replication_adaptive_parse_manifest.v1", "status": "complete_diagnostic_only", "run_identity": config["run_identity"], "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in {"adaptive_config": args.config.resolve(), **inputs}.items()}, "scripts": {path.name: {"path": str(path), "sha256": sha256(path)} for path in scripts}, "outputs": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in generated.items()}, "evidence_tiers_preserved": True, "single_cell_expression_read": False, "accepted_sets_emitted": 0, "importable_artifact_created": False}
    manifest_path = output_root / config["outputs"]["manifest"]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
