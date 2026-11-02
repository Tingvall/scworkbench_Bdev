#!/usr/bin/env python3
"""Refine fold-specific positive sets on training profiles and evaluate holdouts."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import itertools
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
                    writer.writeheader(); writer.writerows(rows)
    else:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    require(config["schema_version"] == "bdev.immgen_positive_refinement_nested_qualification.v1", "Unexpected config schema")
    config_dir = args.config.parent.resolve()
    inputs = {name: resolve(config_dir, spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        observed = sha256(path)
        require(observed == config["inputs"][name]["sha256"], f"Input checksum mismatch for {name}: {observed}")
    contract = yaml.safe_load(inputs["statistical_contract"].read_text())
    output_root = Path(config["outputs"]["root"])
    tables = {name: output_root / rel for name, rel in config["outputs"]["tables"].items()}
    reports = {name: output_root / rel for name, rel in config["outputs"]["reports"].items()}
    require(tables["fold_plan"].is_file() and tables["fold_candidate_pools"].is_file() and tables["fold_model_audit"].is_file(), "Fold-model outputs are incomplete")

    fold_plan = read_tsv(tables["fold_plan"])
    pool_rows = read_tsv(tables["fold_candidate_pools"])
    pools: dict[str, list[str]] = defaultdict(list)
    pool_rank: dict[tuple[str, str], int] = {}
    for row in pool_rows:
        pools[row["fold_id"]].append(row["gene_symbol"])
        pool_rank[(row["fold_id"], row["gene_symbol"])] = int(row["pool_rank"])

    counts = pd.read_csv(inputs["selected_counts"], sep="\t", index_col=0)
    metadata = pd.read_csv(inputs["metadata"], sep="\t", dtype=str).fillna("")
    require(list(counts.columns) == list(metadata["library_id"]), "Counts and metadata are misaligned")
    discovery = metadata[metadata["role"] == "discovery"].copy()
    discovery_ids = list(discovery["library_id"])
    discovery_counts = counts[discovery_ids]
    populations = list(dict.fromkeys(discovery["population"]))
    population_counts = pd.DataFrame({population: discovery_counts[list(discovery.loc[discovery["population"] == population, "library_id"])].sum(axis=1) for population in populations})
    population_ranks = population_counts.rank(axis=0, method="average", pct=True, ascending=True)
    library_ranks = discovery_counts.rank(axis=0, method="average", pct=True, ascending=True)
    population_meta = discovery.drop_duplicates("population").set_index("population").loc[populations]
    library_meta = discovery.set_index("library_id").loc[discovery_ids]

    plan = read_tsv(inputs["contrast_plan"])
    competitor_map: dict[tuple[str, str], list[str]] = {}
    for row in fold_plan:
        key = (row["family"], row["target"])
        if key in competitor_map:
            continue
        competitors = []
        for test in plan:
            if (test["family"], test["target"]) == key:
                competitors.extend(group for group, weight in json.loads(test["weights_json"]).items() if float(weight) < 0)
        competitor_map[key] = list(dict.fromkeys(competitors))

    def family_data(family: str) -> tuple[pd.DataFrame, pd.Series]:
        if family == "b_local_stage":
            ids = list(library_meta.index[library_meta["table_s1_lineage"] == "B"])
            return library_ranks[ids], library_meta.loc[ids, "population"]
        if family == "b_branch":
            pops = list(population_meta.index[population_meta["table_s1_lineage"] == "B"])
            return population_ranks[pops], population_meta.loc[pops, "b_branch"]
        column = "superlineage" if family == "superlineage" else "table_s1_lineage"
        return population_ranks, population_meta[column]

    def training_metrics(score: pd.Series, groups: pd.Series, target: str, competitors: list[str]) -> dict:
        medians = {group: float(np.median(score.loc[groups.index[groups == group]])) for group in [target, *competitors] if (groups == group).any()}
        require(target in medians and len(medians) > 1, f"Missing training groups for {target}")
        comp = {group: val for group, val in medians.items() if group != target}
        strongest = max(comp, key=lambda group: (comp[group], group))
        sd = float(np.std(list(medians.values()), ddof=1))
        gap = (medians[target] - comp[strongest]) / sd if sd > 0 else float("-inf")
        return {"gap": gap, "target_median": medians[target], "strongest_competitor": strongest, "strongest_median": comp[strongest], "competitor_medians": comp}

    production_members: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in read_tsv(inputs["production_validation_sets"]):
        production_members[(row["family"], row["target"])].append(row["gene_symbol"])

    fold_results: list[dict] = []
    membership_rows: list[dict] = []
    trace_rows: list[dict] = []
    final_sets: dict[tuple[str, str], list[set[str]]] = defaultdict(list)
    eligible_rows = [row for row in fold_plan if row["qualification_eligible"].lower() == "true"]
    for fold_number, fold in enumerate(eligible_rows, 1):
        family, target, heldout, fold_id = fold["family"], fold["target"], fold["heldout_unit"], fold["fold_id"]
        matrix, groups = family_data(family)
        require(heldout in matrix.columns, f"Held-out profile absent: {fold_id}")
        training_columns = [column for column in matrix.columns if column != heldout]
        training_matrix = matrix[training_columns]
        training_groups = groups.loc[training_columns]
        competitors = [group for group in competitor_map[(family, target)] if (training_groups == group).any()]
        genes = [gene for gene in pools[fold_id] if gene in matrix.index]
        minimum = int(contract["contrast_families"][family]["minimum_genes"])
        if len(genes) < minimum:
            fold_results.append({"family": family, "target": target, "heldout_unit": heldout, "fold_id": fold_id, "pool_genes": len(genes), "final_genes": 0, "minimum_genes": minimum, "training_initial_gap": "", "training_final_gap": "", "heldout_score": "", "strongest_competitor": "", "strongest_competitor_median": "", "heldout_standardized_gap": "", "heldout_target_rank": "", "heldout_pass": "false", "production_fixed_set_genes": len(production_members[(family, target)]), "production_fixed_set_heldout_pass": "", "fold_status": "pool_below_minimum"})
            continue
        ranks = {gene: pool_rank[(fold_id, gene)] for gene in genes}
        current = list(genes)
        score_sum = training_matrix.loc[current].sum(axis=0)
        current_score = score_sum / len(current)
        initial = training_metrics(current_score, training_groups, target, competitors)
        current_metrics = initial
        iteration = 0
        while len(current) > minimum:
            choices = []
            for gene in current:
                candidate_score = (score_sum - training_matrix.loc[gene]) / (len(current) - 1)
                candidate_metrics = training_metrics(candidate_score, training_groups, target, competitors)
                choices.append((candidate_metrics["gap"], -ranks[gene], gene, candidate_metrics))
            best = max(choices, key=lambda item: (item[0], item[1]))
            if not np.isfinite(best[0]) or best[0] <= current_metrics["gap"] + 1e-12:
                break
            iteration += 1; before = current_metrics; removed = best[2]
            current.remove(removed); score_sum = score_sum - training_matrix.loc[removed]; current_metrics = best[3]
            trace_rows.append({"family": family, "target": target, "heldout_unit": heldout, "fold_id": fold_id, "iteration": iteration, "removed_gene": removed, "genes_before": len(current) + 1, "genes_after": len(current), "gap_before": f"{before['gap']:.10g}", "gap_after": f"{current_metrics['gap']:.10g}", "improvement": f"{current_metrics['gap'] - before['gap']:.10g}"})

        heldout_score = float(matrix.loc[current, heldout].mean())
        strongest_group = max(current_metrics["competitor_medians"], key=lambda group: (current_metrics["competitor_medians"][group], group))
        strongest_value = current_metrics["competitor_medians"][strongest_group]
        comparison_values = [heldout_score, *current_metrics["competitor_medians"].values()]
        heldout_sd = float(np.std(comparison_values, ddof=1))
        heldout_gap = (heldout_score - strongest_value) / heldout_sd if heldout_sd > 0 else float("-inf")
        heldout_rank = 1 + sum(value >= heldout_score for value in current_metrics["competitor_medians"].values())
        heldout_pass = heldout_score > strongest_value

        production = [gene for gene in production_members[(family, target)] if gene in matrix.index]
        if production:
            production_score = matrix.loc[production].mean(axis=0)
            production_training = training_metrics(production_score.loc[training_columns], training_groups, target, competitors)
            production_pass = float(production_score.loc[heldout]) > max(production_training["competitor_medians"].values())
        else:
            production_pass = None
        fold_results.append({"family": family, "target": target, "heldout_unit": heldout, "fold_id": fold_id, "pool_genes": len(genes), "final_genes": len(current), "minimum_genes": minimum, "training_initial_gap": f"{initial['gap']:.10g}", "training_final_gap": f"{current_metrics['gap']:.10g}", "heldout_score": f"{heldout_score:.10g}", "strongest_competitor": strongest_group, "strongest_competitor_median": f"{strongest_value:.10g}", "heldout_standardized_gap": f"{heldout_gap:.10g}", "heldout_target_rank": heldout_rank, "heldout_pass": str(heldout_pass).lower(), "production_fixed_set_genes": len(production), "production_fixed_set_heldout_pass": "" if production_pass is None else str(production_pass).lower(), "fold_status": "evaluated"})
        final_sets[(family, target)].append(set(current))
        for final_rank, gene in enumerate(current, 1):
            membership_rows.append({"family": family, "target": target, "heldout_unit": heldout, "fold_id": fold_id, "gene_symbol": gene, "final_rank": final_rank, "pool_rank": ranks[gene]})
        print(f"evaluated {fold_number}/{len(eligible_rows)} {fold_id}: {heldout_pass}")

    stable_rows: list[dict] = []
    target_rows: list[dict] = []
    pass_cutoff = float(config["qualification"]["minimum_nested_pass_fraction"])
    stable_cutoff = float(config["qualification"]["stable_membership_frequency"])
    target_inventory = list(dict.fromkeys((row["family"], row["target"]) for row in fold_plan))
    for family, target in target_inventory:
        inventory = [row for row in fold_plan if (row["family"], row["target"]) == (family, target)]
        target_units = int(inventory[0]["target_units"])
        minimum = int(contract["contrast_families"][family]["minimum_genes"])
        if inventory[0]["qualification_eligible"].lower() != "true":
            target_rows.append({"family": family, "target": target, "target_units": target_units, "folds": 0, "constructed_folds": 0, "heldout_passes": 0, "heldout_pass_fraction": "", "stable_genes_at_80pct": 0, "minimum_genes": minimum, "median_pairwise_jaccard": "", "production_fixed_set_fold_pass_fraction": "", "method_screen_status": "not_qualified_low_replication", "accepted_gene_set": "false"})
            continue
        results = [row for row in fold_results if (row["family"], row["target"]) == (family, target)]
        folds = len(results); constructed = sum(row["fold_status"] == "evaluated" for row in results); passes = sum(row["heldout_pass"] == "true" for row in results)
        pass_fraction = passes / folds
        sets = final_sets[(family, target)]
        frequencies = Counter(gene for genes in sets for gene in genes)
        stable_genes = sum(count / folds >= stable_cutoff for count in frequencies.values())
        for gene, count in sorted(frequencies.items(), key=lambda item: (-item[1], item[0])):
            stable_rows.append({"family": family, "target": target, "gene_symbol": gene, "selected_folds": count, "folds": folds, "selection_frequency": f"{count / folds:.8f}", "stable_at_80pct": str(count / folds >= stable_cutoff).lower()})
        jaccards = [len(a & b) / len(a | b) if a | b else 1.0 for a, b in itertools.combinations(sets, 2)]
        median_jaccard = float(np.median(jaccards)) if jaccards else 1.0
        production_assessed = [row for row in results if row["production_fixed_set_heldout_pass"]]
        production_fraction = sum(row["production_fixed_set_heldout_pass"] == "true" for row in production_assessed) / len(production_assessed) if production_assessed else None
        if constructed < folds:
            status = "method_screen_fail_pool_construction"
        elif pass_fraction < pass_cutoff:
            status = "method_screen_fail_nested_discrimination"
        elif stable_genes < minimum:
            status = "method_screen_fail_membership_stability"
        else:
            status = "method_screen_pass"
        target_rows.append({"family": family, "target": target, "target_units": target_units, "folds": folds, "constructed_folds": constructed, "heldout_passes": passes, "heldout_pass_fraction": f"{pass_fraction:.8f}", "stable_genes_at_80pct": stable_genes, "minimum_genes": minimum, "median_pairwise_jaccard": f"{median_jaccard:.8f}", "production_fixed_set_fold_pass_fraction": "" if production_fraction is None else f"{production_fraction:.8f}", "method_screen_status": status, "accepted_gene_set": "false"})

    write_tsv(tables["fold_results"], fold_results, list(fold_results[0]))
    write_tsv(tables["fold_final_membership"], membership_rows, list(membership_rows[0]), compressed=True)
    write_tsv(tables["fold_removal_trace"], trace_rows, list(trace_rows[0]) if trace_rows else ["family", "target", "fold_id"], compressed=True)
    write_tsv(tables["target_summary"], target_rows, list(target_rows[0]))
    write_tsv(tables["stable_membership"], stable_rows, list(stable_rows[0]), compressed=True)

    summary = {
        "schema_version": "bdev.immgen_positive_nested_qualification_summary.v1",
        "status": "complete_method_screen_only",
        "run_identity": config["run_identity"],
        "targets": len(target_rows),
        "eligible_targets": sum(row["folds"] > 0 for row in target_rows),
        "low_replication_targets": sum(row["method_screen_status"] == "not_qualified_low_replication" for row in target_rows),
        "folds": len(fold_results),
        "heldout_passes": sum(row["heldout_pass"] == "true" for row in fold_results),
        "method_screen_status_counts": dict(sorted(Counter(row["method_screen_status"] for row in target_rows).items())),
        "positive_only": True,
        "accepted_sets_emitted": 0,
        "importable_artifact_created": False,
        "scworkbench_invoked": False,
    }
    reports["qualification_summary"].parent.mkdir(parents=True, exist_ok=True)
    reports["qualification_summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    generated = {**tables, **reports}
    script_dir = Path(__file__).resolve().parent
    scripts = sorted(script_dir.glob("*.py")) + sorted(script_dir.glob("*.R"))
    manifest = {
        "schema_version": "bdev.immgen_positive_nested_qualification_manifest.v1",
        "status": "complete_method_screen_only",
        "run_identity": config["run_identity"],
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in {"qualification_config": args.config.resolve(), **inputs}.items()},
        "scripts": {path.name: sha256(path) for path in scripts},
        "outputs": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in generated.items()},
        "positive_only": True,
        "accepted_sets_emitted": 0,
        "importable_artifact_created": False,
    }
    manifest_path = output_root / config["outputs"]["manifest"]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
