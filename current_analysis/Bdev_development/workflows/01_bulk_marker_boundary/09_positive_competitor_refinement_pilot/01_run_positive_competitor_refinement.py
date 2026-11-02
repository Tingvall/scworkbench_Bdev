#!/usr/bin/env python3
"""Run a positive-only strongest-competitor refinement diagnostic."""

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
                    writer.writeheader(); writer.writerows(rows)
    else:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)


def as_bool(value: str) -> bool:
    return value.lower() == "true"


def value(value_: str, default: float) -> float:
    return float(value_) if value_ else default


def initialize_state() -> dict:
    return {
        "filter_pass": True,
        "annotation_eligible": True,
        "support": 1.0,
        "expression": 1.0,
        "effect_pass": True,
        "fdr_passes": 0,
        "tests": 0,
        "worst_fdr_ratio": 0.0,
        "worst_effect_margin": float("inf"),
        "worst_auroc": 1.0,
        "ensembl_gene_id": "",
        "chromosome": "",
        "gene_biotype": "",
    }


def thresholds(family: str, test_type: str, contract: dict) -> tuple[float, float]:
    rule = contract["contrast_families"][family]
    if family == "superlineage":
        return float(rule["minimum_global_log2_fold_change"]), float(rule["maximum_global_fdr"])
    if family in {"broad_lineage", "b_branch"}:
        if test_type == "global":
            return float(rule["minimum_global_log2_fold_change"]), float(rule["maximum_global_fdr"])
        return float(rule["minimum_each_competitor_log2_fold_change"]), float(rule["maximum_each_competitor_fdr"])
    return float(rule["minimum_each_neighbor_log2_fold_change"]), float(rule["maximum_each_neighbor_fdr"])


def eligible(family: str, state: dict, contract: dict) -> bool:
    rule = contract["contrast_families"][family]
    support_ok = state["expression"] >= float(rule["minimum_target_library_fraction"]) if family == "b_local_stage" else state["support"] >= float(rule["minimum_target_population_support"])
    return state["filter_pass"] and state["annotation_eligible"] and support_ok and state["effect_pass"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    require(config["schema_version"] == "bdev.immgen_positive_competitor_refinement_pilot.v1", "Unexpected config schema")
    config_dir = args.config.parent.resolve()
    inputs = {name: resolve(config_dir, spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        observed = sha256(path)
        require(observed == config["inputs"][name]["sha256"], f"Checksum mismatch for {name}: {observed}")

    contract = yaml.safe_load(inputs["statistical_contract"].read_text())
    production_config = yaml.safe_load(inputs["production_execution_config"].read_text())
    require(production_config["status"] == "execution_complete_candidates_require_scientific_review", "Production run is not complete")
    output_root = Path(config["outputs"]["root"])
    tables = {name: output_root / rel for name, rel in config["outputs"]["tables"].items()}
    reports = {name: output_root / rel for name, rel in config["outputs"]["reports"].items()}
    pool_max = int(config["method"]["candidate_pool_maximum_genes"])

    # Construct positive-only pools from production statistics. FDR ranks but does not gate.
    pool_rows: list[dict] = []
    pools: dict[tuple[str, str], list[str]] = {}
    current_key: tuple[str, str] | None = None
    genes: dict[str, dict] = {}
    seen: set[tuple[str, str]] = set()

    def flush(key: tuple[str, str], states: dict[str, dict]) -> None:
        family, target = key
        require(key not in seen, f"Noncontiguous target block: {key}")
        ranked = []
        for gene, state in states.items():
            if not eligible(family, state, contract):
                continue
            state["fdr_pass_fraction"] = state["fdr_passes"] / state["tests"]
            ranked.append((gene, state))
        ranked.sort(key=lambda item: (-item[1]["fdr_pass_fraction"], item[1]["worst_fdr_ratio"], -item[1]["support"], -item[1]["worst_effect_margin"], -item[1]["worst_auroc"], item[0]))
        selected = ranked[:pool_max]
        pools[key] = [gene for gene, _ in selected]
        for rank, (gene, state) in enumerate(selected, 1):
            pool_rows.append({
                "family": family, "target": target, "gene_symbol": gene, "pool_rank": rank,
                "required_tests": state["tests"], "fdr_passing_tests": state["fdr_passes"],
                "fdr_pass_fraction": f"{state['fdr_pass_fraction']:.8f}", "worst_fdr_ratio": f"{state['worst_fdr_ratio']:.8g}",
                "target_support_fraction": f"{state['support']:.8g}", "target_expression_fraction": f"{state['expression']:.8g}",
                "worst_effect_margin": f"{state['worst_effect_margin']:.8g}", "worst_auroc": f"{state['worst_auroc']:.8g}",
                "ensembl_gene_id": state["ensembl_gene_id"], "chromosome": state["chromosome"], "gene_biotype": state["gene_biotype"],
            })
        seen.add(key)

    with gzip.open(inputs["complete_statistics"], "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["family"], row["target"])
            if current_key is not None and key != current_key:
                flush(current_key, genes); genes = {}
            current_key = key
            gene = row["gene_symbol"]
            state = genes.setdefault(gene, initialize_state())
            state["filter_pass"] &= as_bool(row["filter_pass"])
            state["annotation_eligible"] &= as_bool(row["core_identity_eligible"])
            state["support"] = min(state["support"], value(row["target_support_fraction"], 0.0))
            state["expression"] = min(state["expression"], value(row["target_expression_fraction"], 0.0))
            lfc_cutoff, fdr_cutoff = thresholds(row["family"], row["test_type"], contract)
            lfc = value(row["logFC"], float("-inf")); fdr = value(row["FDR"], 1.0)
            state["effect_pass"] &= lfc >= lfc_cutoff
            state["tests"] += 1
            state["fdr_passes"] += int(fdr <= fdr_cutoff)
            state["worst_fdr_ratio"] = max(state["worst_fdr_ratio"], fdr / fdr_cutoff)
            state["worst_effect_margin"] = min(state["worst_effect_margin"], lfc - lfc_cutoff)
            state["worst_auroc"] = min(state["worst_auroc"], value(row["target_versus_competitor_auroc"], 0.0))
            state["ensembl_gene_id"] = row["ensembl_gene_id"]
            state["chromosome"] = row["chromosome"]
            state["gene_biotype"] = row["gene_biotype"]
        if current_key is not None:
            flush(current_key, genes)

    set_summary = read_tsv(inputs["production_set_summary"])
    target_order = [(row["family"], row["target"]) for row in set_summary]
    require(set(target_order) == set(pools), "Pilot pools do not cover the production target inventory")
    write_tsv(tables["candidate_pool_membership"], pool_rows, list(pool_rows[0]), compressed=True)

    # Build within-profile expression-percentile ranks from counts.
    counts = pd.read_csv(inputs["selected_counts"], sep="\t", index_col=0)
    metadata = pd.read_csv(inputs["metadata"], sep="\t", dtype=str).fillna("")
    require(list(counts.columns) == list(metadata["library_id"]), "Count columns and metadata are misaligned")
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
    for family, target in target_order:
        competitors = []
        for row in plan:
            if (row["family"], row["target"]) == (family, target):
                weights = json.loads(row["weights_json"])
                competitors.extend(group for group, weight in weights.items() if float(weight) < 0)
        competitor_map[(family, target)] = list(dict.fromkeys(competitors))

    group_columns = {"superlineage": "superlineage", "broad_lineage": "table_s1_lineage", "b_branch": "b_branch", "b_local_stage": "population"}

    def family_data(family: str) -> tuple[pd.DataFrame, pd.Series]:
        if family == "b_local_stage":
            selected_ids = list(library_meta.index[library_meta["table_s1_lineage"] == "B"])
            return library_ranks[selected_ids], library_meta.loc[selected_ids, "population"]
        if family == "b_branch":
            selected_populations = list(population_meta.index[population_meta["table_s1_lineage"] == "B"])
            return population_ranks[selected_populations], population_meta.loc[selected_populations, "b_branch"]
        column = group_columns[family]
        return population_ranks, population_meta[column]

    def metrics(score: pd.Series, groups: pd.Series, target: str, competitors: list[str]) -> dict:
        group_values = {group: float(np.median(score.loc[groups.index[groups == group]])) for group in [target, *competitors] if (groups == group).any()}
        require(target in group_values and len(group_values) > 1, f"Missing score group for {target}")
        competitor_values = {group: val for group, val in group_values.items() if group != target}
        strongest = max(competitor_values, key=lambda group: (competitor_values[group], group))
        sd = float(np.std(list(group_values.values()), ddof=1))
        gap = (group_values[target] - competitor_values[strongest]) / sd if sd > 0 else float("-inf")
        target_scores = score.loc[groups.index[groups == target]]
        unit_wins = float(np.mean(target_scores > competitor_values[strongest]))
        return {"gap": gap, "target_value": group_values[target], "strongest_competitor": strongest, "strongest_value": competitor_values[strongest], "target_units": len(target_scores), "target_unit_win_fraction": unit_wins, "target_rank": 1 + sum(val >= group_values[target] for val in competitor_values.values())}

    production_members: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in read_tsv(inputs["production_validation_sets"]):
        production_members[(row["family"], row["target"])].append(row["gene_symbol"])

    summary_rows: list[dict] = []
    final_rows: list[dict] = []
    trace_rows: list[dict] = []
    profile_rows: list[dict] = []
    comparison_rows: list[dict] = []
    final_by_target: dict[tuple[str, str], list[str]] = {}
    for family, target in target_order:
        matrix, groups = family_data(family)
        competitors = [item for item in competitor_map[(family, target)] if (groups == item).any()]
        pool = [gene for gene in pools[(family, target)] if gene in matrix.index]
        minimum = int(contract["contrast_families"][family]["minimum_genes"])
        initial_size = len(pool)
        if initial_size < minimum:
            final_by_target[(family, target)] = []
            summary_rows.append({"family": family, "target": target, "pool_genes": initial_size, "final_genes": 0, "minimum_genes": minimum, "initial_gap": "", "final_gap": "", "gap_improvement": "", "strongest_competitor": "", "target_units": int((groups == target).sum()), "target_unit_win_fraction": "", "target_rank": "", "screen_status": "failed_pool_below_minimum", "optimistic_full_data_selection": "true"})
            comparison_rows.append({"family": family, "target": target, "production_genes": len(production_members[(family, target)]), "production_gap": "", "pilot_genes": 0, "pilot_gap": "", "pilot_minus_production_gap": "", "pilot_full_data_pass": "false", "pilot_unit_support_pass": "false", "comparison": "pilot_pool_below_minimum"})
            continue

        pool_rank = {gene: rank for rank, gene in enumerate(pool, 1)}
        current = list(pool)
        score_sum = matrix.loc[current].sum(axis=0)
        current_score = score_sum / len(current)
        initial_metrics = metrics(current_score, groups, target, competitors)
        current_metrics = initial_metrics
        iteration = 0
        while len(current) > minimum:
            choices = []
            for gene in current:
                candidate_score = (score_sum - matrix.loc[gene]) / (len(current) - 1)
                candidate_metrics = metrics(candidate_score, groups, target, competitors)
                choices.append((candidate_metrics["gap"], -pool_rank[gene], gene, candidate_metrics, candidate_score))
            best = max(choices, key=lambda item: (item[0], item[1]))
            if not np.isfinite(best[0]) or best[0] <= current_metrics["gap"] + 1e-12:
                break
            iteration += 1
            before = current_metrics
            removed = best[2]
            current.remove(removed)
            score_sum = score_sum - matrix.loc[removed]
            current_metrics = best[3]
            trace_rows.append({"family": family, "target": target, "iteration": iteration, "removed_gene": removed, "genes_before": len(current) + 1, "genes_after": len(current), "gap_before": f"{before['gap']:.10g}", "gap_after": f"{current_metrics['gap']:.10g}", "improvement": f"{current_metrics['gap'] - before['gap']:.10g}", "strongest_competitor_after": current_metrics["strongest_competitor"]})

        full_pass = current_metrics["gap"] > 0
        final_by_target[(family, target)] = list(current)
        support_pass = current_metrics["target_unit_win_fraction"] >= float(config["method"]["diagnostic_support_fraction"])
        if full_pass and support_pass and current_metrics["target_units"] >= 3:
            status = "diagnostic_positive_separation_with_unit_support"
        elif full_pass and support_pass:
            status = "diagnostic_positive_separation_low_replication"
        elif full_pass:
            status = "diagnostic_positive_separation_without_unit_support"
        else:
            status = "diagnostic_no_positive_separation"
        summary_rows.append({"family": family, "target": target, "pool_genes": initial_size, "final_genes": len(current), "minimum_genes": minimum, "initial_gap": f"{initial_metrics['gap']:.10g}", "final_gap": f"{current_metrics['gap']:.10g}", "gap_improvement": f"{current_metrics['gap'] - initial_metrics['gap']:.10g}", "strongest_competitor": current_metrics["strongest_competitor"], "target_units": current_metrics["target_units"], "target_unit_win_fraction": f"{current_metrics['target_unit_win_fraction']:.8f}", "target_rank": current_metrics["target_rank"], "screen_status": status, "optimistic_full_data_selection": "true"})
        for rank, gene in enumerate(current, 1):
            source = next(row for row in pool_rows if row["family"] == family and row["target"] == target and row["gene_symbol"] == gene)
            final_rows.append({"family": family, "target": target, "gene_symbol": gene, "final_rank": rank, "pool_rank": source["pool_rank"], "fdr_pass_fraction": source["fdr_pass_fraction"], "worst_fdr_ratio": source["worst_fdr_ratio"], "target_support_fraction": source["target_support_fraction"], "worst_effect_margin": source["worst_effect_margin"], "status": status})
        final_score = matrix.loc[current].mean(axis=0)
        for profile in final_score.index:
            profile_rows.append({"family": family, "target": target, "profile_id": profile, "profile_group": groups.loc[profile], "is_target": str(groups.loc[profile] == target).lower(), "positive_rank_score": f"{final_score.loc[profile]:.10g}"})

        production = [gene for gene in production_members[(family, target)] if gene in matrix.index]
        if production:
            production_metrics = metrics(matrix.loc[production].mean(axis=0), groups, target, competitors)
            production_gap = production_metrics["gap"]
            delta = current_metrics["gap"] - production_gap
            comparison = "pilot_higher_optimistic_gap" if delta > 1e-12 else "production_equal_or_higher_gap"
        else:
            production_gap = None; delta = None
            comparison = "pilot_adds_diagnostic_set" if full_pass else "neither_has_passing_set"
        comparison_rows.append({"family": family, "target": target, "production_genes": len(production), "production_gap": f"{production_gap:.10g}" if production_gap is not None else "", "pilot_genes": len(current), "pilot_gap": f"{current_metrics['gap']:.10g}", "pilot_minus_production_gap": f"{delta:.10g}" if delta is not None else "", "pilot_full_data_pass": str(full_pass).lower(), "pilot_unit_support_pass": str(support_pass).lower(), "comparison": comparison})

    # Score the three broad validation-only populations that were excluded from discovery.
    heldout_meta = metadata[metadata["role"] == "held_out_validation"].copy()
    heldout_populations = list(dict.fromkeys(heldout_meta["population"]))
    heldout_counts = pd.DataFrame({population: counts[list(heldout_meta.loc[heldout_meta["population"] == population, "library_id"])].sum(axis=1) for population in heldout_populations})
    heldout_ranks = heldout_counts.rank(axis=0, method="average", pct=True, ascending=True)
    expected = heldout_meta.drop_duplicates("population").set_index("population")["table_s1_lineage"].to_dict()
    production_heldout = {row["population"]: row for row in read_tsv(inputs["production_heldout_broad_predictions"])}
    broad_targets = [target for family, target in target_order if family == "broad_lineage" and final_by_target[(family, target)]]
    heldout_rows = []
    for population in heldout_populations:
        scores = {target: float(heldout_ranks.loc[final_by_target[("broad_lineage", target)], population].mean()) for target in broad_targets}
        predicted = max(scores, key=lambda target: (scores[target], target))
        prior = production_heldout[population]
        heldout_rows.append({
            "population": population, "expected_lineage": expected[population],
            "pilot_predicted_lineage": predicted, "pilot_correct": str(predicted == expected[population]).lower(),
            "pilot_winning_score": f"{scores[predicted]:.10g}", "pilot_scores_json": json.dumps(scores, sort_keys=True, separators=(",", ":")),
            "production_predicted_lineage": prior["predicted_lineage"], "production_correct": prior["correct"].lower(),
        })

    write_tsv(tables["refinement_summary"], summary_rows, list(summary_rows[0]))
    write_tsv(tables["final_membership"], final_rows, list(final_rows[0]), compressed=True)
    write_tsv(tables["removal_trace"], trace_rows, list(trace_rows[0]) if trace_rows else ["family", "target", "iteration", "removed_gene"])
    write_tsv(tables["profile_scores"], profile_rows, list(profile_rows[0]), compressed=True)
    write_tsv(tables["production_comparison"], comparison_rows, list(comparison_rows[0]))
    write_tsv(tables["heldout_broad_predictions"], heldout_rows, list(heldout_rows[0]))

    status_counts = Counter(row["screen_status"] for row in summary_rows)
    summary = {
        "schema_version": "bdev.immgen_positive_competitor_refinement_pilot_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "targets": len(summary_rows),
        "positive_only": True,
        "screen_status_counts": dict(sorted(status_counts.items())),
        "pilot_full_data_positive_separation": sum(row["pilot_full_data_pass"] == "true" for row in comparison_rows),
        "pilot_unit_support_pass": sum(row["pilot_unit_support_pass"] == "true" for row in comparison_rows),
        "pilot_additional_diagnostic_sets": sum(row["comparison"] == "pilot_adds_diagnostic_set" for row in comparison_rows),
        "heldout_broad_correct": sum(row["pilot_correct"] == "true" for row in heldout_rows),
        "heldout_broad_total": len(heldout_rows),
        "production_heldout_broad_correct": sum(row["production_correct"] == "true" for row in heldout_rows),
        "production_contract_changed": False,
        "accepted_sets_emitted": 0,
        "importable_artifact_created": False,
        "scworkbench_invoked": False,
    }
    reports["pilot_summary"].parent.mkdir(parents=True, exist_ok=True)
    reports["pilot_summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    generated = {**tables, **reports}
    manifest = {
        "schema_version": "bdev.immgen_positive_competitor_refinement_pilot_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in {"pilot_config": args.config.resolve(), **inputs}.items()},
        "scripts": {Path(__file__).name: sha256(Path(__file__).resolve())},
        "outputs": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in generated.items()},
        "positive_only": True,
        "production_contract_changed": False,
        "accepted_sets_emitted": 0,
        "importable_artifact_created": False,
    }
    manifest_path = output_root / config["outputs"]["manifest"]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
