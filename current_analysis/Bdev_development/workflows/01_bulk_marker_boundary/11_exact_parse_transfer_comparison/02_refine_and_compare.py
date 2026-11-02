#!/usr/bin/env python3
"""Refine the three-fold Parse transfer and compare all ImmGen set methods."""

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
    require(config["schema_version"] == "bdev.immgen_exact_parse_transfer_comparison.v1", "Unexpected config schema")
    config_dir = args.config.parent.resolve()
    inputs = {name: resolve(config_dir, spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        observed = sha256(path)
        require(observed == config["inputs"][name]["sha256"], f"Checksum mismatch for {name}: {observed}")

    output_root = Path(config["outputs"]["root"])
    tables = {name: output_root / rel for name, rel in config["outputs"]["tables"].items()}
    reports = {name: output_root / rel for name, rel in config["outputs"]["reports"].items()}
    for required in ("fold_plan", "fold_rankings", "candidate_membership", "model_audit"):
        require(tables[required].is_file(), f"Missing model output: {tables[required]}")

    counts = pd.read_csv(inputs["selected_counts"], sep="\t", index_col=0)
    metadata = pd.read_csv(inputs["metadata"], sep="\t", dtype=str).fillna("")
    require(list(counts.columns) == list(metadata["library_id"]), "Count columns and metadata are misaligned")
    max_rank = int(config["method"]["scoring"]["max_rank"])
    # Descending within-library ranks, capped exactly as in the Parse UCell cache.
    capped_ranks = counts.rank(axis=0, method="average", ascending=False).clip(upper=max_rank)
    discovery_meta = metadata[metadata["role"] == "discovery"].copy()
    discovery_ids = list(discovery_meta["library_id"])
    heldout_meta = metadata[metadata["role"] == "held_out_validation"].copy()
    approved_targets = defaultdict(list)
    for row in read_tsv(inputs["contrast_plan"]):
        if row["target"] not in approved_targets[row["family"]]:
            approved_targets[row["family"]].append(row["target"])

    group_column = {
        "superlineage": "superlineage",
        "broad_lineage": "table_s1_lineage",
        "b_branch": "b_branch",
        "b_local_stage": "population",
    }

    def family_metadata(family: str, include_heldout: bool = False) -> pd.DataFrame:
        local = metadata.copy() if include_heldout else discovery_meta.copy()
        if family in {"b_branch", "b_local_stage"}:
            local = local[local["table_s1_lineage"] == "B"]
        if family == "b_local_stage":
            local = local[local["population"].isin(approved_targets[family])]
        return local

    def score_from_rank_sum(rank_sum: pd.Series, n_genes: int) -> pd.Series:
        baseline = n_genes * (n_genes + 1) / 2
        denominator = n_genes * max_rank - baseline
        require(denominator > 0, f"Invalid UCell denominator for {n_genes} genes")
        return (1 - (rank_sum - baseline) / denominator).clip(lower=0, upper=1)

    def gene_set_score(genes: list[str], library_ids: list[str]) -> pd.Series:
        matched = [gene for gene in genes if gene in capped_ranks.index]
        require(len(matched) == len(genes), "A membership contains genes absent from the frozen count matrix")
        return score_from_rank_sum(capped_ranks.loc[matched, library_ids].sum(axis=0), len(matched))

    def unit_scores_and_groups(family: str, library_scores: pd.Series, md: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
        md = md.set_index("library_id").loc[list(library_scores.index)]
        if family == "b_local_stage":
            return library_scores, md["population"]
        units = list(dict.fromkeys(md["population"]))
        scores = pd.Series({unit: float(library_scores.loc[md.index[md["population"] == unit]].mean()) for unit in units})
        unit_md = md.drop_duplicates("population").set_index("population").loc[units]
        return scores, unit_md[group_column[family]]

    def metrics_from_library_scores(family: str, target: str, library_scores: pd.Series, md: pd.DataFrame) -> dict:
        unit_scores, groups = unit_scores_and_groups(family, library_scores, md)
        group_order = list(dict.fromkeys(groups))
        require(target in group_order, f"Target {family}/{target} absent from score groups")
        group_values = {group: float(np.median(unit_scores.loc[groups.index[groups == group]])) for group in group_order}
        competitors = {group: score for group, score in group_values.items() if group != target}
        strongest = max(competitors, key=lambda group: (competitors[group], group))
        sd = float(np.std(list(group_values.values()), ddof=1))
        gap = (group_values[target] - competitors[strongest]) / sd if sd > 0 else float("-inf")
        target_unit_scores = unit_scores.loc[groups.index[groups == target]]
        return {
            "gap": gap,
            "target_value": group_values[target],
            "strongest_competitor": strongest,
            "strongest_value": competitors[strongest],
            "target_rank": 1 + sum(value >= group_values[target] for value in competitors.values()),
            "target_units": len(target_unit_scores),
            "target_unit_win_fraction": float(np.mean(target_unit_scores > competitors[strongest])),
        }

    fold_plan = read_tsv(tables["fold_plan"])
    target_order: list[tuple[str, str]] = []
    target_units: dict[tuple[str, str], int] = {}
    exact_eligible: dict[tuple[str, str], bool] = {}
    for row in fold_plan:
        key = (row["family"], row["target"])
        if key not in target_units:
            target_order.append(key)
            target_units[key] = int(row["target_units"])
            exact_eligible[key] = as_bool(row["construction_eligible"])
    require(len(target_order) == 27, f"Expected 27 targets, found {len(target_order)}")

    candidate_rows = read_tsv(tables["candidate_membership"])
    candidate_sets: dict[tuple[str, str], list[str]] = defaultdict(list)
    candidate_ranks: dict[tuple[str, str], dict[str, int]] = defaultdict(dict)
    construction_pass: dict[tuple[str, str], bool] = defaultdict(bool)
    for row in candidate_rows:
        key = (row["family"], row["target"])
        candidate_sets[key].append(row["gene_symbol"])
        candidate_ranks[key][row["gene_symbol"]] = int(row["candidate_rank"])
        construction_pass[key] = as_bool(row["construction_pass"])

    final_sets: dict[tuple[str, str], list[str]] = {}
    final_metrics: dict[tuple[str, str], dict] = {}
    final_rows: list[dict] = []
    trace_rows: list[dict] = []
    refinement_floor = int(config["method"]["refinement"]["floor"])
    for family, target in target_order:
        key = (family, target)
        if not exact_eligible[key] or not construction_pass[key]:
            final_sets[key] = []
            continue
        current = list(candidate_sets[key])
        md = family_metadata(family)
        ids = list(md["library_id"])
        rank_sum = capped_ranks.loc[current, ids].sum(axis=0)
        current_library_scores = score_from_rank_sum(rank_sum, len(current))
        current_metrics = metrics_from_library_scores(family, target, current_library_scores, md)
        iteration = 0
        while len(current) > refinement_floor:
            choices = []
            for gene in current:
                candidate_rank_sum = rank_sum - capped_ranks.loc[gene, ids]
                candidate_scores = score_from_rank_sum(candidate_rank_sum, len(current) - 1)
                candidate_metrics = metrics_from_library_scores(family, target, candidate_scores, md)
                choices.append((candidate_metrics["gap"], -candidate_ranks[key][gene], gene, candidate_metrics, candidate_rank_sum))
            best = max(choices, key=lambda item: (item[0], item[1]))
            improves = np.isfinite(best[0]) and best[0] > current_metrics["gap"]
            iteration += 1
            for gap, _, gene, test_metrics, _ in choices:
                trace_rows.append({
                    "family": family,
                    "target": target,
                    "iteration": iteration,
                    "genes_before": len(current),
                    "tested_deletion": gene,
                    "candidate_rank": candidate_ranks[key][gene],
                    "gap_before": f"{current_metrics['gap']:.12g}",
                    "gap_after": f"{gap:.12g}",
                    "improvement": f"{gap - current_metrics['gap']:.12g}",
                    "strongest_competitor_after": test_metrics["strongest_competitor"],
                    "selected_deletion": str(improves and gene == best[2]).lower(),
                })
            if not improves:
                break
            current.remove(best[2])
            current_metrics = best[3]
            rank_sum = best[4]
        final_sets[key] = current
        final_metrics[key] = current_metrics
        for final_rank, gene in enumerate(current, 1):
            final_rows.append({
                "family": family,
                "target": target,
                "gene_symbol": gene,
                "final_rank": final_rank,
                "candidate_rank": candidate_ranks[key][gene],
                "candidate_genes": len(candidate_sets[key]),
                "final_genes": len(current),
                "final_gap": f"{current_metrics['gap']:.12g}",
                "strongest_competitor": current_metrics["strongest_competitor"],
                "mean_test_pass": str(current_metrics["target_rank"] == 1).lower(),
                "status": "diagnostic_exact_parse_transfer",
            })

    write_tsv(tables["final_membership"], final_rows, list(final_rows[0]), compressed=True)
    write_tsv(tables["refinement_trace"], trace_rows, list(trace_rows[0]), compressed=True)

    memberships: dict[str, dict[tuple[str, str], list[str]]] = {
        "strict_production": defaultdict(list),
        "positive_full_data_pilot": defaultdict(list),
        "positive_nested_stable": defaultdict(list),
        "exact_parse_transfer": final_sets,
    }
    for row in read_tsv(inputs["production_membership"]):
        memberships["strict_production"][(row["family"], row["target"])].append(row["gene_symbol"])
    for row in read_tsv(inputs["positive_pilot_membership"]):
        memberships["positive_full_data_pilot"][(row["family"], row["target"])].append(row["gene_symbol"])
    nested_status = {(row["family"], row["target"]): row["method_screen_status"] for row in read_tsv(inputs["nested_target_summary"])}
    for row in read_tsv(inputs["nested_stable_membership"]):
        key = (row["family"], row["target"])
        if as_bool(row["stable_at_80pct"]) and nested_status.get(key) == "method_screen_pass":
            memberships["positive_nested_stable"][key].append(row["gene_symbol"])

    comparison_rows: list[dict] = []
    for method, method_sets in memberships.items():
        for family, target in target_order:
            key = (family, target)
            genes = method_sets.get(key, [])
            status = "set_available" if genes else "no_set"
            if method == "exact_parse_transfer":
                if not exact_eligible[key]:
                    status = "not_constructible_low_replication"
                elif not construction_pass[key]:
                    status = "threefold_intersection_below_20"
                else:
                    status = "diagnostic_refined_set"
            if genes:
                md = family_metadata(family)
                scores = gene_set_score(genes, list(md["library_id"]))
                result = metrics_from_library_scores(family, target, scores, md)
            else:
                result = None
            comparison_rows.append({
                "method": method,
                "family": family,
                "target": target,
                "independent_target_units": target_units[key],
                "genes": len(genes),
                "status": status,
                "ucell_gap": f"{result['gap']:.12g}" if result else "",
                "strongest_competitor": result["strongest_competitor"] if result else "",
                "target_rank": result["target_rank"] if result else "",
                "target_unit_win_fraction": f"{result['target_unit_win_fraction']:.8f}" if result else "",
                "full_data_mean_pass": str(result["target_rank"] == 1).lower() if result else "false",
                "selection_uses_full_discovery_scores": str(method in {"positive_full_data_pilot", "exact_parse_transfer"}).lower(),
            })
    write_tsv(tables["method_target_comparison"], comparison_rows, list(comparison_rows[0]))

    heldout_rows: list[dict] = []
    heldout_populations = list(dict.fromkeys(heldout_meta["population"]))
    expected = heldout_meta.drop_duplicates("population").set_index("population")["table_s1_lineage"].to_dict()
    for method, method_sets in memberships.items():
        broad_sets = {target: genes for (family, target), genes in method_sets.items() if family == "broad_lineage" and genes}
        for population in heldout_populations:
            local_md = heldout_meta[heldout_meta["population"] == population]
            ids = list(local_md["library_id"])
            scores = {target: float(gene_set_score(genes, ids).mean()) for target, genes in broad_sets.items()}
            require(scores, f"No broad sets for {method}")
            predicted = max(scores, key=lambda target: (scores[target], target))
            heldout_rows.append({
                "method": method,
                "population": population,
                "expected_lineage": expected[population],
                "predicted_lineage": predicted,
                "correct": str(predicted == expected[population]).lower(),
                "winning_score": f"{scores[predicted]:.12g}",
                "broad_panel_sets": len(broad_sets),
                "scores_json": json.dumps(scores, sort_keys=True, separators=(",", ":")),
            })
    write_tsv(tables["heldout_broad_comparison"], heldout_rows, list(heldout_rows[0]))

    method_rows: list[dict] = []
    for method in memberships:
        local = [row for row in comparison_rows if row["method"] == method]
        heldout = [row for row in heldout_rows if row["method"] == method]
        method_rows.append({
            "method": method,
            "targets": len(local),
            "targets_with_set": sum(row["genes"] > 0 for row in local),
            "full_data_mean_passes": sum(row["full_data_mean_pass"] == "true" for row in local),
            "heldout_broad_correct": sum(row["correct"] == "true" for row in heldout),
            "heldout_broad_total": len(heldout),
            "broad_panel_sets": heldout[0]["broad_panel_sets"],
            "single_cell_expression_used": "false",
        })
    write_tsv(tables["method_summary"], method_rows, list(method_rows[0]))

    parse_status_counts = Counter(row["status"] for row in comparison_rows if row["method"] == "exact_parse_transfer")
    parse_summary = next(row for row in method_rows if row["method"] == "exact_parse_transfer")
    summary = {
        "schema_version": "bdev.immgen_exact_parse_transfer_comparison_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "methods_compared": list(memberships),
        "exact_parse_transfer": parse_summary,
        "exact_parse_status_counts": dict(sorted(parse_status_counts.items())),
        "responder_route_run": False,
        "responder_route_reason": config["method"]["unavailable_parse_route"],
        "single_cell_expression_read": False,
        "accepted_sets_emitted": 0,
        "importable_artifact_created": False,
        "scworkbench_invoked": False,
    }
    reports["comparison_summary"].parent.mkdir(parents=True, exist_ok=True)
    reports["comparison_summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    generated = {**tables, **reports}
    script_dir = Path(__file__).resolve().parent
    scripts = [script_dir / "01_fit_threefold_rankings.R", Path(__file__).resolve()]
    manifest = {
        "schema_version": "bdev.immgen_exact_parse_transfer_comparison_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in {"comparison_config": args.config.resolve(), **inputs}.items()},
        "scripts": {path.name: {"path": str(path), "sha256": sha256(path)} for path in scripts},
        "outputs": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in generated.items()},
        "archived_provenance_only": [str(inputs["archived_parse_method"]), str(inputs["archived_parse_implementation"])],
        "single_cell_expression_read": False,
        "accepted_sets_emitted": 0,
        "importable_artifact_created": False,
    }
    manifest_path = output_root / config["outputs"]["manifest"]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
