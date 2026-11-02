#!/usr/bin/env python3
"""Diagnose target-by-target attrition in the completed ImmGen marker run."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

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


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def as_bool(value: str) -> bool:
    return value.lower() == "true"


def number(value: str, default: float) -> float:
    return float(value) if value else default


def gate_spec(family: str, threshold: dict) -> list[tuple[str, str]]:
    gates = [
        ("expression_filter", "all required contrasts pass filterByExpr"),
        ("annotation_identity", "core_identity_eligible is true"),
    ]
    if family == "b_local_stage":
        gates.extend(
            [
                ("target_support", f"target library expression fraction >= {threshold['minimum_target_library_fraction']}"),
                ("competitor_effect", f"every neighbor log2 fold change >= {threshold['minimum_each_neighbor_log2_fold_change']}"),
                ("competitor_fdr", f"every neighbor FDR <= {threshold['maximum_each_neighbor_fdr']}"),
            ]
        )
    else:
        gates.extend(
            [
                ("target_support", f"target population support >= {threshold['minimum_target_population_support']}"),
                ("global_effect", f"global log2 fold change >= {threshold['minimum_global_log2_fold_change']}"),
                ("global_fdr", f"global FDR <= {threshold['maximum_global_fdr']}"),
            ]
        )
        if family in {"broad_lineage", "b_branch"}:
            gates.extend(
                [
                    ("competitor_effect", f"every competitor log2 fold change >= {threshold['minimum_each_competitor_log2_fold_change']}"),
                    ("competitor_fdr", f"every competitor FDR <= {threshold['maximum_each_competitor_fdr']}"),
                ]
            )
    return gates


def gate_values(family: str, state: dict, threshold: dict) -> dict[str, bool]:
    values = {
        "expression_filter": state["filter_pass"],
        "annotation_identity": state["annotation_eligible"],
    }
    if family == "b_local_stage":
        values.update(
            {
                "target_support": state["expression"] >= float(threshold["minimum_target_library_fraction"]),
                "competitor_effect": state["competitor_lfc"] >= float(threshold["minimum_each_neighbor_log2_fold_change"]),
                "competitor_fdr": state["competitor_fdr"] <= float(threshold["maximum_each_neighbor_fdr"]),
            }
        )
    else:
        values.update(
            {
                "target_support": state["support"] >= float(threshold["minimum_target_population_support"]),
                "global_effect": state["global_lfc"] >= float(threshold["minimum_global_log2_fold_change"]),
                "global_fdr": state["global_fdr"] <= float(threshold["maximum_global_fdr"]),
            }
        )
        if family in {"broad_lineage", "b_branch"}:
            values.update(
                {
                    "competitor_effect": state["competitor_lfc"] >= float(threshold["minimum_each_competitor_log2_fold_change"]),
                    "competitor_fdr": state["competitor_fdr"] <= float(threshold["maximum_each_competitor_fdr"]),
                }
            )
    return values


def initialize_state() -> dict:
    return {
        "filter_pass": True,
        "annotation_eligible": True,
        "support": 1.0,
        "expression": 1.0,
        "global_lfc": float("inf"),
        "global_fdr": 0.0,
        "competitor_lfc": float("inf"),
        "competitor_fdr": 0.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    require(config["schema_version"] == "bdev.immgen_marker_attrition_diagnostic.v1", "Unexpected config schema")
    config_dir = args.config.parent.resolve()
    inputs = {name: resolve(config_dir, spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        observed = sha256(path)
        require(observed == config["inputs"][name]["sha256"], f"Input checksum mismatch for {name}: {observed}")

    contract = yaml.safe_load(inputs["statistical_contract"].read_text())
    production_config = yaml.safe_load(inputs["production_execution_config"].read_text())
    require(production_config["status"] == "execution_complete_candidates_require_scientific_review", "Production execution is not complete")
    output_root = Path(config["outputs"]["root"])
    tables = {name: output_root / rel for name, rel in config["outputs"]["tables"].items()}
    reports = {name: output_root / rel for name, rel in config["outputs"]["reports"].items()}

    preliminary = read_tsv(inputs["preliminary_candidates"])
    preliminary_genes: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in preliminary:
        preliminary_genes[(row["family"], row["target"])].add(row["gene_symbol"])
    set_summary = {(row["family"], row["target"]): row for row in read_tsv(inputs["set_summary"])}

    funnel_rows: list[dict] = []
    failure_rows: list[dict] = []
    contrast_counts: dict[tuple[str, str, str, str], Counter] = defaultdict(Counter)
    target_primary_counts: dict[tuple[str, str], int] = {}
    target_largest_loss: dict[tuple[str, str], tuple[str, int]] = {}
    target_largest_specific_loss: dict[tuple[str, str], tuple[str, int]] = {}
    target_first_below_minimum: dict[tuple[str, str], str] = {}
    target_eliminating_gate: dict[tuple[str, str], str] = {}
    target_order: list[tuple[str, str]] = []
    current_key: tuple[str, str] | None = None
    genes: dict[str, dict] = {}
    seen: set[tuple[str, str]] = set()

    def flush_target(key: tuple[str, str], states: dict[str, dict]) -> None:
        family, target = key
        require(key not in seen, f"Noncontiguous target block: {key}")
        require(len(states) == int(production_config["execution"]["expected_gene_rows"]), f"Unexpected gene count for {key}: {len(states)}")
        threshold = contract["contrast_families"][family]
        specs = gate_spec(family, threshold)
        values_by_gene = {gene: gate_values(family, state, threshold) for gene, state in states.items()}
        survivors = set(states)
        largest = ("", -1)
        largest_specific = ("", -1)
        first_below_minimum = ""
        eliminating_gate = ""
        for order, (gate, description) in enumerate(specs, 1):
            input_count = len(survivors)
            survivors = {gene for gene in survivors if values_by_gene[gene][gate]}
            lost = input_count - len(survivors)
            if lost > largest[1]:
                largest = (gate, lost)
            if order >= 3 and lost > largest_specific[1]:
                largest_specific = (gate, lost)
            if not first_below_minimum and len(survivors) < int(threshold["minimum_genes"]):
                first_below_minimum = gate
            if not eliminating_gate and not survivors:
                eliminating_gate = gate
            funnel_rows.append(
                {
                    "family": family,
                    "target": target,
                    "gate_order": order,
                    "gate": gate,
                    "threshold_or_rule": description,
                    "source_genes": len(states),
                    "sequential_input_genes": input_count,
                    "sequential_pass_genes": len(survivors),
                    "lost_at_gate": lost,
                    "independent_pass_genes": sum(value[gate] for value in values_by_gene.values()),
                    "cumulative_pass_fraction": f"{len(survivors) / len(states):.8f}",
                }
            )
        first_failures = Counter()
        for value in values_by_gene.values():
            failure = "primary_eligible"
            for gate, _ in specs:
                if not value[gate]:
                    failure = gate
                    break
            first_failures[failure] += 1
        for order, (gate, _) in enumerate(specs, 1):
            failure_rows.append({"family": family, "target": target, "reason_order": order, "first_failure_reason": gate, "genes": first_failures[gate]})
        failure_rows.append({"family": family, "target": target, "reason_order": len(specs) + 1, "first_failure_reason": "primary_eligible", "genes": first_failures["primary_eligible"]})
        expected_genes = preliminary_genes.get(key, set())
        require(survivors == expected_genes, f"Primary eligibility mismatch for {family}/{target}: diagnostic={len(survivors)}, production={len(expected_genes)}")
        target_primary_counts[key] = len(survivors)
        target_largest_loss[key] = largest
        target_largest_specific_loss[key] = largest_specific
        target_first_below_minimum[key] = first_below_minimum
        target_eliminating_gate[key] = eliminating_gate
        target_order.append(key)
        seen.add(key)

    with gzip.open(inputs["complete_statistics"], "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["family"], row["target"])
            if current_key is not None and key != current_key:
                flush_target(current_key, genes)
                genes = {}
            current_key = key
            gene = row["gene_symbol"]
            state = genes.setdefault(gene, initialize_state())
            state["filter_pass"] &= as_bool(row["filter_pass"])
            state["annotation_eligible"] &= as_bool(row["core_identity_eligible"])
            state["support"] = min(state["support"], number(row["target_support_fraction"], 0.0))
            state["expression"] = min(state["expression"], number(row["target_expression_fraction"], 0.0))
            lfc = number(row["logFC"], float("-inf"))
            fdr = number(row["FDR"], 1.0)
            if row["test_type"] == "global":
                state["global_lfc"] = min(state["global_lfc"], lfc)
                state["global_fdr"] = max(state["global_fdr"], fdr)
            else:
                state["competitor_lfc"] = min(state["competitor_lfc"], lfc)
                state["competitor_fdr"] = max(state["competitor_fdr"], fdr)

            threshold = contract["contrast_families"][row["family"]]
            if row["family"] == "b_local_stage":
                support_ok = number(row["target_expression_fraction"], 0.0) >= float(threshold["minimum_target_library_fraction"])
                effect_cutoff = float(threshold["minimum_each_neighbor_log2_fold_change"])
                fdr_cutoff = float(threshold["maximum_each_neighbor_fdr"])
            else:
                support_ok = number(row["target_support_fraction"], 0.0) >= float(threshold["minimum_target_population_support"])
                if row["test_type"] == "global":
                    effect_cutoff = float(threshold["minimum_global_log2_fold_change"])
                    fdr_cutoff = float(threshold["maximum_global_fdr"])
                else:
                    effect_cutoff = float(threshold["minimum_each_competitor_log2_fold_change"])
                    fdr_cutoff = float(threshold["maximum_each_competitor_fdr"])
            counter = contrast_counts[(row["family"], row["target"], row["test_id"], row["test_type"])]
            counter["source_genes"] += 1
            filter_ok = as_bool(row["filter_pass"])
            annotation_ok = as_bool(row["core_identity_eligible"])
            counter["filter_pass_genes"] += filter_ok
            counter["annotation_eligible_genes"] += annotation_ok
            counter["support_pass_genes"] += support_ok
            universe = filter_ok and annotation_ok and support_ok
            counter["analysis_universe_genes"] += universe
            counter["effect_pass_in_universe"] += universe and lfc >= effect_cutoff
            counter["fdr_pass_in_universe"] += universe and fdr <= fdr_cutoff
            counter["effect_and_fdr_pass_in_universe"] += universe and lfc >= effect_cutoff and fdr <= fdr_cutoff
        if current_key is not None:
            flush_target(current_key, genes)

    contrast_rows = []
    for key, counts in contrast_counts.items():
        family, target, test_id, test_type = key
        universe = counts["analysis_universe_genes"]
        passed = counts["effect_and_fdr_pass_in_universe"]
        contrast_rows.append(
            {
                "family": family,
                "target": target,
                "test_id": test_id,
                "test_type": test_type,
                **counts,
                "effect_or_fdr_lost_from_universe": universe - passed,
                "effect_and_fdr_pass_fraction": f"{passed / universe:.8f}" if universe else "",
            }
        )

    balanced = {(row["family"], row["target"], row["gene_symbol"]): float(row["stability_fraction"]) for row in read_tsv(inputs["balanced_stability"])}
    loo_rows = read_tsv(inputs["population_loo"]) + read_tsv(inputs["library_loo"])
    loo = {(row["family"], row["target"], row["gene_symbol"]): row for row in loo_rows}
    post_rows = []
    diagnosis_rows = []
    for family, target in target_order:
        key = (family, target)
        genes_for_target = preliminary_genes.get(key, set())
        threshold = contract["contrast_families"][family]
        stable = {gene for gene in genes_for_target if balanced[(family, target, gene)] >= float(threshold["minimum_stability_fraction"])}
        assessed = {gene for gene in stable if loo[(family, target, gene)]["status"] == "assessed"}
        nonassessable = stable - assessed
        loo_pass = {gene for gene in assessed if float(loo[(family, target, gene)]["holdout_fraction"]) >= 0.80}
        post_gate = loo_pass | nonassessable
        minimum = int(threshold["minimum_genes"])
        maximum = int(threshold["maximum_genes"])
        expected_emitted = min(len(post_gate), maximum) if len(post_gate) >= minimum else 0
        production = set_summary[key]
        require(int(production["preliminary_eligible_genes"]) == len(genes_for_target), f"Preliminary summary mismatch for {key}")
        require(int(production["emitted_candidate_genes"]) == expected_emitted, f"Emitted summary mismatch for {key}: diagnostic={expected_emitted}, production={production['emitted_candidate_genes']}")
        post_rows.append(
            {
                "family": family,
                "target": target,
                "primary_eligible_genes": len(genes_for_target),
                "balanced_stability_pass_genes": len(stable),
                "lost_at_balanced_stability": len(genes_for_target) - len(stable),
                "loo_assessed_after_stability": len(assessed),
                "loo_not_assessable_low_replication": len(nonassessable),
                "loo_pass_after_stability": len(loo_pass),
                "lost_at_loo_after_stability": len(assessed) - len(loo_pass),
                "post_stability_and_loo_genes": len(post_gate),
                "minimum_genes": minimum,
                "maximum_genes": maximum,
                "lost_to_minimum_size_rule": len(post_gate) if len(post_gate) < minimum else 0,
                "truncated_at_maximum": max(0, len(post_gate) - maximum) if len(post_gate) >= minimum else 0,
                "emitted_candidate_genes": expected_emitted,
                "production_status": production["status"],
            }
        )
        largest_gate, largest_loss = target_largest_loss[key]
        largest_specific_gate, largest_specific_loss = target_largest_specific_loss[key]
        first_below_minimum = target_first_below_minimum[key]
        eliminating_gate = target_eliminating_gate[key]
        if not genes_for_target:
            phase = "no_primary_candidates"
            focus = f"eliminating primary gate: {eliminating_gate}"
        elif len(genes_for_target) < minimum:
            phase = "primary_candidates_below_minimum"
            focus = f"first primary gate below minimum: {first_below_minimum}"
        elif len(stable) < minimum:
            phase = "balanced_stability_below_minimum"
            focus = "balanced stability"
        elif len(post_gate) < minimum:
            phase = "leave_one_out_below_minimum"
            focus = "leave-one-out stability"
        elif nonassessable:
            phase = "size_valid_low_replication_candidate"
            focus = "replication limitation"
        else:
            phase = "size_valid_candidate"
            focus = "scientific review and non-statistical gates"
        diagnosis_rows.append(
            {
                "family": family,
                "target": target,
                "diagnostic_phase": phase,
                "largest_primary_sequential_loss_gate": largest_gate,
                "genes_lost_at_largest_primary_gate": largest_loss,
                "largest_target_specific_primary_loss_gate": largest_specific_gate,
                "genes_lost_at_largest_target_specific_gate": largest_specific_loss,
                "first_primary_gate_below_minimum": first_below_minimum,
                "eliminating_primary_gate": eliminating_gate,
                "review_focus": focus,
                "production_status": production["status"],
                "production_review_reasons": production["review_reasons"],
                "contract_change_authorized": "false",
            }
        )

    write_tsv(tables["attrition_funnel"], funnel_rows, list(funnel_rows[0]))
    write_tsv(tables["first_failure_reasons"], failure_rows, list(failure_rows[0]))
    write_tsv(tables["contrast_bottlenecks"], contrast_rows, list(contrast_rows[0]))
    write_tsv(tables["post_primary_attrition"], post_rows, list(post_rows[0]))
    write_tsv(tables["target_diagnosis"], diagnosis_rows, list(diagnosis_rows[0]))

    summary = {
        "schema_version": "bdev.immgen_marker_attrition_diagnostic_summary.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "targets": len(target_order),
        "contrasts": len(contrast_rows),
        "source_gene_rows_per_target": int(production_config["execution"]["expected_gene_rows"]),
        "production_primary_candidate_rows": len(preliminary),
        "primary_counts_reproduced_exactly": True,
        "emitted_counts_reproduced_exactly": True,
        "diagnostic_phase_counts": dict(sorted(Counter(row["diagnostic_phase"] for row in diagnosis_rows).items())),
        "model_refit": False,
        "thresholds_changed": False,
        "importable_artifact_created": False,
        "scworkbench_invoked": False,
    }
    reports["diagnostic_summary"].parent.mkdir(parents=True, exist_ok=True)
    reports["diagnostic_summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    generated = {**tables, **reports}
    manifest = {
        "schema_version": "bdev.immgen_marker_attrition_diagnostic_manifest.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in {"diagnostic_config": args.config.resolve(), **inputs}.items()},
        "scripts": {Path(__file__).name: sha256(Path(__file__).resolve())},
        "outputs": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in generated.items()},
        "model_refit": False,
        "thresholds_changed": False,
        "importable_artifact_created": False,
    }
    manifest_path = output_root / config["outputs"]["manifest"]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
