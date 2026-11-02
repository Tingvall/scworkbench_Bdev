#!/usr/bin/env python3
"""Apply deterministic selection and validate the synthetic marker workflow."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
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


def as_bool(value: str) -> bool:
    return value.lower() == "true"


def as_float(value: str, default: float = float("nan")) -> float:
    return float(value) if value != "" else default


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--workflow-dir", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    config_dir = args.config.parent
    contract_path = (config_dir / config["inputs"]["statistical_contract"]).resolve()
    hierarchy_path = (config_dir / config["inputs"]["hierarchy"]).resolve()
    software_path = (config_dir / config["inputs"]["software_lock"]).resolve()
    contract = yaml.safe_load(contract_path.read_text())
    hierarchy = yaml.safe_load(hierarchy_path.read_text())
    output_root = Path(config["outputs"]["root"])
    require(args.workflow_dir.resolve() == output_root.resolve(), "Workflow directory does not match the frozen output root")

    fixture_paths = {key: output_root / value for key, value in config["outputs"]["fixture"].items()}
    model_paths = {key: output_root / value for key, value in config["outputs"]["model"].items()}
    selection_paths = {key: output_root / value for key, value in config["outputs"]["selection"].items()}
    manifest_path = output_root / config["outputs"]["manifest"]

    with gzip.open(model_paths["complete_statistics"], "rt", newline="") as handle:
        statistics = list(csv.DictReader(handle, delimiter="\t"))
    with fixture_paths["expected_markers"].open(newline="") as handle:
        expected_rows = list(csv.DictReader(handle, delimiter="\t"))
    with fixture_paths["metadata"].open(newline="") as handle:
        metadata = list(csv.DictReader(handle, delimiter="\t"))
    with model_paths["design_audit"].open(newline="") as handle:
        design_audit = list(csv.DictReader(handle, delimiter="\t"))

    require(len(design_audit) == 4, "Expected four model designs")
    for audit in design_audit:
        require(int(audit["design_rank"]) == int(audit["design_columns"]), f"Rank failure: {audit['family']}")
        require(int(audit["residual_df"]) > 0, f"No residual degrees of freedom: {audit['family']}")
        require(int(audit["filtered_out"]) >= config["synthetic_fixture"]["control_genes"], f"Low-count controls did not filter: {audit['family']}")

    grouped: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in statistics:
        grouped[(row["family"], row["target"], row["gene_symbol"])].append(row)
    expected: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in expected_rows:
        expected[(row["family"], row["target"])].add(row["gene_symbol"])

    family_targets = sorted({(family, target) for family, target, _ in grouped})
    selected_rows: list[dict] = []
    target_summaries: list[dict] = []
    selected_by_target: dict[tuple[str, str], set[str]] = defaultdict(set)

    population_counts = defaultdict(set)
    library_counts = defaultdict(int)
    for row in metadata:
        if row["role"] != "discovery":
            continue
        population_counts[("superlineage", row["superlineage"])].add(row["population"])
        population_counts[("broad_lineage", row["table_s1_lineage"])].add(row["population"])
        if row["b_branch"]:
            population_counts[("b_branch", row["b_branch"])].add(row["population"])
        library_counts[row["population"]] += 1
    cautioned_local = set()
    for edge in hierarchy["b_lineage"]["directed_developmental_edges"]:
        if edge["edge_type"] in {"gap_spanning", "compartment_spanning"}:
            cautioned_local.update([edge["source"], edge["target"]])

    for family, target in family_targets:
        threshold = contract["contrast_families"][family]
        candidates = []
        genes = sorted(gene for fam, tar, gene in grouped if fam == family and tar == target)
        for gene in genes:
            rows = grouped[(family, target, gene)]
            global_rows = [row for row in rows if row["test_type"] == "global"]
            required_rows = [row for row in rows if row["test_type"] in {"competitor", "neighbor"}]
            filter_pass = all(as_bool(row["filter_pass"]) for row in rows)
            annotation_eligible = all(as_bool(row["annotation_eligible"]) for row in rows)
            support = min(as_float(row["target_support_fraction"], 0) for row in rows)
            stability = min(as_float(row["stability_fraction"], 0) for row in rows)
            worst_auroc = min(as_float(row["target_versus_competitor_auroc"], 0) for row in rows)
            global_lfc = min((as_float(row["logFC"], float("-inf")) for row in global_rows), default=float("nan"))
            global_fdr = max((as_float(row["FDR"], 1) for row in global_rows), default=float("nan"))
            competitor_lfc = min((as_float(row["logFC"], float("-inf")) for row in required_rows), default=float("nan"))
            competitor_fdr = max((as_float(row["FDR"], 1) for row in required_rows), default=float("nan"))
            required_lfcs = [as_float(row["logFC"], float("-inf")) for row in rows]
            required_fdrs = [as_float(row["FDR"], 1) for row in rows]
            eligible = filter_pass and annotation_eligible
            if family == "superlineage":
                eligible &= support >= threshold["minimum_target_population_support"]
                eligible &= global_lfc >= threshold["minimum_global_log2_fold_change"] and global_fdr <= threshold["maximum_global_fdr"]
            elif family == "broad_lineage":
                eligible &= support >= threshold["minimum_target_population_support"]
                eligible &= global_lfc >= threshold["minimum_global_log2_fold_change"] and global_fdr <= threshold["maximum_global_fdr"]
                eligible &= competitor_lfc >= threshold["minimum_each_competitor_log2_fold_change"] and competitor_fdr <= threshold["maximum_each_competitor_fdr"]
            elif family == "b_branch":
                eligible &= support >= threshold["minimum_target_population_support"]
                eligible &= global_lfc >= threshold["minimum_global_log2_fold_change"] and global_fdr <= threshold["maximum_global_fdr"]
                eligible &= competitor_lfc >= threshold["minimum_each_competitor_log2_fold_change"] and competitor_fdr <= threshold["maximum_each_competitor_fdr"]
            elif family == "b_local_stage":
                eligible &= min(as_float(row["target_expression_fraction"], 0) for row in rows) >= threshold["minimum_target_library_fraction"]
                eligible &= competitor_lfc >= threshold["minimum_each_neighbor_log2_fold_change"] and competitor_fdr <= threshold["maximum_each_neighbor_fdr"]
            eligible &= stability >= threshold["minimum_stability_fraction"]
            if eligible:
                candidates.append(
                    {
                        "family": family,
                        "target": target,
                        "gene_symbol": gene,
                        "stability_fraction": stability,
                        "target_support_fraction": support,
                        "worst_log2_fold_change": min(required_lfcs),
                        "worst_auroc": worst_auroc,
                        "worst_fdr": max(required_fdrs),
                    }
                )

        candidates.sort(
            key=lambda row: (
                -row["stability_fraction"],
                -row["target_support_fraction"],
                -row["worst_log2_fold_change"],
                -row["worst_auroc"],
                row["worst_fdr"],
                row["gene_symbol"],
            )
        )
        emitted = len(candidates) >= threshold["minimum_genes"]
        selected = candidates[: threshold["maximum_genes"]] if emitted else []
        if family in {"superlineage", "broad_lineage", "b_branch"} and len(population_counts[(family, target)]) < contract["stability"]["leave_one_population_out"]["minimum_target_populations"]:
            predicted_class = "candidate_low_replication"
        elif family == "b_local_stage" and library_counts[target] < contract["stability"]["leave_one_library_out"]["minimum_target_libraries"]:
            predicted_class = "candidate_low_replication"
        elif family == "b_local_stage" and target in cautioned_local:
            predicted_class = "candidate_gap_or_compartment"
        else:
            predicted_class = "accepted_core_if_deferred_gates_pass"
        for rank, row in enumerate(selected, 1):
            row.update(
                {
                    "rank": rank,
                    "qualification_status": "synthetic_selected_not_importable",
                    "predicted_confidence_class": predicted_class,
                }
            )
            selected_rows.append(row)
            selected_by_target[(family, target)].add(row["gene_symbol"])
        target_summaries.append(
            {
                "family": family,
                "target": target,
                "eligible_genes": len(candidates),
                "selected_genes": len(selected),
                "minimum_genes": threshold["minimum_genes"],
                "maximum_genes": threshold["maximum_genes"],
                "set_emitted": emitted,
                "predicted_confidence_class": predicted_class,
            }
        )

    recoveries = []
    no_padding = config["synthetic_fixture"]["no_padding_control"]
    no_padding_key = (no_padding["family"], no_padding["target"])
    for key, expected_genes in sorted(expected.items()):
        recovered = selected_by_target[key] & expected_genes
        recovery = len(recovered) / len(expected_genes)
        if key == no_padding_key:
            require(len(recovered) == 0, "Under-minimum no-padding control unexpectedly emitted genes")
        else:
            require(recovery >= 0.95, f"Injected-marker recovery below 95% for {key}: {recovery:.3f}")
        recoveries.append({"family": key[0], "target": key[1], "expected": len(expected_genes), "recovered": len(recovered), "recovery_fraction": recovery})

    no_padding_summary = next(row for row in target_summaries if (row["family"], row["target"]) == no_padding_key)
    require(no_padding_summary["eligible_genes"] < no_padding_summary["minimum_genes"], "No-padding control unexpectedly met the minimum")
    require(not no_padding_summary["set_emitted"] and no_padding_summary["selected_genes"] == 0, "No-padding control emitted a set")
    truncation_targets = [row for row in target_summaries if row["eligible_genes"] > row["maximum_genes"] and row["selected_genes"] == row["maximum_genes"]]
    require(truncation_targets, "No target exercised deterministic maximum-size truncation")

    selected_symbols = {row["gene_symbol"] for row in selected_rows}
    require(not any(symbol.startswith("SYN_STATE_DECOY") for symbol in selected_symbols), "State decoy selected")
    require(not any(symbol.startswith("SYN_AMBIGUOUS_DECOY") for symbol in selected_symbols), "Ambiguous decoy selected")
    require(not any(symbol.startswith("SYN_FILTER_LOW") for symbol in selected_symbols), "Low-count gene selected")
    require(all(row["qualification_status"] != "accepted_core" for row in selected_rows), "Qualification emitted accepted_core")

    selection_fields = [
        "family", "target", "gene_symbol", "rank", "stability_fraction", "target_support_fraction",
        "worst_log2_fold_change", "worst_auroc", "worst_fdr", "qualification_status",
        "predicted_confidence_class",
    ]
    write_tsv(selection_paths["long_form"], selected_rows, selection_fields)
    summary = {
        "schema_version": "bdev.immgen_marker_workflow_qualification_summary.v1",
        "status": "model_and_selection_core_passed_deferred_gates_required",
        "synthetic_only": True,
        "immgen_expression_values_read": False,
        "libraries": len(metadata),
        "discovery_libraries": sum(row["role"] == "discovery" for row in metadata),
        "held_out_libraries": sum(row["role"] == "held_out_validation" for row in metadata),
        "designs": design_audit,
        "targets": target_summaries,
        "injected_marker_recovery": recoveries,
        "selected_gene_rows": len(selected_rows),
        "state_decoys_selected": 0,
        "ambiguous_decoys_selected": 0,
        "low_count_genes_selected": 0,
        "accepted_core_sets_emitted": 0,
        "no_padding_control_passed": True,
        "maximum_size_truncation_targets": [{"family": row["family"], "target": row["target"]} for row in truncation_targets],
        "importable_gmt_created": False,
        "deferred": config["qualification_scope"]["explicitly_deferred"],
        "scworkbench_invoked": False,
    }
    selection_paths["summary"].parent.mkdir(parents=True, exist_ok=True)
    selection_paths["summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    tracked_inputs = {
        "qualification_config": args.config.resolve(),
        "statistical_contract": contract_path,
        "hierarchy": hierarchy_path,
        "software_lock": software_path,
    }
    script_dir = Path(__file__).resolve().parent
    scripts = sorted(script_dir.glob("*.py")) + sorted(script_dir.glob("*.R"))
    generated = {**fixture_paths, **model_paths, **selection_paths}
    manifest = {
        "schema_version": "bdev.immgen_marker_workflow_qualification_manifest.v1",
        "status": "passed_with_declared_deferred_gates",
        "synthetic_only": True,
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in tracked_inputs.items()},
        "scripts": {path.name: sha256(path) for path in scripts},
        "outputs": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in generated.items()},
        "importable_artifact_created": False,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
