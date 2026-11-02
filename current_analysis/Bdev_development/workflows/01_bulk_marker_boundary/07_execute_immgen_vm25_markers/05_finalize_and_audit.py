#!/usr/bin/env python3
"""Finalize review-only ImmGen candidates, audits, summary, and manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from collections import Counter
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


def resolve(config_dir: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (config_dir / path).resolve()


def read_tsv(path: Path, compressed: bool = False) -> list[dict[str, str]]:
    opener = gzip.open if compressed else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_tsv_gz(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def split_labs(value: str) -> set[str]:
    return {item for item in value.split(";") if item}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    require(config["schema_version"] == "bdev.immgen_vm25_marker_execution.v1", "Wrong config schema")
    require(args.output_root.resolve() == Path(config["outputs"]["root"]).resolve(), "Output root differs from configuration")
    config_dir = args.config.parent
    contract = yaml.safe_load(resolve(config_dir, config["inputs"]["statistical_contract"]["path"]).read_text())
    hierarchy = yaml.safe_load(resolve(config_dir, config["inputs"]["hierarchy"]["path"]).read_text())
    tables = {name: args.output_root / rel for name, rel in config["outputs"]["tables"].items()}
    reports = {name: args.output_root / rel for name, rel in config["outputs"]["reports"].items()}

    metadata = read_tsv(tables["metadata"])
    plan = read_tsv(tables["contrast_plan"])
    preliminary = read_tsv(tables["preliminary_candidates"], compressed=True)
    balanced = read_tsv(tables["balanced_stability"], compressed=True)
    population_loo = read_tsv(tables["population_loo"], compressed=True)
    library_loo = read_tsv(tables["library_loo"])
    validation_sets = read_tsv(tables["validation_gene_sets"])
    broad_predictions = read_tsv(tables["broad_lopo_predictions"])
    heldout_predictions = read_tsv(tables["heldout_predictions"])
    local_auroc = read_tsv(tables["local_neighbor_auroc"])

    # Metadata contingency and contrast separation audits.
    contingency_rows: list[dict] = []
    for population in sorted({row["population"] for row in metadata}):
        population_rows = [row for row in metadata if row["population"] == population]
        for factor, column in (("organ", "organ"), ("laboratory", "laboratories")):
            counts: Counter[str] = Counter()
            for row in population_rows:
                values = split_labs(row[column]) if factor == "laboratory" else {row[column]}
                for value in values:
                    counts[value] += 1
            for level in sorted(counts):
                contingency_rows.append({"population": population, "factor": factor, "level": level, "libraries": counts[level]})
    write_tsv(tables["population_factor_contingency"], contingency_rows, ["population", "factor", "level", "libraries"])

    audit_rows: list[dict] = []
    for test in plan:
        rows = [row for row in metadata if row["role"] == "discovery"]
        if test["sample_scope"] in {"b_discovery", "b_local_discovery"}:
            rows = [row for row in rows if row["table_s1_lineage"] == "B"]
        weights = json.loads(test["weights_json"])
        competitors = {group for group, weight in weights.items() if float(weight) < 0}
        target_rows = [row for row in rows if row[test["group_column"]] == test["target"]]
        competitor_rows = [row for row in rows if row[test["group_column"]] in competitors]
        require(target_rows and competitor_rows, f"Empty audit group: {test['test_id']}")
        target_organs = {row["organ"] for row in target_rows}; competitor_organs = {row["organ"] for row in competitor_rows}
        target_labs = set().union(*(split_labs(row["laboratories"]) for row in target_rows)); competitor_labs = set().union(*(split_labs(row["laboratories"]) for row in competitor_rows))
        organ_separated = target_organs.isdisjoint(competitor_organs); laboratory_separated = target_labs.isdisjoint(competitor_labs)
        audit_rows.append({
            "family": test["family"], "target": test["target"], "test_id": test["test_id"], "test_type": test["test_type"],
            "target_organs": ";".join(sorted(target_organs)), "competitor_organs": ";".join(sorted(competitor_organs)), "organ_completely_separated": str(organ_separated).lower(),
            "target_laboratories": ";".join(sorted(target_labs)), "competitor_laboratories": ";".join(sorted(competitor_labs)), "laboratory_completely_separated": str(laboratory_separated).lower(),
            "requires_scientific_review": str(organ_separated or laboratory_separated).lower(),
        })
    audit_fields = ["family", "target", "test_id", "test_type", "target_organs", "competitor_organs", "organ_completely_separated", "target_laboratories", "competitor_laboratories", "laboratory_completely_separated", "requires_scientific_review"]
    write_tsv(tables["contrast_confounding_audit"], audit_rows, audit_fields)

    prelim_count = Counter((row["family"], row["target"]) for row in preliminary)
    emitted_count = Counter((row["family"], row["target"]) for row in validation_sets)
    balanced_map = {(row["family"], row["target"], row["gene_symbol"]): row["stability_fraction"] for row in balanced}
    pop_map = {(row["family"], row["target"], row["gene_symbol"]): row for row in population_loo}
    lib_map = {(row["family"], row["target"], row["gene_symbol"]): row for row in library_loo}
    confounded_target = {(row["family"], row["target"]) for row in audit_rows if row["requires_scientific_review"] == "true"}
    cautioned_local = set()
    for edge in hierarchy["b_lineage"]["directed_developmental_edges"]:
        if edge["edge_type"] in {"gap_spanning", "compartment_spanning"}:
            cautioned_local.update((edge["source"], edge["target"]))

    broad_overall = sum(row["correct"] == "TRUE" for row in broad_predictions) / len(broad_predictions) if broad_predictions else 0.0
    broad_sensitivity: dict[str, float] = {}
    for lineage in sorted({row["expected_lineage"] for row in broad_predictions}):
        rows = [row for row in broad_predictions if row["expected_lineage"] == lineage]
        broad_sensitivity[lineage] = sum(row["correct"] == "TRUE" for row in rows) / len(rows)
    heldout_correct = {row["expected_lineage"]: row["correct"] == "TRUE" for row in heldout_predictions}
    local_by_target: dict[str, list[dict]] = {}
    for row in local_auroc:
        local_by_target.setdefault(row["target"], []).append(row)

    all_targets = []
    all_targets.extend(("superlineage", target) for target in hierarchy["superlineages"])
    all_targets.extend(("broad_lineage", target) for target in hierarchy["broad_lineages"])
    all_targets.extend(("b_branch", target) for target in hierarchy["b_lineage"]["branch_groups"])
    all_targets.extend(("b_local_stage", target) for target in hierarchy["b_lineage"]["local_stage_neighbors"])
    set_summary: list[dict] = []
    status_by_target: dict[tuple[str, str], tuple[str, str]] = {}
    for family, target in all_targets:
        reasons: list[str] = []
        emitted = emitted_count[(family, target)]
        threshold = contract["contrast_families"][family]
        if emitted < int(threshold["minimum_genes"]):
            status = "rejected_insufficient_eligible_genes"
            reasons.append("fewer_than_family_minimum_after_stability_and_holdout_gates")
        else:
            if family == "b_local_stage":
                rows = [row for key, row in lib_map.items() if key[:2] == (family, target)]
            else:
                rows = [row for key, row in pop_map.items() if key[:2] == (family, target)]
            low_replication = bool(rows) and all(row["status"] == "not_assessable_low_replication" for row in rows)
            if low_replication:
                reasons.append("holdout_not_assessable_low_replication")
            if family == "b_local_stage" and target in cautioned_local:
                reasons.append("gap_or_compartment_spanning")
            if (family, target) in confounded_target:
                reasons.append("one_or_more_required_contrasts_completely_separated_by_organ_or_laboratory")
            validation_failed = False
            if family == "broad_lineage":
                validation_failed = broad_overall < float(contract["internal_validation"]["broad_leave_one_population_out"]["minimum_overall_top1_accuracy"])
                validation_failed |= broad_sensitivity.get(target, 0.0) < float(contract["internal_validation"]["broad_leave_one_population_out"]["minimum_per_lineage_sensitivity"])
                if target in heldout_correct:
                    validation_failed |= not heldout_correct[target]
            if family == "b_local_stage":
                assessed = [float(row["auroc"]) for row in local_by_target.get(target, []) if row["status"] == "assessed" and row["auroc"]]
                validation_failed |= bool(assessed) and min(assessed) < float(contract["internal_validation"]["local_neighbor_validation"]["minimum_pairwise_auroc"])
                if any(row["status"] == "insufficient_replication" for row in local_by_target.get(target, [])):
                    reasons.append("local_validation_insufficient_replication")
            if validation_failed:
                reasons.append("internal_validation_threshold_failed")
            if low_replication:
                status = "candidate_low_replication"
            elif family == "b_local_stage" and target in cautioned_local:
                status = "candidate_gap_or_compartment"
            elif (family, target) in confounded_target:
                status = "candidate_confounding"
            elif validation_failed:
                status = "candidate_internal_validation"
            else:
                status = "candidate_review_ready"
                reasons.append("all_computational_gates_passed_scientific_review_required")
        status_by_target[(family, target)] = (status, ";".join(reasons))
        set_summary.append({
            "family": family, "target": target, "preliminary_eligible_genes": prelim_count[(family, target)], "emitted_candidate_genes": emitted,
            "minimum_genes": threshold["minimum_genes"], "maximum_genes": threshold["maximum_genes"], "status": status, "review_reasons": ";".join(reasons),
        })
    write_tsv(tables["set_summary"], set_summary, ["family", "target", "preliminary_eligible_genes", "emitted_candidate_genes", "minimum_genes", "maximum_genes", "status", "review_reasons"])

    candidate_rows: list[dict] = []
    for row in validation_sets:
        key = (row["family"], row["target"], row["gene_symbol"])
        gate = lib_map.get(key) if row["family"] == "b_local_stage" else pop_map.get(key)
        status, reasons = status_by_target[(row["family"], row["target"])]
        candidate_rows.append({
            **row,
            "balanced_stability_fraction": balanced_map[key],
            "holdout_status": gate["status"] if gate else "",
            "holdout_fraction": gate["holdout_fraction"] if gate else "",
            "candidate_status": status,
            "review_reasons": reasons,
        })
    candidate_fields = list(candidate_rows[0]) if candidate_rows else ["family", "target", "gene_symbol", "candidate_status", "review_reasons"]
    write_tsv_gz(tables["candidate_genesets"], candidate_rows, candidate_fields)

    status_counts = Counter(row["status"] for row in set_summary)
    summary = {
        "schema_version": "bdev.immgen_vm25_execution_summary.v1",
        "status": "execution_complete_candidates_require_scientific_review",
        "run_identity": config["run_identity"],
        "immgen_expression_values_read": True,
        "genes": config["execution"]["expected_gene_rows"],
        "libraries": len(metadata),
        "contrasts": len(plan),
        "preliminary_candidate_rows": len(preliminary),
        "candidate_gene_rows": len(candidate_rows),
        "set_status_counts": dict(sorted(status_counts.items())),
        "broad_lopo_overall_accuracy": broad_overall,
        "broad_lopo_lineage_sensitivity": broad_sensitivity,
        "heldout_predictions": heldout_predictions,
        "local_neighbor_results": len(local_auroc),
        "organ_completely_separated_contrasts": sum(row["organ_completely_separated"] == "true" for row in audit_rows),
        "laboratory_completely_separated_contrasts": sum(row["laboratory_completely_separated"] == "true" for row in audit_rows),
        "accepted_core_sets_emitted": 0,
        "importable_gmt_created": False,
        "scworkbench_invoked": False,
    }
    reports["execution_summary"].parent.mkdir(parents=True, exist_ok=True)
    reports["execution_summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    tracked_inputs = {"execution_config": args.config.resolve()}
    for name, spec in config["inputs"].items():
        if isinstance(spec, dict) and "path" in spec:
            tracked_inputs[name] = resolve(config_dir, spec["path"])
    script_dir = Path(__file__).resolve().parent
    scripts = sorted(script_dir.glob("*.py")) + sorted(script_dir.glob("*.R"))
    generated = {**tables, **reports}
    manifest = {
        "schema_version": "bdev.immgen_vm25_execution_manifest.v1",
        "status": "complete_candidates_require_scientific_review",
        "run_identity": config["run_identity"],
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in tracked_inputs.items()},
        "scripts": {path.name: sha256(path) for path in scripts},
        "outputs": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in generated.items()},
        "accepted_core_sets_emitted": 0,
        "importable_artifact_created": False,
    }
    manifest_path = args.output_root / config["outputs"]["manifest"]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
