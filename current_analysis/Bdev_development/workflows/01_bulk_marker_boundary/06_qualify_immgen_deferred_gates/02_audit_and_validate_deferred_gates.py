#!/usr/bin/env python3
"""Audit confounding and validate the synthetic deferred-gate qualification."""

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


def split_labs(value: str) -> set[str]:
    return {item for item in value.split(";") if item}


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    require(config["schema_version"] == "bdev.immgen_deferred_gate_qualification.v1", "Wrong config schema")
    require(args.output_root.resolve() == Path(config["outputs"]["root"]).resolve(), "Output root differs from frozen configuration")
    config_dir = args.config.parent
    contract_path = (config_dir / config["inputs"]["statistical_contract"]).resolve()
    hierarchy_path = (config_dir / config["inputs"]["hierarchy"]).resolve()
    software_path = (config_dir / config["inputs"]["software_lock"]).resolve()
    core_config_path = (config_dir / config["inputs"]["core_qualification_config"]).resolve()
    contract = yaml.safe_load(contract_path.read_text())
    hierarchy = yaml.safe_load(hierarchy_path.read_text())
    input_root = Path(config["inputs"]["core_qualification_root"])

    input_paths: dict[str, Path] = {}
    for name in (
        "synthetic_counts",
        "synthetic_metadata",
        "expected_markers",
        "contrast_plan",
        "synthetic_selected_genesets",
        "core_qualification_manifest",
    ):
        path = input_root / config["inputs"][name]["path"]
        require(path.is_file(), f"Missing input: {path}")
        require(sha256(path) == config["inputs"][name]["sha256"], f"Checksum mismatch: {name}")
        input_paths[name] = path

    metadata = read_tsv(input_paths["synthetic_metadata"])
    plan = read_tsv(input_paths["contrast_plan"])
    expected = read_tsv(input_paths["expected_markers"])
    selected = read_tsv(input_paths["synthetic_selected_genesets"])
    require(len(metadata) == 126, "Metadata library count changed")
    require(sum(row["role"] == "discovery" for row in metadata) == 121, "Discovery count changed")

    table_paths = {name: args.output_root / rel for name, rel in config["outputs"]["tables"].items()}
    report_paths = {name: args.output_root / rel for name, rel in config["outputs"]["reports"].items()}

    # Exact metadata contingency tables.
    contingency_rows: list[dict] = []
    for population in sorted({row["population"] for row in metadata}):
        population_rows = [row for row in metadata if row["population"] == population]
        for factor, column in (("organ", "organ"), ("laboratory", "laboratories")):
            counts: Counter[str] = Counter()
            for row in population_rows:
                values = split_labs(row[column]) if factor == "laboratory" else {row[column]}
                for value in values:
                    counts[value] += 1
            for value in sorted(counts):
                contingency_rows.append(
                    {
                        "population": population,
                        "factor": factor,
                        "level": value,
                        "libraries": counts[value],
                    }
                )
    write_tsv(table_paths["population_factor_contingency"], contingency_rows, ["population", "factor", "level", "libraries"])

    # Flag complete target/competitor separation for every declared contrast.
    audit_rows: list[dict] = []
    for test in plan:
        rows = [row for row in metadata if row["role"] == "discovery"]
        if test["sample_scope"] in {"b_discovery", "b_local_discovery"}:
            rows = [row for row in rows if row["table_s1_lineage"] == "B"]
        weights = json.loads(test["weights_json"])
        competitor_groups = {group for group, weight in weights.items() if float(weight) < 0}
        group_column = test["group_column"]
        target_rows = [row for row in rows if row[group_column] == test["target"]]
        competitor_rows = [row for row in rows if row[group_column] in competitor_groups]
        require(target_rows and competitor_rows, f"Empty confounding audit group: {test['test_id']}")
        target_organs = {row["organ"] for row in target_rows}
        competitor_organs = {row["organ"] for row in competitor_rows}
        target_labs = set().union(*(split_labs(row["laboratories"]) for row in target_rows))
        competitor_labs = set().union(*(split_labs(row["laboratories"]) for row in competitor_rows))
        audit_rows.append(
            {
                "family": test["family"],
                "target": test["target"],
                "test_id": test["test_id"],
                "test_type": test["test_type"],
                "target_organs": ";".join(sorted(target_organs)),
                "competitor_organs": ";".join(sorted(competitor_organs)),
                "organ_completely_separated": bool_text(target_organs.isdisjoint(competitor_organs)),
                "target_laboratories": ";".join(sorted(target_labs)),
                "competitor_laboratories": ";".join(sorted(competitor_labs)),
                "laboratory_completely_separated": bool_text(target_labs.isdisjoint(competitor_labs)),
            }
        )
    write_tsv(
        table_paths["contrast_confounding_audit"],
        audit_rows,
        [
            "family",
            "target",
            "test_id",
            "test_type",
            "target_organs",
            "competitor_organs",
            "organ_completely_separated",
            "target_laboratories",
            "competitor_laboratories",
            "laboratory_completely_separated",
        ],
    )

    control_rows: list[dict] = []
    tissue = config["synthetic_confounding_controls"]["unrelated_tissue_bias"]
    threshold = float(tissue["rejection_effect_threshold"])
    for case in tissue["cases"]:
        reproduced = float(case["source_effect"]) > 0 and float(case["unrelated_lineage_effect"]) >= threshold
        outcome = "rejected_tissue_bias" if reproduced else "retained_no_reproduced_bias"
        require(outcome == case["expected_outcome"], f"Tissue-bias control failed: {case['control_id']}")
        control_rows.append(
            {
                "control_type": "unrelated_tissue_bias",
                "control_id": case["control_id"],
                "source_effect": case["source_effect"],
                "unrelated_lineage_effect": case["unrelated_lineage_effect"],
                "completely_separated": "",
                "independent_support": "",
                "outcome": outcome,
                "expected_outcome": case["expected_outcome"],
                "passed": "true",
            }
        )
    for case in config["synthetic_confounding_controls"]["unresolved_confounding"]["cases"]:
        outcome = "candidate_confounding" if case["completely_separated"] and not case["independent_support"] else "passes_confounding_gate"
        require(outcome == case["expected_outcome"], f"Unresolved-confounding control failed: {case['control_id']}")
        control_rows.append(
            {
                "control_type": "unresolved_confounding",
                "control_id": case["control_id"],
                "source_effect": "",
                "unrelated_lineage_effect": "",
                "completely_separated": bool_text(case["completely_separated"]),
                "independent_support": bool_text(case["independent_support"]),
                "outcome": outcome,
                "expected_outcome": case["expected_outcome"],
                "passed": "true",
            }
        )
    write_tsv(
        table_paths["confounding_controls"],
        control_rows,
        ["control_type", "control_id", "source_effect", "unrelated_lineage_effect", "completely_separated", "independent_support", "outcome", "expected_outcome", "passed"],
    )

    schedule = read_tsv(table_paths["resample_schedule"])
    iterations = int(config["balanced_resampling"]["iterations"])
    discovery_populations = {row["population"] for row in metadata if row["role"] == "discovery"}
    require(len(schedule) == iterations * len(discovery_populations), "Resample schedule size mismatch")
    valid_library = {row["library_id"]: row["population"] for row in metadata if row["role"] == "discovery"}
    for iteration in range(1, iterations + 1):
        rows = [row for row in schedule if int(row["iteration"]) == iteration]
        require({row["population"] for row in rows} == discovery_populations, f"Population coverage failed in iteration {iteration}")
        require(all(valid_library[row["library_id"]] == row["population"] for row in rows), f"Library mapping failed in iteration {iteration}")
    usage: defaultdict[str, Counter[str]] = defaultdict(Counter)
    for row in schedule:
        usage[row["population"]][row["library_id"]] += 1
    require(all(max(counts.values()) - min(counts.values()) <= 1 for counts in usage.values()), "Library schedule is not balanced")

    balanced = read_tsv(table_paths["balanced_stability"], compressed=True)
    stability = {(row["family"], row["target"], row["gene_symbol"]): float(row["stability_fraction"]) for row in balanced}
    expected_keys = {(row["family"], row["target"], row["gene_symbol"]) for row in expected}
    min_stability = {family: float(contract["contrast_families"][family]["minimum_stability_fraction"]) for family in ("superlineage", "broad_lineage", "b_branch", "b_local_stage")}
    failed_expected = [key for key in expected_keys if stability.get(key, -1) < min_stability[key[0]]]
    require(not failed_expected, f"Injected markers failed balanced stability: {failed_expected[:5]}")
    negative_values = [value for (family, target, gene), value in stability.items() if gene.startswith("SYN_CONTROL_")]
    require(negative_values and max(negative_values) < 0.80, "Negative controls passed stability")

    selected_keys = {(row["family"], row["target"], row["gene_symbol"]) for row in selected}
    population_loo = read_tsv(table_paths["population_loo"])
    library_loo = read_tsv(table_paths["library_loo"])
    for rows, gate_name in ((population_loo, "population LOO"), (library_loo, "library LOO")):
        for row in rows:
            key = (row["family"], row["target"], row["gene_symbol"])
            require(key in selected_keys, f"Unknown selected gene in {gate_name}: {key}")
            if row["status"] == "assessed" and key in expected_keys:
                require(float(row["holdout_fraction"]) >= 0.80, f"Injected marker failed {gate_name}: {key}")

    lopo_predictions = read_tsv(table_paths["broad_lopo_predictions"])
    overall_accuracy = sum(row["correct"] == "TRUE" for row in lopo_predictions) / len(lopo_predictions)
    require(overall_accuracy >= float(contract["internal_validation"]["broad_leave_one_population_out"]["minimum_overall_top1_accuracy"]), "Broad LOPO accuracy failed")
    lineage_sensitivity = {}
    for lineage in sorted({row["expected_lineage"] for row in lopo_predictions}):
        rows = [row for row in lopo_predictions if row["expected_lineage"] == lineage]
        lineage_sensitivity[lineage] = sum(row["correct"] == "TRUE" for row in rows) / len(rows)
    require(min(lineage_sensitivity.values()) >= float(contract["internal_validation"]["broad_leave_one_population_out"]["minimum_per_lineage_sensitivity"]), "Broad LOPO lineage sensitivity failed")

    heldout_predictions = read_tsv(table_paths["heldout_predictions"])
    expected_heldout = set(hierarchy["population_role_overrides"]["broad_validation_only"])
    require({row["population"] for row in heldout_predictions} == expected_heldout, "Held-out population set changed")
    require(all(row["correct"] == "TRUE" for row in heldout_predictions), "Held-out parent-lineage scoring failed")

    local_auroc = read_tsv(table_paths["local_neighbor_auroc"])
    minimum_auroc = float(contract["internal_validation"]["local_neighbor_validation"]["minimum_pairwise_auroc"])
    assessed_aurocs = [float(row["auroc"]) for row in local_auroc if row["status"] == "assessed"]
    require(assessed_aurocs and min(assessed_aurocs) >= minimum_auroc, "Assessed local-neighbor AUROC failed")
    require(any(row["status"] == "insufficient_replication" for row in local_auroc), "Insufficient-replication path was not exercised")

    require(len(audit_rows) == len(plan), "Not every contrast received a confounding audit")
    require(any(row["organ_completely_separated"] == "true" for row in audit_rows), "Organ separation path not exercised")
    require(any(row["laboratory_completely_separated"] == "true" for row in audit_rows), "Laboratory separation path not exercised")

    summary = {
        "schema_version": "bdev.immgen_deferred_gate_qualification_summary.v1",
        "status": "deferred_gate_mechanics_qualified_full_immgen_execution_required",
        "synthetic_only": True,
        "immgen_expression_values_read": False,
        "balanced_resampling_iterations": iterations,
        "balanced_schedule_rows": len(schedule),
        "balanced_schedule_max_library_imbalance": max(max(counts.values()) - min(counts.values()) for counts in usage.values()),
        "injected_markers_passing_balanced_stability": len(expected_keys),
        "injected_markers_failing_balanced_stability": 0,
        "negative_control_max_stability": max(negative_values),
        "population_loo_rows": len(population_loo),
        "library_loo_rows": len(library_loo),
        "broad_lopo_populations": len(lopo_predictions),
        "broad_lopo_overall_accuracy": overall_accuracy,
        "broad_lopo_lineage_sensitivity": lineage_sensitivity,
        "heldout_predictions": heldout_predictions,
        "local_neighbor_tests": len(local_auroc),
        "local_neighbor_assessed_minimum_auroc": min(assessed_aurocs),
        "confounding_contrasts_audited": len(audit_rows),
        "organ_completely_separated_contrasts": sum(row["organ_completely_separated"] == "true" for row in audit_rows),
        "laboratory_completely_separated_contrasts": sum(row["laboratory_completely_separated"] == "true" for row in audit_rows),
        "confounding_controls_passed": len(control_rows),
        "accepted_core_sets_emitted": 0,
        "importable_gmt_created": False,
        "scworkbench_invoked": False,
        "deferred": config["qualification_scope"]["explicitly_deferred"],
    }
    report_paths["qualification_summary"].parent.mkdir(parents=True, exist_ok=True)
    report_paths["qualification_summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    tracked_inputs = {
        "qualification_config": args.config.resolve(),
        "statistical_contract": contract_path,
        "hierarchy": hierarchy_path,
        "software_lock": software_path,
        "core_qualification_config": core_config_path,
        **input_paths,
    }
    script_dir = Path(__file__).resolve().parent
    scripts = sorted(script_dir.glob("*.py")) + sorted(script_dir.glob("*.R"))
    generated = {**table_paths, **report_paths}
    manifest = {
        "schema_version": "bdev.immgen_deferred_gate_qualification_manifest.v1",
        "status": "passed_full_immgen_execution_required",
        "synthetic_only": True,
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in tracked_inputs.items()},
        "scripts": {path.name: sha256(path) for path in scripts},
        "outputs": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in generated.items()},
        "importable_artifact_created": False,
    }
    manifest_path = args.output_root / config["outputs"]["manifest"]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
