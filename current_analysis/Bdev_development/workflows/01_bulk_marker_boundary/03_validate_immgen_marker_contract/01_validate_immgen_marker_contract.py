#!/usr/bin/env python3
"""Validate the approved ImmGen marker statistical contract."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import yaml


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--hierarchy", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    return parser.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = yaml.safe_load(handle)
    require(isinstance(value, dict), f"Expected mapping in {path}")
    return value


def main() -> None:
    args = parse_args()
    with args.manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    hierarchy = load_yaml(args.hierarchy)
    contract = load_yaml(args.contract)

    require(
        hierarchy.get("schema_version") == "bdev.immgen_geneset_hierarchy.v1",
        "Unsupported hierarchy schema",
    )
    require(
        contract.get("schema_version") == "bdev.immgen_marker_statistical_contract.v1",
        "Unsupported contract schema",
    )
    require(contract.get("status") == "approved", "Statistical contract is not approved")
    require(contract["execution_boundary"] == "outside_scworkbench", "Wrong execution boundary")
    require(contract["software"]["primary_engine"] == "edgeR_quasi_likelihood", "Wrong engine")

    selected_rows = [row for row in rows if row["discovery_include"] == "true"]
    require(len(selected_rows) == 126, "Frozen selected-library count changed")
    population_libraries = Counter(
        row["canonical_table_s1_population"] for row in selected_rows
    )
    selected_populations = set(population_libraries)
    require(len(selected_populations) == 64, "Frozen selected-population count changed")

    validation_populations = set(
        hierarchy["population_role_overrides"]["broad_validation_only"]
    )
    contract_validation = set(
        contract["analysis_layers"]["held_out_broad_validation"]["populations"]
    )
    require(contract_validation == validation_populations, "Held-out populations disagree")
    validation_libraries = sum(population_libraries[pop] for pop in validation_populations)
    discovery_populations = selected_populations - validation_populations
    discovery_libraries = sum(population_libraries[pop] for pop in discovery_populations)

    library_layer = contract["analysis_layers"]["library_level"]
    require(library_layer["input_libraries"] == len(selected_rows), "Input-library count mismatch")
    require(
        library_layer["held_out_validation_libraries"] == validation_libraries,
        "Held-out-library count mismatch",
    )
    require(
        library_layer["discovery_model_libraries"] == discovery_libraries,
        "Discovery-library count mismatch",
    )
    require(
        contract["analysis_layers"]["population_level"]["populations"]
        == len(discovery_populations),
        "Discovery-population count mismatch",
    )

    require(contract["count_processing"]["expected_gene_rows"] == 55011, "Gene count changed")
    require(contract["count_processing"]["model"]["robust"] is True, "Robust QL required")
    require(
        contract["count_processing"]["model"]["multiple_testing"]
        == "Benjamini-Hochberg within each declared contrast",
        "Unexpected multiple-testing policy",
    )

    expected_families = {"superlineage", "broad_lineage", "b_branch", "b_local_stage"}
    families = contract["contrast_families"]
    require(set(families) == expected_families, "Contrast-family set changed")
    for family_name, family in families.items():
        require(0 < family["minimum_genes"] <= family["maximum_genes"], f"Bad size for {family_name}")
        if "minimum_stability_fraction" in family:
            require(
                0.5 <= family["minimum_stability_fraction"] <= 1,
                f"Bad stability threshold for {family_name}",
            )
        for key, value in family.items():
            if key.startswith("maximum_") and key.endswith("_fdr"):
                require(0 < value <= 0.10, f"Bad FDR threshold {family_name}.{key}")

    local_targets = set(hierarchy["b_lineage"]["local_stage_neighbors"])
    require(len(local_targets) == 9, "Local-stage target count changed")
    singleton_local_targets = sorted(
        population for population in local_targets if population_libraries[population] == 1
    )

    require(
        contract["stability"]["balanced_library_resampling"]["iterations"] >= 200,
        "Too few stability resamples",
    )
    require(contract["random_seed"] == 20261031, "Unexpected random seed")
    require(contract["deterministic_selection"]["eligibility_first"].startswith("Never pad"), "Padding rule missing")
    require(contract["deterministic_selection"]["positive_only_core_sets"] is True, "Core sets must be positive-only")
    require(contract["scworkbench_boundary"]["core_or_pack_change_required"] is False, "Unexpected framework change")

    summary = {
        "schema_version": contract["schema_version"],
        "status": "structure_passed_contract_approved",
        "selected_input_libraries": len(selected_rows),
        "discovery_model_libraries": discovery_libraries,
        "held_out_validation_libraries": validation_libraries,
        "discovery_populations": len(discovery_populations),
        "held_out_validation_populations": sorted(validation_populations),
        "contrast_families": sorted(families),
        "local_b_targets": len(local_targets),
        "singleton_local_b_targets": singleton_local_targets,
        "balanced_resampling_iterations": contract["stability"]["balanced_library_resampling"]["iterations"],
        "random_seed": contract["random_seed"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
