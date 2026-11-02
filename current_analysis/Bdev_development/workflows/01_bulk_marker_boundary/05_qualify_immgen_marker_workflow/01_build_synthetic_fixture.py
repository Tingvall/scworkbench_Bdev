#!/usr/bin/env python3
"""Build synthetic counts using the exact frozen ImmGen cohort structure."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
from pathlib import Path

import yaml


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").upper()


def stable_noise(gene: str, library: str) -> float:
    value = int(hashlib.sha256(f"{gene}|{library}".encode()).hexdigest()[:8], 16)
    return ((value % 25) - 12) / 100.0


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_counts(path: Path, genes: list[dict], metadata: list[dict], config: dict) -> None:
    baseline = config["synthetic_fixture"]["baseline_count"]
    marker_count = config["synthetic_fixture"]["marker_count"]
    low_count = config["synthetic_fixture"]["low_count"]
    libraries = [row["library_id"] for row in metadata]
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(["GeneSymbol", *libraries])
                for gene in genes:
                    values = []
                    for sample in metadata:
                        if gene["gene_class"] == "low_count":
                            mean = low_count
                        else:
                            target_match = (
                                (gene["family"] == "superlineage" and sample["superlineage"] == gene["target"])
                                or (gene["family"] == "broad_lineage" and sample["table_s1_lineage"] == gene["target"])
                                or (gene["family"] == "b_branch" and sample["b_branch"] == gene["target"])
                                or (gene["family"] == "b_local_stage" and sample["population"] == gene["target"])
                            )
                            mean = marker_count if target_match else baseline
                        count = max(0, int(round(mean * (1 + stable_noise(gene["gene_symbol"], sample["library_id"])))))
                        values.append(count)
                    writer.writerow([gene["gene_symbol"], *values])


def contrast_row(family: str, target: str, test_id: str, test_type: str, group_col: str, scope: str, weights: dict[str, float]) -> dict:
    return {
        "family": family,
        "target": target,
        "test_id": test_id,
        "test_type": test_type,
        "group_column": group_col,
        "sample_scope": scope,
        "weights_json": json.dumps(weights, sort_keys=True, separators=(",", ":")),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    require(config["schema_version"] == "bdev.immgen_marker_workflow_qualification.v1", "Wrong config schema")

    config_dir = args.config.parent
    manifest_path = (config_dir / config["inputs"]["manifest"]).resolve()
    hierarchy_path = (config_dir / config["inputs"]["hierarchy"]).resolve()
    hierarchy = yaml.safe_load(hierarchy_path.read_text())

    with manifest_path.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in manifest if row["discovery_include"] == "true"]
    require(len(selected) == 126, "Selected library count changed")

    superlineage_by_lineage = {}
    for superlineage, definition in hierarchy["superlineages"].items():
        for lineage in definition["table_s1_lineages"]:
            require(lineage not in superlineage_by_lineage, f"Duplicate superlineage mapping: {lineage}")
            superlineage_by_lineage[lineage] = superlineage
    branch_by_population = {}
    for branch, definition in hierarchy["b_lineage"]["branch_groups"].items():
        for population in definition["populations"]:
            require(population not in branch_by_population, f"Duplicate B-branch mapping: {population}")
            branch_by_population[population] = branch
    held_out = set(hierarchy["population_role_overrides"]["broad_validation_only"])

    metadata = []
    for row in selected:
        population = row["canonical_table_s1_population"]
        lineage = row["table_s1_lineage"]
        metadata.append(
            {
                "library_id": row["matrix_library_id"],
                "population": population,
                "replicate": row["matrix_replicate"],
                "table_s1_lineage": lineage,
                "superlineage": superlineage_by_lineage[lineage],
                "b_branch": branch_by_population.get(population, ""),
                "role": "held_out_validation" if population in held_out else "discovery",
                "organ": row["table_s1_organs"],
                "laboratories": row["table_s1_rna_labs"],
            }
        )
    require(sum(row["role"] == "discovery" for row in metadata) == 121, "Discovery library count changed")
    require(len({row["population"] for row in metadata if row["role"] == "discovery"}) == 61, "Discovery population count changed")

    sizes = config["synthetic_fixture"]["marker_genes_per_target"]
    genes: list[dict] = []
    expected: list[dict] = []

    family_targets = {
        "superlineage": list(hierarchy["superlineages"]),
        "broad_lineage": list(hierarchy["broad_lineages"]),
        "b_branch": list(hierarchy["b_lineage"]["branch_groups"]),
        "b_local_stage": list(hierarchy["b_lineage"]["local_stage_neighbors"]),
    }
    for family, targets in family_targets.items():
        for target in targets:
            target_size = sizes[family]
            no_padding = config["synthetic_fixture"]["no_padding_control"]
            if family == no_padding["family"] and target == no_padding["target"]:
                target_size = no_padding["marker_genes"]
            for index in range(1, target_size + 1):
                symbol = f"SYN_{slug(family)}_{slug(target)}_{index:03d}"
                genes.append(
                    {
                        "gene_symbol": symbol,
                        "gene_class": "injected_marker",
                        "family": family,
                        "target": target,
                        "annotation_eligible": True,
                        "state_flag": False,
                        "mapping_status": "unique_reference",
                    }
                )
                expected.append({"family": family, "target": target, "gene_symbol": symbol})

    for index in range(1, config["synthetic_fixture"]["control_genes"] + 1):
        genes.append({"gene_symbol": f"SYN_CONTROL_{index:03d}", "gene_class": "control", "family": "", "target": "", "annotation_eligible": True, "state_flag": False, "mapping_status": "unique_reference"})
        genes.append({"gene_symbol": f"SYN_FILTER_LOW_{index:03d}", "gene_class": "low_count", "family": "", "target": "", "annotation_eligible": True, "state_flag": False, "mapping_status": "unique_reference"})
    for index in range(1, config["synthetic_fixture"]["state_decoy_genes"] + 1):
        target = family_targets["broad_lineage"][(index - 1) % len(family_targets["broad_lineage"])]
        genes.append({"gene_symbol": f"SYN_STATE_DECOY_{index:03d}", "gene_class": "state_decoy", "family": "broad_lineage", "target": target, "annotation_eligible": False, "state_flag": True, "mapping_status": "unique_reference"})
    for index in range(1, config["synthetic_fixture"]["ambiguous_decoy_genes"] + 1):
        target = family_targets["superlineage"][(index - 1) % len(family_targets["superlineage"])]
        genes.append({"gene_symbol": f"SYN_AMBIGUOUS_DECOY_{index:03d}", "gene_class": "ambiguous_decoy", "family": "superlineage", "target": target, "annotation_eligible": False, "state_flag": False, "mapping_status": "ambiguous_reference"})

    contrast_plan: list[dict] = []
    superlevels = list(hierarchy["superlineages"])
    for target in superlevels:
        others = [value for value in superlevels if value != target]
        weights = {value: (-1 / len(others)) for value in others}
        weights[target] = 1.0
        contrast_plan.append(contrast_row("superlineage", target, f"SUPER_{slug(target)}_GLOBAL", "global", "superlineage", "discovery", weights))

    broadlevels = list(hierarchy["broad_lineages"])
    for target, definition in hierarchy["broad_lineages"].items():
        others = [value for value in broadlevels if value != target]
        weights = {value: (-1 / len(others)) for value in others}
        weights[target] = 1.0
        contrast_plan.append(contrast_row("broad_lineage", target, f"BROAD_{slug(target)}_GLOBAL", "global", "table_s1_lineage", "discovery", weights))
        for competitor in definition["closest_competitor_lineages"]:
            contrast_plan.append(contrast_row("broad_lineage", target, f"BROAD_{slug(target)}_VS_{slug(competitor)}", "competitor", "table_s1_lineage", "discovery", {target: 1.0, competitor: -1.0}))

    branchlevels = list(hierarchy["b_lineage"]["branch_groups"])
    for contrast in hierarchy["b_lineage"]["branch_contrasts"]:
        target = contrast["target_group"]
        others = [value for value in branchlevels if value != target]
        weights = {value: (-1 / len(others)) for value in others}
        weights[target] = 1.0
        contrast_plan.append(contrast_row("b_branch", target, f"{contrast['contrast_id']}_GLOBAL", "global", "b_branch", "b_discovery", weights))
        for competitor in contrast["competitor_groups"]:
            contrast_plan.append(contrast_row("b_branch", target, f"{contrast['contrast_id']}_VS_{slug(competitor)}", "competitor", "b_branch", "b_discovery", {target: 1.0, competitor: -1.0}))

    for target, neighbors in hierarchy["b_lineage"]["local_stage_neighbors"].items():
        for neighbor in neighbors:
            contrast_plan.append(contrast_row("b_local_stage", target, f"BLOCAL_{slug(target)}_VS_{slug(neighbor)}", "neighbor", "population", "b_local_discovery", {target: 1.0, neighbor: -1.0}))

    paths = config["outputs"]["fixture"]
    counts_path = args.output_root / paths["counts"]
    counts_path.parent.mkdir(parents=True, exist_ok=True)
    write_counts(counts_path, genes, metadata, config)
    write_tsv(args.output_root / paths["metadata"], metadata, list(metadata[0]))
    annotation_rows = [{key: gene[key] for key in ("gene_symbol", "gene_class", "annotation_eligible", "state_flag", "mapping_status")} for gene in genes]
    write_tsv(args.output_root / paths["annotation"], annotation_rows, list(annotation_rows[0]))
    write_tsv(args.output_root / paths["expected_markers"], expected, ["family", "target", "gene_symbol"])
    write_tsv(args.output_root / paths["contrast_plan"], contrast_plan, list(contrast_plan[0]))
    print(json.dumps({"genes": len(genes), "libraries": len(metadata), "expected_markers": len(expected), "contrasts": len(contrast_plan)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
