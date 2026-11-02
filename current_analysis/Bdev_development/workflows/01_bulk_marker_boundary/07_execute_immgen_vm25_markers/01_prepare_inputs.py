#!/usr/bin/env python3
"""Prepare frozen selected counts, metadata, and contrasts for ImmGen execution."""

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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve(config_dir: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (config_dir / path).resolve()


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").upper()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def contrast_row(family: str, target: str, test_id: str, test_type: str, group_column: str, sample_scope: str, weights: dict[str, float]) -> dict:
    return {
        "family": family,
        "target": target,
        "test_id": test_id,
        "test_type": test_type,
        "group_column": group_column,
        "sample_scope": sample_scope,
        "weights_json": json.dumps(weights, sort_keys=True, separators=(",", ":")),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    require(config["schema_version"] == "bdev.immgen_vm25_marker_execution.v1", "Wrong config schema")
    require(args.output_root.resolve() == Path(config["outputs"]["root"]).resolve(), "Output root differs from configuration")
    config_dir = args.config.parent

    input_paths = {name: resolve(config_dir, spec["path"]) for name, spec in config["inputs"].items() if isinstance(spec, dict) and "path" in spec}
    for name, path in input_paths.items():
        require(path.is_file(), f"Missing input: {path}")
        require(sha256(path) == config["inputs"][name]["sha256"], f"Checksum mismatch: {name}")

    hierarchy = yaml.safe_load(input_paths["hierarchy"].read_text())
    with input_paths["library_manifest"].open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in manifest if row["discovery_include"] == "true"]
    require(len(selected) == config["execution"]["selected_libraries"], "Selected library count changed")

    superlineage_by_lineage: dict[str, str] = {}
    for superlineage, definition in hierarchy["superlineages"].items():
        for lineage in definition["table_s1_lineages"]:
            require(lineage not in superlineage_by_lineage, f"Duplicate lineage mapping: {lineage}")
            superlineage_by_lineage[lineage] = superlineage
    branch_by_population: dict[str, str] = {}
    for branch, definition in hierarchy["b_lineage"]["branch_groups"].items():
        for population in definition["populations"]:
            require(population not in branch_by_population, f"Duplicate branch mapping: {population}")
            branch_by_population[population] = branch
    heldout = set(hierarchy["population_role_overrides"]["broad_validation_only"])

    metadata: list[dict] = []
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
                "role": "held_out_validation" if population in heldout else "discovery",
                "organ": row["table_s1_organs"],
                "laboratories": row["table_s1_rna_labs"],
            }
        )
    require(sum(row["role"] == "discovery" for row in metadata) == config["execution"]["discovery_libraries"], "Discovery library count changed")
    require(sum(row["role"] == "held_out_validation" for row in metadata) == config["execution"]["held_out_libraries"], "Held-out library count changed")

    tables = config["outputs"]["tables"]
    metadata_path = args.output_root / tables["metadata"]
    write_tsv(metadata_path, metadata, list(metadata[0]))

    selected_ids = [row["library_id"] for row in metadata]
    counts_path = args.output_root / tables["selected_counts"]
    counts_path.parent.mkdir(parents=True, exist_ok=True)
    gene_rows = 0
    with gzip.open(input_paths["count_matrix"], "rt", newline="") as source:
        reader = csv.reader(source)
        header = next(reader)
        require(header[0] == "GeneSymbol", "Count matrix identifier column changed")
        column_index = {name: index for index, name in enumerate(header)}
        require(all(name in column_index for name in selected_ids), "Selected library missing from matrix")
        indices = [column_index[name] for name in selected_ids]
        with counts_path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                    writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                    writer.writerow(["GeneSymbol", *selected_ids])
                    for row in reader:
                        require(len(row) == len(header), f"Malformed count row {gene_rows + 2}")
                        values = [row[index] for index in indices]
                        require(all(value.isdigit() for value in values), f"Noninteger count at row {gene_rows + 2}")
                        writer.writerow([row[0], *values])
                        gene_rows += 1
    require(gene_rows == config["execution"]["expected_gene_rows"], "Gene row count changed")

    contrast_plan: list[dict] = []
    superlevels = list(hierarchy["superlineages"])
    for target in superlevels:
        others = [value for value in superlevels if value != target]
        weights = {value: -1 / len(others) for value in others} | {target: 1.0}
        contrast_plan.append(contrast_row("superlineage", target, f"SUPER_{slug(target)}_GLOBAL", "global", "superlineage", "discovery", weights))
    broadlevels = list(hierarchy["broad_lineages"])
    for target, definition in hierarchy["broad_lineages"].items():
        others = [value for value in broadlevels if value != target]
        weights = {value: -1 / len(others) for value in others} | {target: 1.0}
        contrast_plan.append(contrast_row("broad_lineage", target, f"BROAD_{slug(target)}_GLOBAL", "global", "table_s1_lineage", "discovery", weights))
        for competitor in definition["closest_competitor_lineages"]:
            contrast_plan.append(contrast_row("broad_lineage", target, f"BROAD_{slug(target)}_VS_{slug(competitor)}", "competitor", "table_s1_lineage", "discovery", {target: 1.0, competitor: -1.0}))
    branchlevels = list(hierarchy["b_lineage"]["branch_groups"])
    for definition in hierarchy["b_lineage"]["branch_contrasts"]:
        target = definition["target_group"]
        others = [value for value in branchlevels if value != target]
        weights = {value: -1 / len(others) for value in others} | {target: 1.0}
        contrast_plan.append(contrast_row("b_branch", target, f"{definition['contrast_id']}_GLOBAL", "global", "b_branch", "b_discovery", weights))
        for competitor in definition["competitor_groups"]:
            contrast_plan.append(contrast_row("b_branch", target, f"{definition['contrast_id']}_VS_{slug(competitor)}", "competitor", "b_branch", "b_discovery", {target: 1.0, competitor: -1.0}))
    for target, neighbors in hierarchy["b_lineage"]["local_stage_neighbors"].items():
        for neighbor in neighbors:
            contrast_plan.append(contrast_row("b_local_stage", target, f"BLOCAL_{slug(target)}_VS_{slug(neighbor)}", "neighbor", "population", "b_local_discovery", {target: 1.0, neighbor: -1.0}))
    require(len(contrast_plan) == 76, "Contrast count changed")
    plan_path = args.output_root / tables["contrast_plan"]
    write_tsv(plan_path, contrast_plan, list(contrast_plan[0]))

    summary = {
        "schema_version": "bdev.immgen_vm25_preparation_summary.v1",
        "status": "passed",
        "gene_rows": gene_rows,
        "selected_libraries": len(metadata),
        "discovery_libraries": sum(row["role"] == "discovery" for row in metadata),
        "held_out_libraries": sum(row["role"] == "held_out_validation" for row in metadata),
        "discovery_populations": len({row["population"] for row in metadata if row["role"] == "discovery"}),
        "contrasts": len(contrast_plan),
        "selected_counts_sha256": sha256(counts_path),
        "metadata_sha256": sha256(metadata_path),
        "contrast_plan_sha256": sha256(plan_path),
    }
    report_path = args.output_root / config["outputs"]["reports"]["preparation_summary"]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
