#!/usr/bin/env python3
"""Freeze broad and shared-gene-filtered collections for hierarchical scoring."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
from collections import defaultdict
from pathlib import Path

import yaml


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_file(path: Path, expected: str, label: str) -> None:
    observed = sha256_file(path)
    if observed != expected:
        raise RuntimeError(f"{label} checksum mismatch: {observed}")


def read_gmt(path: Path) -> dict[str, tuple[str, list[str]]]:
    output: dict[str, tuple[str, list[str]]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] in output:
                raise RuntimeError(f"Invalid GMT row in {path}")
            output[fields[0]] = (fields[1], list(dict.fromkeys(fields[2:])))
    return output


def write_gmt(path: Path, sets: dict[str, tuple[str, list[str]]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        for set_id, (description, genes) in sets.items():
            handle.write("\t".join([set_id, description, *genes]) + "\n")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def slug(value: str) -> str:
    return "".join(character if character.isalnum() else "_" for character in value).strip("_").upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    config = yaml.safe_load(config_path.read_text())
    inputs_path = Path(config["inputs_manifest"])
    inputs = yaml.safe_load(inputs_path.read_text())
    results_root = Path(config["storage"]["results_root"])
    tables = results_root / "01_tables"
    reports = results_root / "03_reports"
    manifests = results_root / "06_manifests"
    for directory in (tables, reports, manifests):
        directory.mkdir(parents=True, exist_ok=True)

    parent = inputs["parent_execution"]
    broad = inputs["broad_lineage_collection"]
    population = inputs["population_collection"]
    checks = [
        (Path(parent["execution_clone"]) / "manifest.json", parent["manifest_sha256"], "parent manifest"),
        (Path(parent["final_run_manifest"]), parent["final_run_manifest_sha256"], "parent final manifest"),
        (Path(parent["fixture_cells"]), parent["fixture_cells_sha256"], "fixture cells"),
        (Path(parent["fixture_preflight"]), parent["fixture_preflight_sha256"], "fixture preflight"),
        (Path(broad["source_membership"]), broad["source_membership_sha256"], "broad membership"),
        (Path(broad["source_summary"]), broad["source_summary_sha256"], "broad summary"),
        (Path(population["source_gmt"]), population["source_gmt_sha256"], "population GMT"),
        (Path(population["target_summary"]), population["target_summary_sha256"], "population summary"),
    ]
    for path, expected, label in checks:
        check_file(path, expected, label)

    fixture_copy = tables / "01_fixture_cells.tsv"
    shutil.copyfile(parent["fixture_cells"], fixture_copy)

    membership_rows: list[dict[str, str]] = []
    with gzip.open(broad["source_membership"], "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["family"] == broad["family"]:
                membership_rows.append(row)
    broad_by_target: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for row in membership_rows:
        broad_by_target[row["target"]].append((int(row["final_rank"]), row["gene_symbol"]))
    if len(broad_by_target) != int(broad["expected_sets"]):
        raise RuntimeError("Unexpected broad-lineage set count")
    broad_sets: dict[str, tuple[str, list[str]]] = {}
    broad_mapping: list[dict[str, object]] = []
    for target in sorted(broad_by_target):
        set_id = f"IMMGEN_BROAD_{slug(target)}"
        genes = [gene for _, gene in sorted(broad_by_target[target])]
        broad_sets[set_id] = (f"ImmGen diagnostic broad lineage: {target}", genes)
        broad_mapping.append({"broad_set_id": set_id, "lineage": target, "genes": len(genes)})
    broad_gmt = tables / "02_broad_lineage_genesets.gmt"
    write_gmt(broad_gmt, broad_sets)
    write_tsv(tables / "03_broad_lineage_set_map.tsv", ["broad_set_id", "lineage", "genes"], broad_mapping)

    original_sets = read_gmt(Path(population["source_gmt"]))
    summary_rows = list(csv.DictReader(open(population["target_summary"], encoding="utf-8"), delimiter="\t"))
    set_to_lineage = {row["geneset_id"]: row["lineage"] for row in summary_rows}
    if set(original_sets) != set(set_to_lineage) or len(original_sets) != int(population["expected_sets"]):
        raise RuntimeError("Population GMT and target summary disagree")

    gene_sets: dict[str, set[str]] = defaultdict(set)
    gene_lineages: dict[str, set[str]] = defaultdict(set)
    for set_id, (_, genes) in original_sets.items():
        for gene in genes:
            gene_sets[gene].add(set_id)
            gene_lineages[gene].add(set_to_lineage[set_id])
    lineage_threshold = int(config["shared_gene_filter"]["minimum_distinct_lineages"])
    removed = {gene for gene, lineages in gene_lineages.items() if len(lineages) >= lineage_threshold}
    filtered_sets: dict[str, tuple[str, list[str]]] = {}
    set_audit: list[dict[str, object]] = []
    for set_id, (description, genes) in original_sets.items():
        retained = [gene for gene in genes if gene not in removed]
        filtered_sets[set_id] = (description + "; genes shared across >=3 lineages removed", retained)
        set_audit.append({
            "geneset_id": set_id,
            "lineage": set_to_lineage[set_id],
            "original_genes": len(genes),
            "removed_genes": len(genes) - len(retained),
            "retained_genes": len(retained),
            "retained_fraction": len(retained) / len(genes),
        })
    minimum_retained = min(row["retained_genes"] for row in set_audit)
    if minimum_retained < int(config["shared_gene_filter"]["minimum_retained_genes_per_set"]):
        raise RuntimeError("Shared-gene filter leaves an undersized population set")
    filtered_gmt = tables / "04_shared_gene_filtered_population_genesets.gmt"
    write_gmt(filtered_gmt, filtered_sets)
    write_tsv(
        tables / "05_population_filter_set_audit.tsv",
        ["geneset_id", "lineage", "original_genes", "removed_genes", "retained_genes", "retained_fraction"],
        set_audit,
    )
    removed_rows = [
        {
            "gene": gene,
            "population_set_count": len(gene_sets[gene]),
            "broad_lineage_count": len(gene_lineages[gene]),
            "broad_lineages": ";".join(sorted(gene_lineages[gene])),
            "population_sets": ";".join(sorted(gene_sets[gene])),
        }
        for gene in sorted(removed)
    ]
    write_tsv(
        tables / "06_removed_shared_genes.tsv",
        ["gene", "population_set_count", "broad_lineage_count", "broad_lineages", "population_sets"],
        removed_rows,
    )

    preflight = {
        "schema_version": "bdev.immgen_bdev_hierarchy_preflight.v1",
        "status": "passed",
        "run_identity": config["run_identity"],
        "fixture_cells": sum(1 for _ in fixture_copy.open()) - 1,
        "broad_sets": len(broad_sets),
        "population_sets": len(original_sets),
        "removed_shared_genes": len(removed),
        "minimum_retained_genes": minimum_retained,
        "maximum_retained_genes": max(row["retained_genes"] for row in set_audit),
        "checks": {label: True for _, _, label in checks},
    }
    (reports / "01_collection_preflight.json").write_text(json.dumps(preflight, indent=2, sort_keys=True) + "\n")
    input_manifest = {
        "schema_version": "bdev.immgen_bdev_hierarchy_input_manifest.v1",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256_file(config_path)},
        "inputs_manifest": {"path": str(inputs_path), "sha256": sha256_file(inputs_path)},
        "parent_bundle": {
            "path": parent["execution_clone"],
            "fingerprint": parent["bundle_fingerprint"],
            "manifest_sha256": sha256_file(Path(parent["execution_clone"]) / "manifest.json"),
        },
        "generated_collections": {
            "broad_gmt": {"path": str(broad_gmt), "sha256": sha256_file(broad_gmt), "sets": len(broad_sets)},
            "filtered_population_gmt": {"path": str(filtered_gmt), "sha256": sha256_file(filtered_gmt), "sets": len(filtered_sets)},
            "removed_shared_genes": len(removed),
        },
    }
    (manifests / "01_input_manifest.json").write_text(json.dumps(input_manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
