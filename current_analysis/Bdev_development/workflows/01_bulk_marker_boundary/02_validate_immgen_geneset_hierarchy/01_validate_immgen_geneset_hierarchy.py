#!/usr/bin/env python3
"""Validate the approved ImmGen hierarchy against the frozen vM25 manifest."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
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
    return parser.parse_args()


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(bool(rows), "Manifest is empty")
    return rows


def population_table(rows: list[dict[str, str]]) -> dict[str, dict[str, Any]]:
    populations: dict[str, dict[str, Any]] = {}
    for row in rows:
        if row["discovery_include"] != "true":
            continue
        population = row["canonical_table_s1_population"]
        if population not in populations:
            populations[population] = {
                "lineage": row["table_s1_lineage"],
                "libraries": 0,
                "atac_match": row["table_s1_atac_match"] == "true",
            }
        require(
            populations[population]["lineage"] == row["table_s1_lineage"],
            f"Inconsistent lineage for {population}",
        )
        populations[population]["libraries"] += 1
        populations[population]["atac_match"] &= row["table_s1_atac_match"] == "true"
    return populations


def main() -> None:
    args = parse_args()
    rows = load_manifest(args.manifest)
    with args.hierarchy.open(encoding="utf-8") as handle:
        hierarchy = yaml.safe_load(handle)
    require(
        hierarchy.get("schema_version") == "bdev.immgen_geneset_hierarchy.v1",
        "Unsupported hierarchy schema",
    )

    populations = population_table(rows)
    selected = set(populations)
    require(len(selected) == 64, f"Expected 64 selected populations, observed {len(selected)}")
    require(all(value["atac_match"] for value in populations.values()), "Selected ATAC match missing")

    overrides = hierarchy["population_role_overrides"]
    validation_only = set(overrides["broad_validation_only"])
    unique_branch = set(overrides["retained_unique_branch_discovery"])
    require(validation_only <= selected, "A validation-only population is not in the frozen cohort")
    require(unique_branch <= selected, "A retained unique-branch population is not in the frozen cohort")
    require(validation_only.isdisjoint(unique_branch), "Population roles overlap")

    discovery = selected - validation_only
    observed_lineages: dict[str, set[str]] = defaultdict(set)
    for population in discovery:
        observed_lineages[populations[population]["lineage"]].add(population)

    broad_lineages = hierarchy["broad_lineages"]
    require(set(broad_lineages) == set(observed_lineages), "Broad-lineage labels do not match manifest")
    for lineage, definition in broad_lineages.items():
        observed = len(observed_lineages[lineage])
        expected = definition["expected_discovery_populations"]
        require(observed == expected, f"{lineage}: expected {expected} populations, observed {observed}")
        for competitor in definition["closest_competitor_lineages"]:
            require(competitor in broad_lineages, f"Unknown competitor lineage: {competitor}")
            require(competitor != lineage, f"Lineage cannot compete with itself: {lineage}")

    superlineage_members = []
    for definition in hierarchy["superlineages"].values():
        superlineage_members.extend(definition["table_s1_lineages"])
    require(
        Counter(superlineage_members) == Counter(broad_lineages.keys()),
        "Superlineages must partition broad lineages exactly once",
    )

    b_lineage = hierarchy["b_lineage"]
    branch_members: list[str] = []
    for definition in b_lineage["branch_groups"].values():
        branch_members.extend(definition["populations"])
    expected_b_discovery = observed_lineages["B"]
    require(
        Counter(branch_members) == Counter(expected_b_discovery),
        "B branch groups must partition discovery-eligible B populations exactly once",
    )
    require(
        set(b_lineage["broad_validation_populations"]) == validation_only & {
            population for population in selected if populations[population]["lineage"] == "B"
        },
        "B broad-validation declaration does not match role overrides",
    )

    local_neighbors = {
        target: set(neighbors) for target, neighbors in b_lineage["local_stage_neighbors"].items()
    }
    require(set(local_neighbors) <= expected_b_discovery, "Local-stage target is not a B discovery population")
    for target, neighbors in local_neighbors.items():
        require(bool(neighbors), f"Local-stage target lacks neighbors: {target}")
        require(target not in neighbors, f"Local-stage target lists itself: {target}")
        require(neighbors <= expected_b_discovery, f"Unknown local neighbor for {target}")
        for neighbor in neighbors:
            require(
                target in local_neighbors.get(neighbor, set()),
                f"Local-neighbor graph is not symmetric: {target} -> {neighbor}",
            )

    edge_ids: set[str] = set()
    for edge in b_lineage["directed_developmental_edges"]:
        require(edge["edge_id"] not in edge_ids, f"Duplicate edge ID: {edge['edge_id']}")
        edge_ids.add(edge["edge_id"])
        require(edge["source"] in local_neighbors, f"Unknown edge source: {edge['source']}")
        require(edge["target"] in local_neighbors, f"Unknown edge target: {edge['target']}")
        require(
            edge["target"] in local_neighbors[edge["source"]],
            f"Developmental edge is absent from local-neighbor graph: {edge['edge_id']}",
        )

    branch_groups = set(b_lineage["branch_groups"])
    contrast_ids: set[str] = set()
    for contrast in b_lineage["branch_contrasts"]:
        require(contrast["contrast_id"] not in contrast_ids, f"Duplicate contrast ID: {contrast['contrast_id']}")
        contrast_ids.add(contrast["contrast_id"])
        require(contrast["target_group"] in branch_groups, "Unknown branch target")
        require(set(contrast["competitor_groups"]) <= branch_groups, "Unknown branch competitor")
        require(contrast["target_group"] not in contrast["competitor_groups"], "Self branch contrast")

    summary = {
        "schema_version": hierarchy["schema_version"],
        "status": "passed",
        "selected_populations": len(selected),
        "broad_discovery_populations": len(discovery),
        "broad_validation_only_populations": sorted(validation_only),
        "retained_unique_branch_populations": sorted(unique_branch),
        "broad_lineage_population_counts": {
            lineage: len(members) for lineage, members in sorted(observed_lineages.items())
        },
        "b_branch_groups": len(b_lineage["branch_groups"]),
        "b_local_stage_targets": len(local_neighbors),
        "b_developmental_edges": len(edge_ids),
        "b_branch_contrasts": len(contrast_ids),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
