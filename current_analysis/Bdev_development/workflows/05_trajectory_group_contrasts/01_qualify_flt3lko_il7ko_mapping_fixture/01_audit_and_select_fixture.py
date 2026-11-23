#!/usr/bin/env python3
"""Audit immutable KO inputs and freeze a bounded four-library fixture."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED_SHA256 = {
    "qc_metadata": "522c6f2665c192b1d14b20a2f79e0ec3a09f02e76a70c1208c93595eef2019ed",
    "selected_axis": "aeb453c86852d97b6d24a947f3b4e97d3a7954eec3d07b2821f107eab26f126a",
    "forward_chain": "5e7589d8374e299706b5ef0ef06dbbc5818f3bd0aef72003f2d3b767225af9a1",
}

EXPECTED_LIBRARIES = {
    "FLT3LKO_BM_WT": ("Flt3l", "WT"),
    "FLT3LKO_BM_KO": ("Flt3l", "KO"),
    "IL7KO_BM_WT": ("Il7", "WT"),
    "IL7KO_BM_KO": ("Il7", "KO"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def read_selected_axis(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t")
                if row["candidate"] == "hvp_depth20k"]
    if len(rows) != 20000:
        raise RuntimeError(f"Expected 20,000 hvp_depth20k rows, found {len(rows)}")
    order = [int(row["feature_order"]) for row in rows]
    features = [row["feature_id"] for row in rows]
    if order != list(range(1, 20001)) or len(set(features)) != 20000:
        raise RuntimeError("Selected HVP axis is not unique and exactly ordered")
    return rows


def read_fixture(qc_path: Path, cells_per_library: int) -> tuple[list[dict[str, str]], Counter]:
    selected: dict[str, list[dict[str, str]]] = defaultdict(list)
    totals: Counter = Counter()
    with gzip.open(qc_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["qc_use_downstream"].upper() != "TRUE":
                continue
            library = row["library_id"]
            if library not in EXPECTED_LIBRARIES:
                raise RuntimeError(f"Unexpected QC-selected library: {library}")
            target, genotype = EXPECTED_LIBRARIES[library]
            if (row["knockout_target"], row["genotype"]) != (target, genotype):
                raise RuntimeError(f"Identity mismatch for {library}")
            totals[library] += 1
            selected[library].append(row)
    fixture: list[dict[str, str]] = []
    for library in EXPECTED_LIBRARIES:
        rows = sorted(selected[library], key=lambda row: row["cell_id"])
        if len(rows) < cells_per_library:
            raise RuntimeError(f"Insufficient QC-selected cells for {library}")
        fixture.extend(rows[:cells_per_library])
    if len({row["cell_id"] for row in fixture}) != len(fixture):
        raise RuntimeError("Fixture cell IDs are not unique")
    return fixture, totals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc-metadata", type=Path, required=True)
    parser.add_argument("--selected-axis", type=Path, required=True)
    parser.add_argument("--forward-chain", type=Path, required=True)
    parser.add_argument("--legacy-row-names", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cells-per-library", type=int, default=100)
    args = parser.parse_args()

    if args.output.exists():
        raise RuntimeError(f"Refusing to overwrite existing output: {args.output}")
    for name, path in (("qc_metadata", args.qc_metadata),
                       ("selected_axis", args.selected_axis),
                       ("forward_chain", args.forward_chain)):
        observed = sha256(path)
        if observed != EXPECTED_SHA256[name]:
            raise RuntimeError(f"Checksum mismatch for {name}: {observed}")

    axis = read_selected_axis(args.selected_axis)
    fixture, source_counts = read_fixture(args.qc_metadata, args.cells_per_library)
    legacy_features = set(args.legacy_row_names.read_text(encoding="utf-8").splitlines())
    overlap = sum(row["feature_id"] in legacy_features for row in axis)
    if overlap != 0:
        raise RuntimeError(f"Expected zero legacy/current-axis overlap, observed {overlap}")

    tables = args.output / "01_tables"
    reports = args.output / "03_reports"
    manifests = args.output / "06_manifests"
    for directory in (tables, reports, manifests):
        directory.mkdir(parents=True, exist_ok=False)

    fixture_columns = [
        "cell_id", "library_id", "knockout_target", "genotype", "tissue",
        "rna_umi", "detected_genes", "mitochondrial_percent",
        "atac_fragments", "tss_enrichment", "nucleosome_signal",
    ]
    fixture_path = tables / "01_fixture_cells.tsv"
    with fixture_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fixture_columns, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(fixture)

    axis_path = tables / "02_hvp_depth20k_axis.tsv"
    with axis_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle,
                                fieldnames=["feature_order", "feature_id", "selection_rank"],
                                delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(axis)

    report = {
        "status": "passed",
        "dataset_identity": "Flt3l_ligand_KO_and_Il7_KO_with_corresponding_WT",
        "source_qc_selected_cells": dict(source_counts),
        "fixture_cells": len(fixture),
        "fixture_cells_per_library": args.cells_per_library,
        "selected_model_features": len(axis),
        "legacy_axis_features": len(legacy_features),
        "legacy_selected_axis_overlap": overlap,
        "next_action": "reciprocal_liftover_then_fragment_quantification",
        "paper_analysis": False,
    }
    atomic_text(reports / "01_input_audit.json",
                json.dumps(report, indent=2, sort_keys=True) + "\n")

    inventory = {
        "inputs": {
            "qc_metadata": {"path": str(args.qc_metadata), "sha256": EXPECTED_SHA256["qc_metadata"]},
            "selected_axis": {"path": str(args.selected_axis), "sha256": EXPECTED_SHA256["selected_axis"]},
            "forward_chain": {"path": str(args.forward_chain), "sha256": EXPECTED_SHA256["forward_chain"]},
            "legacy_row_names": {"path": str(args.legacy_row_names), "sha256": sha256(args.legacy_row_names)},
        },
        "outputs": {
            "fixture_cells": {"path": str(fixture_path), "sha256": sha256(fixture_path)},
            "selected_axis": {"path": str(axis_path), "sha256": sha256(axis_path)},
        },
    }
    atomic_text(manifests / "01_input_inventory.json",
                json.dumps(inventory, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
