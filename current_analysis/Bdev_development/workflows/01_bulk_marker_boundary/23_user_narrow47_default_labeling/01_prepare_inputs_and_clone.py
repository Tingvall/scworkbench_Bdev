#!/usr/bin/env python3
"""Validate the user-supplied 47x20 matrix, emit GMT, and clone the bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd
import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())

    source_matrix = Path(config["inputs"]["gene_sets"]["path"])
    if sha256(source_matrix) != config["inputs"]["gene_sets"]["sha256"]:
        raise RuntimeError("Gene-set input checksum changed")
    matrix = pd.read_csv(source_matrix, sep="\t", dtype=str)
    expected_shape = (
        int(config["inputs"]["gene_sets"]["genes_per_population"]),
        int(config["inputs"]["gene_sets"]["populations"]),
    )
    if matrix.shape != expected_shape or matrix.isna().any().any():
        raise RuntimeError(f"Expected complete {expected_shape} gene-set matrix; observed {matrix.shape}")
    if matrix.columns.duplicated().any():
        raise RuntimeError("Population names are duplicated")
    if any(matrix[column].duplicated().any() for column in matrix.columns):
        raise RuntimeError("A population contains duplicate genes")

    results_root = Path(config["storage"]["results_root"])
    tables = results_root / "01_tables"
    reports = results_root / "03_reports"
    manifests = results_root / "06_manifests"
    for directory in (tables, reports, manifests):
        directory.mkdir(parents=True, exist_ok=True)

    gmt_path = tables / "01_user_narrow47_genesets.gmt"
    with gmt_path.open("w", encoding="utf-8") as handle:
        for population in matrix.columns:
            handle.write("\t".join([population, "user-supplied narrow population set", *matrix[population]]) + "\n")
    membership = matrix.melt(var_name="population", value_name="gene_symbol")
    membership["rank"] = membership.groupby("population", sort=False).cumcount() + 1
    membership.to_csv(tables / "02_gene_set_membership.tsv", sep="\t", index=False)

    source_bundle = Path(config["inputs"]["source_bundle"]["path"])
    if sha256(source_bundle / "manifest.json") != config["inputs"]["source_bundle"]["manifest_sha256"]:
        raise RuntimeError("Source bundle manifest checksum changed")
    clone = Path(config["storage"]["intermediate_root"]) / "05_objects" / "01_execution_bundle"
    if clone.exists():
        raise RuntimeError(f"Refusing to overwrite execution clone: {clone}")
    clone.parent.mkdir(parents=True, exist_ok=True)
    # Paired storage does not support POSIX permission preservation, so copy
    # bundle content and links without requesting unsupported metadata writes.
    subprocess.run(["cp", "--recursive", "--reflink=auto", "--", str(source_bundle), str(clone)], check=True)
    state = json.loads((clone / "manifest.json").read_text())
    if state["current_state"]["fingerprint"] != config["inputs"]["source_bundle"]["fingerprint"]:
        raise RuntimeError("Execution clone fingerprint changed")

    receipt = {
        "schema_version": "bdev.user_narrow47_preparation.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "source_matrix": str(source_matrix),
        "source_matrix_sha256": sha256(source_matrix),
        "populations": matrix.shape[1],
        "genes_per_population": matrix.shape[0],
        "gmt": str(gmt_path),
        "gmt_sha256": sha256(gmt_path),
        "source_bundle": str(source_bundle),
        "execution_clone": str(clone),
        "initial_bundle_fingerprint": state["current_state"]["fingerprint"],
    }
    (manifests / "01_preparation_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
