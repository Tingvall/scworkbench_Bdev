#!/usr/bin/env python3
"""Prepare the isolated Step 4 bundle clone and freeze the ImmGen ATAC preflight."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())

    results = Path(config["storage"]["results_root"])
    for name in ("01_tables", "02_figures", "03_reports", "04_logs", "05_objects", "06_manifests"):
        (results / name).mkdir(parents=True, exist_ok=True)
    Path(config["storage"]["sessions"]).mkdir(parents=True, exist_ok=True)

    source = Path(config["inputs"]["source_bundle"]["path"])
    source_manifest = source / "manifest.json"
    source_sha = sha256(source_manifest)
    if source_sha != config["inputs"]["source_bundle"]["manifest_sha256"]:
        raise RuntimeError("Active 03c source-bundle manifest checksum changed")

    destination = Path(config["storage"]["execution_clone"])
    if destination.exists():
        raise RuntimeError(f"Refusing to overwrite execution clone: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["cp", "--recursive", "--reflink=auto", "--", str(source), str(destination)],
        check=True,
    )
    if sha256(destination / "manifest.json") != source_sha:
        raise RuntimeError("Execution-clone manifest differs before qualification")

    preflight = config["immgen_preflight"]
    manifest = Path(preflight["manifest"])
    if sha256(manifest) != preflight["manifest_sha256"]:
        raise RuntimeError("ImmGen VM25 library manifest checksum changed")
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = [
        row for row in rows
        if row["discovery_include"].lower() == "true"
        and row["table_s1_lineage"] in set(preflight["selected_lineages"])
        and row["table_s1_atac_match"].lower() == "true"
    ]
    populations: dict[str, dict[str, str]] = {}
    for row in selected:
        population = row["table_s1_atac_population"]
        gsms = sorted(set(filter(None, row["table_s1_atac_gsms"].split(";"))))
        record = {
            "population": population,
            "lineage": row["table_s1_lineage"],
            "cell_family": row["table_s1_cell_family"],
            "organs": row["table_s1_organs"],
            "atac_labs": row["table_s1_atac_labs"],
            "atac_library_count": str(len(gsms)),
            "atac_gsms": ";".join(gsms),
            "geo_series": preflight["geo_series"],
            "genome_build": preflight["genome_build"],
            "input_disposition": preflight["model_input_requirement"],
        }
        if population in populations and populations[population] != record:
            raise RuntimeError(f"Inconsistent repeated ImmGen population record: {population}")
        populations[population] = record
    selected_populations = [populations[key] for key in sorted(populations)]
    selected_libraries = sorted({gsm for row in selected_populations for gsm in row["atac_gsms"].split(";")})
    if len(selected_populations) != int(preflight["expected_populations"]):
        raise RuntimeError("Unexpected selected ImmGen ATAC population count")
    if len(selected_libraries) != int(preflight["expected_libraries"]):
        raise RuntimeError("Unexpected selected ImmGen ATAC library count")
    write_tsv(
        results / "01_tables" / "09_selected_immgen_atac_populations.tsv",
        selected_populations,
        list(selected_populations[0]),
    )
    library_rows = []
    for row in selected_populations:
        for gsm in row["atac_gsms"].split(";"):
            library_rows.append({
                "gsm": gsm,
                "population": row["population"],
                "lineage": row["lineage"],
                "geo_series": preflight["geo_series"],
                "genome_build": preflight["genome_build"],
                "required_source": "raw_SRA",
                "target_feature_universe": "exact_Bdev_peak_universe",
                "status": "accession_frozen_not_downloaded",
            })
    write_tsv(
        results / "01_tables" / "10_selected_immgen_atac_libraries.tsv",
        library_rows,
        list(library_rows[0]),
    )

    receipt = {
        "schema_version": "bdev.reference_mapping_preparation.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "source_bundle": str(source),
        "source_manifest_sha256": source_sha,
        "execution_clone": str(destination),
        "immgen_manifest": str(manifest),
        "immgen_manifest_sha256": sha256(manifest),
        "selected_populations": len(selected_populations),
        "selected_libraries": len(selected_libraries),
        "processed_bigwig_use": preflight["processed_bigwig_use"],
        "model_input_requirement": preflight["model_input_requirement"],
        "archive_accessed": False,
        "sif_built": False,
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
    }
    output = results / "06_manifests" / "01_preparation_and_immgen_preflight.json"
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
