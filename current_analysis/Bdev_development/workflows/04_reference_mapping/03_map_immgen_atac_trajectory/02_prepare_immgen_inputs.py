#!/usr/bin/env python3
"""Resolve the frozen ImmGen GSM panel to exact ENA paired FASTQ objects."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def merge(base: dict, override: dict) -> dict:
    result = dict(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(path: Path) -> dict:
    current = yaml.safe_load(path.read_text())
    base_name = current.pop("base_config", None)
    if base_name is None:
        return current
    return merge(load_config((path.parent / base_name).resolve()), current)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def https_url(value: str) -> str:
    parsed = urlsplit(value if "://" in value else f"ftp://{value}")
    return urlunsplit(("https", parsed.netloc, parsed.path, "", ""))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--ena-report", type=Path, help="Optional cached ENA report")
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = load_config(config_path)
    inputs = config["inputs"]
    results = Path(config["storage"]["results_root"])
    for name in ("01_tables", "02_figures", "03_reports", "04_logs", "05_objects", "06_manifests"):
        (results / name).mkdir(parents=True, exist_ok=True)
    for key in ("raw_fastq", "peak_counts", "sessions"):
        Path(config["storage"][key]).mkdir(parents=True, exist_ok=True)

    population_path = Path(inputs["selected_populations"]["path"])
    library_path = Path(inputs["selected_libraries"]["path"])
    if sha256(population_path) != inputs["selected_populations"]["sha256"]:
        raise RuntimeError("Selected-population table checksum changed")
    if sha256(library_path) != inputs["selected_libraries"]["sha256"]:
        raise RuntimeError("Selected-library table checksum changed")
    populations = read_tsv(population_path)
    libraries = read_tsv(library_path)
    if len(populations) != inputs["selected_populations"]["expected_populations"]:
        raise RuntimeError("Unexpected selected-population count")
    if len(libraries) != inputs["selected_libraries"]["expected_libraries"]:
        raise RuntimeError("Unexpected selected-library count")

    if args.ena_report:
        report_path = args.ena_report.resolve()
    else:
        temporary = tempfile.NamedTemporaryFile(prefix="immgen_ena_", suffix=".tsv", delete=False)
        temporary.close()
        report_path = Path(temporary.name)
        subprocess.run(
            ["curl", "-L", "--fail", "--silent", "--show-error", "--retry", "5",
             "--output", str(report_path), inputs["ena"]["report_url"]],
            check=True,
        )
    ena_rows = read_tsv(report_path)
    by_gsm: dict[str, list[dict[str, str]]] = {}
    for row in ena_rows:
        by_gsm.setdefault(row["sample_alias"], []).append(row)

    output_rows = []
    for library in libraries:
        gsm = library["gsm"]
        matches = by_gsm.get(gsm, [])
        if len(matches) != 1:
            raise RuntimeError(f"{gsm} resolved to {len(matches)} ENA runs")
        ena = matches[0]
        urls = ena["fastq_ftp"].split(";")
        md5s = ena["fastq_md5"].split(";")
        sizes = ena["fastq_bytes"].split(";")
        if ena["library_layout"] != "PAIRED" or ena["library_strategy"] != "ATAC-seq":
            raise RuntimeError(f"{gsm} is not paired ATAC-seq")
        if not (len(urls) == len(md5s) == len(sizes) == 2):
            raise RuntimeError(f"{gsm} does not expose exactly two FASTQ objects")
        output_rows.append({
            **library,
            "run_accession": ena["run_accession"],
            "sample_accession": ena["sample_accession"],
            "experiment_accession": ena["experiment_accession"],
            "fastq_1_url": https_url(urls[0]),
            "fastq_2_url": https_url(urls[1]),
            "fastq_1_md5": md5s[0],
            "fastq_2_md5": md5s[1],
            "fastq_1_bytes": sizes[0],
            "fastq_2_bytes": sizes[1],
            "instrument_model": ena["instrument_model"],
        })

    table_path = results / "01_tables" / "01_immgen_library_run_manifest.tsv"
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    receipt = {
        "schema_version": "bdev.immgen_ena_resolution.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "ena_study": inputs["ena"]["study"],
        "selected_populations": len(populations),
        "selected_libraries": len(libraries),
        "resolved_paired_runs": len(output_rows),
        "compressed_fastq_bytes": sum(int(row["fastq_1_bytes"]) + int(row["fastq_2_bytes"]) for row in output_rows),
        "library_run_manifest": str(table_path),
        "library_run_manifest_sha256": sha256(table_path),
        "archive_modified": False,
    }
    (results / "06_manifests" / "02_ena_resolution.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
