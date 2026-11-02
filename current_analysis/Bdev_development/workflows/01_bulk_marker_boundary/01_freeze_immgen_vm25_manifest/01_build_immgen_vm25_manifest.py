#!/usr/bin/env python3
"""Build and validate the frozen ImmGen GENCODE vM25 library manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import openpyxl


MANIFEST_NAME = "02_immgen_vm25_library_manifest.tsv"
SUMMARY_NAME = "03_immgen_vm25_manifest_summary.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", required=True, type=Path)
    parser.add_argument("--vm25-counts", required=True, type=Path)
    parser.add_argument("--vm25-readme", required=True, type=Path)
    parser.add_argument("--table-s1", required=True, type=Path)
    parser.add_argument("--rna-soft", required=True, type=Path)
    parser.add_argument("--atac-soft", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def clean_cell(value: Any) -> str:
    return "" if value is None else str(value).strip()


def load_policy(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        policy = json.load(handle)
    require(
        policy.get("schema_version") == "bdev.immgen_vm25_manifest_policy.v1",
        "Unsupported manifest policy schema",
    )
    return policy


def validate_source_checksums(args: argparse.Namespace, policy: dict[str, Any]) -> None:
    source_paths = {
        "vm25_count_matrix": args.vm25_counts,
        "vm25_readme": args.vm25_readme,
        "paper_table_s1": args.table_s1,
        "rna_geo_samples_soft": args.rna_soft,
        "atac_geo_samples_soft": args.atac_soft,
    }
    for source_id, path in source_paths.items():
        require(path.is_file(), f"Missing source file: {path}")
        observed = sha256(path)
        expected = policy["sources"][source_id]["sha256"]
        require(
            observed == expected,
            f"Checksum mismatch for {source_id}: expected {expected}, observed {observed}",
        )


def parse_soft_samples(path: Path) -> dict[str, dict[str, str]]:
    samples: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None

    def finish() -> None:
        nonlocal current
        if current is not None:
            samples.append(current)
            current = None

    with path.open(encoding="utf-8", errors="strict") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                finish()
                current = {"record_accession": line.split(" = ", 1)[1], "relations": []}
            elif current is None:
                continue
            elif line.startswith("!Sample_title = "):
                current["title"] = line.split(" = ", 1)[1]
            elif line.startswith("!Sample_geo_accession = "):
                current["gsm"] = line.split(" = ", 1)[1]
            elif line.startswith("!Sample_relation = "):
                current["relations"].append(line.split(" = ", 1)[1])
    finish()

    by_title: dict[str, dict[str, str]] = {}
    for sample in samples:
        require("title" in sample and "gsm" in sample, f"Incomplete SOFT sample: {sample}")
        require(sample["title"] not in by_title, f"Duplicate SOFT title: {sample['title']}")
        biosample = ""
        sra = ""
        for relation in sample["relations"]:
            if relation.startswith("BioSample:"):
                biosample = relation.rsplit("/", 1)[-1]
            elif relation.startswith("SRA:"):
                match = re.search(r"term=([^&]+)", relation)
                sra = match.group(1) if match else relation.rsplit("/", 1)[-1]
        by_title[sample["title"]] = {
            "title": sample["title"],
            "gsm": sample["gsm"],
            "biosample": biosample,
            "sra": sra,
        }
    return by_title


def load_table_s1(
    path: Path, atac_population_aliases: dict[str, str]
) -> tuple[
    dict[str, list[dict[str, Any]]],
    dict[str, list[dict[str, Any]]],
]:
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)

    rna_sheet = workbook["Read Statistics"]
    rna_rows = list(rna_sheet.iter_rows(values_only=True))
    rna_header = {name: index for index, name in enumerate(rna_rows[0])}
    rna_by_population: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rna_rows[1:]:
        population = clean_cell(row[rna_header["population.name"]])
        rna_by_population[population].append(
            {
                "lab": clean_cell(row[rna_header["ImmGen.lab.contributed"]]),
                "sample_name": clean_cell(row[rna_header["sample.name"]]),
            }
        )

    atac_sheet = workbook["SortedPopulations"]
    atac_rows = list(atac_sheet.iter_rows(values_only=True))
    atac_header = {name: index for index, name in enumerate(atac_rows[0])}
    atac_by_population: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in atac_rows[1:]:
        raw_population = clean_cell(row[atac_header["CellType"]])
        canonical = atac_population_aliases.get(raw_population, raw_population)
        atac_by_population[canonical].append(
            {
                "sample_name": clean_cell(row[atac_header["SampleName"]]),
                "raw_population": raw_population,
                "lineage": clean_cell(row[atac_header["Lineage"]]),
                "cell_family": clean_cell(row[atac_header["CellFamily"]]),
                "organ": clean_cell(row[atac_header["Organ"]]),
                "lab": clean_cell(row[atac_header["ImmGenLab"]]),
            }
        )
    return dict(rna_by_population), dict(atac_by_population)


def inspect_matrix(path: Path) -> tuple[list[str], dict[str, int], dict[str, int]]:
    gene_count = 0
    gene_symbols: set[str] = set()
    noninteger_values = 0
    negative_values = 0
    bad_width_rows = 0

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        require(header[0] == "GeneSymbol", "Unexpected first vM25 matrix column")
        libraries = header[1:]
        require(len(libraries) == len(set(libraries)), "Duplicate vM25 library columns")
        for row in reader:
            gene_count += 1
            if len(row) != len(header):
                bad_width_rows += 1
                continue
            gene_symbol = row[0]
            require(bool(gene_symbol), f"Empty gene symbol at matrix row {gene_count + 1}")
            require(gene_symbol not in gene_symbols, f"Duplicate gene symbol: {gene_symbol}")
            gene_symbols.add(gene_symbol)
            for value in row[1:]:
                try:
                    count = int(value)
                except ValueError:
                    noninteger_values += 1
                    continue
                if count < 0:
                    negative_values += 1

    require(bad_width_rows == 0, f"Matrix has {bad_width_rows} malformed rows")
    require(noninteger_values == 0, f"Matrix has {noninteger_values} noninteger values")
    require(negative_values == 0, f"Matrix has {negative_values} negative values")
    qc = {
        "gene_rows": gene_count,
        "unique_gene_symbols": len(gene_symbols),
        "library_columns": len(libraries),
        "bad_width_rows": bad_width_rows,
        "noninteger_values": noninteger_values,
        "negative_values": negative_values,
    }
    population_counts = Counter(library.rsplit("#", 1)[0] for library in libraries)
    return libraries, qc, dict(population_counts)


def unique_join(values: list[str]) -> str:
    return ";".join(sorted(set(values)))


def build_manifest(
    policy: dict[str, Any],
    libraries: list[str],
    matrix_population_counts: dict[str, int],
    rna_by_population: dict[str, list[dict[str, Any]]],
    atac_by_population: dict[str, list[dict[str, Any]]],
    rna_soft: dict[str, dict[str, str]],
    atac_soft: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    population_aliases = policy["vm25_population_aliases"]
    title_aliases = policy["rna_geo_title_aliases"]
    included = set(policy["included_populations"])
    excluded = policy["excluded_table_s1_populations"]
    table_rna_populations = set(rna_by_population)

    require(len(included) == policy["expected"]["selected_population_count"], "Bad allowlist size")
    require(included.isdisjoint(excluded), "Included and excluded populations overlap")
    require(
        included | set(excluded) == table_rna_populations,
        "Population policy does not partition all Table S1 RNA populations",
    )
    require(included <= set(atac_by_population), "A selected population lacks Table S1 ATAC data")

    atac_population_gsms: dict[str, list[str]] = {}
    for population, records in atac_by_population.items():
        gsms = []
        for record in records:
            title = record["sample_name"]
            require(title in atac_soft, f"Table S1 ATAC library lacks exact GEO title: {title}")
            gsms.append(atac_soft[title]["gsm"])
        atac_population_gsms[population] = gsms

    canonical_matrix_counts: Counter[str] = Counter()
    for raw_population, count in matrix_population_counts.items():
        canonical_matrix_counts[population_aliases.get(raw_population, raw_population)] += count

    rows: list[dict[str, str]] = []
    seen_rna_gsms: set[str] = set()
    for column_index, library in enumerate(libraries, start=1):
        require("#" in library, f"Matrix library lacks replicate delimiter: {library}")
        raw_population, replicate = library.rsplit("#", 1)
        require(replicate.isdigit(), f"Nonnumeric matrix replicate: {library}")
        canonical = population_aliases.get(raw_population, raw_population)
        population_mapping = "approved_alias" if raw_population in population_aliases else "exact"

        geo_stem = title_aliases.get(library, library)
        geo_title = f"{geo_stem}_RNA-seq"
        require(geo_title in rna_soft, f"vM25 library lacks GEO mapping: {library} -> {geo_title}")
        rna_sample = rna_soft[geo_title]
        require(rna_sample["gsm"] not in seen_rna_gsms, f"GSM reused by matrix: {rna_sample['gsm']}")
        seen_rna_gsms.add(rna_sample["gsm"])
        geo_mapping = "approved_title_alias" if library in title_aliases else "exact_title"

        in_table = canonical in table_rna_populations
        if canonical in included:
            role = "discovery"
            include = "true"
            exclusion_code = ""
        elif canonical in excluded:
            role = "excluded_table_s1"
            include = "false"
            exclusion_code = excluded[canonical]
        else:
            require(not in_table, f"Unclassified Table S1 population: {canonical}")
            role = "excluded_later_population"
            include = "false"
            exclusion_code = "not_in_paper_table_s1"

        atac_records = atac_by_population.get(canonical, [])
        atac_raw_populations = [record["raw_population"] for record in atac_records]
        lineages = [record["lineage"] for record in atac_records]
        cell_families = [record["cell_family"] for record in atac_records]
        organs = [record["organ"] for record in atac_records]
        atac_labs = [record["lab"] for record in atac_records]
        rna_labs = [record["lab"] for record in rna_by_population.get(canonical, [])]

        rows.append(
            {
                "matrix_column_index": str(column_index),
                "matrix_library_id": library,
                "matrix_population_label": raw_population,
                "matrix_replicate": replicate,
                "canonical_table_s1_population": canonical if in_table else "",
                "population_mapping_method": population_mapping if in_table else "not_in_table_s1",
                "rna_gsm": rna_sample["gsm"],
                "rna_geo_title": rna_sample["title"],
                "rna_geo_mapping_method": geo_mapping,
                "rna_biosample": rna_sample["biosample"],
                "rna_sra_experiment": rna_sample["sra"],
                "table_s1_rna_match": str(in_table).lower(),
                "table_s1_rna_labs": unique_join(rna_labs),
                "table_s1_rna_library_count": str(len(rna_by_population.get(canonical, []))) if in_table else "0",
                "vm25_population_library_count": str(canonical_matrix_counts[canonical]),
                "table_s1_lineage": unique_join(lineages),
                "table_s1_cell_family": unique_join(cell_families),
                "table_s1_organs": unique_join(organs),
                "table_s1_atac_population": unique_join(atac_raw_populations),
                "table_s1_atac_match": str(bool(atac_records)).lower(),
                "table_s1_atac_labs": unique_join(atac_labs),
                "table_s1_atac_library_count": str(len(atac_records)),
                "table_s1_atac_gsms": unique_join(atac_population_gsms.get(canonical, [])),
                "discovery_include": include,
                "discovery_role": role,
                "exclusion_code": exclusion_code,
            }
        )

    require(len(seen_rna_gsms) == len(libraries), "Matrix-to-GSM mapping is not one-to-one")
    return rows


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    require(bool(rows), "Cannot write an empty manifest")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            dialect="excel-tab",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def build_summary(
    policy: dict[str, Any],
    matrix_qc: dict[str, int],
    rows: list[dict[str, str]],
    manifest_path: Path,
) -> dict[str, Any]:
    selected = [row for row in rows if row["discovery_include"] == "true"]
    excluded = [row for row in rows if row["discovery_include"] == "false"]
    selected_populations = sorted({row["canonical_table_s1_population"] for row in selected})
    lineage_libraries = Counter(row["table_s1_lineage"] for row in selected)
    lineage_populations: dict[str, set[str]] = defaultdict(set)
    for row in selected:
        lineage_populations[row["table_s1_lineage"]].add(row["canonical_table_s1_population"])

    paper_counts = {
        row["canonical_table_s1_population"]: int(row["table_s1_rna_library_count"])
        for row in selected
    }
    vm25_counts = {
        row["canonical_table_s1_population"]: int(row["vm25_population_library_count"])
        for row in selected
    }
    replicate_differences = {
        population: {
            "paper_table_s1": paper_counts[population],
            "vm25": vm25_counts[population],
        }
        for population in selected_populations
        if paper_counts[population] != vm25_counts[population]
    }

    summary = {
        "schema_version": "bdev.immgen_vm25_manifest_summary.v1",
        "approved_on": policy["approved_on"],
        "source_checksums_validated": True,
        "matrix_annotation": "GENCODE vM25",
        "matrix_qc": matrix_qc,
        "manifest": {
            "file": MANIFEST_NAME,
            "sha256": sha256(manifest_path),
            "total_libraries": len(rows),
            "selected_libraries": len(selected),
            "excluded_libraries": len(excluded),
            "selected_populations": len(selected_populations),
            "selected_lineages": {
                lineage: {
                    "populations": len(lineage_populations[lineage]),
                    "libraries": lineage_libraries[lineage],
                }
                for lineage in sorted(lineage_libraries)
            },
            "excluded_roles": dict(sorted(Counter(row["discovery_role"] for row in excluded).items())),
            "excluded_reasons": dict(sorted(Counter(row["exclusion_code"] for row in excluded).items())),
            "rna_gsm_mapping_methods": dict(
                sorted(Counter(row["rna_geo_mapping_method"] for row in rows).items())
            ),
            "selected_population_replicate_differences": replicate_differences,
        },
        "acceptance_checks": {
            "all_matrix_values_are_nonnegative_integers": True,
            "all_matrix_columns_map_one_to_one_to_rna_gsm": True,
            "all_selected_populations_are_in_paper_table_s1": True,
            "all_selected_populations_have_population_matched_atac": all(
                row["table_s1_atac_match"] == "true" for row in selected
            ),
            "no_unapproved_population_selected": set(selected_populations)
            == set(policy["included_populations"]),
        },
        "status": "passed",
    }
    return summary


def validate_expected(
    policy: dict[str, Any],
    matrix_qc: dict[str, int],
    matrix_population_counts: dict[str, int],
    rows: list[dict[str, str]],
    rna_by_population: dict[str, list[dict[str, Any]]],
    atac_by_population: dict[str, list[dict[str, Any]]],
) -> None:
    expected = policy["expected"]
    selected = [row for row in rows if row["discovery_include"] == "true"]
    require(matrix_qc["gene_rows"] == expected["matrix_gene_rows"], "Unexpected gene-row count")
    require(matrix_qc["library_columns"] == expected["matrix_library_columns"], "Unexpected library count")
    require(len(matrix_population_counts) == expected["matrix_population_labels"], "Unexpected population-label count")
    require(len(rna_by_population) == expected["paper_rna_population_count"], "Unexpected Table S1 RNA population count")
    require(sum(map(len, rna_by_population.values())) == expected["paper_rna_library_count"], "Unexpected Table S1 RNA library count")
    require(len(atac_by_population) == expected["paper_atac_population_count"], "Unexpected Table S1 ATAC population count")
    require(sum(map(len, atac_by_population.values())) == expected["paper_atac_library_count"], "Unexpected Table S1 ATAC library count")
    require(len(selected) == expected["selected_library_count"], "Unexpected selected-library count")
    require(len(rows) - len(selected) == expected["excluded_library_count"], "Unexpected excluded-library count")
    require(
        len({row["canonical_table_s1_population"] for row in selected})
        == expected["selected_population_count"],
        "Unexpected selected-population count",
    )


def main() -> None:
    args = parse_args()
    policy = load_policy(args.policy)
    validate_source_checksums(args, policy)
    rna_soft = parse_soft_samples(args.rna_soft)
    atac_soft = parse_soft_samples(args.atac_soft)
    rna_by_population, atac_by_population = load_table_s1(
        args.table_s1, policy["atac_population_aliases"]
    )
    libraries, matrix_qc, matrix_population_counts = inspect_matrix(args.vm25_counts)
    rows = build_manifest(
        policy,
        libraries,
        matrix_population_counts,
        rna_by_population,
        atac_by_population,
        rna_soft,
        atac_soft,
    )
    validate_expected(
        policy,
        matrix_qc,
        matrix_population_counts,
        rows,
        rna_by_population,
        atac_by_population,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / MANIFEST_NAME
    summary_path = args.output_dir / SUMMARY_NAME
    write_manifest(manifest_path, rows)
    summary = build_summary(policy, matrix_qc, rows, manifest_path)
    require(all(summary["acceptance_checks"].values()), "An acceptance check failed")
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
