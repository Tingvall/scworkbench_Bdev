#!/usr/bin/env python3
"""Freeze and validate ImmGen annotation, state-gene, and software inputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import sqlite3
import tarfile
import tempfile
from collections import defaultdict
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


def check_source(entry: dict, label: str) -> Path:
    path = Path(entry["path"])
    require(path.is_file(), f"Missing {label}: {path}")
    actual = sha256(path)
    require(actual == entry["sha256"], f"Checksum mismatch for {label}: {actual}")
    return path


def read_matrix_symbols(path: Path, expected_rows: int, expected_libraries: int) -> list[str]:
    symbols: list[str] = []
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        require(header[0] == "GeneSymbol", "Count matrix does not begin with GeneSymbol")
        require(len(header) - 1 == expected_libraries, "Unexpected matrix library count")
        for row in reader:
            require(len(row) == len(header), f"Bad count-matrix row width for {row[0]}")
            symbols.append(row[0])
    require(len(symbols) == expected_rows, "Unexpected matrix gene-row count")
    require(len(set(symbols)) == len(symbols), "Matrix gene symbols are not unique")
    return symbols


ATTR_RE = re.compile(r'(\w+) "([^"]*)";')


def read_gtf_genes(path: Path) -> dict[str, list[dict[str, str]]]:
    genes: dict[str, list[dict[str, str]]] = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(ATTR_RE.findall(fields[8]))
            genes[attrs["gene_name"]].append(
                {
                    "ensembl_gene_id": attrs["gene_id"],
                    "gene_biotype": attrs["gene_type"],
                    "chromosome": fields[0],
                }
            )
    return genes


def choose_mapping(records: list[dict[str, str]], reference: set[str]) -> tuple[str, dict[str, str] | None]:
    reference_records = [record for record in records if record["chromosome"] in reference]
    if len(reference_records) == 1:
        return "unique_reference", reference_records[0]
    if len(reference_records) > 1:
        return "ambiguous_reference", None
    if len(records) == 1:
        return "unique_nonreference", records[0]
    if len(records) > 1:
        return "ambiguous_nonreference", None
    return "unmapped", None


def extract_orgdb_sqlite(package: Path, destination: Path) -> None:
    member_name = "org.Mm.eg.db/inst/extdata/org.Mm.eg.sqlite"
    with tarfile.open(package, "r:gz") as archive:
        member = archive.getmember(member_name)
        source = archive.extractfile(member)
        require(source is not None, "Cannot read OrgDb SQLite member")
        destination.write_bytes(source.read())


def query_go_symbols(connection: sqlite3.Connection, ontology: str, roots: list[dict]) -> dict[str, set[str]]:
    require(ontology in {"bp", "cc", "mf"}, f"Unsupported GO ontology: {ontology}")
    table = f"go_{ontology}_all"
    result: dict[str, set[str]] = {}
    for root in roots:
        rows = connection.execute(
            f"SELECT DISTINCT gene_info.symbol FROM {table} "
            "JOIN gene_info USING (_id) WHERE go_id = ? ORDER BY gene_info.symbol",
            (root["id"],),
        )
        result[root["id"]] = {row[0] for row in rows}
    return result


def add_flag(flags: dict[str, list[dict[str, str]]], symbol: str, category: str, source_type: str, source_id: str) -> None:
    flags[symbol].append(
        {
            "category": category,
            "source_type": source_type,
            "source_id": source_id,
        }
    )


def write_tsv_gz(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as gzip_handle:
            with io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--software-lock", required=True, type=Path)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    software = yaml.safe_load(args.software_lock.read_text())
    require(config["schema_version"] == "bdev.immgen_analysis_preflight.v1", "Wrong preflight schema")
    require(software["schema_version"] == "bdev.immgen_bulk_analysis_software_lock.v1", "Wrong software schema")
    require(software["qualification"]["immgen_model_run"] is False, "Software lock claims an ImmGen model run")

    checked_inputs: dict[str, dict[str, str | int]] = {}
    for label, entry in config["sources"].items():
        path = check_source(entry, label)
        checked_inputs[label] = {"path": str(path), "sha256": entry["sha256"], "bytes": path.stat().st_size}
    for label, entry in software["source_packages"].items():
        path = check_source(entry, f"software source {label}")
        checked_inputs[f"software_{label}"] = {"path": str(path), "sha256": entry["sha256"], "bytes": path.stat().st_size}
    sif = check_source(software["base_sif"], "base SIF")
    checked_inputs["base_sif"] = {"path": str(sif), "sha256": software["base_sif"]["sha256"], "bytes": sif.stat().st_size}
    library_archive = Path(software["installed_library"]["archive"])
    require(library_archive.is_file(), f"Missing installed library archive: {library_archive}")
    require(sha256(library_archive) == software["installed_library"]["archive_sha256"], "Installed library archive checksum mismatch")
    checked_inputs["installed_library_archive"] = {
        "path": str(library_archive),
        "sha256": software["installed_library"]["archive_sha256"],
        "bytes": library_archive.stat().st_size,
    }

    matrix = Path(config["sources"]["count_matrix"]["path"])
    symbols = read_matrix_symbols(
        matrix,
        config["sources"]["count_matrix"]["expected_gene_rows"],
        config["sources"]["count_matrix"]["expected_library_columns"],
    )
    symbol_set = set(symbols)
    gtf_records = read_gtf_genes(Path(config["sources"]["gencode_all_annotation"]["path"]))
    reference = set(config["annotation_mapping"]["reference_chromosomes"])

    annotation_rows: list[dict] = []
    mapping_by_symbol: dict[str, dict] = {}
    for symbol in symbols:
        status, selected = choose_mapping(gtf_records.get(symbol, []), reference)
        row = {
            "gene_symbol": symbol,
            "mapping_status": status,
            "ensembl_gene_id": selected["ensembl_gene_id"] if selected else "",
            "chromosome": selected["chromosome"] if selected else "",
            "gene_biotype": selected["gene_biotype"] if selected else "",
            "gtf_record_count": len(gtf_records.get(symbol, [])),
            "state_flag": False,
            "mandatory_identity_exclusion": False,
            "core_identity_eligible": False,
            "exclusion_reasons": "",
        }
        annotation_rows.append(row)
        mapping_by_symbol[symbol] = row

    flags: dict[str, list[dict[str, str]]] = defaultdict(list)
    identity_flags: dict[str, list[dict[str, str]]] = defaultdict(list)
    with tempfile.TemporaryDirectory(prefix="bdev_orgdb_") as temp_dir:
        sqlite_path = Path(temp_dir) / "org.Mm.eg.sqlite"
        extract_orgdb_sqlite(Path(config["sources"]["mouse_orgdb"]["path"]), sqlite_path)
        connection = sqlite3.connect(sqlite_path)
        for category, definition in config["state_gene_policy"]["go_sets"].items():
            for go_id, go_symbols in query_go_symbols(connection, definition["ontology"], definition["roots"]).items():
                for symbol in sorted(go_symbols & symbol_set):
                    add_flag(flags, symbol, category, "GO_ALL", go_id)
        for category, definition in config["state_gene_policy"]["curated_sets"].items():
            for symbol in definition["symbols"]:
                require(symbol in symbol_set, f"Curated state gene is absent from matrix: {symbol}")
                add_flag(flags, symbol, category, "curated", definition["doi"])
        for category, definition in config["mandatory_identity_exclusions"]["go_sets"].items():
            for go_id, go_symbols in query_go_symbols(connection, definition["ontology"], definition["roots"]).items():
                for symbol in sorted(go_symbols & symbol_set):
                    add_flag(identity_flags, symbol, category, "GO_ALL", go_id)
        metadata = dict(connection.execute("SELECT name, value FROM metadata"))
        connection.close()

    require(metadata["GOSOURCEDATE"] == config["sources"]["mouse_orgdb"]["go_source_date"], "OrgDb GO source date changed")
    require(metadata["EGSOURCEDATE"] == config["sources"]["mouse_orgdb"]["entrez_source_date"], "OrgDb Entrez source date changed")

    excluded_biotypes = set(config["mandatory_identity_exclusions"]["gencode_biotypes"])
    for row in annotation_rows:
        symbol = row["gene_symbol"]
        if row["chromosome"] == "chrY":
            add_flag(identity_flags, symbol, "chromosome_Y_genes", "GENCODE", "chrY")
        if symbol == "Xist":
            add_flag(identity_flags, symbol, "Xist", "GENCODE", "symbol")
        if row["chromosome"] == "chrM":
            add_flag(identity_flags, symbol, "mitochondrial_genes", "GENCODE", "chrM")
        if row["gene_biotype"] in excluded_biotypes:
            add_flag(identity_flags, symbol, "immunoglobulin_or_tcr_vdj", "GENCODE", row["gene_biotype"])

        state_reasons = sorted({flag["category"] for flag in flags.get(symbol, [])})
        identity_reasons = sorted({flag["category"] for flag in identity_flags.get(symbol, [])})
        reasons: list[str] = []
        if row["mapping_status"].startswith("ambiguous") or row["mapping_status"] == "unmapped":
            reasons.append(row["mapping_status"])
        if row["gene_biotype"] != config["annotation_mapping"]["accepted_export_biotype"]:
            reasons.append("non_protein_coding")
        reasons.extend(f"state:{reason}" for reason in state_reasons)
        reasons.extend(f"identity:{reason}" for reason in identity_reasons)
        row["state_flag"] = bool(state_reasons)
        row["mandatory_identity_exclusion"] = bool(identity_reasons)
        row["core_identity_eligible"] = not reasons
        row["exclusion_reasons"] = ";".join(reasons)

    flag_rows: list[dict] = []
    for symbol in sorted(set(flags) | set(identity_flags)):
        for flag_class, collection in (("state", flags), ("identity", identity_flags)):
            for flag in sorted(collection.get(symbol, []), key=lambda item: (item["category"], item["source_type"], item["source_id"])):
                flag_rows.append(
                    {
                        "gene_symbol": symbol,
                        "flag_class": flag_class,
                        "category": flag["category"],
                        "source_type": flag["source_type"],
                        "source_id": flag["source_id"],
                    }
                )

    output_root = Path(config["outputs"]["root"])
    output_paths = {name: output_root / relative for name, relative in config["outputs"].items() if name != "root"}
    for path in output_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)

    annotation_fields = [
        "gene_symbol", "mapping_status", "ensembl_gene_id", "chromosome", "gene_biotype",
        "gtf_record_count", "state_flag", "mandatory_identity_exclusion", "core_identity_eligible",
        "exclusion_reasons",
    ]
    write_tsv_gz(output_paths["gene_annotation"], annotation_rows, annotation_fields)
    write_tsv(output_paths["state_gene_flags"], flag_rows, ["gene_symbol", "flag_class", "category", "source_type", "source_id"])

    mapping_counts: dict[str, int] = defaultdict(int)
    for row in annotation_rows:
        mapping_counts[row["mapping_status"]] += 1
    state_counts = {
        category: len({symbol for symbol, entries in flags.items() if any(entry["category"] == category for entry in entries)})
        for category in sorted({entry["category"] for entries in flags.values() for entry in entries})
    }
    identity_counts = {
        category: len({symbol for symbol, entries in identity_flags.items() if any(entry["category"] == category for entry in entries)})
        for category in sorted({entry["category"] for entries in identity_flags.values() for entry in entries})
    }
    summary = {
        "schema_version": "bdev.immgen_analysis_preflight_summary.v1",
        "status": "passed",
        "matrix_gene_rows": len(symbols),
        "mapping_counts": dict(sorted(mapping_counts.items())),
        "state_flagged_unique_genes": len(flags),
        "state_category_counts": state_counts,
        "identity_excluded_unique_genes": len(identity_flags),
        "identity_category_counts": identity_counts,
        "core_identity_eligible_genes": sum(bool(row["core_identity_eligible"]) for row in annotation_rows),
        "software": {
            "r_version": software["base_sif"]["r_version"],
            "bioconductor_version": software["base_sif"]["bioconductor_version"],
            "edgeR": software["installed_library"]["packages"]["edgeR"],
            "limma": software["base_sif"]["inherited_packages"]["limma"],
            "locfit": software["installed_library"]["packages"]["locfit"],
            "qualification": software["qualification"],
        },
        "immgen_model_run": False,
        "scworkbench_change_required": False,
    }
    output_paths["preflight_summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    generated = {}
    for name in ("gene_annotation", "state_gene_flags", "preflight_summary"):
        path = output_paths[name]
        generated[name] = {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}
    manifest = {
        "schema_version": "bdev.immgen_analysis_preflight_run_manifest.v1",
        "status": "passed",
        "configuration": {"path": str(args.config), "sha256": sha256(args.config)},
        "software_lock": {"path": str(args.software_lock), "sha256": sha256(args.software_lock)},
        "inputs": checked_inputs,
        "outputs": generated,
        "immgen_model_run": False,
    }
    output_paths["run_manifest"].write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
