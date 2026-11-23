#!/usr/bin/env python3
"""Export a pinned GMT as one row per set with comma-separated genes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_GMT_SHA256 = "22f041528a1e882c0bec7bd920b7638646b20589020669bf0f42b59880e0bdb7"
EXPECTED_SETS = 38
EXPECTED_MEMBERSHIPS = 934


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmt", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run-identity", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    actual_gmt_sha256 = sha256(args.gmt)
    if actual_gmt_sha256 != EXPECTED_GMT_SHA256:
        raise ValueError(
            f"GMT checksum mismatch: expected {EXPECTED_GMT_SHA256}, "
            f"observed {actual_gmt_sha256}"
        )

    records: list[tuple[str, str, list[str]]] = []
    seen_names: set[str] = set()
    with args.gmt.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"GMT line {line_number} has fewer than three fields")
            name, description, *genes = fields
            if not name or name in seen_names:
                raise ValueError(f"GMT line {line_number} has a missing or duplicate set name")
            if any(not gene for gene in genes):
                raise ValueError(f"GMT line {line_number} contains an empty gene")
            if any("," in gene for gene in genes):
                raise ValueError(f"GMT line {line_number} contains a gene with a comma")
            if len(genes) != len(set(genes)):
                raise ValueError(f"GMT line {line_number} contains duplicate genes")
            seen_names.add(name)
            records.append((name, description, genes))

    membership_count = sum(len(genes) for _, _, genes in records)
    if len(records) != EXPECTED_SETS or membership_count != EXPECTED_MEMBERSHIPS:
        raise ValueError(
            f"Unexpected collection dimensions: {len(records)} sets and "
            f"{membership_count} memberships"
        )

    tables = args.output / "01_tables"
    manifests = args.output / "06_manifests"
    tables.mkdir(parents=True, exist_ok=False)
    manifests.mkdir(parents=True, exist_ok=False)

    export_path = tables / "01_parse38_gene_sets_comma_separated.tsv"
    with export_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["set_name", "description", "n_genes", "genes_comma_separated"])
        for name, description, genes in records:
            writer.writerow([name, description, len(genes), ",".join(genes)])

    manifest_path = manifests / "01_export_manifest.json"
    manifest = {
        "schema_version": "bdev.parse38_comma_separated_export.v1",
        "run_identity": args.run_identity,
        "created_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source": {
            "path": str(args.gmt),
            "sha256": actual_gmt_sha256,
            "format": "GMT",
        },
        "export": {
            "path": str(export_path),
            "sha256": sha256(export_path),
            "format": "TSV",
            "gene_delimiter": ",",
            "set_count": len(records),
            "membership_count": membership_count,
            "set_and_gene_order": "preserved_from_source_gmt",
        },
    }
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")

    print(export_path)
    print(f"sets={len(records)} memberships={membership_count}")
    print(f"sha256={manifest['export']['sha256']}")


if __name__ == "__main__":
    main()
