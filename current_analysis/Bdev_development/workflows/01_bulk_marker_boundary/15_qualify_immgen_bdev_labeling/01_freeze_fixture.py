#!/usr/bin/env python3
"""Freeze a deterministic Bdev acceptance-cell fixture without mutating it."""

from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
import struct
from pathlib import Path

import yaml


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pt_digest(rows: list[tuple[str, float | None]]) -> str:
    digest = hashlib.sha256()
    for cell_id, value in rows:
        encoded = cell_id.encode("utf-8")
        digest.update(struct.pack(">I", len(encoded)))
        digest.update(encoded)
        if value is None:
            digest.update(b"\x00")
        else:
            digest.update(b"\x01")
            digest.update(struct.pack(">d", float(value)))
    return digest.hexdigest()


def write_tsv(path: Path, header: list[str], rows: list[tuple[object, ...]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config = yaml.safe_load(config_path.read_text())
    inputs_path = Path(config["inputs_manifest"])
    inputs = yaml.safe_load(inputs_path.read_text())
    results_root = Path(config["storage"]["results_root"])
    table_dir = results_root / "01_tables"
    report_dir = results_root / "03_reports"
    manifest_dir = results_root / "06_manifests"
    for directory in (table_dir, report_dir, manifest_dir):
        directory.mkdir(parents=True, exist_ok=True)

    source_bundle = Path(inputs["bdev_reference"]["bundle_path"])
    manifest_path = source_bundle / "manifest.json"
    database_path = source_bundle / "db" / "app.sqlite"
    gmt_path = Path(inputs["immgen_collection"]["candidate_gmt"])
    target_summary_path = Path(inputs["immgen_collection"]["target_summary"])
    checks = {
        "source_bundle_exists": source_bundle.is_dir(),
        "source_manifest_sha256": sha256_file(manifest_path)
        == inputs["bdev_reference"]["manifest_sha256"],
        "candidate_gmt_sha256": sha256_file(gmt_path)
        == inputs["immgen_collection"]["candidate_gmt_sha256"],
        "target_summary_sha256": sha256_file(target_summary_path)
        == inputs["immgen_collection"]["target_summary_sha256"],
    }
    if not all(checks.values()):
        raise RuntimeError(f"Input identity check failed: {checks}")

    connection = sqlite3.connect(
        f"file:{database_path}?mode=ro&immutable=1", uri=True
    )
    cluster_field = config["fixture"]["cluster_field"]
    pt_field = config["fixture"]["pseudotime_field"]
    query = (
        f'SELECT cell_id, "{cluster_field}", "{pt_field}" '
        "FROM cells ORDER BY cell_id"
    )
    source_rows = connection.execute(query).fetchall()
    schema_version = connection.execute(
        "SELECT MAX(version) FROM bundle_schema_migrations"
    ).fetchone()[0]
    connection.close()

    expected_cells = int(inputs["bdev_reference"]["expected_cells"])
    expected_defined = int(
        inputs["bdev_reference"]["authoritative_pseudotime"]["defined_cells"]
    )
    expected_missing = int(
        inputs["bdev_reference"]["authoritative_pseudotime"]["missing_cells"]
    )
    observed_defined = sum(value is not None for _, _, value in source_rows)
    observed_missing = sum(value is None for _, _, value in source_rows)
    if (len(source_rows), observed_defined, observed_missing) != (
        expected_cells,
        expected_defined,
        expected_missing,
    ):
        raise RuntimeError("Bdev cell or pseudotime counts do not match the frozen input")

    seed = str(config["random_seed"])
    cap = int(config["fixture"]["maximum_cells_per_stratum"])
    strata: dict[tuple[str, str], list[tuple[str, float | None]]] = {}
    for cell_id, cluster, pt in source_rows:
        key = (str(cluster), "defined" if pt is not None else "missing")
        strata.setdefault(key, []).append((cell_id, pt))

    selected: list[tuple[str, str, str, int, str]] = []
    stratum_summary: list[tuple[object, ...]] = []
    for (cluster, pt_status), members in sorted(strata.items()):
        ranked = sorted(
            members,
            key=lambda item: hashlib.sha256(
                f"{seed}\0{item[0]}".encode("utf-8")
            ).hexdigest(),
        )
        keep = ranked[:cap]
        for rank, (cell_id, _) in enumerate(keep, start=1):
            selected.append(
                (
                    cell_id,
                    cluster,
                    pt_status,
                    rank,
                    hashlib.sha256(f"{seed}\0{cell_id}".encode()).hexdigest(),
                )
            )
        stratum_summary.append(
            (cluster, pt_status, len(members), len(keep), cap)
        )
    selected.sort(key=lambda row: row[0])
    expected_fixture = int(config["fixture"]["expected_cells"])
    if len(selected) != expected_fixture:
        raise RuntimeError(
            f"Expected {expected_fixture} fixture cells, observed {len(selected)}"
        )

    fixture_path = table_dir / "01_fixture_cells.tsv"
    strata_path = table_dir / "02_fixture_strata.tsv"
    write_tsv(
        fixture_path,
        ["cell_id", "cluster", "pt_status", "selection_rank", "selection_sha256"],
        selected,
    )
    write_tsv(
        strata_path,
        ["cluster", "pt_status", "source_cells", "selected_cells", "cap"],
        stratum_summary,
    )

    source_pt_sha256 = pt_digest([(cell, pt) for cell, _, pt in source_rows])
    report = {
        "schema_version": "bdev.immgen_bdev_fixture_preflight.v1",
        "status": "passed",
        "run_identity": config["run_identity"],
        "source_schema_version": schema_version,
        "source_cells": len(source_rows),
        "source_pt_defined": observed_defined,
        "source_pt_missing": observed_missing,
        "source_pt_binary_sha256": source_pt_sha256,
        "fixture_cells": len(selected),
        "fixture_strata": len(strata),
        "checks": checks,
    }
    (report_dir / "01_fixture_preflight.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    manifest = {
        "schema_version": "bdev.immgen_bdev_labeling_input_manifest.v1",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256_file(config_path)},
        "inputs_manifest": {
            "path": str(inputs_path),
            "sha256": sha256_file(inputs_path),
        },
        "source_bundle": {
            "path": str(source_bundle),
            "accepted_state_fingerprint": inputs["bdev_reference"][
                "accepted_state_fingerprint"
            ],
            "manifest_sha256": sha256_file(manifest_path),
            "pt_binary_sha256": source_pt_sha256,
        },
        "candidate_gmt": {
            "path": str(gmt_path),
            "sha256": sha256_file(gmt_path),
            "sets": sum(1 for _ in gmt_path.open(encoding="utf-8")),
        },
        "fixture": {
            "path": str(fixture_path),
            "sha256": sha256_file(fixture_path),
            "cells": len(selected),
        },
    }
    (manifest_dir / "01_input_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()

