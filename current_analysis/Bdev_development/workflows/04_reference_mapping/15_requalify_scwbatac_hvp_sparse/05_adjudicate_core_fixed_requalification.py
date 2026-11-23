#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_ledger(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: str) -> float:
    if value in {"NA", "NaN", ""}:
        return math.nan
    return float(value)


def within(left: str, right: str, tolerance: float) -> bool:
    a, b = number(left), number(right)
    return (math.isnan(a) and math.isnan(b)) or abs(a - b) <= tolerance


def max_difference(current: list[dict[str, str]], prior_by_id: dict[str, dict[str, str]],
                   column: str) -> float:
    values = []
    for row in current:
        a, b = number(row[column]), number(prior_by_id[row["feature_id"]][column])
        if math.isnan(a) and math.isnan(b):
            continue
        values.append(abs(a - b))
    return max(values, default=0.0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    config_path = parser.parse_args().config.resolve()
    config = yaml.safe_load(config_path.read_text())
    acceptance = config["acceptance"]
    tolerance = float(acceptance["numerical_tolerance"])
    execution = config["execution"]
    execution_report_path = Path(execution["report"])
    execution_report = json.loads(execution_report_path.read_text())

    ledgers: dict[str, tuple[list[dict[str, str]], list[dict[str, str]]]] = {}
    for mode in ("hvp_depth20k", "hvp_raw20k"):
        ledgers[mode] = (
            read_ledger(Path(execution["current_ledgers"][mode])),
            read_ledger(Path(execution["prior_ledgers"][mode])),
        )

    checks: dict[str, object] = {}
    expected_features = int(acceptance["expected_features"])
    expected_selected = int(acceptance["expected_selected_per_mode"])
    common_columns = ["mean_accessibility", "detected_cells", "detection_pct",
                      "raw_variance", "residual_variance", "depth_slope",
                      "depth_r_squared", "score"]

    for mode, (current, prior) in ledgers.items():
        prior_by_id = {row["feature_id"]: row for row in prior}
        feature_axis_identical = ([row["feature_id"] for row in current] ==
                                  [row["feature_id"] for row in prior])
        feature_set_identical = (
            len(prior_by_id) == len(prior) == len(current) and
            set(prior_by_id) == {row["feature_id"] for row in current})
        selected_identical = all(
            row["selected"] == prior_by_id[row["feature_id"]]["selected"]
            for row in current)
        selected_count = sum(row["selected"] == "TRUE" for row in current)
        base = {
            "rows": len(current),
            "feature_axis_identical": feature_axis_identical,
            "feature_set_identical": feature_set_identical,
            "selected_membership_identical": selected_identical,
            "selected_count": selected_count,
            "method_identical": all(
                row["method"] == prior_by_id[row["feature_id"]]["method"]
                for row in current),
            "depth_adjustment_status_identical": all(
                row["depth_adjustment_status"] ==
                prior_by_id[row["feature_id"]]["depth_adjustment_status"]
                for row in current),
        }
        if mode == "hvp_depth20k":
            differences = {column: max_difference(current, prior_by_id, column)
                           for column in common_columns}
            base.update({
                "rank_identical": all(
                    row["rank"] == prior_by_id[row["feature_id"]]["rank"]
                    for row in current),
                "maximum_absolute_differences": differences,
                "complete_numerical_equivalence":
                    all(value <= tolerance for value in differences.values()),
            })
        else:
            score_columns = list(acceptance["raw_mode_score_columns"])
            differences = {column: max_difference(current, prior_by_id, column)
                           for column in score_columns}
            current_by_rank = {int(row["rank"]): row for row in current}
            prior_by_rank = {int(row["rank"]): row for row in prior}
            rank_mismatches = []
            tied_only = True
            for row in current:
                old = prior_by_id[row["feature_id"]]
                new_rank, old_rank = int(row["rank"]), int(old["rank"])
                if new_rank == old_rank:
                    continue
                tied_current = within(row["score"], current_by_rank[old_rank]["score"], tolerance)
                tied_prior = within(old["score"], prior_by_rank[new_rank]["score"], tolerance)
                tied_only = tied_only and tied_current and tied_prior
                rank_mismatches.append({
                    "feature_id": row["feature_id"],
                    "current_rank": new_rank,
                    "prior_rank": old_rank,
                    "rank_displacement": abs(new_rank - old_rank),
                })
            sentinels = {
                "residual_variance_equals_raw_variance": all(
                    within(row["residual_variance"], row["raw_variance"], tolerance)
                    for row in current),
                "depth_slope_zero": all(abs(number(row["depth_slope"])) <= tolerance
                                        for row in current),
                "depth_r_squared_zero": all(
                    abs(number(row["depth_r_squared"])) <= tolerance for row in current),
                "status_not_requested": all(
                    row["depth_adjustment_status"] == "not_requested" for row in current),
            }
            base.update({
                "maximum_absolute_differences_score_columns": differences,
                "score_columns_equivalent":
                    all(value <= tolerance for value in differences.values()),
                "rank_mismatch_count": len(rank_mismatches),
                "maximum_rank_displacement": max(
                    (item["rank_displacement"] for item in rank_mismatches), default=0),
                "rank_permutations_tied_only": tied_only,
                "selection_changes_among_rank_mismatches": sum(
                    row["selected"] != prior_by_id[row["feature_id"]]["selected"]
                    for row in current if row["rank"] !=
                    prior_by_id[row["feature_id"]]["rank"]),
                "canonical_not_requested_sentinels": sentinels,
            })
        checks[mode] = base

    depth = checks["hvp_depth20k"]
    raw = checks["hvp_raw20k"]
    source_manifest = Path(yaml.safe_load(
        (Path(config_path.parent) / Path(execution["config"]).name).read_text()
    )["inputs"]["source_bundle"]["path"]) / "manifest.json"
    source_checks = {
        "scientific_state_fingerprint_unchanged": bool(
            execution_report["source_scientific_state_unchanged"]),
        "publication_timestamp_only_change": bool(
            execution_report["source_manifest_timestamp_only_change"]),
        "current_manifest_matches_execution_after":
            sha256(source_manifest) == execution_report["source_manifest_after"],
    }

    gate = all([
        depth["rows"] == expected_features,
        depth["feature_axis_identical"], depth["selected_membership_identical"],
        depth["selected_count"] == expected_selected, depth["method_identical"],
        depth["depth_adjustment_status_identical"], depth["rank_identical"],
        depth["complete_numerical_equivalence"],
        raw["rows"] == expected_features,
        raw["feature_set_identical"], raw["selected_membership_identical"],
        raw["selected_count"] == expected_selected, raw["method_identical"],
        raw["depth_adjustment_status_identical"], raw["score_columns_equivalent"],
        raw["rank_permutations_tied_only"],
        raw["selection_changes_among_rank_mismatches"] == 0,
        all(raw["canonical_not_requested_sentinels"].values()),
        all(source_checks.values()),
    ])

    root = Path(execution["result_root"])
    outputs = config.get("outputs", {})
    report_path = root / "03_reports" / outputs.get(
        "report", "02_hvp_sparse_contract_adjudication.json")
    receipt_path = root / "06_manifests" / outputs.get(
        "receipt", "03_hvp_sparse_contract_adjudication_receipt.json")
    report = {
        "schema_version": "bdev.scwbatac_hvp_sparse_adjudication_report.v1",
        "status": "passed" if gate else "failed",
        "adjudication_identity": config["adjudication_identity"],
        "execution_run_identity": execution["run_identity"],
        "numerical_tolerance": tolerance,
        "checks": checks,
        "source_checks": source_checks,
        "implementation_unlocked": False,
        "sif_built": False,
        "gate_passed": gate,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    receipt = {
        "schema_version": "bdev.scwbatac_hvp_sparse_adjudication_receipt.v1",
        "status": "complete_passed" if gate else "complete_failed",
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "adjudication_script_sha256": sha256(Path(__file__).resolve()),
        "execution_report_sha256": sha256(execution_report_path),
        "adjudication_report_sha256": sha256(report_path),
        "gate_passed": gate,
    }
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    if not gate:
        raise RuntimeError("HVP sparse-contract adjudication failed")


if __name__ == "__main__":
    main()
