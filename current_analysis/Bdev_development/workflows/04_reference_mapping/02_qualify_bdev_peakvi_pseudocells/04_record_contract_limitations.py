#!/usr/bin/env python3
"""Record the contract boundary around the completed 2D scientific gate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

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


def load_config(path: Path) -> tuple[dict, list[Path]]:
    current = yaml.safe_load(path.read_text())
    base_name = current.pop("base_config", None)
    if base_name is None:
        return current, [path]
    base, dependencies = load_config((path.parent / base_name).resolve())
    return merge(base, current), dependencies + [path]


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config, dependencies = load_config(config_path)
    root = Path(config["storage"]["results_root"])
    summary_path = root / "03_reports" / "01_peakvi_pseudocell_qualification_summary.json"
    receipt_path = root / "06_manifests" / "02_execution_receipt.json"
    summary = json.loads(summary_path.read_text())
    receipt = json.loads(receipt_path.read_text())
    dependency_records = [
        {"path": str(path), "sha256": sha256(path)} for path in dependencies
    ]
    limitations = {
        "higher_dimensional_projection": {
            "status": "unqualified",
            "failed_run": "bdev_peakvi_pseudocell_20261110T054149Z",
            "fitted_dimensions": 10,
            "granted_projection_dimensions": 2,
            "failure": "Output space dimensions did not match the current peakvi_project grant.",
        },
        "model_bound_latent_union": {
            "status": "unqualified",
            "failure": "Complete typed-space union source changed after authentication.",
            "fallback": "Exact row-wise attachment of the retained fitted and projected spaces.",
            "exact_coordinate_identity_verified": True,
            "rows": 1938,
            "dimensions": 2,
        },
        "duplicate_fit_artifacts": {
            "status": "preserved_not_deleted",
            "reason": "A continuation reran the deterministic 2D fit before registry-based recovery was implemented.",
            "projection_model_selection": "The model consumed by the retained projection was resolved from authenticated provenance.",
        },
    }
    summary.update({
        "status": "passed_bounded_2d_scientific_gate_with_contract_limitations",
        "scientific_gate_passed": True,
        "full_step4_qualified": False,
        "contract_limitations": limitations,
        "config_dependencies": dependency_records,
        "decision": (
            "The two-dimensional Bdev pseudobulk/pseudo-cell scientific gate passed, "
            "but real ImmGen projection remains blocked until the higher-dimensional "
            "projection and model-bound union contracts are corrected and requalified."
        ),
    })
    receipt.update({
        "status": "complete_with_contract_limitations",
        "scientific_gate_passed": True,
        "full_step4_qualified": False,
        "contract_limitations": limitations,
        "config_dependencies": dependency_records,
    })
    write_json(summary_path, summary)
    write_json(receipt_path, receipt)


if __name__ == "__main__":
    main()
