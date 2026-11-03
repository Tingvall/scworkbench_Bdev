#!/usr/bin/env python3
"""Freeze checksums and completion claims for the Bdev trajectory recreation."""

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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    config = yaml.safe_load(config_path.read_text())
    results = Path(config["storage"]["results_root"])
    report_path = results / "03_reports" / "01_trajectory_recreation_summary.json"
    report = json.loads(report_path.read_text())
    if report["status"] != "complete_unpromoted_development_validation":
        raise RuntimeError("Trajectory report is not complete")
    source_manifest = Path(config["inputs"]["source_bundle"]["path"]) / "manifest.json"
    if sha256(source_manifest) != config["inputs"]["source_bundle"]["manifest_sha256"]:
        raise RuntimeError("Archived source manifest checksum changed")

    workflow = Path(__file__).resolve().parent
    scripts = sorted(workflow.glob("0[1-3]_*"))
    output_files = sorted(
        path for path in results.glob("*/*")
        if path.is_file() and path.name != "04_run_manifest.json"
    )
    manifest = {
        "schema_version": "bdev.slingshot_trajectory_run_manifest.v1",
        "status": "complete_unpromoted_development_validation",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "scripts": {path.name: sha256(path) for path in scripts},
        "runtime": config["runtime"],
        "source_bundle": {
            "path": config["inputs"]["source_bundle"]["path"],
            "manifest_sha256": sha256(source_manifest),
            "mutated": False,
        },
        "result_summary": report,
        "outputs": {
            str(path.relative_to(results)): {"sha256": sha256(path), "bytes": path.stat().st_size}
            for path in output_files
        },
        "scworkbench_policy": {
            "operation": "scwbtrajectory::slingshot_fit",
            "outputs_promoted": False,
            "staged_session_retained": True,
            "preserved_pt_overwritten": False,
        },
    }
    destination = results / "06_manifests" / "04_run_manifest.json"
    destination.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": manifest["status"],
        "run_identity": manifest["run_identity"],
        "output_count": len(output_files),
        "manifest": str(destination),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
