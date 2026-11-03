#!/usr/bin/env python3
"""Finalize the hierarchical population-calibration diagnostic manifest."""

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
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    root = Path(config["storage"]["results_root"])
    output = root / "06_manifests" / "03_final_run_manifest.json"
    if output.exists():
        raise RuntimeError(f"Refusing to overwrite manifest: {output}")
    for name, spec in config["inputs"].items():
        path = Path(spec["path"])
        if sha256(path) != spec["sha256"]:
            raise RuntimeError(f"Input checksum changed: {name}")
    summary = json.loads((root / "03_reports" / "02_combined_summary.json").read_text())
    bundle_manifest = Path(summary["execution_bundle"]) / "manifest.json"
    bundle_state = json.loads(bundle_manifest.read_text())
    if bundle_state["current_state"]["fingerprint"] != config["bdev_test"]["bundle_fingerprint"]:
        raise RuntimeError("Session bundle registry fingerprint changed")
    if not summary["session_bundle_registry_unchanged"] or summary["accepted_labels_emitted"]:
        raise RuntimeError("Acceptance boundary violated")
    files = [
        {
            "path": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(root.rglob("*"))
        if path.is_file() and path != output
    ]
    workflow_dir = Path(__file__).resolve().parent
    workflows = [
        {"path": str(path), "sha256": sha256(path)}
        for path in sorted(workflow_dir.iterdir())
        if path.is_file()
    ]
    manifest = {
        "schema_version": "bdev.immgen_hierarchical_population_calibration_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "inputs": {
            name: {"path": spec["path"], "sha256": sha256(Path(spec["path"]))}
            for name, spec in config["inputs"].items()
        },
        "execution": {
            "session_bundle": summary["execution_bundle"],
            "bundle_fingerprint": summary["session_bundle_fingerprint"],
            "bundle_manifest_sha256_after": sha256(bundle_manifest),
            "registry_unchanged": True,
            "artifacts_promoted": False,
        },
        "result": {
            "reference": summary["reference"],
            "bdev_policy_summary": summary["bdev_policy_summary"],
            "accepted_labels_emitted": False,
        },
        "workflow_files": workflows,
        "files": files,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"manifest_sha256={sha256(output)}")


if __name__ == "__main__":
    main()
