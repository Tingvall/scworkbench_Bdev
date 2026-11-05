#!/usr/bin/env python3
"""Finalize the bounded 03b qualification manifest."""

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
    summary = json.loads((root / "03_reports" / "01_qualification_summary.json").read_text())
    if summary["status"] != "complete_development_qualification":
        raise RuntimeError("03b qualification summary is not complete")
    if summary["methods"]["intervals"]["permutation_mode"] != "unblocked":
        raise RuntimeError("Signal intervals did not use unblocked permutations")
    if summary["methods"]["gsea"]["permutation_mode"] != "unblocked":
        raise RuntimeError("GSEA did not use unblocked permutations")
    if summary["gene_sets"]["rna_included"] != 38 or summary["gene_sets"]["sct_included"] != 38:
        raise RuntimeError("Not all 38 signatures passed overlap qualification")

    output = root / "06_manifests" / "03_final_run_manifest.json"
    if output.exists():
        raise RuntimeError(f"Refusing to overwrite manifest: {output}")
    files = [
        {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(root.rglob("*"))
        if path.is_file() and path != output
    ]
    workflow_dir = Path(__file__).resolve().parent
    workflow_files = [
        {"path": str(path), "sha256": sha256(path)}
        for path in sorted(workflow_dir.iterdir()) if path.is_file()
    ]
    source_manifest = Path(config["inputs"]["source_bundle"]["path"]) / "manifest.json"
    if sha256(source_manifest) != config["inputs"]["source_bundle"]["manifest_sha256"]:
        raise RuntimeError("Archived source manifest changed")
    manifest = {
        "schema_version": "bdev.parse38_trajectory_final_manifest.v1",
        "status": "complete_development_qualification",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "summary": summary,
        "runtime": config["runtime"],
        "source_bundle": {
            "path": str(source_manifest.parent),
            "manifest_sha256": sha256(source_manifest),
            "mutated": False,
        },
        "workflow_files": workflow_files,
        "files": files,
        "implementation_unlocked": False,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"manifest_sha256={sha256(output)}")


if __name__ == "__main__":
    main()
