#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git(path: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(path), *args], text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--resume-existing", action="store_true")
    arguments = parser.parse_args()
    config_path = arguments.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    runtime = config["runtime"]
    inputs = config["inputs"]
    storage = config["storage"]

    qualified = Path(inputs["qualified_bundle"]["path"])
    bulk = Path(inputs["bulk_query_bundle"]["path"])
    checks = [
        (qualified / "manifest.json",
         inputs["qualified_bundle"]["manifest_sha256"]),
        (bulk / "manifest.json",
         inputs["bulk_query_bundle"]["manifest_sha256"]),
        (Path(inputs["fixed_umap_bundle"]["path"]) / "manifest.json",
         inputs["fixed_umap_bundle"]["manifest_sha256"]),
        (Path(inputs["library_manifest"]["path"]),
         inputs["library_manifest"]["sha256"]),
        (Path(inputs["qualification_report"]["path"]),
         inputs["qualification_report"]["sha256"]),
        (Path(inputs["authenticated_l4"]["model_path"]),
         inputs["authenticated_l4"]["model_sha256"]),
        (Path(inputs["authenticated_l4"]["pseudotime_path"]),
         inputs["authenticated_l4"]["pseudotime_sha256"]),
        (Path(runtime["sif"]), runtime["sif_sha256"]),
    ]
    for path, expected in checks:
        if not path.exists() or sha256(path) != expected:
            raise RuntimeError(f"input checksum mismatch: {path}")

    packages = []
    for key in ("scworkbench", "scwbLatent", "scwbRepresentation",
                "scwbReferenceMap"):
        record = runtime[key]
        source = Path(record["source_path"])
        if git(source, "rev-parse", "HEAD") != record["source_commit"]:
            raise RuntimeError(f"{key} source commit changed")
        if git(source, "status", "--porcelain"):
            raise RuntimeError(f"{key} source is dirty")
        manifest = source / "inst" / "scworkbench" / "operation-pack.json"
        if "manifest_sha256" in record and sha256(manifest) != record["manifest_sha256"]:
            raise RuntimeError(f"{key} operation-pack manifest changed")
        packages.append({"package": key, "version": record["version"],
                         "commit": record["source_commit"]})

    results = Path(storage["results_root"])
    for name in ("01_tables", "02_figures", "03_reports", "04_logs",
                 "05_objects", "06_manifests"):
        (results / name).mkdir(parents=True, exist_ok=True)
    Path(storage["sessions"]).mkdir(parents=True, exist_ok=True)

    clone = Path(storage["mapping_bundle"])
    if clone.exists():
        if not arguments.resume_existing:
            raise RuntimeError(f"refusing to overwrite mapping clone: {clone}")
        previous_receipt = results / "06_manifests" / "01_prepare_receipt.json"
        if not previous_receipt.exists() or not (clone / "manifest.json").exists():
            raise RuntimeError("resume requested without an authenticated prior clone")
        previous = json.loads(previous_receipt.read_text())
        if previous.get("run_identity") != config["run_identity"]:
            raise RuntimeError("existing clone belongs to a different run")
        initial_manifest_sha256 = previous["clone_initial_manifest_sha256"]
    else:
        if arguments.resume_existing:
            raise RuntimeError("resume requested but mapping clone is absent")
        clone.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([
            "cp", "--recursive", "--reflink=auto", "--", str(qualified),
            str(clone),
        ], check=True)
        if sha256(clone / "manifest.json") != inputs["qualified_bundle"]["manifest_sha256"]:
            raise RuntimeError("mapping clone is not initially identical")
        initial_manifest_sha256 = sha256(clone / "manifest.json")

    receipt = {
        "schema_version": "bdev.immgen_bulk_gated_20d_prepare.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "qualified_bundle_manifest_sha256": sha256(qualified / "manifest.json"),
        "bulk_bundle_manifest_sha256": sha256(bulk / "manifest.json"),
        "clone_initial_manifest_sha256": initial_manifest_sha256,
        "clone_resume_manifest_sha256": sha256(clone / "manifest.json"),
        "copy_mode": "cp_recursive_reflink_auto",
        "mapping_bundle": str(clone),
        "runtime_library": runtime["r_library"],
        "sif_sha256": sha256(Path(runtime["sif"])),
        "packages": packages,
        "source_mutation": False,
        "model_refit": False,
        "umap_refit": False,
        "paper_analysis": False,
        "sif_built": False,
        "resumed_existing_clone": arguments.resume_existing,
    }
    (results / "06_manifests" / "01_prepare_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
