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
    config_path = parser.parse_args().config.resolve()
    config = yaml.safe_load(config_path.read_text())
    runtime = config["runtime"]
    storage = config["storage"]
    inputs = config["inputs"]

    source = Path(inputs["qualified_bundle"]["path"])
    source_manifest = source / "manifest.json"
    if sha256(source_manifest) != inputs["qualified_bundle"]["manifest_sha256"]:
        raise RuntimeError("qualified bundle manifest changed")
    for key in ("qualification_report", "qualification_inventory",
                "prepared_cells", "quantification_audit"):
        path = Path(inputs[key]["path"])
        if sha256(path) != inputs[key]["sha256"]:
            raise RuntimeError(f"input changed: {key}")
    sif = Path(runtime["sif"])
    if sha256(sif) != runtime["sif_sha256"]:
        raise RuntimeError("baseline SIF changed")

    packages = []
    for key in ("scworkbench", "scwbLatent", "scwbRepresentation",
                "scwbReferenceMap"):
        record = runtime[key]
        path = Path(record["source_path"])
        if git(path, "rev-parse", "HEAD") != record["source_commit"]:
            raise RuntimeError(f"{key} commit changed")
        if git(path, "status", "--porcelain"):
            raise RuntimeError(f"{key} source is dirty")
        manifest = path / "inst" / "scworkbench" / "operation-pack.json"
        if "manifest_sha256" in record and sha256(manifest) != record["manifest_sha256"]:
            raise RuntimeError(f"{key} operation-pack manifest changed")
        packages.append((key, path))

    results = Path(storage["results_root"])
    for name in ("01_tables", "02_figures", "03_reports", "04_logs",
                 "05_objects", "06_manifests"):
        (results / name).mkdir(parents=True, exist_ok=True)
    run_library = Path(storage["run_library"])
    sessions = Path(storage["sessions"])
    run_library.mkdir(parents=True, exist_ok=True)
    sessions.mkdir(parents=True, exist_ok=True)

    clone = Path(storage["mapping_bundle"])
    if clone.exists():
        raise RuntimeError(f"refusing to overwrite mapping clone: {clone}")
    clone.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([
        "cp", "--archive", "--reflink=auto", "--", str(source), str(clone),
    ], check=True)
    if sha256(clone / "manifest.json") != inputs["qualified_bundle"]["manifest_sha256"]:
        raise RuntimeError("mapping clone is not initially identical")

    r_libraries = f"{run_library}:{runtime['shared_development_library']}"
    for _, path in packages:
        subprocess.run([
            runtime["apptainer"], "exec", "--cleanenv",
            "--bind", "/mnt/sata2:/mnt/sata2",
            "--bind", "/home/jtg/msi_projects/scworkbench_code:/scworkbench_code",
            "--env", f"R_LIBS_USER={r_libraries}", runtime["sif"],
            "R", "CMD", "INSTALL", f"--library={run_library}",
            f"/scworkbench_code/{path.name}",
        ], check=True)

    receipt = {
        "schema_version": "bdev.full_gated_20d_runtime.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "qualified_bundle_manifest_sha256": sha256(source_manifest),
        "mapping_bundle": str(clone),
        "clone_initial_manifest_sha256": sha256(clone / "manifest.json"),
        "copy_mode": "cp_archive_reflink_auto",
        "sif_sha256": sha256(sif),
        "packages": [
            {"package": key, "version": runtime[key]["version"],
             "commit": git(path, "rev-parse", "HEAD")}
            for key, path in packages
        ],
        "source_mutation": False,
        "paper_analysis": False,
        "sif_built": False,
    }
    (results / "06_manifests" / "01_runtime_and_clone_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
