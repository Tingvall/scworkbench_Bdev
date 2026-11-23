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
    source = Path(config["inputs"]["source_bundle"]["path"])
    source_manifest = source / "manifest.json"
    expected_manifest = config["inputs"]["source_bundle"]["manifest_sha256"]
    if sha256(source_manifest) != expected_manifest:
        raise RuntimeError("qualification-parent bundle manifest changed")
    original = config["inputs"]["source_bundle"].get("original_source")
    original_manifest = None
    if original:
        original_manifest = Path(original["path"]) / "manifest.json"
        if sha256(original_manifest) != original["manifest_sha256"]:
            raise RuntimeError("original source bundle manifest changed")
    if sha256(Path(runtime["sif"])) != runtime["sif_sha256"]:
        raise RuntimeError("baseline SIF changed")

    results = Path(config["storage"]["results_root"])
    for name in ("01_tables", "02_figures", "03_reports", "04_logs",
                 "05_objects", "06_manifests"):
        (results / name).mkdir(parents=True, exist_ok=True)
    Path(config["storage"]["sessions"]).mkdir(parents=True, exist_ok=True)
    run_library = Path(runtime["run_library"])
    run_library.mkdir(parents=True, exist_ok=True)

    packages = []
    for key in ("scworkbench", "scwbATAC"):
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

    r_libraries = f"{run_library}:{runtime['shared_development_library']}"
    for _, path in packages:
        subprocess.run([
            runtime["apptainer"], "exec", "--cleanenv",
            "--bind", "/mnt/sata2:/mnt/sata2",
            "--bind", "/home/jtg/msi_projects/scworkbench_code:/scworkbench_code",
            "--env", f"R_LIBS_USER={r_libraries}", runtime["sif"],
            "R", "CMD", "INSTALL", f"--library={run_library}",
            f"/scworkbench_code/{path.name}"
        ], check=True)

    receipt = {
        "schema_version": "bdev.scwbatac_hvp_sparse_runtime.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "source_manifest_sha256": sha256(source_manifest),
        "original_source_manifest_sha256": (
            sha256(original_manifest) if original_manifest else None),
        "sif_sha256": sha256(Path(runtime["sif"])),
        "packages": [
            {"package": key, "commit": git(path, "rev-parse", "HEAD")}
            for key, path in packages
        ],
        "source_access": "writable_qualification_parent_no_promotion",
        "original_source_access": "checksum_only_no_open"
    }
    (results / "06_manifests" / "01_runtime_and_source_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
