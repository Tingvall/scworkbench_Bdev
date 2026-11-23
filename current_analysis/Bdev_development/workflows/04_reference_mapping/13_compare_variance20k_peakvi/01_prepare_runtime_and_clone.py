#!/usr/bin/env python3
"""Install pinned packages and reflink-clone the source for variance selection."""

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


def git(source: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(source), *args], text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    runtime = config["runtime"]

    source_bundle = Path(config["inputs"]["source_bundle"]["path"])
    if sha256(source_bundle / "manifest.json") != config["inputs"]["source_bundle"]["manifest_sha256"]:
        raise RuntimeError("source bundle manifest changed")
    if sha256(Path(runtime["sif"])) != runtime["sif_sha256"]:
        raise RuntimeError("SIF checksum changed")

    run_library = Path(runtime["run_library"])
    run_library.mkdir(parents=True, exist_ok=True)
    results = Path(config["storage"]["results_root"])
    for name in ("01_tables", "02_figures", "03_reports", "04_logs", "05_objects", "06_manifests"):
        (results / name).mkdir(parents=True, exist_ok=True)
    Path(config["storage"]["sessions"]).mkdir(parents=True, exist_ok=True)

    packages = []
    for key in ("scworkbench", "scwbRNA", "scwbLatent"):
        record = runtime[key]
        source = Path(record["source_path"])
        if git(source, "rev-parse", "HEAD") != record["source_commit"]:
            raise RuntimeError(f"{key} source commit differs from configuration")
        if git(source, "status", "--porcelain"):
            raise RuntimeError(f"{key} source tree is not clean")
        manifest = source / "inst" / "scworkbench" / "operation-pack.json"
        if "manifest_sha256" in record and sha256(manifest) != record["manifest_sha256"]:
            raise RuntimeError(f"{key} operation-pack manifest changed")
        packages.append((key, source))

    r_libs = ":".join((str(run_library), runtime["shared_development_library"]))
    for _, source in packages:
        subprocess.run([
            runtime["apptainer"], "exec", "--cleanenv", "--bind", "/mnt/sata2:/mnt/sata2",
            "--bind", "/home/jtg/msi_projects/scworkbench_code:/scworkbench_code",
            "--env", f"R_LIBS_USER={r_libs}", runtime["sif"], "R", "CMD", "INSTALL",
            f"--library={run_library}", f"/scworkbench_code/{source.name}",
        ], check=True)

    clone = Path(config["storage"]["selection_clone"])
    if clone.exists():
        if sha256(clone / "manifest.json") != config["inputs"]["source_bundle"]["manifest_sha256"]:
            raise RuntimeError("existing selection clone is not the pristine source clone")
    else:
        clone.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["cp", "--recursive", "--reflink=auto", "--", str(source_bundle), str(clone)], check=True)
    if sha256(clone / "manifest.json") != config["inputs"]["source_bundle"]["manifest_sha256"]:
        raise RuntimeError("selection clone differs before operation execution")

    receipt = {
        "schema_version": "bdev.variance20k_runtime_and_clone.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "sif_sha256": sha256(Path(runtime["sif"])),
        "run_library": str(run_library),
        "source_bundle_manifest_sha256": sha256(source_bundle / "manifest.json"),
        "selection_clone": str(clone),
        "packages": [
            {"package": key, "source": str(source), "commit": git(source, "rev-parse", "HEAD")}
            for key, source in packages
        ],
    }
    (results / "06_manifests" / "01_runtime_and_clone_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
