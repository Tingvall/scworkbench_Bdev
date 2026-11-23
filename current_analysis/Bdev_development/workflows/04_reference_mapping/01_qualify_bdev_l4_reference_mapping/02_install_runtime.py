#!/usr/bin/env python3
"""Install the exact Step 4 development core and operation-pack commits."""

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
    config = yaml.safe_load(args.config.resolve().read_text())
    runtime = config["runtime"]
    run_library = Path(runtime["run_library"])
    run_library.mkdir(parents=True, exist_ok=True)

    packages = []
    for key in ("scworkbench", "scwbRepresentation", "scwbReferenceMap"):
        record = runtime[key]
        source = Path(record["source_path"])
        if git(source, "rev-parse", "HEAD") != record["source_commit"]:
            raise RuntimeError(f"{key} source commit differs from configuration")
        if git(source, "status", "--porcelain"):
            raise RuntimeError(f"{key} source tree is not clean")
        if key != "scworkbench":
            manifest = source / "inst" / "scworkbench" / "operation-pack.json"
            if sha256(manifest) != record["manifest_sha256"]:
                raise RuntimeError(f"{key} operation-pack manifest checksum changed")
        packages.append((key, source))

    r_libs = ":".join([str(run_library), runtime["shared_development_library"]])
    for key, source in packages:
        command = [
            runtime["apptainer"], "exec", "--cleanenv",
            "--bind", "/mnt/sata2:/mnt/sata2",
            "--bind", "/home/jtg/msi_projects/scworkbench_code:/scworkbench_code",
            "--env", f"R_LIBS_USER={r_libs}", runtime["sif"],
            "R", "CMD", "INSTALL", f"--library={run_library}",
            f"/scworkbench_code/{source.name}",
        ]
        subprocess.run(command, check=True)

    receipt = {
        "schema_version": "bdev.reference_mapping_runtime_install.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "run_library": str(run_library),
        "packages": [
            {"package": key, "source": str(source), "commit": runtime[key]["source_commit"]}
            for key, source in packages
        ],
        "baseline_sif_unchanged": True,
        "new_sif_built": False,
    }
    output = Path(config["storage"]["results_root"]) / "06_manifests" / "02_runtime_install.json"
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
