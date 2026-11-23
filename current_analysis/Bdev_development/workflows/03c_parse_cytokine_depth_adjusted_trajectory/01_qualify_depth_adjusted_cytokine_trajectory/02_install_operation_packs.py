#!/usr/bin/env python3
"""Install the pinned core, scwbRNA, and scwbTrajectory sources for 03c."""

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
    run_library = Path(runtime["run_library"])
    run_library.mkdir(parents=True, exist_ok=True)

    packages = []
    for key in ("scworkbench", "scwbRNA", "scwbTrajectory"):
        record = runtime[key]
        source = Path(record["source_path"])
        if git(source, "rev-parse", "HEAD") != record["source_commit"]:
            raise RuntimeError(f"{key} source commit differs from configuration")
        if git(source, "status", "--porcelain"):
            raise RuntimeError(f"{key} source tree is not clean")
        if key != "scworkbench":
            manifest = source / "inst" / "scworkbench" / "operation-pack.json"
            if sha256(manifest) != record["source_manifest_file_sha256"]:
                raise RuntimeError(f"{key} operation-pack manifest checksum changed")
        packages.append((key, source))

    r_libs = ":".join(
        [str(run_library), runtime["dependency_library"], runtime["shared_development_library"]]
    )
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
        "schema_version": "bdev.operation_pack_install.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "run_library": str(run_library),
        "packages": [
            {"package": key, "source": str(source), "commit": runtime[key]["source_commit"]}
            for key, source in packages
        ],
        "r_libs_user": r_libs,
    }
    output = Path(config["storage"]["results_root"]) / "06_manifests" / "02_operation_pack_install.json"
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
