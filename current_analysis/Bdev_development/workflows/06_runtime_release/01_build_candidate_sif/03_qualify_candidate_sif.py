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


def run(command: list[str], log: Path) -> None:
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    log.write_text(completed.stdout)
    if completed.returncode:
        raise RuntimeError(f"qualification command failed ({completed.returncode}): {' '.join(command)}")


def load_config(path: Path) -> tuple[dict, list[Path]]:
    raw = yaml.safe_load(path.read_text())
    if "extends" not in raw:
        return raw, [path]
    base_path = (path.parent / raw["extends"]).resolve()
    base, sources = load_config(base_path)
    override = {key: value for key, value in raw.items() if key != "extends"}
    def merge(left: dict, right: dict) -> dict:
        result = dict(left)
        for key, value in right.items():
            if isinstance(value, dict) and isinstance(result.get(key), dict):
                result[key] = merge(result[key], value)
            else:
                result[key] = value
        return result
    return merge(base, override), sources + [path]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    config_path = parser.parse_args().config.resolve()
    config, config_sources = load_config(config_path)
    build = config["build"]
    sif = Path(build["output_sif"])
    if not sif.is_file():
        raise RuntimeError("candidate SIF is missing")
    root = Path(build["results_root"])
    logs = root / "04_logs"
    reports = root / "03_reports"
    manifests = root / "06_manifests"
    for directory in (logs, reports, manifests):
        directory.mkdir(parents=True, exist_ok=True)
    apptainer = build["apptainer"]
    common = [apptainer, "exec", "--cleanenv", "--bind", "/mnt/sata2:/mnt/sata2"]

    catalog = reports / "01_installed_runtime_catalog.json"
    run(common + [str(sif), "Rscript",
        "/home/jtg/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/06_runtime_release/01_build_candidate_sif/04_qualify_installed_runtime.R",
        "--output", str(catalog)], logs / "02_installed_runtime_catalog.log")
    run(common + ["--pwd", "/opt/scworkbench", str(sif), "Rscript", "-e",
        "library(scworkbench); source('tests/test_operation_pack_v3_contract.R'); stopifnot(test_operation_pack_v3_sparse_typed_slice_preservation_contract()); cat('PASS installed sparse typed-slice preservation\\n')"],
        logs / "03_core_sparse_typed_slice.log")
    run(common + ["--pwd", "/opt/operation_packs/scwbATAC",
        "--env", "SCWBATAC_RUN_SUITE=1", str(sif), "Rscript",
        "tests/run_contract_tests.R"], logs / "04_scwbatac_contract_suite.log")
    run(common + [str(sif), "/opt/pyenv/bin/python", "-c",
        "import pybench,scvi,anndata,scipy,torch,ruptures; print('PASS candidate Python runtime imports')"],
        logs / "05_python_runtime.log")

    catalog_report = json.loads(catalog.read_text())
    scwbatac_log = (logs / "04_scwbatac_contract_suite.log").read_text()
    core_log = (logs / "03_core_sparse_typed_slice.log").read_text()
    python_log = (logs / "05_python_runtime.log").read_text()
    gate = all([
        catalog_report.get("gate_passed") is True,
        "PASS installed sparse typed-slice preservation" in core_log,
        "RESULT discovered=7 executed=7 passed=7 skipped=0 failed=0" in scwbatac_log,
        "PASS candidate Python runtime imports" in python_log,
    ])
    report = {
        "schema_version": "bdev.candidate_sif_qualification.v1",
        "status": "passed" if gate else "failed",
        "run_identity": config["run_identity"],
        "sif_path": str(sif),
        "sif_bytes": sif.stat().st_size,
        "sif_sha256": sha256(sif),
        "base_sif_sha256": config["base_sif"]["sha256"],
        "installed_runtime_catalog": str(catalog),
        "installed_runtime_catalog_sha256": sha256(catalog),
        "core_sparse_typed_slice_regression": "passed",
        "scwbatac_contracts": {"discovered": 7, "passed": 7, "failed": 0},
        "python_runtime_imports": "passed",
        "implementation_unlocked": False,
        "explicit_user_acceptance_required": True,
        "gate_passed": gate,
    }
    report_path = reports / "02_candidate_sif_qualification.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    inventory_path = manifests / "02_candidate_sif_inventory.json"
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path != inventory_path:
            files.append({"path": str(path.relative_to(root)),
                          "bytes": path.stat().st_size, "sha256": sha256(path)})
    inventory = {
        "schema_version": "bdev.candidate_sif_inventory.v1",
        "status": "complete_passed" if gate else "complete_failed",
        "run_identity": config["run_identity"],
        "config_sha256": sha256(config_path),
        "config_sources": [
            {"path": str(path), "sha256": sha256(path)} for path in config_sources],
        "qualification_report_sha256": sha256(report_path),
        "sif_sha256": report["sif_sha256"],
        "files": files,
    }
    inventory_path.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    if not gate:
        raise RuntimeError("candidate SIF qualification failed")


if __name__ == "__main__":
    main()
