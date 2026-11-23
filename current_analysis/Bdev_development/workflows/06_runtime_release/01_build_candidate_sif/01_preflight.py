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
    return subprocess.check_output(["git", "-C", str(path), *args], text=True).strip()


def tree_sha256(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        digest.update(str(path.relative_to(root)).encode())
        digest.update(b"\0")
        digest.update(bytes.fromhex(sha256(path)))
    return digest.hexdigest()


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
    base = Path(config["base_sif"]["path"])
    if sha256(base) != config["base_sif"]["sha256"]:
        raise RuntimeError("base SIF checksum changed")
    output = Path(config["build"]["output_sif"])
    if output.exists():
        raise RuntimeError(f"candidate SIF already exists: {output}")

    resolved = []
    for component in config["components"]:
        path = Path(component["source_path"])
        observed = git(path, "rev-parse", "HEAD")
        if observed != component["commit"]:
            raise RuntimeError(f"{component['package']} commit changed")
        if git(path, "status", "--porcelain"):
            raise RuntimeError(f"{component['package']} checkout is dirty")
        record = dict(component)
        if component["role"] == "operation_pack":
            manifest = path / "inst" / "scworkbench" / "operation-pack.json"
            observed_manifest = sha256(manifest)
            if observed_manifest != component["manifest_sha256"]:
                raise RuntimeError(f"{component['package']} manifest changed")
            record["observed_manifest_sha256"] = observed_manifest
        resolved.append(record)

    dependencies = []
    for dependency in config.get("runtime_dependencies", []):
        record = dict(dependency)
        if "source_path" in dependency:
            path = Path(dependency["source_path"])
            observed_tree = tree_sha256(path)
            if observed_tree != dependency["tree_sha256"]:
                raise RuntimeError(f"{dependency['package']} installed tree changed")
            record["observed_tree_sha256"] = observed_tree
        dependencies.append(record)

    results = Path(config["build"]["results_root"])
    for name in ("01_tables", "02_figures", "03_reports", "04_logs",
                 "05_objects", "06_manifests"):
        (results / name).mkdir(parents=True, exist_ok=True)
    Path(config["build"]["intermediate_root"]).mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    receipt = {
        "schema_version": "bdev.candidate_sif_preflight.v1",
        "status": "passed",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "config_sources": [
            {"path": str(path), "sha256": sha256(path)} for path in config_sources],
        "base_sif_sha256": sha256(base),
        "output_sif": str(output),
        "components": resolved,
        "runtime_dependencies": dependencies,
    }
    (results / "06_manifests" / "01_build_preflight.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
