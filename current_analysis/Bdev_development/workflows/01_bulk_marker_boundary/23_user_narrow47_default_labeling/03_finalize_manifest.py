#!/usr/bin/env python3
"""Create a checksummed manifest for the exploratory narrow47 run."""

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
    summary = json.loads((root / "03_reports" / "01_labeling_summary.json").read_text())
    files = [
        {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(root.rglob("*")) if path.is_file() and path != output
    ]
    workflow_dir = Path(__file__).resolve().parent
    manifest = {
        "schema_version": "bdev.user_narrow47_default_labeling_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "summary": summary,
        "workflow_files": [
            {"path": str(path), "sha256": sha256(path)}
            for path in sorted(workflow_dir.iterdir()) if path.is_file()
        ],
        "files": files,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"manifest_sha256={sha256(output)}")


if __name__ == "__main__":
    main()

