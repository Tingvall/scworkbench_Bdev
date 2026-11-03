#!/usr/bin/env python3
"""Finalize the lineage-anchored narrow-label diagnostic run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
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
    internal = json.loads((root / "03_reports" / "01_immgen_internal_summary.json").read_text())
    bdev = json.loads((root / "03_reports" / "02_bdev_summary.json").read_text())
    lines = [
        "# Lineage-anchored narrow-label exploration",
        "",
        "## ImmGen-only internal result",
        "",
        f"- All {internal['populations_competing']} populations competed directly; no broad gate was used.",
        f"- Unthresholded exact population: {internal['unthresholded_exact_correct']}/{internal['heldout_libraries']} ({100 * internal['unthresholded_exact_accuracy']:.1f}%).",
        f"- Unthresholded correct lineage: {internal['unthresholded_lineage_correct']}/{internal['heldout_libraries']} ({100 * internal['unthresholded_lineage_accuracy']:.1f}%).",
        f"- Accepted: {internal['accepted']}/{internal['heldout_libraries']}; exact {internal['accepted_exact_correct']}/{internal['accepted']}, correct lineage {internal['accepted_lineage_correct']}/{internal['accepted']}.",
        "",
        "## Bdev result by sampleId",
        "",
    ]
    for row in bdev["sample_ids"]:
        lines.append(f"- sampleId {row['sampleId']}: cells={row['cells']}, accepted={row['accepted']}, borderline={row['borderline']}, unknown={row['unknown']}, accepted B={row['accepted_B']}")
    lines.extend(["", "All outputs remain diagnostic. Bdev did not affect gene selection, and nothing was promoted."])
    report_path = root / "03_reports" / "03_combined_report.md"
    report_path.write_text("\n".join(lines) + "\n")

    output = root / "06_manifests" / "04_final_run_manifest.json"
    if output.exists():
        raise RuntimeError(f"Refusing to overwrite manifest: {output}")
    bundle_manifest = Path(config["storage"]["intermediate_root"]) / "05_objects" / "01_session_bundle" / "manifest.json"
    state = json.loads(bundle_manifest.read_text())
    if state["current_state"]["fingerprint"] != config["inputs"]["source_bundle"]["fingerprint"]:
        raise RuntimeError("Bundle registry fingerprint changed")
    files = [
        {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(root.rglob("*")) if path.is_file() and path != output
    ]
    workflow_dir = Path(__file__).resolve().parent
    workflows = [
        {"path": str(path), "sha256": sha256(path)}
        for path in sorted(workflow_dir.iterdir()) if path.is_file()
    ]
    manifest = {
        "schema_version": "bdev.lineage_anchored_narrow_manifest.v1",
        "status": "complete_diagnostic_only", "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "result": {"internal": internal, "bdev": bdev},
        "execution": {"bundle": str(bundle_manifest.parent), "bundle_manifest_sha256": sha256(bundle_manifest),
          "bundle_fingerprint": state["current_state"]["fingerprint"], "registry_unchanged": True, "artifacts_promoted": False},
        "workflow_files": workflows, "files": files,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"manifest_sha256={sha256(output)}")


if __name__ == "__main__":
    main()
