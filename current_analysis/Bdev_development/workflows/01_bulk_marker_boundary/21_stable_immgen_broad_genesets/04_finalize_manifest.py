#!/usr/bin/env python3
"""Write the immutable diagnostic manifest and concise combined report."""

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
    sample_labels = pd.read_csv(root / "01_tables" / "09_bdev_sampleid_label_summary.tsv", sep="\t")
    report_lines = [
        "# Stable ImmGen broad-lineage exploration",
        "",
        "## ImmGen-only internal result",
        "",
        f"- Nested leave-one-population-out winners: {internal['nested_leave_population_out_winner_correct']}/{internal['libraries']} correct ({100 * internal['nested_leave_population_out_winner_accuracy']:.1f}%).",
        f"- Fixed score/margin gate: {internal['fixed_gate_accepted']}/{internal['libraries']} accepted; {internal['fixed_gate_correct']} correct ({100 * internal['fixed_gate_precision']:.1f}% precision).",
        "- Bdev was not used for marker discovery, ranking, or internal threshold evaluation.",
        "",
        "## Bdev diagnostic summary by sampleId",
        "",
    ]
    for sample_id in sorted(sample_labels["sampleId"].astype(str).unique()):
        rows = sample_labels[sample_labels["sampleId"].astype(str) == sample_id]
        text = ", ".join(f"{row.reported_label}={int(row.cells)}" for row in rows.itertuples(index=False))
        report_lines.append(f"- sampleId {sample_id}: {text}")
    report_lines.extend([
        "",
        "All labels are provisional. No feature set, label, or score artifact was promoted.",
    ])
    report_path = root / "03_reports" / "03_combined_report.md"
    report_path.write_text("\n".join(report_lines) + "\n")

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
    workflow_files = [
        {"path": str(path), "sha256": sha256(path)}
        for path in sorted(workflow_dir.iterdir()) if path.is_file()
    ]
    manifest = {
        "schema_version": "bdev.stable_immgen_broad_manifest.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "result": {"internal": internal, "bdev": bdev},
        "execution": {
            "bundle": str(bundle_manifest.parent),
            "bundle_manifest_sha256": sha256(bundle_manifest),
            "bundle_fingerprint": state["current_state"]["fingerprint"],
            "registry_unchanged": True,
            "artifacts_promoted": False,
        },
        "workflow_files": workflow_files,
        "files": files,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"manifest_sha256={sha256(output)}")


if __name__ == "__main__":
    main()
