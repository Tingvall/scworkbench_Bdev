#!/usr/bin/env python3
"""Validate and inventory the label-blind ImmGen sparsity-adapter diagnostic."""

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
    output = root / "06_manifests" / "02_run_file_inventory.json"

    required = {
        "01_tables/01_global_reference_structure.tsv",
        "01_tables/02_sparsity_pool_distribution.tsv",
        "01_tables/03_pseudocell_design.tsv",
        "01_tables/04_bulk_rank_preservation_audit.tsv",
        "01_tables/05_pseudocell_peakvi_latent.tsv",
        "01_tables/06_pseudocell_neighbor_ledger.tsv",
        "01_tables/07_pseudocell_fixed_umap_positions.tsv",
        "01_tables/08_pseudocell_number_convergence.tsv",
        "01_tables/09_library_consensus_neighbors.tsv",
        "01_tables/10_consensus_neighbor_reuse.tsv",
        "01_tables/11_library_adapter_summary.tsv",
        "01_tables/12_full_vs_adapter_latent.tsv",
        "01_tables/13_library_neighbor_technical_composition.tsv",
        "01_tables/14_recurrent_reference_hubs.tsv",
        "02_figures/01_sparsity_matching.png",
        "02_figures/01_sparsity_matching.pdf",
        "02_figures/02_fixed_umap_adapter_positions.png",
        "02_figures/02_fixed_umap_adapter_positions.pdf",
        "02_figures/03_library_l4_distributions.png",
        "02_figures/03_library_l4_distributions.pdf",
        "02_figures/04_pseudocell_number_convergence.png",
        "02_figures/04_pseudocell_number_convergence.pdf",
        "02_figures/05_reference_hub_concentration.png",
        "02_figures/05_reference_hub_concentration.pdf",
        "03_reports/01_sparsity_adapter_summary.json",
        "03_reports/02_hub_concentration_summary.json",
        "04_logs/01_sparsity_adapter.log",
        "04_logs/02_hub_concentration_audit.log",
        "06_manifests/01_execution_receipt.json",
    }
    files = [
        {
            "path": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(root.rglob("*"))
        if path.is_file() and path != output
    ]
    present = {item["path"] for item in files}
    missing = sorted(required - present)
    if missing:
        raise RuntimeError(f"Required outputs missing: {missing}")

    summary = json.loads(
        (root / "03_reports" / "01_sparsity_adapter_summary.json").read_text()
    )
    expected_pseudocells = int(config["acceptance"]["expected_pseudocells"])
    if summary["status"] != "complete_diagnostic":
        raise RuntimeError("Diagnostic did not complete")
    if int(summary["pseudocells"]) != expected_pseudocells:
        raise RuntimeError("Unexpected pseudo-cell count")
    construction = summary["adapter_construction"]
    if construction["labels_used"] or construction["expected_positions_used"]:
        raise RuntimeError("Adapter construction was not label blind")
    if construction["metadata_loaded_for_construction"]:
        raise RuntimeError("Reference metadata entered adapter construction")
    if not summary["source_manifests_unchanged"]:
        raise RuntimeError("A source manifest changed")
    if summary["mapping_accepted"]:
        raise RuntimeError("Diagnostic mapping was incorrectly accepted")
    hub_audit = json.loads(
        (root / "03_reports" / "02_hub_concentration_summary.json").read_text()
    )
    if not hub_audit["audit_ran_after_projection"]:
        raise RuntimeError("Hub audit was not post-projection")
    if hub_audit["labels_used_for_adapter"]:
        raise RuntimeError("Post-projection audit changed adapter construction")

    inventory = {
        "schema_version": "bdev.immgen_label_blind_sparsity_adapter_inventory.v1",
        "status": "complete_diagnostic",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "file_count": len(files),
        "files": files,
        "labels_used_for_adapter": False,
        "expected_positions_used_for_adapter": False,
        "source_mutated": False,
        "mapping_accepted": False,
        "paper_analysis": False,
        "implementation_unlocked": False,
        "sif_built": False,
    }
    output.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {"status": inventory["status"], "files": len(files), "output": str(output)},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
