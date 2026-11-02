#!/usr/bin/env python3
"""Stream primary statistics and retain pre-stability eligible candidates."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
from pathlib import Path

import yaml


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def as_bool(value: str) -> bool:
    return value.lower() == "true"


def value_or(value: str, default: float) -> float:
    return float(value) if value else default


def finalize(family: str, target: str, genes: dict[str, dict], contract: dict) -> list[dict]:
    threshold = contract["contrast_families"][family]
    candidates: list[dict] = []
    for gene, state in genes.items():
        eligible = state["filter_pass"] and state["annotation_eligible"]
        if family == "superlineage":
            eligible &= state["support"] >= threshold["minimum_target_population_support"]
            eligible &= state["global_lfc"] >= threshold["minimum_global_log2_fold_change"]
            eligible &= state["global_fdr"] <= threshold["maximum_global_fdr"]
        elif family in {"broad_lineage", "b_branch"}:
            eligible &= state["support"] >= threshold["minimum_target_population_support"]
            eligible &= state["global_lfc"] >= threshold["minimum_global_log2_fold_change"]
            eligible &= state["global_fdr"] <= threshold["maximum_global_fdr"]
            eligible &= state["competitor_lfc"] >= threshold["minimum_each_competitor_log2_fold_change"]
            eligible &= state["competitor_fdr"] <= threshold["maximum_each_competitor_fdr"]
        else:
            eligible &= state["expression"] >= threshold["minimum_target_library_fraction"]
            eligible &= state["competitor_lfc"] >= threshold["minimum_each_neighbor_log2_fold_change"]
            eligible &= state["competitor_fdr"] <= threshold["maximum_each_neighbor_fdr"]
        if eligible:
            candidates.append(
                {
                    "family": family,
                    "target": target,
                    "gene_symbol": gene,
                    "target_support_fraction": state["support"],
                    "worst_log2_fold_change": state["worst_lfc"],
                    "worst_auroc": state["auroc"],
                    "worst_fdr": state["worst_fdr"],
                    "ensembl_gene_id": state["ensembl_gene_id"],
                    "chromosome": state["chromosome"],
                    "gene_biotype": state["gene_biotype"],
                }
            )
    candidates.sort(key=lambda row: (-row["target_support_fraction"], -row["worst_log2_fold_change"], -row["worst_auroc"], row["worst_fdr"], row["gene_symbol"]))
    for rank, row in enumerate(candidates, 1):
        row["preliminary_rank"] = rank
    return candidates


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    require(config["schema_version"] == "bdev.immgen_vm25_marker_execution.v1", "Wrong config schema")
    config_dir = args.config.parent
    contract_path = (config_dir / config["inputs"]["statistical_contract"]["path"]).resolve()
    contract = yaml.safe_load(contract_path.read_text())
    tables = config["outputs"]["tables"]
    statistics_path = args.output_root / tables["complete_statistics"]
    output_path = args.output_root / tables["preliminary_candidates"]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["family", "target", "gene_symbol", "preliminary_rank", "target_support_fraction", "worst_log2_fold_change", "worst_auroc", "worst_fdr", "ensembl_gene_id", "chromosome", "gene_biotype"]

    target_counts: dict[str, int] = {}
    seen_targets: set[tuple[str, str]] = set()
    current_key: tuple[str, str] | None = None
    genes: dict[str, dict] = {}
    total_candidates = 0

    with output_path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                with gzip.open(statistics_path, "rt", newline="") as handle:
                    for row in csv.DictReader(handle, delimiter="\t"):
                        key = (row["family"], row["target"])
                        if current_key is not None and key != current_key:
                            require(key not in seen_targets, f"Noncontiguous target block: {key}")
                            candidates = finalize(current_key[0], current_key[1], genes, contract)
                            writer.writerows(candidates)
                            target_counts["|".join(current_key)] = len(candidates)
                            total_candidates += len(candidates)
                            seen_targets.add(current_key)
                            genes = {}
                        current_key = key
                        gene = row["gene_symbol"]
                        state = genes.setdefault(
                            gene,
                            {
                                "filter_pass": True,
                                "annotation_eligible": True,
                                "support": 1.0,
                                "expression": 1.0,
                                "auroc": 1.0,
                                "worst_lfc": float("inf"),
                                "worst_fdr": 0.0,
                                "global_lfc": float("inf"),
                                "global_fdr": 0.0,
                                "competitor_lfc": float("inf"),
                                "competitor_fdr": 0.0,
                                "ensembl_gene_id": row["ensembl_gene_id"],
                                "chromosome": row["chromosome"],
                                "gene_biotype": row["gene_biotype"],
                            },
                        )
                        state["filter_pass"] &= as_bool(row["filter_pass"])
                        state["annotation_eligible"] &= as_bool(row["core_identity_eligible"])
                        state["support"] = min(state["support"], value_or(row["target_support_fraction"], 0.0))
                        state["expression"] = min(state["expression"], value_or(row["target_expression_fraction"], 0.0))
                        state["auroc"] = min(state["auroc"], value_or(row["target_versus_competitor_auroc"], 0.0))
                        lfc = value_or(row["logFC"], float("-inf"))
                        fdr = value_or(row["FDR"], 1.0)
                        state["worst_lfc"] = min(state["worst_lfc"], lfc)
                        state["worst_fdr"] = max(state["worst_fdr"], fdr)
                        if row["test_type"] == "global":
                            state["global_lfc"] = min(state["global_lfc"], lfc)
                            state["global_fdr"] = max(state["global_fdr"], fdr)
                        else:
                            state["competitor_lfc"] = min(state["competitor_lfc"], lfc)
                            state["competitor_fdr"] = max(state["competitor_fdr"], fdr)
                if current_key is not None:
                    candidates = finalize(current_key[0], current_key[1], genes, contract)
                    writer.writerows(candidates)
                    target_counts["|".join(current_key)] = len(candidates)
                    total_candidates += len(candidates)

    summary = {"schema_version": "bdev.immgen_preliminary_candidate_summary.v1", "total_candidates": total_candidates, "targets": target_counts}
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
