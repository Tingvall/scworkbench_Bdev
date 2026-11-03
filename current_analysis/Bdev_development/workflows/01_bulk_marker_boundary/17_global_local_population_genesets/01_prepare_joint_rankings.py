#!/usr/bin/env python3
"""Combine pinned local crossfit rankings with population-balanced global effects."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve(base: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (base / path).resolve()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str], compressed: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if compressed:
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
                    writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n")
                    writer.writeheader(); writer.writerows(rows)
    else:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    if config["schema_version"] != "bdev.immgen_global_local_population_genesets.v1":
        raise RuntimeError("Unexpected config schema")
    inputs = {name: resolve(config_path.parent, spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        if sha256(path) != config["inputs"][name]["sha256"]:
            raise RuntimeError(f"Checksum mismatch for {name}")

    counts = pd.read_csv(inputs["selected_counts"], sep="\t", index_col=0)
    metadata = pd.read_csv(inputs["metadata"], sep="\t", dtype=str).fillna("")
    if list(counts.columns) != list(metadata["library_id"]):
        raise RuntimeError("Count columns and metadata are misaligned")
    plan_rows = read_tsv(inputs["local_fold_plan"])
    local_rankings = read_tsv(inputs["local_fold_rankings"])
    local_audit = read_tsv(inputs["local_model_audit"])
    population_info = metadata.drop_duplicates("population").set_index("population")
    discovery_role = config["target_scope"]["discovery_role"]
    heldout_role = config["target_scope"]["held_out_role"]

    library_sizes = counts.sum(axis=0).replace(0, np.nan)
    log_cpm = np.log2(counts.divide(library_sizes, axis=1) * 1_000_000 + 0.5)
    rankings_by_fold: dict[tuple[str, int], list[dict[str, str]]] = defaultdict(list)
    for row in local_rankings:
        rankings_by_fold[(row["target"], int(row["fold"]))].append(row)
    audit_by_fold = {(row["target"], int(row["fold"])): row for row in local_audit}

    joint_rows: list[dict] = []
    audit_rows: list[dict] = []
    threshold = float(config["method"]["minimum_cross_lineage_log2_cpm_difference"])
    local_threshold = float(config["method"]["minimum_within_lineage_log2_fold_change"])
    top_n = int(config["method"]["fold_top_genes"])
    for index, plan in enumerate(plan_rows, 1):
        target = plan["target"]
        fold = int(plan["fold"])
        selected = metadata["role"] == discovery_role
        if population_info.loc[target, "role"] == heldout_role:
            selected = selected | (metadata["population"] == target)
        heldout = [value for value in plan["heldout_units"].split(";") if value]
        if heldout:
            selected = selected & ~metadata["library_id"].isin(heldout)
        md = metadata[selected].copy()
        target_ids = list(md.loc[md["population"] == target, "library_id"])
        if not target_ids:
            raise RuntimeError(f"No training target libraries for {target} fold {fold}")
        target_profile = log_cpm[target_ids].mean(axis=1)
        other = md[md["table_s1_lineage"] != plan["lineage"]]
        lineage_profiles = []
        for lineage in sorted(other["table_s1_lineage"].unique()):
            population_profiles = []
            lineage_md = other[other["table_s1_lineage"] == lineage]
            for population in sorted(lineage_md["population"].unique()):
                ids = list(lineage_md.loc[lineage_md["population"] == population, "library_id"])
                population_profiles.append(log_cpm[ids].mean(axis=1))
            lineage_profiles.append(pd.concat(population_profiles, axis=1).mean(axis=1))
        global_reference = pd.concat(lineage_profiles, axis=1).mean(axis=1)
        cross_effect = target_profile - global_reference

        local = rankings_by_fold[(target, fold)]
        eligible = []
        for row in local:
            gene = row["gene_symbol"]
            global_effect = float(cross_effect.loc[gene])
            local_effect = float(row["logFC"])
            if global_effect <= threshold or local_effect <= local_threshold:
                continue
            eligible.append({
                "target": target, "source_role": plan["source_role"], "lineage": plan["lineage"],
                "independent_units": plan["independent_units"], "evidence_tier": plan["evidence_tier"],
                "fold": fold, "heldout_units": plan["heldout_units"], "gene_symbol": gene,
                "cross_lineage_log2_cpm_difference": global_effect,
                "within_lineage_logFC": local_effect, "within_lineage_F": float(row["F"]),
                "within_lineage_FDR": float(row["FDR"]),
                "within_lineage_strict_significant": str(row["strict_significant"]).lower(),
                "minimum_effect": min(global_effect, local_effect),
            })
        eligible.sort(key=lambda row: (
            -int(row["within_lineage_strict_significant"] == "true"),
            -row["minimum_effect"], -row["within_lineage_F"], row["within_lineage_FDR"], row["gene_symbol"]))
        eligible = eligible[:top_n]
        for rank, row in enumerate(eligible, 1):
            row["fold_rank"] = rank
            joint_rows.append(row)
        prior = audit_by_fold[(target, fold)]
        audit_rows.append({
            "target": target, "source_role": plan["source_role"], "lineage": plan["lineage"],
            "independent_units": plan["independent_units"], "evidence_tier": plan["evidence_tier"],
            "fold": fold, "heldout_units": plan["heldout_units"],
            "local_top_genes": prior["top_genes"], "joint_positive_eligible_genes": len(eligible),
            "cross_lineage_reference_lineages": len(lineage_profiles),
            "cross_lineage_method": "population_balanced_mean_log2_cpm",
        })
        print(f"prepared {index} of {len(plan_rows)} {target} fold {fold}: {len(local)} local -> {len(eligible)} joint", flush=True)

    intermediate_root = Path(config["outputs"]["intermediate_root"])
    output_specs = config["outputs"]["intermediate_tables"]
    write_tsv(intermediate_root / output_specs["fold_plan"], plan_rows, list(plan_rows[0]))
    joint_fields = [
        "target", "source_role", "lineage", "independent_units", "evidence_tier", "fold", "heldout_units",
        "gene_symbol", "cross_lineage_log2_cpm_difference", "within_lineage_logFC", "within_lineage_F",
        "within_lineage_FDR", "within_lineage_strict_significant", "minimum_effect", "fold_rank",
    ]
    write_tsv(intermediate_root / output_specs["fold_rankings"], joint_rows, joint_fields, compressed=True)
    write_tsv(intermediate_root / output_specs["model_audit"], audit_rows, list(audit_rows[0]))
    summary = {
        "schema_version": "bdev.immgen_global_local_model_summary.v2",
        "status": "complete_adaptive_global_screen",
        "run_identity": config["run_identity"],
        "targets": len({row["target"] for row in plan_rows}),
        "folds": len(plan_rows),
        "folds_with_joint_candidates": len({(row["target"], row["fold"]) for row in joint_rows}),
        "single_cell_expression_read": False,
        "scworkbench_invoked": False,
        "reused_local_crossfit_rankings": True,
    }
    report = Path(config["outputs"]["results_root"]) / config["outputs"]["reports"]["model_summary"]
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
