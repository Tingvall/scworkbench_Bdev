#!/usr/bin/env python3
"""Add user-authorized small two-fold panels to the adaptive ImmGen collection."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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
                    writer.writeheader()
                    writer.writerows(rows)
    else:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    require(config["schema_version"] == "bdev.immgen_small_panel_parse_exception.v1", "Unexpected config schema")
    config_dir = args.config.parent.resolve()
    inputs = {name: resolve(config_dir, spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in inputs.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        require(sha256(path) == config["inputs"][name]["sha256"], f"Checksum mismatch: {name}")

    output_root = Path(config["outputs"]["root"])
    tables = {name: output_root / rel for name, rel in config["outputs"]["tables"].items()}
    reports = {name: output_root / rel for name, rel in config["outputs"]["reports"].items()}
    counts = pd.read_csv(inputs["selected_counts"], sep="\t", index_col=0)
    metadata = pd.read_csv(inputs["metadata"], sep="\t", dtype=str).fillna("")
    require(list(counts.columns) == list(metadata["library_id"]), "Counts and metadata are misaligned")
    discovery = metadata[metadata["role"] == "discovery"].copy()
    heldout = metadata[metadata["role"] == "held_out_validation"].copy()
    max_rank = int(config["scoring"]["max_rank"])
    capped_ranks = counts.rank(axis=0, method="average", ascending=False).clip(upper=max_rank)
    approved_targets: dict[str, list[str]] = defaultdict(list)
    for row in read_tsv(inputs["contrast_plan"]):
        if row["target"] not in approved_targets[row["family"]]:
            approved_targets[row["family"]].append(row["target"])
    group_column = {"superlineage": "superlineage", "broad_lineage": "table_s1_lineage", "b_branch": "b_branch", "b_local_stage": "population"}

    def family_metadata(family: str) -> pd.DataFrame:
        local = discovery.copy()
        if family in {"b_branch", "b_local_stage"}:
            local = local[local["table_s1_lineage"] == "B"]
        if family == "b_local_stage":
            local = local[local["population"].isin(approved_targets[family])]
        return local

    def score(genes: list[str], library_ids: list[str]) -> pd.Series:
        n = len(genes)
        baseline = n * (n + 1) / 2
        denominator = n * max_rank - baseline
        require(n >= int(config["policy"]["small_panel_minimum_genes"]) and denominator > 0, "Invalid small-panel score")
        rank_sum = capped_ranks.loc[genes, library_ids].sum(axis=0)
        return (1 - (rank_sum - baseline) / denominator).clip(lower=0, upper=1)

    def unit_scores_and_groups(family: str, library_scores: pd.Series, md: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
        local = md.set_index("library_id").loc[list(library_scores.index)]
        if family == "b_local_stage":
            return library_scores, local["population"]
        units = list(dict.fromkeys(local["population"]))
        scores = pd.Series({unit: float(library_scores.loc[local.index[local["population"] == unit]].mean()) for unit in units})
        unit_md = local.drop_duplicates("population").set_index("population").loc[units]
        return scores, unit_md[group_column[family]]

    def metrics(family: str, target: str, library_scores: pd.Series, md: pd.DataFrame) -> dict:
        unit_scores, groups = unit_scores_and_groups(family, library_scores, md)
        values = {group: float(np.median(unit_scores.loc[groups.index[groups == group]])) for group in dict.fromkeys(groups)}
        competitors = {group: value for group, value in values.items() if group != target}
        strongest = max(competitors, key=lambda group: (competitors[group], group))
        sd = float(np.std(list(values.values()), ddof=1))
        gap = (values[target] - competitors[strongest]) / sd if sd > 0 else float("-inf")
        target_scores = unit_scores.loc[groups.index[groups == target]]
        return {"gap": gap, "strongest_competitor": strongest, "target_rank": 1 + sum(value >= values[target] for value in competitors.values()), "target_unit_win_fraction": float(np.mean(target_scores > competitors[strongest]))}

    prior_summary = read_tsv(inputs["adaptive_target_summary"])
    failed = {(row["family"], row["target"]): row for row in prior_summary if row["status"] == "candidate_below_20"}
    require(len(failed) == 2, f"Expected two eligible failures, found {len(failed)}")
    candidate_rows = read_tsv(inputs["low_replication_candidates"])
    small_sets: dict[tuple[str, str], list[tuple[int, str]]] = defaultdict(list)
    for row in candidate_rows:
        key = (row["family"], row["target"])
        if key in failed:
            small_sets[key].append((int(row["candidate_rank"]), row["gene_symbol"]))
    minimum = int(config["policy"]["small_panel_minimum_genes"])
    standard = int(config["policy"]["standard_minimum_genes"])
    for key, members in small_sets.items():
        require(minimum <= len(members) < standard, f"Small panel outside policy: {key} has {len(members)} genes")
        require(failed[key]["evidence_tier"] == config["policy"]["eligible_evidence_tier"], f"Wrong evidence tier: {key}")

    complete_sets: dict[tuple[str, str], list[tuple[int, str, dict]]] = defaultdict(list)
    for row in read_tsv(inputs["adaptive_final_membership"]):
        key = (row["family"], row["target"])
        complete_sets[key].append((int(row["final_rank"]), row["gene_symbol"], row))
    for key, members in small_sets.items():
        family, target = key
        for rank, gene in sorted(members):
            complete_sets[key].append((rank, gene, {"family": family, "target": target, "independent_units": failed[key]["independent_units"], "evidence_tier": config["policy"]["evidence_label"], "candidate_genes": len(members), "final_genes": len(members), "status": "small_panel_exception_retained_complete_intersection"}))
    require(len(complete_sets) == 27, f"Expected 27 complete sets, found {len(complete_sets)}")

    small_rows: list[dict] = []
    complete_rows: list[dict] = []
    target_rows: list[dict] = []
    for key in [(family, target) for family in ("superlineage", "broad_lineage", "b_branch", "b_local_stage") for target in approved_targets[family]]:
        family, target = key
        members = sorted(complete_sets[key], key=lambda item: item[0])
        genes = [gene for _, gene, _ in members]
        md = family_metadata(family)
        result = metrics(family, target, score(genes, list(md["library_id"])), md)
        meta = members[0][2]
        evidence = meta["evidence_tier"]
        status = meta["status"]
        independent_units = int(meta["independent_units"])
        candidate_genes = int(meta["candidate_genes"])
        for final_rank, (_, gene, _) in enumerate(members, 1):
            row = {"family": family, "target": target, "independent_units": independent_units, "evidence_tier": evidence, "gene_symbol": gene, "final_rank": final_rank, "candidate_genes": candidate_genes, "final_genes": len(genes), "status": status, "small_panel": str(key in small_sets).lower()}
            complete_rows.append(row)
            if key in small_sets:
                small_rows.append(row)
        target_rows.append({"family": family, "target": target, "independent_units": independent_units, "evidence_tier": evidence, "genes": len(genes), "status": status, "small_panel": str(key in small_sets).lower(), "ucell_gap": f"{result['gap']:.12g}", "strongest_competitor": result["strongest_competitor"], "target_rank": result["target_rank"], "target_unit_win_fraction": f"{result['target_unit_win_fraction']:.8f}", "full_data_mean_pass": str(result["target_rank"] == 1).lower()})

    write_tsv(tables["small_panel_membership"], small_rows, list(small_rows[0]))
    write_tsv(tables["complete_membership"], complete_rows, list(complete_rows[0]), compressed=True)
    write_tsv(tables["target_summary"], target_rows, list(target_rows[0]))

    prior_comparison = read_tsv(inputs["prior_target_comparison"])
    comparison_rows = list(prior_comparison)
    for row in target_rows:
        comparison_rows.append({"method": "replication_adaptive_small_panel", "family": row["family"], "target": row["target"], "independent_target_units": row["independent_units"], "genes": row["genes"], "status": row["status"], "ucell_gap": row["ucell_gap"], "strongest_competitor": row["strongest_competitor"], "target_rank": row["target_rank"], "target_unit_win_fraction": row["target_unit_win_fraction"], "full_data_mean_pass": row["full_data_mean_pass"], "selection_uses_full_discovery_scores": "true", "evidence_tier": row["evidence_tier"]})
    write_tsv(tables["method_target_comparison"], comparison_rows, list(comparison_rows[0]))

    prior_heldout = read_tsv(inputs["prior_heldout_comparison"])
    heldout_rows = list(prior_heldout)
    broad_sets = {target: [gene for _, gene, _ in sorted(rows)] for (family, target), rows in complete_sets.items() if family == "broad_lineage"}
    expected = heldout.drop_duplicates("population").set_index("population")["table_s1_lineage"].to_dict()
    for population in dict.fromkeys(heldout["population"]):
        ids = list(heldout.loc[heldout["population"] == population, "library_id"])
        scores = {target: float(score(genes, ids).mean()) for target, genes in broad_sets.items()}
        predicted = max(scores, key=lambda target: (scores[target], target))
        heldout_rows.append({"method": "replication_adaptive_small_panel", "population": population, "expected_lineage": expected[population], "predicted_lineage": predicted, "correct": str(predicted == expected[population]).lower(), "winning_score": f"{scores[predicted]:.12g}", "broad_panel_sets": len(broad_sets), "scores_json": json.dumps(scores, sort_keys=True, separators=(",", ":"))})
    write_tsv(tables["heldout_broad_comparison"], heldout_rows, list(heldout_rows[0]))

    prior_methods = read_tsv(inputs["prior_method_summary"])
    method_fields = list(prior_methods[0]) + ["small_twofold_sets"]
    method_rows = [{**row, "small_twofold_sets": ""} for row in prior_methods]
    new_heldout = [row for row in heldout_rows if row["method"] == "replication_adaptive_small_panel"]
    method_rows.append({"method": "replication_adaptive_small_panel", "targets": 27, "targets_with_set": 27, "full_data_mean_passes": sum(row["full_data_mean_pass"] == "true" for row in target_rows), "heldout_broad_correct": sum(row["correct"] == "true" for row in new_heldout), "heldout_broad_total": len(new_heldout), "broad_panel_sets": len(broad_sets), "single_cell_expression_used": "false", "high_threefold_sets": 14, "moderate_twofold_sets": 8, "provisional_singleton_sets": 3, "small_twofold_sets": len(small_sets)})
    write_tsv(tables["method_summary"], method_rows, method_fields)

    summary = {"schema_version": "bdev.immgen_small_panel_parse_exception_summary.v1", "status": "complete_diagnostic_only", "run_identity": config["run_identity"], "targets_with_set": 27, "full_data_mean_passes": sum(row["full_data_mean_pass"] == "true" for row in target_rows), "small_panel_targets": [{"family": key[0], "target": key[1], "genes": len(value)} for key, value in sorted(small_sets.items())], "heldout_broad_correct": sum(row["correct"] == "true" for row in new_heldout), "heldout_broad_total": len(new_heldout), "single_cell_expression_read": False, "accepted_sets_emitted": 0, "importable_artifact_created": False, "scworkbench_invoked": False}
    reports["summary"].parent.mkdir(parents=True, exist_ok=True)
    reports["summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    generated = {**tables, **reports}
    script = Path(__file__).resolve()
    manifest = {"schema_version": "bdev.immgen_small_panel_parse_exception_manifest.v1", "status": "complete_diagnostic_only", "run_identity": config["run_identity"], "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in {"small_panel_config": args.config.resolve(), **inputs}.items()}, "scripts": {script.name: {"path": str(script), "sha256": sha256(script)}}, "outputs": {name: {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in generated.items()}, "small_panel_policy_preserved": True, "single_cell_expression_read": False, "accepted_sets_emitted": 0, "importable_artifact_created": False}
    manifest_path = output_root / config["outputs"]["manifest"]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
