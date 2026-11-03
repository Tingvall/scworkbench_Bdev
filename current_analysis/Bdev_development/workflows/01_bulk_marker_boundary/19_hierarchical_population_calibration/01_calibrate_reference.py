#!/usr/bin/env python3
"""Calibrate diagnostic population margins from hierarchical ImmGen holdouts."""

from __future__ import annotations

import argparse
import hashlib
import json
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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_gmt(path: Path) -> dict[str, list[str]]:
    output = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        require(len(fields) >= 3 and fields[0] not in output, f"Invalid GMT: {path}")
        output[fields[0]] = list(dict.fromkeys(fields[2:]))
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    require(config["schema_version"] == "bdev.immgen_hierarchical_population_calibration.v1", "Unexpected schema")
    paths = {name: Path(spec["path"]) for name, spec in config["inputs"].items()}
    for name, path in paths.items():
        require(path.is_file(), f"Missing input {name}: {path}")
        require(sha256(path) == config["inputs"][name]["sha256"], f"Checksum mismatch: {name}")

    results = Path(config["storage"]["results_root"])
    tables = results / "01_tables"
    reports = results / "03_reports"
    manifests = results / "06_manifests"
    require(not results.exists(), f"Refusing to overwrite run: {results}")
    for directory in (tables, reports, manifests):
        directory.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(paths["selected_counts"], sep="\t", index_col=0)
    metadata = pd.read_csv(paths["metadata"], sep="\t", dtype=str).fillna("")
    require(list(counts.columns) == list(metadata["library_id"]), "Counts and metadata are misaligned")
    ids = list(metadata["library_id"])
    max_rank = int(config["reference_test"]["ucell_max_rank"])
    ranks = counts.rank(axis=0, method="average", ascending=False).clip(upper=max_rank)

    def score(genes: list[str]) -> pd.Series:
        n_genes = len(genes)
        baseline = n_genes * (n_genes + 1) / 2
        denominator = n_genes * max_rank - baseline
        raw = (1 - (ranks.loc[genes, ids].sum(axis=0) - baseline) / denominator).clip(0, 1)
        return raw.rank(method="average", pct=True)

    membership = pd.read_csv(paths["population_membership"], sep="\t")
    population_summary = pd.read_csv(paths["population_summary"], sep="\t", dtype=str).fillna("")
    population_order = list(dict.fromkeys(population_summary["target"]))
    lineage = dict(zip(population_summary["target"], population_summary["lineage"]))
    population_sets = {
        target: list(rows.sort_values("final_rank")["gene_symbol"])
        for target, rows in membership.groupby("target", sort=False)
    }
    require(set(population_order) == set(population_sets), "Population membership is incomplete")
    population_scores = pd.DataFrame({target: score(population_sets[target]) for target in population_order}, index=ids)

    broad_sets = read_gmt(paths["broad_gmt"])
    broad_map = pd.read_csv(paths["broad_map"], sep="\t")
    broad_lineage = dict(zip(broad_map["broad_set_id"], broad_map["lineage"]))
    require(set(broad_sets) == set(broad_lineage), "Broad mapping is incomplete")
    broad_scores = pd.DataFrame({set_id: score(genes) for set_id, genes in broad_sets.items()}, index=ids)
    broad_minimum = float(config["reference_test"]["broad_minimum_score"])
    broad_delta = float(config["reference_test"]["broad_maximum_delta_from_top"])
    eligible_states = {}
    for library_id in ids:
        values = broad_scores.loc[library_id]
        top = float(values.max())
        eligible_states[library_id] = sorted(
            broad_lineage[set_id]
            for set_id, value in values.items()
            if value >= broad_minimum and top - value <= broad_delta + 1e-12
        )

    plans = pd.read_csv(paths["fold_plan"], sep="\t", keep_default_na=False)
    rankings = pd.read_csv(paths["fold_rankings"], sep="\t")
    minimum_score = float(config["reference_test"]["population_minimum_score"])
    tie_tolerance = float(config["reference_test"]["tie_tolerance"])
    decision_rows = []
    for plan in plans.itertuples(index=False):
        heldout = [value for value in plan.heldout_units.split(";") if value]
        if not heldout:
            continue
        fold_genes = list(
            rankings[(rankings["target"] == plan.target) & (rankings["fold"] == plan.fold)]
            .sort_values("fold_rank")["gene_symbol"]
        )
        require(bool(fold_genes), f"Missing fold program: {plan.target} fold {plan.fold}")
        fold_scores = score(fold_genes)
        for library_id in heldout:
            values = population_scores.loc[library_id].copy()
            values.loc[plan.target] = float(fold_scores.loc[library_id])
            states = eligible_states[library_id]
            local_targets = [target for target in population_order if lineage[target] in states]
            if not local_targets:
                winner = ""
                top_score = runner = margin = np.nan
                reason = "broad_unassigned"
            else:
                ordered = values.loc[local_targets].sort_values(ascending=False, kind="mergesort")
                winner = str(ordered.index[0])
                top_score = float(ordered.iloc[0])
                runner = float(ordered.iloc[1]) if len(ordered) > 1 else np.nan
                margin = top_score - runner if np.isfinite(runner) else np.inf
                if top_score < minimum_score:
                    reason = "min_score"
                elif margin <= tie_tolerance:
                    reason = "tie"
                else:
                    reason = "none"
            decision_rows.append({
                "target": plan.target,
                "fold": plan.fold,
                "heldout_library": library_id,
                "true_lineage": lineage[plan.target],
                "eligible_broad_states": ";".join(states) if states else "unassigned",
                "true_lineage_eligible": lineage[plan.target] in states,
                "winner_population": winner,
                "top_score": top_score,
                "runner_up_score": runner,
                "score_margin": margin,
                "base_reason": reason,
                "winner_correct": winner == plan.target,
            })
    decisions = pd.DataFrame(decision_rows)
    decisions.to_csv(tables / "01_reference_hierarchical_holdout_decisions.tsv", sep="\t", index=False)

    calibration = config["reference_test"]["calibration"]
    minimum_predictions = int(calibration["minimum_predictions"])
    required_precision = float(calibration["required_observed_precision"])
    threshold_rows = []
    for target in population_order:
        predicted = decisions[
            (decisions["winner_population"] == target)
            & (decisions["base_reason"] == "none")
        ].copy()
        selected = None
        for threshold in sorted(predicted["score_margin"].unique()) if len(predicted) else []:
            accepted = predicted[predicted["score_margin"] >= threshold - 1e-12]
            precision = float(accepted["winner_correct"].mean())
            if len(accepted) >= minimum_predictions and precision >= required_precision:
                selected = (float(threshold), len(accepted), int(accepted["winner_correct"].sum()), precision)
                break
        threshold_rows.append({
            "target": target,
            "lineage": lineage[target],
            "heldout_predictions": len(predicted),
            "heldout_correct": int(predicted["winner_correct"].sum()),
            "observed_precision": float(predicted["winner_correct"].mean()) if len(predicted) else np.nan,
            "threshold_available": selected is not None,
            "calibrated_minimum_margin": selected[0] if selected else np.nan,
            "accepted_calibration_predictions": selected[1] if selected else 0,
            "accepted_calibration_correct": selected[2] if selected else 0,
            "accepted_observed_precision": selected[3] if selected else np.nan,
            "calibration_limitation": calibration["limitation"],
        })
    thresholds = pd.DataFrame(threshold_rows)
    thresholds.to_csv(tables / "02_population_margin_calibration.tsv", sep="\t", index=False)

    usable = decisions[decisions["base_reason"] == "none"]
    summary = {
        "schema_version": "bdev.immgen_hierarchical_reference_calibration_summary.v1",
        "status": "complete_diagnostic_only",
        "run_identity": config["run_identity"],
        "heldout_libraries": len(decisions),
        "true_lineage_eligible": int(decisions["true_lineage_eligible"].sum()),
        "true_lineage_eligible_fraction": float(decisions["true_lineage_eligible"].mean()),
        "no_margin_nontie_predictions": len(usable),
        "no_margin_correct": int(usable["winner_correct"].sum()),
        "no_margin_precision": float(usable["winner_correct"].mean()),
        "populations_with_calibrated_threshold": int(thresholds["threshold_available"].sum()),
        "calibration_minimum_predictions": minimum_predictions,
        "calibration_required_observed_precision": required_precision,
        "calibration_is_independent_validation": False,
        "accepted_labels_emitted": False,
    }
    (reports / "01_reference_calibration_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    manifest = {
        "schema_version": "bdev.immgen_hierarchical_reference_input_manifest.v1",
        "run_identity": config["run_identity"],
        "config": {"path": str(config_path), "sha256": sha256(config_path)},
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in paths.items()},
    }
    (manifests / "01_reference_input_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
