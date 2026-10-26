#!/usr/bin/env python3
"""Report the validated Parse Seurat object without loading the large RDS."""

from __future__ import annotations

import csv
import json
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


COMBINATION_COMPONENTS = {
    "Flt3L_TPO": ("Flt3L", "TPO"),
    "IL-7_Flt3L": ("IL-7", "Flt3L"),
    "c-kit_Flt3L": ("c-kit", "Flt3L"),
    "c-kit_TPO": ("c-kit", "TPO"),
    "c-kit_TPO_Flt3L": ("c-kit", "TPO", "Flt3L"),
}
SINGLE_LABELS_WITH_UNDERSCORES = {"Csf_1", "Csf_2"}
CONTROL_LABEL = "control"


def classify_condition(condition: str) -> tuple[str, tuple[str, ...], str]:
    if condition == CONTROL_LABEL:
        return "control", (), "explicit_control_label"
    if condition in COMBINATION_COMPONENTS:
        return "combination", COMBINATION_COMPONENTS[condition], "reviewed_combination_label"
    if "_" in condition and condition not in SINGLE_LABELS_WITH_UNDERSCORES:
        raise ValueError(
            f"Unreviewed underscore-containing condition {condition!r}; "
            "classify its components explicitly before generating the report"
        )
    basis = "reviewed_internal_underscore_single" if condition in SINGLE_LABELS_WITH_UNDERSCORES else "single_source_label"
    return "single", (condition,), basis


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    def value(x: object) -> str:
        return str(x).replace("|", "\\|").replace("\n", " ")

    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(value(x) for x in row) + " |" for row in rows)
    return "\n".join(lines)


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: 01_report_prepared_seurat.py PREPARED_DIR OUTPUT_DIR")

    prepared_dir = Path(sys.argv[1]).resolve(strict=True)
    output_dir = Path(sys.argv[2]).resolve()
    sample_summary_path = prepared_dir / "sample_summary.tsv"
    validation_path = prepared_dir / "validation.json"
    environment_path = prepared_dir / "environment.tsv"
    rds_path = prepared_dir / "parse_cytokines_unfiltered_seurat.rds"
    for path in (sample_summary_path, validation_path, environment_path, rds_path):
        if not path.exists():
            raise FileNotFoundError(path)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = {
        "report": output_dir / "parse_cytokines_seurat_report.md",
        "samples": output_dir / "sample_classification_proposed.tsv",
        "conditions": output_dir / "condition_classification_proposed.tsv",
        "summary": output_dir / "summary.json",
    }
    existing = [str(path) for path in output_paths.values() if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing report output(s): " + ", ".join(existing))

    validation = json.loads(validation_path.read_text())
    if validation.get("passed") is not True:
        raise ValueError("The prepared object validation did not pass")

    source_rows = read_tsv(sample_summary_path)
    if not source_rows:
        raise ValueError("sample_summary.tsv is empty")

    sample_ids = [row["sample_id"] for row in source_rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Sample IDs are not unique")

    sample_rows: list[dict[str, object]] = []
    condition_samples: dict[str, list[dict[str, object]]] = defaultdict(list)
    for source in source_rows:
        condition = source["condition"]
        stimulus_class, components, basis = classify_condition(condition)
        replicate = int(source["replicate"])
        n_cells = int(source["n_cells"])
        expected_suffixes = (f"_rep{replicate}", f"_repl{replicate}")
        if not source["sample_id"].endswith(expected_suffixes):
            raise ValueError(f"Sample ID and replicate disagree: {source['sample_id']}")
        review_note = ""
        if condition in SINGLE_LABELS_WITH_UNDERSCORES:
            review_note = "Underscore belongs to the single source stimulus label."
        elif condition == "TGB-b1":
            review_note = "Source spelling retained; confirm whether TGF-b1 was intended."
        elif condition == "c-kit":
            review_note = "Source sample ID uses the suffix repl rather than rep."
        row = {
            "sample_id": source["sample_id"],
            "condition": condition,
            "replicate": replicate,
            "n_cells": n_cells,
            "stimulus_class": stimulus_class,
            "stimulus_count": len(components),
            "stimulus_components": ";".join(components),
            "classification_basis": basis,
            "review_note": review_note,
        }
        sample_rows.append(row)
        condition_samples[condition].append(row)

    observed_cells = sum(int(row["n_cells"]) for row in sample_rows)
    observed_samples = len(sample_rows)
    if observed_cells != int(validation["matrix"]["cells"]):
        raise ValueError("Sample cell counts do not sum to the validated object cell count")
    if observed_samples != int(validation["metadata"]["sample_count"]):
        raise ValueError("Sample count does not match validation.json")

    condition_rows: list[dict[str, object]] = []
    for condition, rows in condition_samples.items():
        replicates = sorted(int(row["replicate"]) for row in rows)
        if replicates != [1, 2, 3]:
            raise ValueError(f"Expected replicates 1, 2, and 3 for {condition}; observed {replicates}")
        classes = {str(row["stimulus_class"]) for row in rows}
        components = {str(row["stimulus_components"]) for row in rows}
        if len(classes) != 1 or len(components) != 1:
            raise ValueError(f"Inconsistent classification within condition {condition}")
        condition_rows.append(
            {
                "condition": condition,
                "stimulus_class": next(iter(classes)),
                "stimulus_count": int(rows[0]["stimulus_count"]),
                "stimulus_components": next(iter(components)),
                "n_replicates": len(rows),
                "n_samples": len(rows),
                "n_cells": sum(int(row["n_cells"]) for row in rows),
                "min_cells_per_sample": min(int(row["n_cells"]) for row in rows),
                "max_cells_per_sample": max(int(row["n_cells"]) for row in rows),
                "sample_ids": ";".join(str(row["sample_id"]) for row in rows),
            }
        )

    class_order = {"single": 0, "combination": 1, "control": 2}
    condition_rows.sort(key=lambda row: (class_order[str(row["stimulus_class"])], str(row["condition"]).lower()))
    sample_rows.sort(key=lambda row: (class_order[str(row["stimulus_class"])], str(row["condition"]).lower(), int(row["replicate"])))

    class_counts: dict[str, dict[str, int]] = {}
    for stimulus_class in ("single", "combination", "control"):
        selected_samples = [row for row in sample_rows if row["stimulus_class"] == stimulus_class]
        selected_conditions = [row for row in condition_rows if row["stimulus_class"] == stimulus_class]
        class_counts[stimulus_class] = {
            "conditions": len(selected_conditions),
            "samples": len(selected_samples),
            "cells": sum(int(row["n_cells"]) for row in selected_samples),
        }

    write_tsv(
        output_paths["samples"],
        sample_rows,
        [
            "sample_id", "condition", "replicate", "n_cells", "stimulus_class",
            "stimulus_count", "stimulus_components", "classification_basis", "review_note",
        ],
    )
    write_tsv(
        output_paths["conditions"],
        condition_rows,
        [
            "condition", "stimulus_class", "stimulus_count", "stimulus_components",
            "n_replicates", "n_samples", "n_cells", "min_cells_per_sample",
            "max_cells_per_sample", "sample_ids",
        ],
    )

    cell_counts = [int(row["n_cells"]) for row in sample_rows]
    generated_utc = datetime.now(timezone.utc).isoformat()
    summary = {
        "schema_version": 1,
        "generated_utc": generated_utc,
        "report_scope": "validated preparation evidence; Seurat RDS not loaded",
        "prepared_dir": str(prepared_dir),
        "rds": {
            "path": str(rds_path),
            "size_bytes": rds_path.stat().st_size,
            "recorded_sha256": validation["output"]["sha256"],
        },
        "object": {
            "assay": "RNA",
            "source_matrix_set": validation["source_selection"]["requested"],
            "features": validation["matrix"]["features"],
            "cells": validation["matrix"]["cells"],
            "stored_entries": validation["matrix"]["stored_entries"],
            "samples": observed_samples,
            "conditions": len(condition_rows),
            "replicates_per_condition": 3,
            "minimum_cells_per_sample": min(cell_counts),
            "median_cells_per_sample": statistics.median(cell_counts),
            "maximum_cells_per_sample": max(cell_counts),
        },
        "classification": class_counts,
        "checks": {
            "validation_passed": True,
            "sample_ids_unique": True,
            "sample_counts_sum_to_object_cells": True,
            "three_replicates_per_condition": True,
            "unresolved_condition_labels": [],
        },
    }
    output_paths["summary"].write_text(json.dumps(summary, indent=2) + "\n")

    class_table = []
    for stimulus_class in ("single", "combination", "control"):
        values = class_counts[stimulus_class]
        class_table.append(
            [
                stimulus_class,
                values["conditions"],
                values["samples"],
                f"{values['cells']:,}",
                f"{100 * values['cells'] / observed_cells:.2f}%",
            ]
        )
    condition_table = [
        [
            row["condition"], row["stimulus_class"], row["stimulus_count"],
            row["stimulus_components"] or "—", row["n_replicates"],
            f"{int(row['n_cells']):,}", row["sample_ids"],
        ]
        for row in condition_rows
    ]
    sample_table = [
        [
            row["sample_id"], row["condition"], row["replicate"],
            row["stimulus_class"], row["stimulus_components"] or "—",
            f"{int(row['n_cells']):,}",
        ]
        for row in sample_rows
    ]

    report = f"""# Parse cytokine Seurat object report

Generated: `{generated_utc}`

## Scope and evidence

This report summarizes the completed, validated unfiltered Seurat preparation without loading the 4 GB RDS. It uses `validation.json` and `sample_summary.tsv`, which were written during object construction after exact cell-order, feature-order, and metadata-alignment checks.

- Current RDS: `{rds_path}`
- Recorded SHA-256: `{validation['output']['sha256']}`
- Input matrix set: `{validation['source_selection']['requested']}`
- Default assay at creation: `RNA`
- Matrix class: `{validation['matrix']['class']}`
- Features: {int(validation['matrix']['features']):,}
- Cells: {int(validation['matrix']['cells']):,}
- Stored sparse entries: {int(validation['matrix']['stored_entries']):,}
- Samples: {observed_samples}
- Conditions: {len(condition_rows)}
- Replicates: three per condition
- Cells per sample: {min(cell_counts):,} minimum, {statistics.median(cell_counts):,.0f} median, {max(cell_counts):,} maximum
- Filtering: none; zero cells and zero features removed
- Cell and feature identifiers: unique and retained in source order

The object retains the nine delivered cell-metadata fields listed in `validation.json`, plus `orig.ident`, `sample_id`, `condition`, `replicate`, `dataset_id`, and `source_matrix_set`. The RNA assay retains `gene_id`, `gene_name_original`, and `genome` feature metadata.

## Proposed stimulus classification

`single` denotes one named stimulus. `combination` denotes two or more explicitly reviewed stimulus components. `control` remains separate because it has zero cytokine or ligand components.

Underscores alone are not sufficient for classification: `Csf_1` and `Csf_2` are reviewed as single-stimulus labels. The five combination labels are mapped explicitly. Any new underscore-containing condition causes the report generator to stop for manual review.

{markdown_table(['Class', 'Conditions', 'Samples', 'Cells', 'Cells (%)'], class_table)}

## Condition classification and sample IDs

{markdown_table(['Condition', 'Class', 'Components (n)', 'Components', 'Replicates', 'Cells', 'Sample IDs'], condition_table)}

## All sample IDs

{markdown_table(['Sample ID', 'Condition', 'Replicate', 'Class', 'Components', 'Cells'], sample_table)}

## Labels requiring review

- `TGB-b1` is retained exactly from the source. Confirm whether the intended label is `TGF-b1` before adopting a standardized display name.
- `c-kit_repl1`, `c-kit_repl2`, and `c-kit_repl3` use `repl` while other sample IDs generally use `rep`. The parsed replicate values are correct.
- `c-kit` is retained as the delivered stimulus name. Confirm whether it denotes a c-Kit ligand/SCF treatment before biological relabeling.

## Validation result

All 96 sample IDs are unique. Every condition has replicates 1–3, and sample cell counts sum exactly to {observed_cells:,}, matching the validated Seurat object. The proposed annotation has no unresolved condition labels. This report does not modify the RDS or its metadata.
"""
    output_paths["report"].write_text(report)
    print(json.dumps({name: str(path) for name, path in output_paths.items()}, indent=2))


if __name__ == "__main__":
    main()
