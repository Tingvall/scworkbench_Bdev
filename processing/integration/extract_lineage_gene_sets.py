#!/usr/bin/env python3
"""Extract the exact Table S1B B/T lineage gene-list rows from the source XLSX."""

import argparse
import csv
import os
import tempfile

import openpyxl


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("workbook")
    parser.add_argument("output_tsv")
    args = parser.parse_args()

    workbook = openpyxl.load_workbook(
        args.workbook, read_only=True, data_only=True
    )
    sheet_name = "Table S1B"
    if sheet_name not in workbook.sheetnames:
        raise RuntimeError(f"Missing expected worksheet: {sheet_name}")
    sheet = workbook[sheet_name]
    if sheet.cell(4, 1).value != "B lineage" or sheet.cell(4, 2).value != "T lineage":
        raise RuntimeError("Unexpected Table S1B lineage headers")

    records = []
    for source_row, row in enumerate(
        sheet.iter_rows(min_row=5, max_col=2, values_only=True), start=5
    ):
        for lineage, value in zip(("B lineage", "T lineage"), row):
            if value is None or not str(value).strip():
                continue
            records.append(
                {
                    "source_sheet": sheet_name,
                    "source_row": source_row,
                    "lineage": lineage,
                    "feature_id": str(value).strip(),
                }
            )
    counts = {
        lineage: sum(record["lineage"] == lineage for record in records)
        for lineage in ("B lineage", "T lineage")
    }
    if counts != {"B lineage": 254, "T lineage": 234}:
        raise RuntimeError(f"Unexpected source row counts: {counts}")

    output_dir = os.path.dirname(os.path.abspath(args.output_tsv))
    os.makedirs(output_dir, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix="lineage_genes_", suffix=".tsv", dir=output_dir)
    try:
        with os.fdopen(fd, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("source_sheet", "source_row", "lineage", "feature_id"),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(records)
        os.replace(temporary, args.output_tsv)
    except Exception:
        if os.path.exists(temporary):
            os.unlink(temporary)
        raise


if __name__ == "__main__":
    main()
