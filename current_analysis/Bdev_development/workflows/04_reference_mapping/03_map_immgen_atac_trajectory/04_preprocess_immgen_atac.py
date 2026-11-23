#!/usr/bin/env python3
"""Download and count ImmGen paired-end ATAC fragments on the Bdev peak axis."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import yaml


def checksum(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def merge(base: dict, override: dict) -> dict:
    result = dict(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(path: Path) -> dict:
    current = yaml.safe_load(path.read_text())
    base_name = current.pop("base_config", None)
    if base_name is None:
        return current
    return merge(load_config((path.parent / base_name).resolve()), current)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(command: list[str], *, stdout=None, stderr=None, env=None) -> None:
    subprocess.run(command, check=True, stdout=stdout, stderr=stderr, env=env)


def count_lines(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def download(url: str, output: Path, expected_md5: str) -> None:
    if output.exists() and checksum(output, "md5") == expected_md5:
        return
    run([
        "curl", "-L", "--fail", "--silent", "--show-error",
        "--retry", "50", "--retry-all-errors", "--retry-delay", "10",
        "--retry-max-time", "1800", "--connect-timeout", "60",
        "--continue-at", "-", "--output", str(output), url,
    ])
    observed = checksum(output, "md5")
    if observed != expected_md5:
        raise RuntimeError(f"MD5 mismatch for {output}: {observed} != {expected_md5}")


def tool_version(command: list[str], env=None) -> str:
    result = subprocess.run(
        command, check=True, text=True, encoding="utf-8", errors="replace",
        capture_output=True, env=env,
    )
    return (result.stdout or result.stderr).splitlines()[0]


def process_library(row: dict[str, str], config: dict, peak_bed: Path, env: dict,
                    temporary_root: Path) -> dict:
    storage = config["storage"]
    prep = config["preprocessing"]
    gsm = row["gsm"]
    run_accession = row["run_accession"]
    raw_root = Path(storage["raw_fastq"]) / gsm
    raw_root.mkdir(parents=True, exist_ok=True)
    fastq_1 = raw_root / f"{run_accession}_1.fastq.gz"
    fastq_2 = raw_root / f"{run_accession}_2.fastq.gz"
    download(row["fastq_1_url"], fastq_1, row["fastq_1_md5"])
    download(row["fastq_2_url"], fastq_2, row["fastq_2_md5"])

    counts_root = Path(storage["peak_counts"])
    counts_root.mkdir(parents=True, exist_ok=True)
    output = counts_root / f"{gsm}.bdev_peak_counts.tsv.gz"
    qc_output = counts_root / f"{gsm}.qc.json"
    if output.exists() and qc_output.exists():
        qc = json.loads(qc_output.read_text())
        if qc.get("status") == "complete" and qc.get("output_sha256") == checksum(output, "sha256"):
            return qc

    started = time.time()
    work = Path(tempfile.mkdtemp(prefix=f"{gsm}_", dir=temporary_root))
    log_path = Path(config["storage"]["results_root"]) / "04_logs" / f"04_{gsm}_alignment.log"
    raw_bam = work / "01_raw.bam"
    name_bam = work / "02_name.bam"
    fixmate_bam = work / "03_fixmate.bam"
    coord_bam = work / "04_coord.bam"
    dedup_bam = work / "05_dedup.bam"
    proper_bam = work / "06_proper_primary.bam"
    proper_name_bam = work / "07_proper_primary_name.bam"
    shifted = work / "08_tn5_fragments.bed"
    clean_fragments = work / "09_clean_fragments.bed"
    coverage = work / "10_peak_coverage.tsv"
    try:
        with log_path.open("wb") as log:
            bowtie = subprocess.Popen([
                prep["bowtie2"], "--very-sensitive", "--no-mixed", "--no-discordant",
                "-X", str(prep["maximum_fragment_length"]), "-p", str(prep["threads"]),
                "-x", config["inputs"]["archived_mm10"]["bowtie2_prefix"],
                "-1", str(fastq_1), "-2", str(fastq_2),
            ], stdout=subprocess.PIPE, stderr=log, env=env)
            view = subprocess.Popen([
                prep["samtools"], "view", "-@", str(prep["threads"]), "-b", "-o", str(raw_bam), "-"
            ], stdin=bowtie.stdout, stderr=log)
            assert bowtie.stdout is not None
            bowtie.stdout.close()
            view_code = view.wait()
            bowtie_code = bowtie.wait()
            if bowtie_code or view_code:
                raise RuntimeError(f"Alignment pipeline failed for {gsm}: bowtie2={bowtie_code}, samtools={view_code}")
            run([prep["samtools"], "sort", "-@", str(prep["threads"]), "-n", "-o", str(name_bam), str(raw_bam)], stderr=log)
            run([prep["samtools"], "fixmate", "-@", str(prep["threads"]), "-m", str(name_bam), str(fixmate_bam)], stderr=log)
            run([prep["samtools"], "sort", "-@", str(prep["threads"]), "-o", str(coord_bam), str(fixmate_bam)], stderr=log)
            run([prep["samtools"], "markdup", "-@", str(prep["threads"]), "-r", str(coord_bam), str(dedup_bam)], stderr=log)
            run([prep["samtools"], "view", "-@", str(prep["threads"]), "-b", "-f", "2", "-F", "3852", "-o", str(proper_bam), str(dedup_bam)], stderr=log)
            run([prep["samtools"], "sort", "-@", str(prep["threads"]), "-n", "-o", str(proper_name_bam), str(proper_bam)], stderr=log)

            bedpe = subprocess.Popen([prep["bedtools"], "bamtobed", "-bedpe", "-i", str(proper_name_bam)], stdout=subprocess.PIPE, stderr=log)
            # BEDTools requires tab-delimited BED records; awk otherwise uses spaces.
            awk_program = (
                'BEGIN {OFS="\\t"} '
                '$1==$4 && $1!="chrM" && $8>=q {'
                'p1=($9=="+"?$2+4:$3-5); p2=($10=="+"?$5+4:$6-5); '
                's=(p1<p2?p1:p2); e=(p1>p2?p1:p2)+1; '
                'if(s>=0 && e>s) print $1,s,e,$7}'
            )
            with shifted.open("wb") as handle:
                awk = subprocess.Popen(["awk", "-v", f"q={prep['minimum_mapping_quality']}", awk_program], stdin=bedpe.stdout, stdout=handle, stderr=log)
                assert bedpe.stdout is not None
                bedpe.stdout.close()
                awk_code = awk.wait()
                bedpe_code = bedpe.wait()
                if bedpe_code or awk_code:
                    raise RuntimeError(f"Fragment extraction failed for {gsm}")
            if shifted.stat().st_size == 0:
                raise RuntimeError(f"Fragment extraction produced no retained fragments for {gsm}")
            with shifted.open() as handle:
                first_fragment = handle.readline().rstrip("\n").split("\t")
            if len(first_fragment) != 4:
                raise RuntimeError(f"Fragment extraction produced malformed BED output for {gsm}")
            with clean_fragments.open("wb") as handle:
                run([prep["bedtools"], "intersect", "-v", "-a", str(shifted),
                     "-b", config["inputs"]["archived_mm10"]["blacklist"]], stdout=handle, stderr=log)
            with coverage.open("wb") as handle:
                run([prep["bedtools"], "coverage", "-counts", "-a", str(peak_bed),
                     "-b", str(clean_fragments)], stdout=handle, stderr=log)

        feature_count = 0
        total_peak_overlaps = 0
        with coverage.open() as source, gzip.open(output, "wt", newline="") as target:
            target.write("feature_id\tcount\n")
            for line in source:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 5:
                    raise RuntimeError(f"Malformed coverage row for {gsm}")
                count = int(fields[4])
                target.write(f"{fields[3]}\t{count}\n")
                feature_count += 1
                total_peak_overlaps += count
        if feature_count != config["inputs"]["expected_features"]:
            raise RuntimeError(f"{gsm} produced {feature_count} peak rows")
        mapped_fragments = count_lines(clean_fragments)
        qc = {
            "schema_version": "bdev.immgen_atac_library_qc.v1",
            "status": "complete",
            "gsm": gsm,
            "population": row["population"],
            "run_accession": run_accession,
            "fastq_1_bytes": fastq_1.stat().st_size,
            "fastq_2_bytes": fastq_2.stat().st_size,
            "fastq_md5_valid": True,
            "mapped_deduplicated_blacklist_filtered_fragments": mapped_fragments,
            "total_peak_overlaps": total_peak_overlaps,
            "feature_count": feature_count,
            "elapsed_seconds": round(time.time() - started, 3),
            "output": str(output),
            "output_sha256": checksum(output, "sha256"),
            "alignment_log": str(log_path),
        }
        if mapped_fragments < config["acceptance"]["minimum_mapped_fragments"]:
            raise RuntimeError(f"{gsm} has only {mapped_fragments} retained fragments")
        qc_output.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")
        return qc
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--temp-root", type=Path, default=Path("/tmp"))
    parser.add_argument("--only-gsm")
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    config = load_config(args.config.resolve())
    results = Path(config["storage"]["results_root"])
    manifest = results / "01_tables" / "01_immgen_library_run_manifest.tsv"
    peak_bed = results / "01_tables" / "02_bdev_exact_peak_axis.bed"
    if not manifest.exists() or not peak_bed.exists():
        raise RuntimeError("Run input preparation and peak-axis export first")

    for filename, expected in config["inputs"]["archived_mm10"]["bowtie2_index_sha256"].items():
        path = Path(config["inputs"]["archived_mm10"]["bowtie2_prefix"]).parent / filename
        if checksum(path, "sha256") != expected:
            raise RuntimeError(f"Bowtie2 index checksum changed: {path}")
    blacklist = Path(config["inputs"]["archived_mm10"]["blacklist"])
    if checksum(blacklist, "sha256") != config["inputs"]["archived_mm10"]["blacklist_sha256"]:
        raise RuntimeError("Blacklist checksum changed")

    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = config["preprocessing"]["bowtie2_library_path"] + (
        ":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else ""
    )
    versions = {
        "bowtie2": tool_version([config["preprocessing"]["bowtie2"], "--version"], env=env),
        "samtools": tool_version([config["preprocessing"]["samtools"], "--version"]),
        "bedtools": tool_version([config["preprocessing"]["bedtools"], "--version"]),
    }
    rows = read_tsv(manifest)
    if args.only_gsm:
        rows = [row for row in rows if row["gsm"] == args.only_gsm]
        if not rows:
            raise RuntimeError(f"Unknown GSM: {args.only_gsm}")
    if args.workers < 1:
        raise RuntimeError("--workers must be positive")
    qc_rows = []
    errors = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {}
        for index, row in enumerate(rows, 1):
            print(f"[submit {index}/{len(rows)}] {row['gsm']} {row['population']}", flush=True)
            future = executor.submit(
                process_library, row, config, peak_bed, env, args.temp_root.resolve()
            )
            futures[future] = row
        for future in as_completed(futures):
            row = futures[future]
            try:
                qc = future.result()
                qc_rows.append(qc)
                print(f"[complete {len(qc_rows)}/{len(rows)}] {row['gsm']}", flush=True)
            except Exception as error:
                errors.append((row["gsm"], repr(error)))
                print(f"[failed] {row['gsm']}: {error}", flush=True)

    all_qc = []
    for row in read_tsv(manifest):
        path = Path(config["storage"]["peak_counts"]) / f"{row['gsm']}.qc.json"
        if path.exists():
            all_qc.append(json.loads(path.read_text()))
    qc_table = results / "01_tables" / "03_immgen_preprocessing_qc.tsv"
    if all_qc:
        with qc_table.open("w", newline="") as handle:
            fields = list(all_qc[0])
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(all_qc)
    complete = len(all_qc) == config["inputs"]["selected_libraries"]["expected_libraries"]
    receipt = {
        "schema_version": "bdev.immgen_atac_preprocessing.v1",
        "status": "complete" if complete else "partial",
        "run_identity": config["run_identity"],
        "libraries_complete": len(all_qc),
        "libraries_expected": config["inputs"]["selected_libraries"]["expected_libraries"],
        "tool_versions": versions,
        "parameters": config["preprocessing"],
        "peak_axis_sha256": checksum(peak_bed, "sha256"),
        "archive_modified": False,
    }
    (results / "06_manifests" / "04_preprocessing_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if errors:
        raise RuntimeError(f"Libraries failed after retries: {errors}")


if __name__ == "__main__":
    main()
