#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: 01_prepare_atac_counts.sh SAMPLE_MANIFEST RAW_ROOT WORK_ROOT PREPARED_ROOT RUNTIME THREADS" >&2
  exit 2
fi

manifest="$(realpath "$1")"
raw_root="$2"
work_root="$3"
prepared_root="$4"
runtime="$(realpath "$5")"
threads="$6"
export PATH="${runtime}/bin:${PATH}"

if command -v macs3 >/dev/null 2>&1; then
  peak_caller="macs3"
elif command -v macs2 >/dev/null 2>&1; then
  peak_caller="macs2"
else
  echo "Neither macs3 nor macs2 is available in ${runtime}" >&2
  exit 1
fi

reference_root="${raw_root}/reference/mm10_ucsc"
sra_root="${raw_root}/sra"
fastq_root="${work_root}/fastq"
trim_root="${work_root}/trimmed_fastq"
bam_root="${prepared_root}/bam"
peak_root="${prepared_root}/sample_peaks"
qc_root="${prepared_root}/qc"
count_root="${prepared_root}/counts"
mkdir -p "${reference_root}" "${sra_root}" "${fastq_root}" "${trim_root}" \
  "${bam_root}" "${peak_root}" "${qc_root}" "${count_root}"

reference_gz="${reference_root}/mm10.fa.gz"
reference_fa="${reference_root}/mm10.fa"
blacklist_gz="${reference_root}/mm10-blacklist.v2.bed.gz"
blacklist_bed="${reference_root}/mm10-blacklist.v2.bed"
index_prefix="${reference_root}/bowtie2/mm10"
mkdir -p "$(dirname "${index_prefix}")"

if [[ ! -s "${reference_gz}" ]]; then
  curl -L --fail --retry 5 --retry-delay 5 \
    https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz \
    -o "${reference_gz}.partial"
  mv "${reference_gz}.partial" "${reference_gz}"
fi
if [[ ! -s "${reference_fa}" ]]; then
  gzip -dc "${reference_gz}" > "${reference_fa}.partial"
  mv "${reference_fa}.partial" "${reference_fa}"
fi
if [[ ! -s "${blacklist_gz}" ]]; then
  curl -L --fail --retry 5 --retry-delay 5 \
    https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz \
    -o "${blacklist_gz}.partial"
  mv "${blacklist_gz}.partial" "${blacklist_gz}"
fi
if [[ ! -s "${blacklist_bed}" ]]; then
  gzip -dc "${blacklist_gz}" > "${blacklist_bed}.partial"
  mv "${blacklist_bed}.partial" "${blacklist_bed}"
fi
if [[ ! -s "${index_prefix}.1.bt2" && ! -s "${index_prefix}.1.bt2l" ]]; then
  bowtie2-build --threads "${threads}" "${reference_fa}" "${index_prefix}"
fi

primary_contigs=()
for chromosome in $(seq 1 19) X Y; do
  primary_contigs+=("chr${chromosome}")
done

tail -n +2 "${manifest}" | while IFS=$'\t' read -r sample_id condition treatment_time_h geo_accession sra_run sra_experiment spots bases sra_size_mb source_label; do
  if [[ -n "${ONLY_SAMPLE:-}" && "${sample_id}" != "${ONLY_SAMPLE}" ]]; then
    continue
  fi
  final_bam="${bam_root}/${sample_id}.filtered.bam"
  narrow_peak="${peak_root}/${sample_id}_peaks.narrowPeak"
  if [[ -s "${final_bam}" && -s "${final_bam}.bai" && -s "${narrow_peak}" ]]; then
    echo "Skipping completed sample ${sample_id}"
    continue
  fi

  run_dir="${sra_root}/${sra_run}"
  sra_file="${run_dir}/${sra_run}.sra"
  sample_fastq_dir="${fastq_root}/${sample_id}"
  fasterq_tmp_dir="${work_root}/fasterq_tmp_${sample_id}"
  read1="${sample_fastq_dir}/${sra_run}_1.fastq"
  read2="${sample_fastq_dir}/${sra_run}_2.fastq"
  trim1="${trim_root}/${sample_id}_R1.fastq.gz"
  trim2="${trim_root}/${sample_id}_R2.fastq.gz"
  name_bam="${work_root}/${sample_id}.name.bam"
  fixmate_bam="${work_root}/${sample_id}.fixmate.bam"
  position_bam="${work_root}/${sample_id}.position.bam"
  dedup_bam="${work_root}/${sample_id}.dedup.bam"
  primary_bam="${work_root}/${sample_id}.primary.bam"

  echo "Preparing ${sample_id} (${sra_run})"
  if [[ ! -s "${dedup_bam}" ]]; then
    if [[ ! -s "${sra_file}" ]]; then
      prefetch --max-size 10G --output-directory "${sra_root}" "${sra_run}"
    fi
    vdb-validate "${sra_file}" > "${qc_root}/${sample_id}.vdb_validate.txt" 2>&1

    mkdir -p "${sample_fastq_dir}" "${fasterq_tmp_dir}"
    if [[ ! -s "${read1}" || ! -s "${read2}" ]]; then
      fasterq-dump "${sra_file}" --split-files --threads "${threads}" \
        --temp "${fasterq_tmp_dir}" --outdir "${sample_fastq_dir}"
    fi
    fastp --in1 "${read1}" --in2 "${read2}" --out1 "${trim1}" --out2 "${trim2}" \
      --detect_adapter_for_pe --thread "${threads}" \
      --json "${qc_root}/${sample_id}.fastp.json" \
      --html "${qc_root}/${sample_id}.fastp.html"
    rm -f -- "${read1}" "${read2}"

    bowtie2 --very-sensitive --no-mixed --no-discordant --dovetail -X 2000 \
      -p "${threads}" -x "${index_prefix}" -1 "${trim1}" -2 "${trim2}" \
      2> "${qc_root}/${sample_id}.bowtie2.txt" \
      | samtools view -@ "${threads}" -b -f 2 -q 30 - \
      | samtools sort -@ "${threads}" -n -o "${name_bam}" -
    samtools fixmate -@ "${threads}" -m "${name_bam}" "${fixmate_bam}"
    samtools sort -@ "${threads}" -o "${position_bam}" "${fixmate_bam}"
    samtools markdup -@ "${threads}" -r "${position_bam}" "${dedup_bam}"
  else
    echo "Resuming ${sample_id} from completed deduplicated BAM"
  fi
  samtools index -@ "${threads}" "${dedup_bam}"
  samtools view -@ "${threads}" -b "${dedup_bam}" "${primary_contigs[@]}" > "${primary_bam}"
  bedtools intersect -v -abam "${primary_bam}" -b "${blacklist_bed}" > "${final_bam}.partial"
  mv "${final_bam}.partial" "${final_bam}"
  samtools index -@ "${threads}" "${final_bam}"
  samtools flagstat -@ "${threads}" "${final_bam}" > "${qc_root}/${sample_id}.flagstat.txt"

  "${peak_caller}" callpeak -t "${final_bam}" -f BAMPE -g mm --keep-dup all -q 0.01 \
    -n "${sample_id}" --outdir "${peak_root}" \
    > "${qc_root}/${sample_id}.macs2.txt" 2>&1

  rm -f -- "${trim1}" "${trim2}" "${name_bam}" "${fixmate_bam}" \
    "${position_bam}" "${dedup_bam}" "${dedup_bam}.bai" "${primary_bam}"
  rmdir "${sample_fastq_dir}" 2>/dev/null || true
  rmdir "${fasterq_tmp_dir}" 2>/dev/null || true
done

if [[ -n "${ONLY_SAMPLE:-}" ]]; then
  echo "Completed requested sample ${ONLY_SAMPLE}"
  exit 0
fi

consensus_bed="${count_root}/gse259355_tcf3_prob_consensus_peaks.mm10.bed"
if [[ ! -s "${consensus_bed}" ]]; then
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "${peak_root}"/*_peaks.narrowPeak \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "${consensus_bed}.partial"
  mv "${consensus_bed}.partial" "${consensus_bed}"
fi

consensus_saf="${count_root}/gse259355_tcf3_prob_consensus_peaks.mm10.saf"
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"}
  {print $1"-"$2"-"$3,$1,$2+1,$3,"."}' "${consensus_bed}" > "${consensus_saf}"

bams=()
while IFS=$'\t' read -r sample_id condition treatment_time_h geo_accession sra_run sra_experiment spots bases sra_size_mb source_label; do
  [[ "${sample_id}" == "sample_id" ]] && continue
  bams+=("${bam_root}/${sample_id}.filtered.bam")
done < "${manifest}"

featurecounts_out="${count_root}/gse259355_tcf3_prob_featurecounts.tsv"
featureCounts -F SAF -a "${consensus_saf}" -o "${featurecounts_out}.partial" \
  -T "${threads}" -p --countReadPairs -B -C -s 0 "${bams[@]}"
mv "${featurecounts_out}.partial" "${featurecounts_out}"
if [[ -s "${featurecounts_out}.partial.summary" ]]; then
  mv "${featurecounts_out}.partial.summary" "${featurecounts_out}.summary"
fi

cp "${manifest}" "${prepared_root}/sample_manifest.tsv"
sha256sum "${manifest}" "${reference_gz}" "${blacklist_gz}" "${consensus_bed}" \
  "${featurecounts_out}" > "${prepared_root}/input_checksums.sha256"
{
  bowtie2 --version | head -n 1
  samtools --version | head -n 1
  bedtools --version
  "${peak_caller}" --version
  featureCounts -v 2>&1 | head -n 1
  fastp --version 2>&1 | head -n 1
  fasterq-dump --version | head -n 2
} > "${prepared_root}/software_versions.txt"

echo "Prepared counts: ${featurecounts_out}"
