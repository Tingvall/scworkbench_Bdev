#!/usr/bin/env Rscript
source("src/integration_workflow.R")
main <- function() {
  stopifnot(read_json("creation.json")$status=="complete")
  check_parents()
  b <- wb_open_bundle(bundle_dir,validate=FALSE); on.exit(wb_close_bundle(b),add=TRUE)
  accepted <- read.delim(file.path(report_dir,"accepted_source_gene_sets.tsv"),check.names=FALSE)
  coverage <- read.delim(file.path(report_dir,"feature_coverage.tsv"),check.names=FALSE)
  csv <- file.path(report_dir,"accepted_source_gene_sets.csv")
  write.csv(accepted,csv,row.names=FALSE)
  stopifnot(identical(read.csv(csv,check.names=FALSE),accepted))
  imported <- run_operation(b,"signature_import","feature_set_import",list(
    matrix="RNA_counts",cell_set="all",source_path=csv,source_format="table",
    namespace="parse_cytokine_signatures",collection="individual_cytokine",
    species="Mus musculus",organism="Mus musculus",gene_id_type="gene_symbol",
    set_id_col="feature_set",set_name_col="set_name",gene_col="feature_id",description_col="description",
    split_by_collection=FALSE,overwrite=FALSE,strict=FALSE,created_by="jtg",
    source_id="parse_cytokine_signature_workflow_v1",source_version="2026-09-21",origin="derived",is_readonly=TRUE),
    c("import_summary","import_manifest","rejected_sets"),
    "Transfer all accepted Parse memberships to Bdev RNA; record missing source genes explicitly; apply coverage eligibility separately.")
  write_tsv(load_output(b,imported,"import_manifest")$table,"signature_import_manifest.tsv")
  qc <- run_operation(b,"signature_overlap_qc","feature_set_overlap_qc",list(
    matrix="RNA_counts",cell_set="all",namespace="parse_cytokine_signatures",
    min_set_size=5L,min_detected_fraction=0.8,deduplicate_sets=FALSE),
    c("overlap_qc","namespace_summary"),"Prespecified minimum five genes and 80 percent source membership; detected means matrix feature presence.")
  qc_table <- load_output(b,qc,"overlap_qc")$table
  write_tsv(qc_table,"signature_overlap_qc.tsv")
  selected <- coverage$feature_set[coverage$eligible]
  stopifnot(setequal(selected,qc_table$feature_set[qc_table$screening_status=="include"]))
  signatures <- setNames(lapply(selected,function(id) accepted$feature_id[accepted$feature_set==id]),selected)
  score <- run_operation(b,"cytokine_ucell","ucell",list(matrix="RNA_counts",cell_set="all",min_cells=10L,
    features=signatures,name="parse_cytokine",max_rank=1500L,ncores=4L,parallel_backend="auto",
    chunk_size=400L,missing_features_policy="drop"),c("signature_scores","signature_summary"),
    "Bdev RNA scoring of 22 coverage-eligible accepted Parse directional signatures. UCell max_rank 1500; no graph smoothing or cross-signature calibration. Missing source genes dropped and reported.")
  scores <- load_output(b,score,"signature_scores")
  cells <- wb_cells(b,"all")
  stopifnot(length(scores$value)==length(cells)*length(selected), all(is.finite(scores$value)),
            all(scores$value>=0 & scores$value<=1),setequal(unique(scores$key),selected))
  for(id in selected) stopifnot(identical(scores$cell_id[scores$key==id],cells))
  saveRDS(scores,file.path(report_dir,"ucell_scores.rds"))
  write_tsv(load_output(b,score,"signature_summary")$table,"ucell_summary.tsv")
  write_json(list(status="complete",score_artifact_id=artifact_id(score,"signature_scores"),
    summary_artifact_id=artifact_id(score,"signature_summary"),signatures=length(selected),cells=length(cells),
    values=length(scores$value),all_cells_aligned=TRUE,all_values_finite=TRUE,
    excluded=coverage$feature_set[!coverage$eligible],exclusion_reason="less than 80 percent source gene coverage"),"scoring.json")
  log_event("Scoring complete",length(scores$value),"retained per-cell scores")
}
main()
