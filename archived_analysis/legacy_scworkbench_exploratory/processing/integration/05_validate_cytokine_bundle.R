#!/usr/bin/env Rscript
source("src/integration_workflow.R")
main <- function() {
  stopifnot(read_json("dynamics.json")$status=="complete")
  parents <- check_parents()
  bdev <- open_parent("bdev_reference"); on.exit(wb_close_bundle(bdev),add=TRUE)
  b <- wb_open_bundle(bundle_dir,validate=FALSE)
  closed <- FALSE; on.exit(if(!closed) try(wb_close_bundle(b),silent=TRUE),add=TRUE)
  cells <- wb_cells(bdev,"all")
  stopifnot(identical(wb_cells(b,"all"),cells))
  source_md <- wb_meta(bdev,cells=cells); md <- wb_meta(b,cells=cells)
  fields <- setdiff(names(source_md),"dataset_id")
  stopifnot(all(vapply(fields,function(nm) reference_values_equal(source_md[[nm]],md[[nm]]),logical(1))),
    identical(md$reference_source_cell_id,cells),all(md$reference_source_bundle=="bdev_reference"),
    identical(is.na(source_md$pt),is.na(md$pt)),identical(as.numeric(source_md$pt),as.numeric(md$pt)),
    identical(pt_family(b),readRDS(file.path(report_dir,"exact_pt_family.rds"))))
  for(id in ref_cfg$bundle$spaces) stopifnot(integration_space_values_equal(wb_space(b,id,cells=cells),wb_space(bdev,id,cells=cells)))
  scoring <- read_json("scoring.json")
  scores <- wb_art_load(b,scoring$score_artifact_id,expected_kind="cell_vector_family")
  stopifnot(identical(scores,readRDS(file.path(report_dir,"ucell_scores.rds"))),
    all(is.finite(scores$value)),length(scores$value)==22L*26723L)
  for(key in unique(scores$key)) stopifnot(identical(scores$cell_id[scores$key==key],cells))
  sets <- read.delim(file.path(report_dir,"accepted_source_gene_sets.tsv"))
  features <- wb_feature_ids(b,"RNA_counts")
  for(id in unique(sets$feature_set)) {
    expected <- intersect(sets$feature_id[sets$feature_set==id],features)
    stopifnot(identical(as.character(wb_feature_set_get(b,id,"RNA_counts")),as.character(expected)),
      identical(as.character(wb_feature_set_get(b,id,"RNA_log1p_cp10k")),as.character(expected)))
  }
  populations <- c("pooled_30","pooled_20","pooled_50",paste0("sample_",1:4),paste0("phase_",c("G1","G2M","S")))
  bin_checks <- do.call(rbind,lapply(populations,function(label) {
    curve <- read.delim(file.path(report_dir,paste0(label,"_signal_table.tsv")))
    fit <- read.delim(file.path(report_dir,paste0(label,"_fit_summary.tsv")))
    selected <- is.finite(md$pt)
    if(startsWith(label,"sample_")) selected <- selected & md$orig.ident==sub("sample_","",label)
    if(startsWith(label,"phase_")) selected <- selected & md$Phase==sub("phase_","",label)
    minimum <- if(startsWith(label,"pooled_")) 30L else 20L
    stopifnot(nrow(fit)==1L,fit$n_cells_used==sum(selected),
      setequal(unique(curve$signal_key),unique(scores$key)),all(is.finite(curve$value)),
      all(curve$n_cells>=minimum),all(curve$pseudotime>=min(md$pt[selected])),
      all(curve$pseudotime<=max(md$pt[selected])),!anyDuplicated(curve[c("signal_key","bin")]))
    retained <- curve[curve$signal_key==curve$signal_key[[1]],]
    stopifnot(sum(retained$n_cells)<=sum(selected))
    data.frame(population=label,pt_cells=sum(selected),bins_retained=nrow(retained),
      cells_in_retained_bins=sum(retained$n_cells),cells_in_omitted_bins=sum(selected)-sum(retained$n_cells),
      min_cells_per_bin=minimum)
  }))
  write_tsv(bin_checks,"trajectory_bin_coverage.tsv")
  scientific <- list(cells_exact=TRUE,source_metadata_exact=TRUE,pseudotime_doubles_exact=TRUE,
    pseudotime_missingness_exact=TRUE,all_coordinates_exact=TRUE,all_rna_values_exact_at_ingest=TRUE,
    signature_memberships_exact_intersections=TRUE,score_cells_exact=TRUE,score_values_finite=TRUE,
    score_export_exact=TRUE,trajectory_cohorts_exact=TRUE,bin_thresholds_verified=TRUE,
    parent_fingerprints_unchanged=TRUE,parents=parents)
  write_json(scientific,"scientific_validation.json")
  workflow_files <- c("src/integration_workflow.R","src/bundle_config.R","src/bundle_validation.R",
    "src/reference_validation.R","processing/integration/run_stage.sh",
    list.files("processing/integration",pattern="[.]R$",full.names=TRUE),
    "figure_scripts/bdev_integrated/01_cytokine_pseudotime_report.R")
  write_json(list(recorded_at="final scientific validation",runtime=integration_runtime,
    files=lapply(workflow_files,function(path) list(path=path,sha256=wb_source_sha256(path)))),
    "workflow_source_identity.json")
  file.copy("config/bdev_integrated.yml",file.path(report_dir,"analysis_configuration.yml"),overwrite=TRUE)
  file.copy("config/bundles.yml",file.path(report_dir,"registry_at_validation.yml"),overwrite=TRUE)
  readiness <- wb_finalize_recipient_bundle(b)
  write_json(readiness,"recipient_readiness.json")
  stopifnot(isTRUE(readiness$ready))
  manifest <- wb_manifest(b)
  write_json(manifest,"final_manifest.json")
  wb_close_bundle(b); closed <- TRUE
  log_event("Running full recipient validation")
  validation <- wb_bundle_validation_report(bundle_dir,level="recipient",progress=TRUE)
  write_json(validation,"validation.json")
  stopifnot(isTRUE(validation$ok))
  write_json(list(status="complete",bundle_id="bdev_integrated",bundle_dir=bundle_dir,
    active_phase="parse_cytokine_pseudotime",completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    bundle_state_fingerprint=manifest$current_state$fingerprint,
    recipient_validation_fingerprint=validation$post_validation_fingerprint,
    validation_status=validation$status,parents=parents,signatures=22L,cytokines=14L,
    cells=26723L,trajectory_cells=12642L,report_dir=report_dir),"completion.json")
  log_event("Recipient validation passed; integration ready for acceptance")
}
main()
