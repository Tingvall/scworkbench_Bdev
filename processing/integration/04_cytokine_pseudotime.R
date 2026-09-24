#!/usr/bin/env Rscript
source("src/integration_workflow.R")
main <- function() {
  if(identical(cfg$cytokine_analysis$population,"pooled"))
    stop("This is the historical stratified promotion workflow. Use the figures stage for the current pooled-only review; do not rerun or promote historical sample/phase analyses.")
  scoring <- read_json("scoring.json"); stopifnot(scoring$status=="complete")
  b <- wb_open_bundle(bundle_dir,validate=FALSE); on.exit(wb_close_bundle(b),add=TRUE)
  pt <- pt_family(b)
  stopifnot(identical(pt,readRDS(file.path(report_dir,"exact_pt_family.rds"))))
  signal_ref <- wb_ref_artifact(scoring$score_artifact_id,"signature_scores")
  dynamics_outputs <- c("signal_table","dynamics_table","top_signals","fit_summary")
  batch_receipt_path <- file.path(report_dir,"dynamics_batch_promotion.rds")
  batch_receipt <- if(file.exists(batch_receipt_path)) readRDS(batch_receipt_path) else NULL
  # Reuse the first run's interrupted next checkpoint, if present.
  workspace <- file.path(session_root,"pooled_50")
  if(!dir.exists(workspace)) workspace <- file.path(session_root,"trajectory_batch")
  session <- NULL
  if(is.null(batch_receipt)) {
    log_event("Opening one durable workspace for remaining trajectory operations")
    session <- if(dir.exists(workspace)) wb_session_resume(b,workspace) else
      wb_session_create(b,workspace_root=session_root,session_id=basename(workspace))
    if(!is.null(session$durable$promotion_receipt))
      stop("Batch already promoted but project receipt missing; inspect native receipt before restoring, never repeat promotion.")
  }
  nodes <- list()
  queue <- function(label,op,inputs,outputs,notes) {
    if(!is.null(batch_receipt)) {
      node <- batch_receipt$nodes[[label]]
      stopifnot(!is.null(node))
      nodes[[label]] <<- node
      return(node)
    }
    old <- file.path(report_dir,paste0(label,"_promotion.rds"))
    if(file.exists(old)) {
      node <- list(outputs=outputs,promotion=readRDS(old),run_id=NULL)
      log_event("Reusing retained",label)
    } else {
      write_json(wb_op_interface(op,bundle=wb_session_view(session)),paste0(label,"_interface.json"))
      marker <- paste0("[",label,"] ")
      existing <- names(session$runs)[vapply(session$runs,function(r) {
        marked <- startsWith(as.character(r$notes %||% ""),marker)
        legacy <- identical(label,"pooled_50") && identical(r$op_name,"trajectory_signal_dynamics") &&
          identical(as.integer(r$params$n_bins),50L)
        isTRUE(marked) || isTRUE(legacy)
      },logical(1))]
      stopifnot(length(existing)<=1L)
      if(length(existing)) {
        run_id <- existing[[1]]
        log_event("Reusing completed checkpoint",label)
      } else {
        log_event("Running",label,op)
        rr <- wb_session_run_op(session,b,op,inputs=inputs,notes=paste0(marker,notes))
        session <<- rr$session
        run_id <- rr$run$run_id
      }
      stopifnot(all(outputs %in% names(session$runs[[run_id]]$outputs)))
      node <- list(outputs=outputs,promotion=NULL,run_id=run_id)
    }
    nodes[[label]] <<- node
    node
  }
  node_ref <- function(node,output) {
    if(!is.null(node$promotion)) wb_ref_artifact(artifact_id(node$promotion,output),output) else
      wb_ref_run(node$run_id,output)
  }
  run_dynamics <- function(label,cell_set,bins=30L,min_bin=30L,signal=signal_ref) {
    queue(label,"trajectory_signal_dynamics",list(signal=signal,pseudotime=pt,
      pseudotime_key="bdev_pt",pseudotime_direction="as_is",cell_set=cell_set,min_cells=10L,
      n_bins=as.integer(bins),summary_fun="median",min_cells_per_bin=as.integer(min_bin),
      smooth_method="none",bootstrap_n=0L,bootstrap_seed=1729L,top_n=50L,min_bins_observed=5L),
      dynamics_outputs,paste("Descriptive median dynamics on exact Bdev metadata pt;",cell_set,
        "equal-width bins in the selected population's original pseudotime range; no refit, smoothing or inferential confidence intervals."))
  }
  for(bins in c(30L,20L,50L)) run_dynamics(paste0("pooled_",bins),"set:bdev_pt_defined",bins)
  md <- wb_meta(b,cells=wb_cells(b,"all"))
  for(id in sort(unique(md$orig.ident))) run_dynamics(paste0("sample_",id),paste0("set:bdev_sample_",id),min_bin=20L)
  for(phase in sort(unique(md$Phase))) run_dynamics(paste0("phase_",phase),paste0("set:bdev_phase_",phase),min_bin=20L)
  numeric_fields <- c("nCount_RNA","nFeature_RNA","S.Score","G2M.Score")
  covariates <- wb_prim_cell_vector_family(rep(numeric_fields,each=nrow(md)),
    rep(md$cell_id,times=length(numeric_fields)),unlist(md[numeric_fields],use.names=FALSE),
    meta=list(source="exact Bdev reference metadata",fields=numeric_fields))
  run_dynamics("rna_cycle_covariates","set:bdev_pt_defined",signal=covariates)
  grouping <- queue("source_sample_vector","metadata_cell_vector",list(
    source_column="orig.ident",value_type="categorical",missing_policy="error",vector_name="bdev_source_sample",cell_set="all"),
    "vector","Source sample identity; biological independence has not been established.")
  queue("score_relationships","cell_vector_family_qc",list(signal=signal_ref,
    grouping=node_ref(grouping,"vector"),cell_set="set:bdev_pt_defined",correlation_methods=c("pearson","spearman")),
    c("key_summary","pairwise_summary","group_summary"),
    "Signature correlations and sample composition. Ignore winner-frequency fields: raw scores are not calibrated across signatures.")
  normal <- queue("rna_log_normalize","matrix_log_normalize",list(matrix="RNA_counts",cell_set="all",
    output_matrix_id="RNA_log1p_cp10k",scale_factor=10000,overwrite=FALSE),"summary",
    "RNA log1p CP10K matrix for secondary member-gene expression; primary scores remain raw-count UCell.")
  normal_matrix <- if(!is.null(normal$promotion)) "RNA_log1p_cp10k" else wb_ref_run(normal$run_id,"matrix")
  queue("normalized_signature_transfer","feature_set_namespace_transfer",list(
    source_matrix="RNA_counts",target_matrix=normal_matrix,cell_set="all",namespace="parse_cytokine_signatures",
    mode="mirror",overwrite=FALSE,strict=TRUE,created_by="jtg"),
    c("transfer_summary","transfer_manifest","transfer_rejected"),"Preserve memberships on the lossless RNA_log1p_cp10k gene universe.")
  coverage <- read.delim(file.path(report_dir,"feature_coverage.tsv"))
  eligible <- coverage$feature_set[coverage$eligible]
  expression_outputs <- c("signature_table","dynamics_table","top_sets","driver_table","driver_curve_table","fit_summary")
  queue("member_expression","trajectory_geneset_dynamics",list(matrix=normal_matrix,pseudotime=pt,
    pseudotime_key="bdev_pt",pseudotime_direction="as_is",cell_set="set:bdev_pt_defined",
    namespace="parse_cytokine_signatures",feature_sets=eligible,score_method="mean",n_bins=30L,
    summary_fun="mean",min_cells_per_bin=30L,min_set_size=5L,top_n=50L,top_driver_genes=5L,
    smooth_method="none",bootstrap_n=0L,n_permutations=0L,bootstrap_seed=1729L,permutation_seed=1729L),
    expression_outputs,"Secondary mean log-normalized member-gene expression to retain sparse RNA support. These are expression summaries, not UCell; no hypothesis tests or bootstrap confidence bands.")

  if(is.null(batch_receipt)) {
    selected <- unlist(lapply(nodes,function(node) {
      if(!is.null(node$promotion)) return(list())
      lapply(node$outputs,function(output) wb_ref_run(node$run_id,output))
    }),recursive=FALSE)
    log_event("Promoting completed trajectory comparisons together")
    promoted <- wb_session_promote_current(b,session,outputs=selected,name_prefix="cytokine_pt",
      description="Preserved-pt cytokine dynamics, sample/phase diagnostics and secondary gene expression",created_by="jtg",fail_fast=TRUE)
    for(label in names(nodes)) if(is.null(nodes[[label]]$promotion)) {
      rows <- promoted$promoted[promoted$promoted$run_id==nodes[[label]]$run_id,,drop=FALSE]
      nodes[[label]]$promotion <- list(promoted=rows)
    }
    batch_receipt <- list(nodes=nodes,promoted=promoted$promoted,workspace=workspace)
    saveRDS(batch_receipt,batch_receipt_path)
    write_json(promoted$promoted,"dynamics_batch_promotion.json")
    wb_close_bundle(wb_session_view(promoted$session))
  }
  receipts <- list()
  for(label in names(nodes)) {
    node <- nodes[[label]]
    ids <- setNames(vapply(node$outputs,function(output) artifact_id(node$promotion,output),character(1)),node$outputs)
    receipts[[label]] <- as.list(ids)
    for(output in node$outputs) {
      primitive <- wb_art_load(b,ids[[output]])
      if(!is.null(primitive$table)) write_tsv(primitive$table,paste0(label,"_",output,".tsv"))
    }
  }
  genes <- unique(unlist(lapply(eligible,function(id) wb_feature_set_get(b,id,"RNA_counts"))))
  counts <- as.matrix(wb_matrix_read_features(b,"RNA_counts",genes,pt$cell_id,missing="NA"))
  support <- data.frame(feature_id=genes,cells_detected=rowSums(counts>0),detection_fraction=rowMeans(counts>0),total_counts=rowSums(counts))
  write_tsv(support,"member_gene_expression_support.tsv")
  write_json(list(status="complete",pseudotime_exact=TRUE,defined_pt=length(pt$value),
    score_artifact_id=scoring$score_artifact_id,dynamics=receipts,
    expression=receipts$member_expression,inference="descriptive; source sample identity is not evidence of biological independence"),"dynamics.json")
  log_event("Pooled, sample, phase, QC and member-gene dynamics complete")
}
main()
