#!/usr/bin/env Rscript
source("src/integration_workflow.R")
suppressPackageStartupMessages(library(ggplot2))
report_inputs <- character()
read_table <- function(name) {
  report_inputs <<- unique(c(report_inputs,name))
  read.delim(file.path(report_dir,name),check.names=FALSE)
}
fig_dir <- file.path(report_dir,"figures")
dir.create(fig_dir,showWarnings=FALSE)
theme_set(theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(),plot.title.position="plot"))
save_plot <- function(p,name,width=11,height=7) {
  ggsave(file.path(fig_dir,paste0(name,".png")),p,width=width,height=height,dpi=180,bg="white")
  ggsave(file.path(fig_dir,paste0(name,".pdf")),p,width=width,height=height,bg="white")
}
# Never connect across bins excluded for insufficient cell coverage.
segment_lines <- function(d,keys) {
  grouping <- interaction(d[keys],drop=TRUE,lex.order=TRUE)
  d$segment <- NA_character_
  for(g in levels(grouping)) {
    ix <- which(grouping==g); ix <- ix[order(d$bin[ix])]
    d$segment[ix] <- paste(g,cumsum(c(TRUE,diff(d$bin[ix])>1)),sep=":")
  }
  d
}
main <- function() {
  stopifnot(read_json("dynamics.json")$status=="complete",read_json("scoring.json")$status=="complete")
  log_event("Rendering review report; all bundles are read-only")
  b <- wb_open_bundle(bundle_dir,read_only=TRUE,validate=FALSE); on.exit(wb_close_bundle(b),add=TRUE)
  bdev <- open_parent("bdev_reference"); on.exit(wb_close_bundle(bdev),add=TRUE)
  manifest_before <- wb_manifest(b)
  cells <- wb_cells(b,"all"); md <- wb_meta(b,cells=cells)
  source_cells <- wb_cells(bdev,"all"); source_md <- wb_meta(bdev,cells=source_cells)
  scores <- readRDS(file.path(report_dir,"ucell_scores.rds"))
  checks <- list()
  check <- function(name,expr) {
    checks[[name]] <<- tryCatch(list(ok=isTRUE(force(expr))),
      error=function(e) list(ok=FALSE,error=conditionMessage(e)))
  }
  check("parent_fingerprints",{check_parents(); TRUE})
  check("cell_ids_and_order",identical(cells,source_cells))
  check("display_universe_all_cells",identical(wb_cells(b,"set:bdev_reference"),cells))
  check("source_metadata",all(vapply(setdiff(names(source_md),"dataset_id"),
    function(nm) reference_values_equal(source_md[[nm]],md[[nm]]),logical(1))))
  check("pseudotime_exact",identical(as.numeric(md$pt),as.numeric(source_md$pt)))
  check("pseudotime_missingness",identical(is.na(md$pt),is.na(source_md$pt)))
  check("analysis_pseudotime_exact",identical(pt_family(b),readRDS(file.path(report_dir,"exact_pt_family.rds"))))
  for(id in ref_cfg$bundle$spaces) check(paste0("coordinates_",id),
    integration_space_values_equal(wb_space(b,id,cells=cells),wb_space(bdev,id,cells=cells)))
  check("scores_match_retained_artifact",identical(scores,
    wb_art_load(b,read_json("scoring.json")$score_artifact_id,expected_kind="cell_vector_family")))
  check("score_cells_and_values",length(scores$value)==22L*26723L && length(unique(scores$key))==22L && all(is.finite(scores$value)) &&
    all(scores$value>=0 & scores$value<=1) && all(vapply(unique(scores$key),
      function(key) identical(scores$cell_id[scores$key==key],cells),logical(1))))
  coverage <- read_table("feature_coverage.tsv")
  samples <- read_table("sample_coverage.tsv")
  sets <- read_table("accepted_source_gene_sets.tsv")
  features <- wb_feature_ids(b,"RNA_counts")
  check("mapped_signature_memberships",all(vapply(unique(sets$feature_set),function(id) {
    expected <- intersect(sets$feature_id[sets$feature_set==id],features)
    identical(as.character(wb_feature_set_get(b,id,"RNA_counts")),as.character(expected)) &&
      identical(as.character(wb_feature_set_get(b,id,"RNA_log1p_cp10k")),as.character(expected))
  },logical(1))))
  populations <- c("pooled_30","pooled_20","pooled_50",paste0("sample_",1:4),paste0("phase_",c("G1","G2M","S")))
  bin_coverage <- do.call(rbind,lapply(populations,function(label) {
    curve <- read_table(paste0(label,"_signal_table.tsv")); fit <- read_table(paste0(label,"_fit_summary.tsv"))
    selected <- is.finite(md$pt)
    if(startsWith(label,"sample_")) selected <- selected & md$orig.ident==sub("sample_","",label)
    if(startsWith(label,"phase_")) selected <- selected & md$Phase==sub("phase_","",label)
    minimum <- if(startsWith(label,"pooled_")) 30L else 20L
    retained <- curve[curve$signal_key==curve$signal_key[[1]],]
    check(paste0("trajectory_",label),nrow(fit)==1L && fit$n_cells_used==sum(selected) &&
      setequal(unique(curve$signal_key),unique(scores$key)) && all(is.finite(curve$value)) &&
      all(curve$n_cells>=minimum) && all(curve$pseudotime>=min(md$pt[selected])) &&
      all(curve$pseudotime<=max(md$pt[selected])) && !anyDuplicated(curve[c("signal_key","bin")]) &&
      sum(retained$n_cells)<=sum(selected))
    data.frame(population=label,pt_cells=sum(selected),bins_retained=nrow(retained),
      cells_in_retained_bins=sum(retained$n_cells),cells_in_omitted_bins=sum(selected)-sum(retained$n_cells),
      min_cells_per_bin=minimum)
  }))
  write_tsv(bin_coverage,"review_trajectory_bin_coverage.tsv")
  checks_ok <- all(vapply(checks,function(x) isTRUE(x$ok),logical(1)))
  write_json(list(status=if(checks_ok) "passed" else "issues_found",checks=checks,
    scope="read-only scientific report checks; does not replace full recipient validation",
    coordinate_comparison="All values, dimensions and IDs must match exactly; only bundle-specific scworkbench_payload_read evidence is excluded",
    all_rna_values_exact_at_ingest=read_json("creation.json")$all_rna_values_exact),"review_checks.json")
  log_event("Read-only report checks:",if(checks_ok) "passed" else "issues found; see review_checks.json")
  add_labels <- function(d) {
    ix <- match(d$feature_set,coverage$feature_set)
    d$cytokine <- coverage$cytokine[ix]; d$direction <- coverage$direction[ix]
    d
  }
  pooled <- add_labels(read_table("pooled_30_signal_table.tsv"))
  pooled$population <- "Pooled"
  sample_curves <- do.call(rbind,lapply(samples$sample_id,function(id) {
    d <- add_labels(read_table(paste0("sample_",id,"_signal_table.tsv")))
    d$population <- paste("Sample",id); d
  }))
  curves <- rbind(pooled,sample_curves)
  curves <- segment_lines(curves,c("feature_set","population"))
  colors <- c(Pooled="#111827","Sample 1"="#0072B2","Sample 2"="#D55E00","Sample 3"="#009E73","Sample 4"="#CC79A7")
  line_plot <- function(d,title) ggplot(d,aes(pseudotime,value,color=population,group=segment)) +
    geom_line(linewidth=0.65) + geom_point(size=0.7) +
    scale_color_manual(values=colors) + labs(title=title,x="Existing Bdev pseudotime (pt)",y="Median UCell score",color=NULL,
      caption="Raw directional scores; lower down-set expression can support a response-like pattern. Gaps indicate insufficient cells.\nSamples have different pseudotime ranges; sample endpoints and bin numbers are not directly comparable.") +
    theme(legend.position="bottom")
  grDevices::pdf(file.path(fig_dir,"cytokine_curves_all.pdf"),width=11,height=6.5)
  for(cytokine in unique(coverage$cytokine)) {
    d <- curves[curves$cytokine==cytokine,,drop=FALSE]
    missing <- coverage$direction[coverage$cytokine==cytokine & !coverage$eligible]
    title <- paste("Parse response signature:",cytokine)
    if(length(missing)) title <- paste0(title," (",paste(missing,collapse=",")," excluded: insufficient gene coverage)")
    p <- line_plot(d,title) + facet_wrap(~direction,scales="free_y",nrow=1)
    print(p)
    ggsave(file.path(fig_dir,paste0("cytokine_",cytokine,".png")),p,width=11,height=6.5,dpi=180,bg="white")
  }
  grDevices::dev.off()
  primary <- segment_lines(pooled,"feature_set")
  p <- ggplot(primary,aes(pseudotime,value,group=segment,color=direction)) +
    geom_line(linewidth=0.7) + facet_wrap(~feature_set,scales="free_y",ncol=4) +
    scale_color_manual(values=c(up="#B2182B",down="#2166AC")) +
    labs(title="Cytokine-response signatures along Bdev pseudotime",subtitle="22 eligible signatures; unsmoothed pooled bin medians",
      x="Existing Bdev pseudotime (pt)",y="Median UCell score",color="Parse direction",
      caption="Each panel has its own y scale. Raw values are not calibrated for comparison between gene sets.")
  save_plot(p,"pooled_curves",14,13)
  sensitivity <- do.call(rbind,lapply(c(20,30,50),function(n) {
    d <- read_table(paste0("pooled_",n,"_signal_table.tsv")); d$bins <- factor(n); d
  }))
  sensitivity <- segment_lines(sensitivity,c("feature_set","bins"))
  p <- ggplot(sensitivity,aes(pseudotime,value,color=bins,group=segment)) + geom_line(linewidth=.6) +
    facet_wrap(~feature_set,scales="free_y",ncol=4) +
    labs(title="Sensitivity to the number of pseudotime bins",x="Existing Bdev pseudotime (pt)",y="Median UCell score",color="Bins",
      caption="All settings retain the same scores and exact pt coordinate; bins with fewer than 30 cells are omitted.")
  save_plot(p,"binning_sensitivity",14,13)
  log_event("Pooled curves and binning-sensitivity figures rendered")

  # Within-signature scaling is a display transformation of engine curve tables.
  heat <- pooled
  heat$display_z <- ave(heat$value,heat$feature_set,FUN=function(x) if(sd(x)>0) (x-mean(x))/sd(x) else rep(0,length(x)))
  ordered <- coverage$feature_set[order(coverage$cytokine,match(coverage$direction,c("up","down")))]
  heat$row <- match(heat$feature_set,ordered)
  empty <- data.frame(row=seq_along(ordered),xmin=min(pooled$pseudotime_min),xmax=max(pooled$pseudotime_max))
  labels <- ordered
  labels[!coverage$eligible[match(ordered,coverage$feature_set)]] <- paste0(labels[!coverage$eligible[match(ordered,coverage$feature_set)]]," [excluded]")
  p <- ggplot() + geom_rect(data=empty,aes(xmin=xmin,xmax=xmax,ymin=row-.46,ymax=row+.46),fill="#D1D5DB") +
    geom_rect(data=heat,aes(xmin=pseudotime_min,xmax=pseudotime_max,ymin=row-.46,ymax=row+.46,fill=display_z)) +
    scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,name="Within-set\ncurve z score") +
    scale_y_reverse(breaks=seq_along(ordered),labels=labels) +
    labs(title="Location of cytokine-response programs along Bdev pseudotime",
      subtitle="Scaled separately within each signature; gray marks excluded or unsupported regions",
      x="Existing Bdev pseudotime (pt)",y=NULL,
      caption="Up and down directions remain separate. Colors compare positions within a signature, not cytokine activity between signatures.") +
    theme(panel.grid=element_blank())
  save_plot(p,"signature_heatmap",13,9)

  # Counts are taken from the native engine summaries, once per population/bin.
  support <- curves[!duplicated(curves[c("population","bin")]),]
  p <- ggplot(support,aes(pseudotime,n_cells,color=population)) + geom_point() +
    scale_color_manual(values=colors) + labs(title="Cell support for retained pseudotime bins",x="Existing Bdev pseudotime (pt)",y="Cells per bin",color=NULL,
      caption="Pooled bins require at least 30 cells; sample bins require at least 20. Missing bins are not interpolated.")
  save_plot(p,"bin_coverage",11,5)
  covariates <- segment_lines(read_table("rna_cycle_covariates_signal_table.tsv"),"signal_key")
  p <- ggplot(covariates,aes(pseudotime,value,group=segment)) + geom_line(color="#374151") +
    facet_wrap(~signal_key,scales="free_y",ncol=2) +
    labs(title="RNA depth and cell cycle along the same coordinate",x="Existing Bdev pseudotime (pt)",y="Median source metadata value")
  save_plot(p,"rna_cycle_covariates",11,7)
  phases <- do.call(rbind,lapply(c("G1","S","G2M"),function(phase) {
    d <- add_labels(read_table(paste0("phase_",phase,"_signal_table.tsv"))); d$phase <- phase; d
  }))
  phases <- segment_lines(phases,c("feature_set","phase"))
  p <- ggplot(phases,aes(pseudotime,value,color=phase,group=segment)) + geom_line() +
    facet_wrap(~feature_set,scales="free_y",ncol=4) +
    labs(title="Signature curves within source cell-cycle phases",x="Existing Bdev pseudotime (pt)",y="Median UCell score",color="Phase")
  save_plot(p,"phase_stratified_curves",14,13)

  correlations <- read_table("score_relationships_pairwise_summary.tsv")
  # Display both orientations of the native unordered pairs, without computing
  # new correlations. Self-correlations are omitted and visibly marked gray.
  pair_values <- correlations[,c("key_a","key_b","spearman_cor")]
  mirrored <- pair_values; mirrored$key_a <- pair_values$key_b; mirrored$key_b <- pair_values$key_a
  correlation_keys <- sort(unique(c(pair_values$key_a,pair_values$key_b)))
  correlation_display <- rbind(pair_values,mirrored,
    data.frame(key_a=correlation_keys,key_b=correlation_keys,spearman_cor=NA_real_))
  stopifnot(length(correlation_keys)==22L,nrow(correlation_display)==22L^2,
    !anyDuplicated(correlation_display[c("key_a","key_b")]))
  p <- ggplot(correlation_display,aes(key_a,key_b,fill=spearman_cor)) + geom_tile() +
    scale_x_discrete(limits=correlation_keys,labels=function(x) sub("^parse_cytokine[.]","",x)) +
    scale_y_discrete(limits=rev(correlation_keys),labels=function(x) sub("^parse_cytokine[.]","",x)) +
    scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",limits=c(-1,1),na.value="#D1D5DB",name="Spearman\ncorrelation") +
    coord_equal() + labs(title="Relationships between directional signatures",x=NULL,y=NULL,
      caption="Native scWorkbench coefficients across cells with defined pt; each pair is mirrored. Gray marks omitted self-correlations or unavailable values.\nCorrelation can reflect shared genes, developmental state or technical effects.") +
    theme(axis.text.x=element_text(angle=60,hjust=1),panel.grid=element_blank())
  save_plot(p,"signature_correlations",13,11)

  genes <- read_table("member_expression_driver_curve_table.tsv")
  genes <- segment_lines(genes,c("feature_set","feature_id"))
  grDevices::pdf(file.path(fig_dir,"member_gene_curves.pdf"),width=11,height=6)
  for(id in unique(genes$feature_set)) {
    d <- genes[genes$feature_set==id,]
    p <- ggplot(d,aes(pseudotime,value,color=feature_id,group=segment)) + geom_line() +
      labs(title=paste("Member-gene expression:",id),x="Existing Bdev pseudotime (pt)",y="Mean log1p(CP10K) expression",color="Gene",
        caption="Descriptive member genes selected by scWorkbench dynamics; these curves use expression, not UCell.")
    print(p)
  }
  grDevices::dev.off()
  pt_values <- setNames(md$pt,cells)
  # Select one exhaustive set explicitly: cell_set="all" overlays every stored
  # subset separately, which would draw cells repeatedly in this derived bundle.
  p <- wb_plot_space_2d(b,"umap.atac",cell_set="set:bdev_reference",color_by=pt_values,color_label="pt",overlay_alpha=0.7,overlay_size=0.3) +
    labs(title="Preserved Bdev pseudotime on the original ATAC UMAP",subtitle="14,081 cells without pt remain gray; embedding coordinates unchanged")
  save_plot(p,"reference_pseudotime",9,7)
  grDevices::pdf(file.path(fig_dir,"signature_umaps.pdf"),width=9,height=7)
  for(id in unique(scores$key)) {
    ix <- scores$key==id; values <- setNames(scores$value[ix],scores$cell_id[ix])
    p <- wb_plot_space_2d(b,"umap.atac",cell_set="set:bdev_reference",color_by=values,color_label="UCell",overlay_alpha=0.65,overlay_size=0.3) +
      labs(title=id,subtitle="All Bdev cells; original umap.atac coordinates; raw UCell score")
    print(p)
    if(id %in% c("parse_cytokine.flt3l.down","parse_cytokine.c_kit.up","parse_cytokine.il_6.up"))
      ggsave(file.path(fig_dir,paste0("umap_",id,".png")),p,width=9,height=7,dpi=180,bg="white")
  }
  grDevices::dev.off()
  log_event("Sample, phase, member-gene and UMAP figures rendered")

  summary <- read_table("pooled_30_dynamics_table.tsv")
  summary <- merge(coverage,summary,by="feature_set",all.x=TRUE,sort=FALSE)
  for(bins in c(20,50)) {
    other <- read_table(paste0("pooled_",bins,"_dynamics_table.tsv"))
    ix <- match(summary$feature_set,other$feature_set)
    summary[[paste0("late_minus_early_",bins)]] <- other$late_minus_early[ix]
    summary[[paste0("peak_pt_",bins)]] <- other$peak_pseudotime[ix]
  }
  summary$binning_direction_consistent <- with(summary,sign(late_minus_early)==sign(late_minus_early_20) & sign(late_minus_early)==sign(late_minus_early_50))
  write_tsv(summary,"coverage_and_binning_robustness.tsv")
  n_consistent <- sum(summary$binning_direction_consistent,na.rm=TRUE)
  check_table <- data.frame(check=names(checks),result=vapply(checks,function(x) if(isTRUE(x$ok)) "passed" else
    paste("issue",x$error %||% "comparison did not match"),character(1)))
  validation_path <- file.path(report_dir,"validation.json")
  recipient_ok <- file.exists(validation_path) && isTRUE(read_json("validation.json")$ok)
  validation_note <- if(recipient_ok) "Full recipient validation passed; bundle acceptance remains pending user review." else
    "Full recipient validation has not completed. This is a review draft, and the integration bundle remains working and unaccepted."
  report_sources <- c("figure_scripts/bdev_integrated/01_cytokine_pseudotime_report.R",
    "processing/integration/run_stage.sh","src/integration_workflow.R","config/bdev_integrated.yml","config/bundles.yml")
  write_json(list(run_id=integration_run,bundle_dir=bundle_dir,review_status="awaiting_user_review",
    bundle_state_fingerprint=manifest_before$current_state$fingerprint,
    parents=read_json("creation.json")$parents,runtime=integration_runtime,
    generated_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    sources=lapply(report_sources,function(path) list(path=path,sha256=wb_source_sha256(path))),
    inputs=lapply(c(report_inputs,"ucell_scores.rds","exact_pt_family.rds","creation.json","scoring.json","dynamics.json"),
      function(name) list(file=name,sha256=wb_source_sha256(file.path(report_dir,name))))),"review_provenance.json")
  file.copy("config/bdev_integrated.yml",file.path(report_dir,"review_configuration.yml"),overwrite=TRUE)
  capture.output(sessionInfo(),file=file.path(report_dir,"review_session_info.txt"))
  display <- summary[,c("feature_set","coverage","eligible","curve_class","amplitude","late_minus_early","peak_pseudotime","binning_direction_consistent")]
  html_table <- function(d) as.character(knitr::kable(d,format="html",digits=3,escape=TRUE,row.names=FALSE))
  image_tag <- function(name) paste0('<img src="figures/',name,'.png" alt="',name,'">')
  html <- c('<!doctype html><html><head><meta charset="utf-8"><title>Bdev cytokine pseudotime</title>',
    '<style>body{font:16px system-ui,sans-serif;max-width:1200px;margin:40px auto;padding:0 24px;color:#172033}img{max-width:100%;height:auto}table{border-collapse:collapse;font-size:13px}th,td{border:1px solid #ddd;padding:6px;text-align:left}p{line-height:1.5}a{color:#1261a0}</style></head><body>',
    '<h1>Parse cytokine-response programs along Bdev pseudotime</h1>',
    '<p><strong>Review draft — awaiting user review.</strong> Analysis artifacts had already been promoted into the working derived bundle before the request for this report. This report performs no further promotion or bundle finalization. Parent bundles and the derived bundle are mounted read-only during rendering.</p>',
    paste0('<p>',validation_note,'</p>'),
    '<p>22 of 23 accepted Parse signatures meet the prespecified Bdev feature-coverage threshold. They cover 14 cytokines. IFN-gamma up is excluded because only 35/50 genes (70%) map, below the required 80%. Missing genes are retained in the coverage table; no replacement genes were selected.</p>',
    '<p>UCell scores were computed on RNA counts for all 26,723 Bdev cells, with max_rank=1500 and no graph smoothing. Trajectory summaries use the 12,642 cells with existing pt. All source pt values, 14,081 missing entries, metadata, RNA counts and embedding coordinates were preserved. Summaries are descriptive and do not establish cytokine exposure or causality.</p>',
    '<p>Higher up-set expression together with lower down-set expression can support a concordant response-like pattern. Directional sets remain separate. Raw scores and heatmap colors must not be used to compare absolute activity between cytokines. Four source sample IDs are present; biological independence is unconfirmed, and their pseudotime ranges differ. No replicate-level confidence intervals or p-values are claimed.</p>',
    '<h2>Findings to review</h2><ul>',
    '<li><strong>c-Kit induced signature:</strong> pooled expression decreases toward higher pt (late-minus-early UCell summary -0.0258). The induced set also decreases in each cell-cycle stratum. Its repressed set is heterogeneous across samples and phases, so the paired signatures do not support a uniform response interpretation.</li>',
    '<li><strong>IL-1β induced signature:</strong> the pooled curve increases (+0.0428), with similar direction in all three cell-cycle strata. Sample 2 shows the clearest increase; the other sample curves are flat or much weaker. The repressed set also increases in the pooled early/late comparison. This is not a concordant activation call.</li>',
    '<li><strong>FLT3L repressed signature:</strong> the pooled net change is small (+0.0039), its native curve classification changes between 30 and 50 bins, and sample directions differ. It does not provide a robust directional FLT3L activity result in this analysis. An accepted FLT3L induced set is unavailable.</li>',
    '<li><strong>Context matters:</strong> RNA depth and cell-cycle scores vary substantially with pt; samples cover different parts of the trajectory. Some signatures are strongly correlated (for example CSF1-down and IL3-down, Spearman 0.7781), limiting cytokine-specific interpretation.</li></ul>',
    '<p>These values are native early/late summaries of retained bin medians, not matched biological time points or replicate effect estimates. Increasing pt follows the original numerical direction; no biological start/end labels were imposed. Source identifier <code>tgb_b1</code> is preserved as supplied.</p>',
    '<h2>Signature patterns</h2>',image_tag("signature_heatmap"),image_tag("pooled_curves"),
    '<h2>Sample coverage and potential confounding</h2>',html_table(samples),
    html_table(bin_coverage),image_tag("bin_coverage"),image_tag("rna_cycle_covariates"),
    '<p>Inspect sample and cell-cycle curves before interpreting a pooled pattern. Each sample/phase is binned across its own observed range; compare actual pt coordinates within shared support. Lines do not bridge unsupported bins.</p>',image_tag("phase_stratified_curves"),
    '<h2>Coverage and binning sensitivity</h2>',html_table(display),
    paste0('<p>',n_consistent,' of 22 scored signatures retain the same early/late change sign across 20, 30 and 50 bins, including unchanged zero summaries. This does not establish consistent curve shape.</p>'),image_tag("binning_sensitivity"),
    '<p>Curve classes and peak locations describe the retained bins. Small changes can receive non-flat classes under engine defaults; inspect amplitudes and raw curves. Direction agreement across 20, 30 and 50 bins is a sensitivity check, not a significance test.</p>',
    '<h2>Signature relationships</h2>',image_tag("signature_correlations"),
    '<p>Correlated signatures may describe shared programs rather than distinct cytokine effects. The complete accepted source memberships and native correlation table are retained for inspection.</p>',
    '<h2>Per-cytokine sample curves</h2>',vapply(unique(coverage$cytokine),function(x) image_tag(paste0("cytokine_",x)),character(1)),
    '<h2>Reference geometry</h2>',image_tag("reference_pseudotime"),
    '<h2>Methods and validation state</h2>',
    '<p>Accepted Parse control/single-cytokine gene sets were transferred by exact gene ID to Bdev RNA_counts. Eligibility requires at least 5 mapped genes and 80% source membership. UCell uses max_rank=1500 and scores all cells without graph smoothing. Primary trajectory curves are unsmoothed medians in 30 bins with at least 30 cells; sensitivity uses 20 and 50 bins. Four source samples and G1/S/G2M phases use 30 bins with at least 20 cells. No pseudotime imputation, reversal, rescaling or fitting was performed. Native tabular summaries round values to four decimals; source metadata pt and per-cell scores retain full precision.</p>',
    '<p>Member-gene figures show the five driver genes selected by the native dynamics operation for each eligible set, using mean log1p(CP10K) RNA expression. No bootstrap, permutation, confidence interval or significance claim is included (bootstrap_n=0, permutations=0; recorded seed 1729). Cell-cycle stratification and depth curves are diagnostics, not covariate adjustment. The KO multiome dataset is not included in this phase.</p>',
    '<p>Coordinate checks require identical values, dimensions and cell identifiers; bundle-specific storage-read evidence is excluded from the comparison. Score counts use the values in the scWorkbench vector-family object. These checks do not modify the stored objects.</p>',
    paste0('<p>',validation_note,' Source ingestion verified all RNA count values. The following checks were repeated without modifying either bundle:</p>'),html_table(check_table),
    '<h2>Downloads and provenance</h2><ul>',
    '<li><a href="figures/cytokine_curves_all.pdf">All cytokine curves (PDF)</a></li><li><a href="figures/member_gene_curves.pdf">Member-gene expression curves (PDF)</a></li>',
    '<li><a href="figures/signature_umaps.pdf">All signature UMAP overlays (PDF)</a></li><li><a href="coverage_and_binning_robustness.tsv">Coverage and robustness table</a></li>',
    '<li><a href="member_gene_expression_support.tsv">Actual RNA detection of member genes</a></li><li><a href="feature_coverage.tsv">Missing genes and coverage</a></li>',
    '<li><a href="accepted_source_gene_sets.tsv">Complete source memberships</a></li><li><a href="score_relationships_pairwise_summary.tsv">Signature correlation table</a></li>',
    '<li><a href="review_checks.json">Read-only report checks</a></li><li><a href="creation.json">Source preservation at ingestion</a></li><li><a href="review_trajectory_bin_coverage.tsv">Trajectory cell support</a></li></ul>',
    if(recipient_ok) '<p><a href="validation.json">Full recipient validation</a></p>' else character(),
    '<p><a href="review_provenance.json">Accepted parent identities, input checksums and runtime</a> | <a href="review_configuration.yml">Analysis configuration</a> | <a href="review_session_info.txt">R session information</a></p>',
    paste0('<p>Run: ',integration_run,'. scWorkbench ',as.character(packageVersion("scworkbench")),'. Output bundle: ',bundle_dir,'.</p>'),
    '</body></html>')
  writeLines(html,file.path(report_dir,"REPORT.html"))
  writeLines(c("# Bdev cytokine pseudotime integration — review draft", "", validation_note,
    "Analysis artifacts were already promoted before the report request. Report rendering makes no further promotions; all bundles are read-only.",
    "", "- 22 eligible signatures across 14 cytokines; IFN-gamma up excluded at 70% feature coverage.",
    "- 26,723 scored cells; 12,642 defined pt values; 14,081 missing values preserved.",
    "- Raw UCell scores, 20/30/50-bin pooled curves, four sample curves, three cell-cycle strata, and member-gene expression retained.",
    "- Descriptive response-like transcriptional programs; no cytokine-exposure, causal, or biological-replicate inference.",
    "", "[Open the complete report](REPORT.html)","",paste("Bundle fingerprint:",manifest_before$current_state$fingerprint)),file.path(report_dir,"REPORT.md"))
  stopifnot(identical(manifest_before,wb_manifest(b)))
  write_json(list(status="complete",review_status="awaiting_user_review",bundle_acceptance="pending",
    read_only_checks_passed=checks_ok,full_recipient_validation_passed=recipient_ok,
    html=file.path(report_dir,"REPORT.html"),figures=list.files(fig_dir),
    report_sha256=wb_source_sha256(file.path(report_dir,"REPORT.html"))),"figures.json")
  log_event("Report and figures complete")
}
main()
