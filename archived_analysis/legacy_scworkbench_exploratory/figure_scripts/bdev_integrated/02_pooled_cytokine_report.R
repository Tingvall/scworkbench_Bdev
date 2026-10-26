#!/usr/bin/env Rscript
source("src/integration_workflow.R")
source("src/pooled_signature_summary.R")
suppressPackageStartupMessages(library(ggplot2))
input_dir <- report_dir
report_dir <- file.path(input_dir,"pooled_review")
fig_dir <- file.path(report_dir,"figures")
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
inputs_used <- character()
read_input <- function(name) {
  inputs_used <<- unique(c(inputs_used,name))
  read.delim(file.path(input_dir,name),check.names=FALSE)
}
input_json <- function(name) jsonlite::read_json(file.path(input_dir,name),simplifyVector=TRUE)
save_plot <- function(p,name,width=11,height=7) {
  ggsave(file.path(fig_dir,paste0(name,".png")),p,width=width,height=height,dpi=180,bg="white")
  ggsave(file.path(fig_dir,paste0(name,".pdf")),p,width=width,height=height,bg="white")
}
short_id <- function(x) sub("^parse_cytokine[.]","",x)
theme_set(theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title.position="plot"))
main <- function() {
  stopifnot(identical(cfg$cytokine_analysis$population,"pooled"),input_json("dynamics.json")$status=="complete")
  log_event("Pooled-only review: native median, quartiles and SD; all bundles read-only")
  b <- wb_open_bundle(bundle_dir,read_only=TRUE,validate=FALSE); on.exit(wb_close_bundle(b),add=TRUE)
  ref <- open_parent("bdev_reference"); on.exit(wb_close_bundle(ref),add=TRUE)
  manifest <- wb_manifest(b)
  parents <- check_parents()
  cells <- wb_cells(b,"all"); md <- wb_meta(b,cells=cells)
  ref_md <- wb_meta(ref,cells=cells)
  stopifnot(identical(cells,wb_cells(ref,"all")),
    identical(as.numeric(md$pt),as.numeric(ref_md$pt)),identical(is.na(md$pt),is.na(ref_md$pt)))
  pt <- pt_family(b)
  stopifnot(identical(pt,readRDS(file.path(input_dir,"exact_pt_family.rds"))))
  scores <- wb_art_load(b,input_json("scoring.json")$score_artifact_id,expected_kind="cell_vector_family")
  stopifnot(identical(scores,readRDS(file.path(input_dir,"ucell_scores.rds"))))
  coverage <- read_input("feature_coverage.tsv")
  old_curve <- read_input("pooled_30_signal_table.tsv")
  write_json(wb_op_interface("cell_vector_family_qc",bundle=b),"native_spread_interface.json")
  result <- pooled_signature_summary(scores,pt,old_curve,n_bins=30L,min_cells=30L)
  saveRDS(result[c("native","native_pt","grouping","parameters")],file.path(report_dir,"native_pooled_spread.rds"))
  write_tsv(result$spread,"pooled_bin_spread.tsv")
  write_tsv(result$support,"pooled_bin_support.tsv")
  write_tsv(data.frame(cell_id=pt$cell_id,pt=pt$value,bin=as.integer(sub("bin_","",result$grouping$value))),"pooled_cell_bins.tsv")
  write_json(result$parameters,"parameters.json")
  d <- result$primary
  match_coverage <- match(d$feature_set,coverage$feature_set)
  stopifnot(!anyNA(match_coverage),setequal(unique(d$feature_set),coverage$feature_set[coverage$eligible]))
  d$cytokine <- coverage$cytokine[match_coverage]
  d$direction <- factor(coverage$direction[match_coverage],levels=c("up","down"))
  low_support <- result$support[!result$support$eligible_bin,]
  palette <- c(up="#B2182B",down="#2166AC")
  gap_layer <- function() geom_rect(data=low_support,
    aes(xmin=pt_left,xmax=pt_right,ymin=-Inf,ymax=Inf),inherit.aes=FALSE,fill="#64748B",alpha=.10)
  caption <- "Line: pooled median; band: 25th–75th percentiles (IQR); dots: retained bin centers.\nGray intervals have fewer than 30 cells. Lines and band edges connect across these intervals for display."
  curve_plot <- function(dat,title) ggplot(dat,aes(pseudotime,score_median,group=feature_set,color=direction,fill=direction))+
    gap_layer()+geom_ribbon(aes(ymin=score_q25,ymax=score_q75),alpha=.19,colour=NA)+
    geom_line(linewidth=.8)+geom_point(size=1.2)+scale_color_manual(values=palette)+scale_fill_manual(values=palette)+
    labs(title=title,x="Existing Bdev pseudotime (pt)",y="UCell score",color="Parse direction",fill="Parse direction",caption=caption)+
    theme(legend.position="bottom")
  p <- curve_plot(d,"Pooled cytokine-response signatures: median and cell-to-cell spread")+
    facet_wrap(~feature_set,scales="free_y",ncol=4,labeller=labeller(feature_set=short_id))+
    labs(subtitle="All cells with defined source pt pooled together; 22 directional signatures")
  # A single finite line group per signature means retained bins are connected.
  built <- ggplot_build(p)$data[[3]]
  stopifnot(length(unique(built$group))==22L,all(is.finite(built$x)),all(is.finite(built$y)),nrow(built)==nrow(d))
  save_plot(p,"pooled_curves_iqr",14,13)
  grDevices::pdf(file.path(fig_dir,"cytokine_curves_pooled_iqr.pdf"),width=11,height=6.5)
  for(cytokine in unique(coverage$cytokine)) {
    z <- d[d$cytokine==cytokine,]
    title <- paste("Pooled Parse response signature:",cytokine)
    if(any(!coverage$eligible[coverage$cytokine==cytokine])) title <- paste(title,"(up excluded: 70% gene coverage)")
    pp <- curve_plot(z,title)+facet_wrap(~direction,scales="free_y",nrow=1)
    print(pp)
    ggsave(file.path(fig_dir,paste0("cytokine_",cytokine,".png")),pp,width=11,height=6.5,dpi=180,bg="white")
  }
  grDevices::dev.off()
  log_event("Pooled median/IQR panels rendered; generating spread and coverage figures")
  p <- ggplot(d,aes(pseudotime,score_sd,group=feature_set,color=direction))+gap_layer()+
    geom_line(linewidth=.75)+geom_point(size=.8)+scale_color_manual(values=palette)+
    facet_wrap(~feature_set,scales="free_y",ncol=4,labeller=labeller(feature_set=short_id))+
    labs(title="Within-bin standard deviation of pooled cell scores",x="Existing Bdev pseudotime (pt)",y="Cell-score SD",color="Parse direction",
      caption="SD measures cell-to-cell spread, including tails beyond the IQR. Gray marks low cell support; connecting segments are for display.")
  save_plot(p,"pooled_spread_sd",14,13)
  p <- ggplot(result$support,aes((pt_left+pt_right)/2,n_cells))+gap_layer()+geom_line(color="#334155")+
    geom_point(aes(shape=eligible_bin),size=2)+geom_hline(yintercept=30,linetype="dotted")+
    scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1),labels=c(`TRUE`="Retained",`FALSE`="Below 30 cells"))+
    labs(title="Cell support across the full pooled trajectory",x="Existing Bdev pseudotime (pt)",y="Cells per bin",shape=NULL,
      caption="All 12,642 defined-pt cells are accounted for. Primary curve anchors retain bins with at least 30 cells.")
  save_plot(p,"pooled_bin_support",11,5)
  sensitivity <- do.call(rbind,lapply(c(20L,30L,50L),function(n) {
    x <- read_input(paste0("pooled_",n,"_signal_table.tsv")); x$bins <- factor(n); x
  }))
  p <- ggplot(sensitivity,aes(pseudotime,value,color=bins,group=interaction(feature_set,bins)))+
    geom_line(linewidth=.65)+facet_wrap(~feature_set,scales="free_y",ncol=4,labeller=labeller(feature_set=short_id))+
    labs(title="Pooled median sensitivity to binning",x="Existing Bdev pseudotime (pt)",y="Median UCell score",color="Bins",
      caption="Each curve connects its retained bin centers. All settings use the same pooled cohort and a minimum of 30 cells per bin.")
  save_plot(p,"pooled_binning_sensitivity",14,13)
  covariates <- read_input("rna_cycle_covariates_signal_table.tsv")
  p <- ggplot(covariates,aes(pseudotime,value,group=signal_key))+gap_layer()+geom_line(color="#334155")+
    facet_wrap(~signal_key,scales="free_y",ncol=2)+labs(title="Pooled RNA depth and cell-cycle covariates",
      x="Existing Bdev pseudotime (pt)",y="Pooled median metadata value",caption="Pooled diagnostics; no sample or cell-cycle stratification. Connecting segments across gray intervals are for display.")
  save_plot(p,"pooled_covariates",11,7)
  summary <- read_input("pooled_30_dynamics_table.tsv")
  summary <- merge(coverage,summary,by="feature_set",all.x=TRUE,sort=FALSE)
  for(n in c(20,50)) {
    other <- read_input(paste0("pooled_",n,"_dynamics_table.tsv"))
    summary[[paste0("delta_",n)]] <- other$late_minus_early[match(summary$feature_set,other$feature_set)]
  }
  summary$binning_sign_consistent <- with(summary,sign(late_minus_early)==sign(delta_20)&sign(late_minus_early)==sign(delta_50))
  write_tsv(summary,"pooled_coverage_and_robustness.tsv")
  stopifnot(identical(manifest,wb_manifest(b)))
  checks <- c(result$checks,list(pseudotime_exact=TRUE,pseudotime_missingness_exact=TRUE,
    score_export_matches_bundle=TRUE,parent_fingerprints_unchanged=TRUE,
    single_finite_line_per_signature=TRUE,bundle_manifest_unchanged=TRUE))
  write_json(list(status="passed",checks=checks,scope="pooled review checks; full recipient validation remains pending"),"review_checks.json")
  source_files <- c("figure_scripts/bdev_integrated/02_pooled_cytokine_report.R","src/pooled_signature_summary.R",
    "src/integration_workflow.R","processing/integration/run_stage.sh","config/bdev_integrated.yml","config/bundles.yml")
  write_json(list(run_id=integration_run,revision="pooled_review",bundle_dir=bundle_dir,
    bundle_state_fingerprint=manifest$current_state$fingerprint,parents=parents,runtime=integration_runtime,
    score_artifact_id=input_json("scoring.json")$score_artifact_id,parameters=result$parameters,
    execution="exported scWorkbench wb_op_cell_vector_family_qc on in-memory primitives; no promotion",
    generated_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    sources=lapply(source_files,function(path) list(path=path,sha256=wb_source_sha256(path))),
    inputs=lapply(unique(c(inputs_used,"ucell_scores.rds","exact_pt_family.rds","scoring.json","creation.json")),
      function(name) list(path=file.path(input_dir,name),sha256=wb_source_sha256(file.path(input_dir,name))))),"review_provenance.json")
  file.copy("config/bdev_integrated.yml",file.path(report_dir,"review_configuration.yml"),overwrite=TRUE)
  capture.output(sessionInfo(),file=file.path(report_dir,"session_info.txt"))
  image_tag <- function(name) paste0('<img src="figures/',name,'.png" alt="',name,'">')
  html_table <- function(x) as.character(knitr::kable(x,format="html",digits=4,escape=TRUE,row.names=FALSE))
  html <- c('<!doctype html><html><head><meta charset="utf-8"><title>Pooled Bdev cytokine trajectories</title>',
    '<style>body{font:16px system-ui,sans-serif;max-width:1200px;margin:36px auto;padding:0 24px;color:#172033}img{max-width:100%;height:auto}p,li{line-height:1.55}table{border-collapse:collapse;font-size:13px}th,td{border:1px solid #ddd;padding:6px;text-align:left}a{color:#1261a0}</style></head><body>',
    '<h1>Pooled cytokine-response signatures along Bdev pseudotime</h1>',
    '<p><strong>Pooled-only review revision.</strong> All source samples and cell-cycle states are combined. The line is the pooled median UCell score; shading shows the 25th–75th percentiles (interquartile range, IQR) of cell scores in each retained bin. SD is also plotted and exported. These are measures of cell-to-cell spread, not confidence intervals.</p>',
    '<p>The full cohort has 26,723 scored cells. Trajectory summaries pool all 12,642 cells with defined original <code>pt</code>; 14,081 missing entries are preserved and cannot be positioned on the trajectory. Twenty-two signatures across 14 cytokines pass coverage filtering. IFN-gamma up is excluded at 35/50 mapped genes (70%, below 80%).</p>',
    '<p>Every curve connects its retained bin centers with a continuous solid line. Gray intervals have fewer than 30 cells: connecting line segments and band edges across them are display interpolation, not measurements within those intervals. No new pseudotime values or cell scores are imputed. The primary 30-bin setting retains 26 bins containing 12,576 cells; 66 cells in four low-coverage bins remain in the exported summaries.</p>',
    '<h2>Pooled median and spread</h2>',image_tag("pooled_curves_iqr"),
    '<p>The IQR may collapse to zero where most scores are zero; that does not imply every cell has the same score. The SD figure below captures spread in the tails. Panels use separate y scales. Up and down sets remain separate; higher up-set scores and lower down-set scores would be a concordant response-like pattern.</p>',
    '<h2>Per-cytokine curves</h2>',vapply(unique(coverage$cytokine),function(x) image_tag(paste0("cytokine_",x)),character(1)),
    '<h2>Within-bin standard deviation</h2>',image_tag("pooled_spread_sd"),
    '<h2>Coverage and binning sensitivity</h2>',image_tag("pooled_bin_support"),image_tag("pooled_binning_sensitivity"),
    html_table(summary[,c("feature_set","coverage","eligible","late_minus_early","binning_sign_consistent")]),
    '<p>All 22 scored signatures retain their early/late change sign across 20, 30 and 50 bins, including unchanged zero summaries. This is a sensitivity check, not a significance test. In the pooled summaries, c-Kit up declines (-0.0258), IL-1β up increases (+0.0428), and FLT3L down has a small net change (+0.0039) with a binning-sensitive curve classification. These values are descriptive summaries of bin medians and do not establish cytokine exposure.</p>',
    '<h2>Pooled covariates</h2>',image_tag("pooled_covariates"),
    '<h2>Methods and review state</h2>',
    '<p>Existing raw-count UCell scores (max_rank=1500, no graph smoothing) are reused. Native scWorkbench <code>cell_vector_family_qc</code> computes pooled quartiles (type 7), medians, means and sample SD within the same equal-width bins used by the retained trajectory operation. Cell identifiers, bin counts and medians were checked against retained outputs. Full precision is retained in the new native summary exports; the earlier trajectory tables have four-decimal rounding. No bootstrap or sampling is used.</p>',
    '<p>The original numerical pseudotime direction and missingness are unchanged. No sample-specific or phase-specific analysis is run or included in this revision. Earlier reports remain as historical records. All bundles are read-only; the working integration remains unaccepted and full recipient validation is pending. No new artifacts are promoted.</p>',
    '<h2>Downloads and reproducibility</h2><ul>',
    '<li><a href="figures/cytokine_curves_pooled_iqr.pdf">All 14 cytokine panels (PDF)</a></li>',
    '<li><a href="figures/pooled_curves_iqr.pdf">Pooled median/IQR overview (PDF)</a></li>',
    '<li><a href="figures/pooled_spread_sd.pdf">Pooled SD overview (PDF)</a></li>',
    '<li><a href="pooled_bin_spread.tsv">All bin summaries: quartiles, median, mean, SD and cell counts</a></li>',
    '<li><a href="pooled_bin_support.tsv">Cell support in every bin</a></li><li><a href="pooled_cell_bins.tsv">Exact cell-to-bin assignment and original pt</a></li>',
    '<li><a href="pooled_coverage_and_robustness.tsv">Gene coverage and binning sensitivity</a></li>',
    '<li><a href="review_checks.json">Validation checks</a></li><li><a href="review_provenance.json">Inputs, runtime and source checksums</a></li>',
    '<li><a href="parameters.json">Parameters</a></li><li><a href="review_configuration.yml">Project configuration</a></li>',
    '<li><a href="session_info.txt">R session information</a></li><li><a href="native_pooled_spread.rds">Native scWorkbench results and grouping (RDS)</a></li></ul>',
    paste0('<p>Source run: ',integration_run,'. Report revision: pooled_review. Source bundle: ',bundle_dir,'.</p>'),'</body></html>')
  writeLines(html,file.path(report_dir,"REPORT.html"))
  writeLines(c("# Pooled cytokine pseudotime review","","[Open the report](REPORT.html)","",
    "Only pooled data are analyzed. Continuous median lines carry 25th–75th percentile (IQR) bands; SD and bin counts are exported.",
    "Gray low-coverage intervals are crossed by display connectors. Missing source pt values are preserved.",
    "No new promotion or bundle modification. Review and full recipient validation remain pending."),file.path(report_dir,"REPORT.md"))
  write_json(list(status="complete",population="pooled",review_status="awaiting_user_review",promoted=FALSE,
    n_trajectory_cells=length(pt$value),signatures=length(unique(d$feature_set)),
    html=file.path(report_dir,"REPORT.html"),figures=list.files(fig_dir),
    report_sha256=wb_source_sha256(file.path(report_dir,"REPORT.html"))),"figures.json")
  log_event("Pooled-only report complete:",file.path(report_dir,"REPORT.html"))
}
main()
