#!/usr/bin/env Rscript
source("src/integration_workflow.R")
source("src/sliding_window_signature_summary.R")
source("src/sliding_window_iqr_plot.R")
suppressPackageStartupMessages(library(ggplot2))
input_dir <- report_dir
report_dir <- file.path(input_dir,"sliding_window_review")
fig_dir <- file.path(report_dir,"figures")
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
input_json <- function(name) jsonlite::read_json(file.path(input_dir,name),simplifyVector=TRUE)
short_id <- function(x) sub("^parse_cytokine[.]","",x)
save_plot <- function(p,name,width=11,height=7) {
  ggsave(file.path(fig_dir,paste0(name,".png")),p,width=width,height=height,dpi=180,bg="white")
  ggsave(file.path(fig_dir,paste0(name,".pdf")),p,width=width,height=height,bg="white")
}
theme_set(theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title.position="plot"))
main <- function() {
  stopifnot(identical(cfg$cytokine_analysis$population,"pooled"))
  settings <- cfg$cytokine_analysis$sliding_window
  log_event("Pooled sliding-window median/IQR analysis; all bundles read-only")
  b <- wb_open_bundle(bundle_dir,read_only=TRUE,validate=FALSE); on.exit(wb_close_bundle(b),add=TRUE)
  ref <- open_parent("bdev_reference"); on.exit(wb_close_bundle(ref),add=TRUE)
  manifest <- wb_manifest(b); parents <- check_parents()
  cells <- wb_cells(b,"all"); md <- wb_meta(b,cells=cells); ref_md <- wb_meta(ref,cells=cells)
  stopifnot(identical(cells,wb_cells(ref,"all")),identical(as.numeric(md$pt),as.numeric(ref_md$pt)),
    identical(is.na(md$pt),is.na(ref_md$pt)))
  pt <- pt_family(b)
  stopifnot(identical(pt,readRDS(file.path(input_dir,"exact_pt_family.rds"))))
  scores <- wb_art_load(b,input_json("scoring.json")$score_artifact_id,expected_kind="cell_vector_family")
  stopifnot(identical(scores,readRDS(file.path(input_dir,"ucell_scores.rds")))); rm(ref_md)
  coverage <- read.delim(file.path(input_dir,"feature_coverage.tsv"))
  stopifnot(setequal(unique(scores$key),coverage$feature_set[coverage$eligible]))
  write_json(wb_op_interface("cell_vector_family_qc",bundle=b),"native_spread_interface.json")
  result <- sliding_window_signature_summary(scores,pt,unlist(settings$width_fractions),
    settings$step_fraction,settings$min_cells)
  result$parameters$primary_width_fraction <- settings$primary_width_fraction
  write_json(result$parameters,"parameters.json")
  write_tsv(result$summaries,"sliding_window_summary.tsv")
  write_tsv(result$windows,"window_support.tsv")
  saveRDS(result$membership,file.path(report_dir,"window_membership.rds"))
  saveRDS(result$native_summaries,file.path(report_dir,"native_window_summaries.rds"))
  write_tsv(data.frame(cell_id=pt$cell_id,pt=pt$value),"source_pseudotime.tsv")
  all <- result$summaries
  primary_fraction <- settings$primary_width_fraction
  d <- all[all$width_fraction==primary_fraction & all$eligible_window,]
  windows <- result$windows[result$windows$width_fraction==primary_fraction,]
  stopifnot(nrow(d)>0,setequal(unique(d$feature_set),unique(scores$key)))
  ix <- match(d$feature_set,coverage$feature_set)
  d$cytokine <- coverage$cytokine[ix]
  d$direction <- factor(coverage$direction[ix],levels=c("up","down"))
  low <- windows[!windows$eligible_window,]
  half_step <- result$parameters$step_pt_units/2
  low$left <- pmax(min(pt$value),low$pseudotime-half_step)
  low$right <- pmin(max(pt$value),low$pseudotime+half_step)
  gap_layer <- function() geom_rect(data=low,aes(xmin=left,xmax=right,ymin=-Inf,ymax=Inf),
    inherit.aes=FALSE,fill="#64748B",alpha=.13)
  palette <- c(up="#B2182B",down="#2166AC")
  primary_width <- diff(range(pt$value))*primary_fraction
  caption <- sprintf("Window: %.0f%% of pt range (%.3f units); step: %.0f%% (%.3f units). Line: median; band: 25th-75th percentiles.\nAt least %d cells per window; gray marks unsupported centers crossed by display connectors. Windows are clipped at the endpoints.",
    primary_fraction*100,primary_width,settings$step_fraction*100,result$parameters$step_pt_units,settings$min_cells)
  curve_plot <- function(dat,title) ggplot(dat,aes(pseudotime,score_median,group=feature_set,color=direction,fill=direction))+
    gap_layer()+geom_ribbon(aes(ymin=score_q25,ymax=score_q75),alpha=.20,colour=NA)+
    geom_line(linewidth=.85)+scale_color_manual(values=palette)+scale_fill_manual(values=palette)+
    labs(title=title,x="Existing Bdev pseudotime (pt)",y="UCell score",color="Parse direction",fill="Parse direction",caption=caption)+
    theme(legend.position="bottom")
  p <- curve_plot(d,"Pooled cytokine signatures: sliding-window median and IQR")+
    facet_wrap(~feature_set,scales="free_y",ncol=4,labeller=labeller(feature_set=short_id))+
    labs(subtitle="All 12,642 cells with defined source pseudotime pooled together")
  line <- ggplot_build(p)$data[[3]]
  stopifnot(length(unique(line$group))==22L,nrow(line)==nrow(d),all(is.finite(line$x)),all(is.finite(line$y)))
  save_plot(p,"sliding_window_iqr",14,13)
  p5 <- sliding_window_iqr_overview(all,result$windows,coverage,.05,result$parameters)
  save_plot(p5,"sliding_window_iqr_5pct",14,13)
  grDevices::pdf(file.path(fig_dir,"cytokine_sliding_window_iqr.pdf"),width=11,height=6.5)
  for(cytokine in unique(coverage$cytokine)) {
    title <- paste("Pooled sliding-window response signature:",cytokine)
    if(any(!coverage$eligible[coverage$cytokine==cytokine])) title <- paste(title,"(up excluded: 70% gene coverage)")
    p <- curve_plot(d[d$cytokine==cytokine,],title)+facet_wrap(~direction,scales="free_y",nrow=1)
    print(p)
    ggsave(file.path(fig_dir,paste0("cytokine_",cytokine,".png")),p,width=11,height=6.5,dpi=180,bg="white")
  }
  grDevices::dev.off()
  log_event("Sliding-window IQR panels rendered; preparing width comparison")
  sensitivity <- all[all$eligible_window,]
  sensitivity$width <- factor(paste0(100*sensitivity$width_fraction,"%"),levels=paste0(100*unlist(settings$width_fractions),"%"))
  p <- ggplot(sensitivity,aes(pseudotime,score_median,color=width,group=interaction(feature_set,width)))+
    geom_line(linewidth=.65)+facet_wrap(~feature_set,scales="free_y",ncol=4,labeller=labeller(feature_set=short_id))+
    labs(title="Sensitivity to sliding-window width",x="Existing Bdev pseudotime (pt)",y="Median UCell score",color="Window width",
      caption="All widths use the same pooled cells, 1% steps and at least 30 cells per window. Retained centers are connected; no extra smoother is fitted.")
  save_plot(p,"window_width_sensitivity",14,13)
  support <- result$windows
  support$width <- factor(paste0(support$width_fraction*100,"%"),levels=paste0(100*unlist(settings$width_fractions),"%"))
  p <- ggplot(support,aes(pseudotime,n_cells,color=width,group=width))+geom_line()+
    geom_hline(yintercept=settings$min_cells,linetype="dotted")+
    labs(title="Number of pooled cells in each overlapping window",x="Existing Bdev pseudotime (pt)",y="Cells per window",color="Window width",
      caption="Counts overlap and must not be summed to obtain the number of unique cells. Endpoint windows are truncated to the observed range.")
  save_plot(p,"window_cell_support",11,5)
  support_summary <- do.call(rbind,lapply(unlist(settings$width_fractions),function(w) {
    x <- support[support$width_fraction==w,]
    data.frame(width_percent=100*w,width_pt_units=w*diff(range(pt$value)),centers=nrow(x),
      retained_centers=sum(x$eligible_window),unsupported_centers=sum(!x$eligible_window),
      min_cells=min(x$n_cells),max_cells=max(x$n_cells))
  }))
  write_tsv(support_summary,"width_support_summary.tsv")
  stopifnot(identical(manifest,wb_manifest(b)))
  checks <- c(result$checks,list(pseudotime_exact=TRUE,pseudotime_missingness_exact=TRUE,
    score_export_matches_bundle=TRUE,parent_fingerprints_unchanged=TRUE,
    single_finite_line_per_signature=TRUE,bundle_manifest_unchanged=TRUE))
  write_json(list(status="passed",checks=checks,scope="sliding-window review checks; recipient validation remains pending"),"review_checks.json")
  source_files <- c("figure_scripts/bdev_integrated/03_sliding_window_cytokine_report.R","src/sliding_window_signature_summary.R",
    "src/sliding_window_iqr_plot.R",
    "src/integration_workflow.R","processing/integration/run_stage.sh","config/bdev_integrated.yml","config/bundles.yml")
  write_json(list(run_id=integration_run,revision="sliding_window_review",bundle_dir=bundle_dir,
    bundle_state_fingerprint=manifest$current_state$fingerprint,parents=parents,runtime=integration_runtime,
    score_artifact_id=input_json("scoring.json")$score_artifact_id,parameters=result$parameters,
    execution="native scWorkbench distribution operation applied separately to overlapping pooled cell subsets; no promotion",
    generated_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    sources=lapply(source_files,function(path) list(path=path,sha256=wb_source_sha256(path))),
    inputs=lapply(c("ucell_scores.rds","exact_pt_family.rds","scoring.json","feature_coverage.tsv"),
      function(name) list(path=file.path(input_dir,name),sha256=wb_source_sha256(file.path(input_dir,name))))),"review_provenance.json")
  file.copy("config/bdev_integrated.yml",file.path(report_dir,"review_configuration.yml"),overwrite=TRUE)
  capture.output(sessionInfo(),file=file.path(report_dir,"session_info.txt"))
  image_tag <- function(name) paste0('<img src="figures/',name,'.png" alt="',name,'">')
  html_table <- function(x) as.character(knitr::kable(x,format="html",digits=4,escape=TRUE,row.names=FALSE))
  html <- c('<!doctype html><html><head><meta charset="utf-8"><title>Pooled sliding-window cytokine signatures</title>',
    '<style>body{font:16px system-ui,sans-serif;max-width:1200px;margin:36px auto;padding:0 24px;color:#172033}img{max-width:100%;height:auto}p,li{line-height:1.55}table{border-collapse:collapse;font-size:13px}th,td{border:1px solid #ddd;padding:6px;text-align:left}a{color:#1261a0}</style></head><body>',
    '<h1>Pooled sliding-window cytokine signatures</h1>',
    sprintf('<p><strong>Median curves with interquartile-range bands.</strong> The primary window spans %.0f%% of the original pt range (%.3f pseudotime units), advances by %.0f%% (%.3f units), and is evaluated at %d centers. Shading covers the 25th–75th percentiles of cell scores within each window.</p>',
      100*primary_fraction,primary_width,100*settings$step_fraction,result$parameters$step_pt_units,result$parameters$centers_per_width),
    '<p>All samples and cell-cycle states are pooled. All 12,642 cells with defined original pt enter the window analysis; the 14,081 missing pt values remain missing. The existing 22 UCell signatures across 14 cytokines are reused. IFN-gamma up remains excluded because only 70% of its source genes map to Bdev.</p>',
    sprintf('<p>Windows require at least %d cells; %d of %d primary centers are retained. Curves connect all retained centers. Gray marks unsupported centers; any segment across them is display interpolation. At the trajectory endpoints, windows are clipped to observed pt limits. Every defined-pt cell belongs to at least one window for each tested width.</p>',settings$min_cells,sum(windows$eligible_window),nrow(windows)),
    '<p>The IQR shows cell-to-cell spread, not confidence in the median. Neighboring windows share cells and are not independent. No additional smoothing, pseudotime refitting or score imputation is performed. A zero IQR can coexist with variation in the distribution tails; within-window SD is included in the table.</p>',
    '<h2>Sliding-window median and IQR</h2>',image_tag("sliding_window_iqr"),
    '<h2 id="five-percent-window">All gene sets: 5% sliding window with IQR</h2>',image_tag("sliding_window_iqr_5pct"),
    '<p>The 5% window uses the same pooled cells, 1% steps and 30-cell threshold. All 22 scored gene sets are shown; IFN-gamma up remains excluded by gene coverage. Shading shows cell-score IQR, not a confidence interval. <a href="figures/sliding_window_iqr_5pct.pdf">Download PDF</a> | <a href="figures/sliding_window_iqr_5pct.png">Download PNG</a>.</p>',
    '<h2>Per-cytokine curves</h2>',vapply(unique(coverage$cytokine),function(x) image_tag(paste0("cytokine_",x)),character(1)),
    '<h2>Window-width sensitivity and cell support</h2>',image_tag("window_width_sensitivity"),html_table(support_summary),image_tag("window_cell_support"),
    '<p>Widths of 5%, 10% and 15% expose the tradeoff between local detail and pooling over broader portions of the trajectory. Each panel uses its own y scale. Up and down gene sets remain separate; raw scores are not calibrated for comparing absolute cytokine activity.</p>',
    '<h2>Methods and provenance</h2>',
    '<p>Window membership uses inclusive pt bounds on the exact source metadata values. The exported scWorkbench operation <code>wb_op_cell_vector_family_qc</code> calculates medians, type-7 quartiles, means and sample SD separately in each window. Tables and saved memberships retain exact cell identifiers. The procedure is deterministic, with no random sampling or bootstrap. Original RNA-count UCell scores use max_rank=1500 and no graph smoothing.</p>',
    '<p>All bundles are mounted read-only. No new artifacts are promoted; acceptance and full recipient validation remain pending. The earlier <a href="../pooled_review/REPORT.html">fixed-bin pooled report</a> is preserved for comparison.</p>',
    '<h2>Downloads</h2><ul>',
    '<li><a href="figures/cytokine_sliding_window_iqr.pdf">All 14 cytokine panels (PDF)</a></li>',
    '<li><a href="figures/sliding_window_iqr.pdf">Sliding-window overview (PDF)</a></li>',
    '<li><a href="figures/sliding_window_iqr_5pct.pdf">All gene sets, 5% window with IQR (PDF)</a></li>',
    '<li><a href="figures/window_width_sensitivity.pdf">Window-width sensitivity (PDF)</a></li>',
    '<li><a href="sliding_window_summary.tsv">All window medians, quartiles, SD and cell counts</a></li>',
    '<li><a href="window_support.tsv">Window bounds and support</a></li><li><a href="window_membership.rds">Exact window membership (RDS)</a></li>',
    '<li><a href="source_pseudotime.tsv">Original cell pseudotime values</a></li><li><a href="native_window_summaries.rds">Native scWorkbench summaries (RDS)</a></li>',
    '<li><a href="review_checks.json">Validation checks</a></li><li><a href="review_provenance.json">Input and source checksums</a></li>',
    '<li><a href="parameters.json">Parameters</a></li><li><a href="review_configuration.yml">Configuration</a></li>',
    '<li><a href="session_info.txt">Runtime information</a></li></ul>',
    paste0('<p>Source run: ',integration_run,'. Source bundle: ',bundle_dir,'.</p>'),'</body></html>')
  writeLines(html,file.path(report_dir,"REPORT.html"))
  writeLines(c("# Pooled sliding-window cytokine report","","[Open report](REPORT.html)","",
    sprintf("Median and IQR: %.0f%% window, %.0f%% steps; 5%%/15%% width sensitivity included.",100*primary_fraction,100*settings$step_fraction),
    "All source samples and phases pooled; exact source pt preserved. No new promotion."),file.path(report_dir,"REPORT.md"))
  write_json(list(status="complete",population="pooled",method="sliding_window",review_status="awaiting_user_review",promoted=FALSE,
    signatures=22L,trajectory_cells=length(pt$value),primary_retained_centers=sum(windows$eligible_window),
    html=file.path(report_dir,"REPORT.html"),figures=list.files(fig_dir),
    report_sha256=wb_source_sha256(file.path(report_dir,"REPORT.html"))),"figures.json")
  log_event("Sliding-window report complete:",file.path(report_dir,"REPORT.html"))
}
main()
