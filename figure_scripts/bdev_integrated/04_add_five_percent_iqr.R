#!/usr/bin/env Rscript
# Add a figure from existing native results without rerunning window analysis.
source("src/integration_workflow.R")
source("src/sliding_window_iqr_plot.R")
suppressPackageStartupMessages(library(ggplot2))
input_dir <- report_dir
report_dir <- file.path(input_dir,"sliding_window_review")
fig_dir <- file.path(report_dir,"figures")
theme_set(theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title.position="plot"))
main <- function() {
  log_event("Adding 5% IQR overview from saved native window summaries")
  receipt <- read_json("figures.json"); stopifnot(receipt$status=="complete",identical(receipt$promoted,FALSE))
  provenance <- read_json("review_provenance.json")
  for(i in seq_len(nrow(provenance$inputs))) stopifnot(identical(
    wb_source_sha256(provenance$inputs$path[[i]]),provenance$inputs$sha256[[i]]))
  summaries <- read.delim(file.path(report_dir,"sliding_window_summary.tsv"))
  support <- read.delim(file.path(report_dir,"window_support.tsv"))
  coverage <- read.delim(file.path(input_dir,"feature_coverage.tsv"))
  parameters <- read_json("parameters.json")
  native <- readRDS(file.path(report_dir,"native_window_summaries.rds"))
  # Check the plotted values against retained native summaries, allowing only
  # floating-point export rounding; no window statistics are recomputed.
  selected <- summaries[summaries$width_fraction==.05 & summaries$eligible_window,]
  columns <- c("score_q25","score_median","score_q75")
  for(id in unique(selected$window_id)) {
    d <- selected[selected$window_id==id,]
    original <- native[[id]]$table
    ix <- match(d$feature_set,original$key); stopifnot(!anyNA(ix))
    stopifnot(all(d$n_cells==original$n_finite[ix]),
      max(abs(as.matrix(d[columns])-as.matrix(original[ix,columns])))<1e-12)
  }
  p <- sliding_window_iqr_overview(summaries,support,coverage,.05,parameters)
  stem <- file.path(fig_dir,"sliding_window_iqr_5pct")
  ggsave(paste0(stem,".png"),p,width=14,height=13,dpi=180,bg="white")
  ggsave(paste0(stem,".pdf"),p,width=14,height=13,bg="white")
  html <- readLines(file.path(report_dir,"REPORT.html"),warn=FALSE)
  if(!any(grepl('id="five-percent-window"',html,fixed=TRUE))) {
    anchor <- which(grepl('<h2>Per-cytokine curves</h2>',html,fixed=TRUE)); stopifnot(length(anchor)==1L)
    html <- append(html,c('<h2 id="five-percent-window">All gene sets: 5% sliding window with IQR</h2>',
      '<img src="figures/sliding_window_iqr_5pct.png" alt="5% sliding-window median and IQR for all scored gene sets">',
      '<p>The 5% window uses the same pooled cells, 1% steps and 30-cell threshold. All 22 scored gene sets are shown; IFN-gamma up remains excluded by gene coverage. Shading shows cell-score IQR, not a confidence interval. <a href="figures/sliding_window_iqr_5pct.pdf">Download PDF</a> | <a href="figures/sliding_window_iqr_5pct.png">Download PNG</a>.</p>'),after=anchor-1L)
  }
  writeLines(html,file.path(report_dir,"REPORT.html"))
  write_json(list(status="complete",width_fraction=.05,step_fraction=parameters$step_fraction,
    signatures=length(unique(selected$feature_set)),retained_centers=length(unique(selected$window_id)),
    native_summary_values_verified=TRUE,analysis_recomputed=FALSE,promoted=FALSE,
    sources=lapply(c("figure_scripts/bdev_integrated/04_add_five_percent_iqr.R","src/sliding_window_iqr_plot.R"),
      function(path) list(path=path,sha256=wb_source_sha256(path))),
    inputs=lapply(c("sliding_window_summary.tsv","window_support.tsv","parameters.json","native_window_summaries.rds"),
      function(name) list(file=name,sha256=wb_source_sha256(file.path(report_dir,name))))),"five_percent_figure_provenance.json")
  receipt$figures <- list.files(fig_dir)
  receipt$five_percent_overview <- list(signatures=22L,retained_centers=length(unique(selected$window_id)),
    provenance="five_percent_figure_provenance.json")
  receipt$report_sha256 <- wb_source_sha256(file.path(report_dir,"REPORT.html"))
  write_json(receipt,"figures.json")
  log_event("5% window IQR overview added to report")
}
main()
