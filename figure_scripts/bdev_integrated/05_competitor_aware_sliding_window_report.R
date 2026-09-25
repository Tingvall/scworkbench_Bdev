#!/usr/bin/env Rscript

# Reproduce the previous pooled Bdev trajectory display for the newly accepted
# competitor-aware full signatures. Scientific summaries use the same native
# scWorkbench distribution operation and window definitions as the prior run.

source("src/competitor_aware_integration_workflow.R")
source("src/sliding_window_signature_summary.R")
suppressPackageStartupMessages(library(ggplot2))

input_dir <- report_dir
report_dir <- file.path(input_dir, "trajectory_review")
figure_dir <- file.path(report_dir, "figures")
if (file.exists(report_dir)) stop("Report path already exists: ", report_dir)
stage_dir <- paste0(report_dir, ".incomplete")
if (file.exists(stage_dir)) stop("Incomplete report path already exists: ", stage_dir)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
report_dir <- stage_dir
figure_dir <- file.path(report_dir, "figures")

save_plot <- function(plot, name, width = 12, height = 8) {
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot,
    width = width, height = height, dpi = 180, bg = "white")
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot,
    width = width, height = height, bg = "white")
}

slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
  plot.title.position = "plot"))

main <- function() {
  score_state <- jsonlite::read_json(file.path(input_dir, "scoring.json"),
    simplifyVector = TRUE)
  stopifnot(identical(score_state$status, "complete"),
    identical(score_state$promoted, FALSE))
  parents <- check_competitor_aware_parents()
  settings <- competitor_cfg$sliding_window

  bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
  on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  ref <- open_parent("bdev_reference")
  on.exit(try(wb_close_bundle(ref), silent = TRUE), add = TRUE)
  manifest <- wb_manifest(bundle)
  cells <- wb_cells(bundle, "all")
  metadata <- wb_meta(bundle, cells = cells)
  ref_metadata <- wb_meta(ref, cells = cells)
  stopifnot(identical(as.numeric(metadata$pt), as.numeric(ref_metadata$pt)),
    identical(is.na(metadata$pt), is.na(ref_metadata$pt)))
  pt <- pt_family(bundle)
  stopifnot(identical(pt, readRDS(file.path(input_dir, "exact_pt_family.rds"))))
  scores <- readRDS(file.path(input_dir, "ucell_scores.rds"))
  coverage <- read.delim(file.path(input_dir, "feature_coverage.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(setequal(unique(scores$key), coverage$feature_set[coverage$included]),
    length(pt$value) == 12642L, all(is.finite(scores$value)))

  log_event("Calculating pooled sliding-window native median/IQR summaries")
  result <- sliding_window_signature_summary(scores, pt,
    unlist(settings$width_fractions), settings$step_fraction,
    as.integer(settings$min_cells))
  result$parameters$primary_width_fraction <- settings$primary_width_fraction
  write_json(result$parameters, "parameters.json")
  write_tsv(result$summaries, "sliding_window_summary.tsv")
  write_tsv(result$windows, "window_support.tsv")
  saveRDS(result$membership, file.path(report_dir, "window_membership.rds"))
  saveRDS(result$native_summaries,
    file.path(report_dir, "native_window_summaries.rds"))
  write_tsv(data.frame(cell_id = pt$cell_id, pt = pt$value),
    "source_pseudotime.tsv")

  all_summaries <- result$summaries
  primary_fraction <- settings$primary_width_fraction
  primary <- all_summaries[all_summaries$width_fraction == primary_fraction &
    all_summaries$eligible_window, , drop = FALSE]
  primary_windows <- result$windows[
    result$windows$width_fraction == primary_fraction, , drop = FALSE]
  stopifnot(nrow(primary) > 0L,
    setequal(unique(primary$feature_set), coverage$feature_set[coverage$included]))
  primary$cytokine <- coverage$cytokine[match(primary$feature_set,
    coverage$feature_set)]
  primary$direction <- factor(coverage$direction[match(primary$feature_set,
    coverage$feature_set)], levels = c("up", "down"))

  label_map <- setNames(paste(coverage$cytokine,
    toupper(coverage$direction)), coverage$feature_set)
  label_feature <- function(x) unname(label_map[x])
  low <- primary_windows[!primary_windows$eligible_window, , drop = FALSE]
  half_step <- result$parameters$step_pt_units / 2
  low$left <- pmax(min(pt$value), low$pseudotime - half_step)
  low$right <- pmin(max(pt$value), low$pseudotime + half_step)
  palette <- c(up = "#B2182B", down = "#2166AC")
  primary_width <- diff(range(pt$value)) * primary_fraction
  caption <- sprintf(paste0(
    "Window: %.0f%% of pt range (%.3f units); step: %.0f%% (%.3f units). ",
    "Line: median; band: 25th–75th percentiles.\nAt least %d cells per window; ",
    "gray marks unsupported centers. Endpoint windows are clipped."),
    primary_fraction * 100, primary_width, settings$step_fraction * 100,
    result$parameters$step_pt_units, settings$min_cells)
  gap_layer <- function() geom_rect(data = low,
    aes(xmin = left, xmax = right, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#64748B", alpha = 0.13)
  curve_plot <- function(data, title) {
    ggplot(data, aes(pseudotime, score_median, group = feature_set,
      color = direction, fill = direction)) + gap_layer() +
      geom_ribbon(aes(ymin = score_q25, ymax = score_q75),
        alpha = 0.20, colour = NA) + geom_line(linewidth = 0.85) +
      scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
      labs(title = title, x = "Existing Bdev pseudotime (pt)",
        y = "UCell score", color = "Parse direction",
        fill = "Parse direction", caption = caption) +
      theme(legend.position = "bottom")
  }

  overview <- curve_plot(primary,
    "Passing competitor-aware cytokine signatures along Bdev pseudotime") +
    facet_wrap(~feature_set, scales = "free_y", ncol = 4,
      labeller = labeller(feature_set = label_feature)) +
    labs(subtitle = paste("All", length(pt$value),
      "cells with defined source pseudotime pooled together"))
  built <- ggplot_build(overview)
  stopifnot(any(vapply(built$data, function(x) nrow(x) == nrow(primary),
    logical(1))))
  save_plot(overview, "sliding_window_iqr", 14, 18)

  # Render the 5% overview from the already-calculated native summaries.
  five <- all_summaries[all_summaries$width_fraction == 0.05 &
    all_summaries$eligible_window, , drop = FALSE]
  five$direction <- factor(coverage$direction[match(five$feature_set,
    coverage$feature_set)], levels = c("up", "down"))
  five_windows <- result$windows[result$windows$width_fraction == 0.05, ]
  five_low <- five_windows[!five_windows$eligible_window, , drop = FALSE]
  five_low$left <- pmax(min(pt$value), five_low$pseudotime - half_step)
  five_low$right <- pmin(max(pt$value), five_low$pseudotime + half_step)
  five_caption <- sprintf(paste0(
    "Window: 5%% of pt range (%.3f units); step: %.0f%% (%.3f units). ",
    "Line: median; band: 25th–75th percentiles.\nGray marks windows below %d cells; endpoint windows are clipped."),
    0.05 * diff(range(pt$value)), settings$step_fraction * 100,
    result$parameters$step_pt_units, settings$min_cells)
  five_plot <- ggplot(five, aes(pseudotime, score_median,
      group = feature_set, color = direction, fill = direction)) +
    geom_rect(data = five_low,
      aes(xmin = left, xmax = right, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "#64748B", alpha = 0.13) +
    geom_ribbon(aes(ymin = score_q25, ymax = score_q75),
      alpha = 0.20, colour = NA) + geom_line(linewidth = 0.85) +
    scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
    facet_wrap(~feature_set, scales = "free_y", ncol = 4,
      labeller = labeller(feature_set = label_feature)) +
    labs(title = "Passing signatures: 5% sliding-window median and IQR",
      subtitle = paste(length(unique(five$feature_set)),
        "passing competitor-aware gene sets; Bdev coverage is reported separately"),
      x = "Existing Bdev pseudotime (pt)", y = "UCell score",
      color = "Parse direction", fill = "Parse direction",
      caption = five_caption) + theme(legend.position = "bottom")
  save_plot(five_plot, "sliding_window_iqr_5pct", 14, 18)

  # One page per cytokine, preserving separate up/down curves when both pass.
  cytokines <- unique(coverage$cytokine[coverage$included])
  multipage <- file.path(figure_dir,
    "cytokine_competitor_aware_sliding_window_iqr.pdf")
  grDevices::pdf(multipage, width = 11, height = 6.5)
  for (cytokine in cytokines) {
    data <- primary[primary$cytokine == cytokine, , drop = FALSE]
    plot <- curve_plot(data,
      paste("Competitor-aware response signature:", cytokine)) +
      facet_wrap(~direction, scales = "free_y", nrow = 1,
        drop = TRUE)
    print(plot)
    ggsave(file.path(figure_dir,
      paste0("cytokine_", slug(cytokine), ".png")), plot,
      width = 11, height = 6.5, dpi = 180, bg = "white")
  }
  grDevices::dev.off()

  sensitivity <- all_summaries[all_summaries$eligible_window, , drop = FALSE]
  sensitivity$width <- factor(paste0(100 * sensitivity$width_fraction, "%"),
    levels = paste0(100 * unlist(settings$width_fractions), "%"))
  sensitivity$label <- factor(label_feature(sensitivity$feature_set),
    levels = unique(label_feature(coverage$feature_set[coverage$included])))
  width_plot <- ggplot(sensitivity,
      aes(pseudotime, score_median, color = width,
        group = interaction(feature_set, width))) +
    geom_line(linewidth = 0.65) +
    facet_wrap(~label, scales = "free_y", ncol = 4) +
    labs(title = "Sensitivity to sliding-window width",
      x = "Existing Bdev pseudotime (pt)", y = "Median UCell score",
      color = "Window width",
      caption = paste("All widths use the same pooled cells, 1% steps and at least",
        settings$min_cells, "cells per window. No extra smoother is fitted."))
  save_plot(width_plot, "window_width_sensitivity", 14, 18)

  support <- result$windows
  support$width <- factor(paste0(support$width_fraction * 100, "%"),
    levels = paste0(100 * unlist(settings$width_fractions), "%"))
  support_plot <- ggplot(support,
      aes(pseudotime, n_cells, color = width, group = width)) +
    geom_line() + geom_hline(yintercept = settings$min_cells,
      linetype = "dotted") +
    labs(title = "Pooled cells in each overlapping window",
      x = "Existing Bdev pseudotime (pt)", y = "Cells per window",
      color = "Window width",
      caption = "Endpoint windows are clipped; overlapping counts must not be summed as unique cells.")
  save_plot(support_plot, "window_cell_support", 11, 5)

  support_summary <- do.call(rbind, lapply(unlist(settings$width_fractions),
    function(width) {
      x <- support[support$width_fraction == width, ]
      data.frame(width_percent = 100 * width,
        width_pt_units = width * diff(range(pt$value)), centers = nrow(x),
        retained_centers = sum(x$eligible_window),
        unsupported_centers = sum(!x$eligible_window),
        min_cells = min(x$n_cells), max_cells = max(x$n_cells))
    }))
  write_tsv(support_summary, "width_support_summary.tsv")

  checks <- c(result$checks, list(pseudotime_exact = TRUE,
    pseudotime_missingness_exact = TRUE, source_memberships_accepted = TRUE,
    scoring_session_unpromoted = TRUE, bundle_manifest_unchanged = TRUE,
    all_passing_signatures_plotted =
      setequal(unique(primary$feature_set), coverage$feature_set[coverage$included])))
  stopifnot(identical(wb_manifest(bundle), manifest),
    all(unlist(checks[vapply(checks, is.logical, logical(1))])))
  write_json(list(status = "passed", checks = checks,
    scope = "competitor-aware UCell trajectory review; no promotion"),
    "review_checks.json")

  source_files <- c(
    "processing/integration/06_score_competitor_aware_signatures.R",
    "figure_scripts/bdev_integrated/05_competitor_aware_sliding_window_report.R",
    "src/competitor_aware_integration_workflow.R",
    "src/sliding_window_signature_summary.R",
    "config/bdev_integrated.yml", "config/bundles.yml")
  input_files <- c("scoring.json", "feature_coverage.tsv", "ucell_scores.rds",
    "exact_pt_family.rds", "accepted_source_gene_sets.tsv")
  write_json(list(run_id = integration_run,
    bundle_dir = bundle_dir,
    bundle_state_fingerprint = manifest$current_state$fingerprint,
    parents = parents, runtime = integration_runtime,
    source_bundle = "parse_cytokines_rest_only",
    source_namespace = competitor_cfg$source_namespace,
    source_membership_sha256 = wb_source_sha256(source_membership_tsv),
    parameters = result$parameters,
    execution = paste("native scWorkbench UCell in a durable unpromoted session;",
      "native scWorkbench distribution summaries on overlapping pooled cell subsets"),
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    sources = lapply(source_files, function(path) list(path = path,
      sha256 = wb_source_sha256(path))),
    inputs = lapply(input_files, function(name) list(
      path = file.path(input_dir, name),
      sha256 = wb_source_sha256(file.path(input_dir, name))))),
    "review_provenance.json")
  file.copy("config/bdev_integrated.yml",
    file.path(report_dir, "review_configuration.yml"), overwrite = TRUE)
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))

  scored_n <- sum(coverage$included)
  below_gate <- coverage[!coverage$coverage_gate_pass, , drop = FALSE]
  coverage_text <- if (!nrow(below_gate)) {
    "All 33 passing Parse signatures also met the previous Bdev coverage gate."
  } else {
    paste0(nrow(below_gate), " signatures were retained but fall below the previous Bdev coverage gate: ",
      paste(paste0(below_gate$cytokine, " ", toupper(below_gate$direction),
        " (", below_gate$mapped_genes, "/", below_gate$original_genes, ")"),
        collapse = "; "), ".")
  }
  image <- function(name) paste0('<img src="figures/', name,
    '.png" alt="', name, '">')
  html_table <- function(x) as.character(knitr::kable(x, format = "html",
    digits = 4, escape = TRUE, row.names = FALSE))
  coverage_display <- coverage[, c("cytokine", "direction", "selected_size",
    "original_genes", "mapped_genes", "coverage", "coverage_gate_pass",
    "included")]
  html <- c(
    '<!doctype html><html><head><meta charset="utf-8"><title>Competitor-aware cytokine trajectories</title>',
    '<style>body{font:16px system-ui,sans-serif;max-width:1200px;margin:36px auto;padding:0 24px;color:#172033}img{max-width:100%;height:auto}p,li{line-height:1.55}table{border-collapse:collapse;font-size:13px}th,td{border:1px solid #ddd;padding:6px;text-align:left}a{color:#1261a0}</style></head><body>',
    '<h1>Passing competitor-aware cytokine signatures along Bdev pseudotime</h1>',
    sprintf(paste0('<p><strong>%d of 33</strong> passing Parse full signatures were scored on all 26,723 Bdev cells and profiled over the 12,642 cells with defined original <code>pt</code>. ',
      '%s Exact reference pseudotime values and all 14,081 missing entries were preserved.</p>'), scored_n, coverage_text),
    sprintf(paste0('<p>The primary window spans %.0f%% of the observed pt range (%.3f units), advances by %.0f%% (%.3f units), and requires at least %d cells. ',
      'Lines show the within-window median UCell score and bands show the 25th–75th percentiles. Windows are clipped at trajectory endpoints.</p>'),
      100 * primary_fraction, primary_width, 100 * settings$step_fraction,
      result$parameters$step_pt_units, settings$min_cells),
    '<p>All reference samples and cell-cycle states are pooled. The bands describe cell-to-cell spread, not confidence in the median. Neighboring windows overlap and are not independent. Raw UCell values are not calibrated for comparing absolute activity between different signatures.</p>',
    '<h2>Primary 10% windows</h2>', image("sliding_window_iqr"),
    '<p><a href="figures/sliding_window_iqr.pdf">Download overview PDF</a> · <a href="figures/cytokine_competitor_aware_sliding_window_iqr.pdf">Download per-cytokine PDF</a></p>',
    '<h2>5% windows</h2>', image("sliding_window_iqr_5pct"),
    '<h2>Window-width sensitivity</h2>', image("window_width_sensitivity"),
    '<h2>Window support</h2>', image("window_cell_support"),
    html_table(support_summary),
    '<h2>Bdev feature coverage</h2>', html_table(coverage_display),
    '<h2>Methods and provenance</h2>',
    '<p>The accepted full memberships from <code>parse_cytokines_rest_only</code> were verified gene-for-gene against its read-only feature-set registry. Bdev scoring used raw RNA counts, native scWorkbench UCell, max_rank=1500, missing-feature dropping, and no graph smoothing. The previous five-gene/80% Bdev coverage rule is retained as a diagnostic column, while all 33 passing Parse signatures are included as requested. Curves below that threshold rely on fewer mapped genes and should be interpreted with added caution. Stable core sets were not substituted for the full sets.</p>',
    '<p>Window membership uses inclusive bounds on exact <code>Bdev@meta.data$pt</code> values. Native <code>cell_vector_family_qc</code> summaries provide medians, type-7 quartiles, means and sample SD. No pseudotime refitting, imputation, extra smoothing, bootstrap, or hypothesis test was applied. Scores and imports remain in a durable review session and were not promoted.</p>',
    '<ul><li><a href="../feature_coverage.tsv">Feature coverage</a></li>',
    '<li><a href="sliding_window_summary.tsv">Window-level scores</a></li>',
    '<li><a href="window_support.tsv">Window support</a></li>',
    '<li><a href="review_checks.json">Review checks</a></li>',
    '<li><a href="review_provenance.json">Provenance</a></li>',
    '<li><a href="review_configuration.yml">Configuration</a></li></ul>',
    paste0('<p>Run: <code>', integration_run, '</code>. Integration bundle fingerprint: <code>',
      manifest$current_state$fingerprint, '</code>.</p></body></html>'))
  writeLines(html, file.path(report_dir, "REPORT.html"))
  writeLines(c("# Competitor-aware cytokine trajectory report", "",
    "[Open report](REPORT.html)", "",
    sprintf("All %d passing full signatures were profiled; Bdev gene coverage is reported as a diagnostic.", scored_n),
    "Pooled median and IQR use 10% windows with 1% steps; 5% and 15% sensitivity is included.",
    "Exact source pt is preserved. No integration result was promoted."),
    file.path(report_dir, "REPORT.md"))
  write_json(list(status = "complete", promoted = FALSE,
    review_status = "awaiting_user_review", population = "pooled",
    method = "sliding_window_median_iqr", source_signatures = 33L,
    scored_signatures = scored_n,
    below_prior_coverage_gate = coverage$feature_set[!coverage$coverage_gate_pass],
    coverage_policy = "diagnostic_only_all_33_passing_sets_retained",
    trajectory_cells = length(pt$value),
    primary_retained_centers = sum(primary_windows$eligible_window),
    figures = list.files(figure_dir),
    report_sha256 = wb_source_sha256(file.path(report_dir, "REPORT.html"))),
    "figures.json")

  final_report_dir <- sub("[.]incomplete$", "", report_dir)
  if (!file.rename(report_dir, final_report_dir)) {
    stop("Failed to atomically finalize report directory")
  }
  log_event("Competitor-aware trajectory report complete:",
    file.path(final_report_dir, "REPORT.html"))
}

main()
