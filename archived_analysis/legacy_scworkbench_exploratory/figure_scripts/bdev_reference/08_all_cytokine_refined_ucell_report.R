#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 08_all_cytokine_refined_ucell_report.R",
  "<bdev_bundle> <parse_refinement_report> <output_dir>",
  "<expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
parse_report <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
expected_fingerprint <- args[[4L]]
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(jsonlite)
})
source("src/sliding_window_signature_summary.R")

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  na.strings = c("NA", ""), stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}
slug <- function(x) gsub("(^_+|_+$)", "",
  tolower(gsub("[^A-Za-z0-9]+", "_", x)))
as_flag <- function(x) !is.na(x) & tolower(as.character(x)) %in%
  c("true", "t", "1", "yes")
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
ncores <- as.integer(Sys.getenv("BDEV_THREADS", "4"))
stopifnot(is.finite(ncores), ncores >= 1L)

source_completion_path <- file.path(parse_report, "completion.json")
summary_path <- file.path(parse_report, "signature_summary.tsv")
membership_path <- file.path(parse_report, "final_gene_membership.tsv")
stopifnot(file.exists(source_completion_path), file.exists(summary_path),
  file.exists(membership_path))
source_completion <- read_json(source_completion_path, simplifyVector = TRUE)
stopifnot(identical(source_completion$status, "complete"),
  identical(source_completion$analysis_id,
    "all_cytokine_competitor_refinement_top100_v1"),
  identical(source_completion$promoted, FALSE),
  identical(as.integer(source_completion$candidate_top_n), 100L),
  identical(as.integer(source_completion$minimum_final_size), 20L))

inventory <- read_tsv(summary_path)
required_inventory <- c("signature_id", "cytokine", "direction",
  "construction_ok", "final_size", "final_mean_pass",
  "final_responder_pass", "final_overall_pass")
stopifnot(all(required_inventory %in% names(inventory)), nrow(inventory) == 52L,
  !anyDuplicated(inventory$signature_id),
  setequal(inventory$direction, c("up", "down")))
passing <- inventory[as_flag(inventory$final_overall_pass), , drop = FALSE]
if (!nrow(passing)) stop("No Parse signatures passed the final tests")
stopifnot(all(as_flag(passing$construction_ok)),
  all(as.integer(passing$final_size) >= 20L))
passing <- passing[order(passing$cytokine,
  match(passing$direction, c("up", "down"))), ]

members <- read_tsv(membership_path)
required_members <- c("signature_id", "cytokine", "direction", "final_rank",
  "feature_id")
stopifnot(all(required_members %in% names(members)))
members <- members[members$signature_id %in% passing$signature_id,
  required_members, drop = FALSE]
members$final_rank <- as.integer(members$final_rank)
members <- members[order(match(members$signature_id, passing$signature_id),
  members$final_rank), ]
stopifnot(setequal(unique(members$signature_id), passing$signature_id),
  !anyDuplicated(members[c("signature_id", "feature_id")]))
member_counts <- table(factor(members$signature_id,
  levels = passing$signature_id))
stopifnot(identical(as.integer(member_counts), as.integer(passing$final_size)))
signatures <- setNames(lapply(passing$signature_id, function(id)
  members$feature_id[members$signature_id == id]), passing$signature_id)
write_tsv(passing, file.path(stage_dir, "passing_signature_inventory.tsv"))
write_tsv(members, file.path(stage_dir, "passing_gene_membership.tsv"))

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint),
  identical(source_completion$bdev_bundle_fingerprint, fingerprint_before))
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt", "bb_up"))
stopifnot(length(cells) == 26723L,
  identical(as.character(metadata$cell_id), cells),
  sum(is.finite(metadata$pt)) == 12642L, sum(is.na(metadata$pt)) == 14081L)
rna_features <- as.character(wb_feature_ids(bundle, "RNA_counts"))
missing_features <- setdiff(unique(members$feature_id), rna_features)
if (length(missing_features)) stop("Passing signatures contain genes absent from Bdev RNA_counts: ",
  paste(missing_features, collapse = ", "))

# Native UCell scores are staged in a durable, deliberately unpromoted session.
session_id <- "all_cytokine_refined_top100_bdev_ucell_v1"
session_root <- file.path(dirname(bundle_dir), "sessions",
  "all_cytokine_refined_top100_v1")
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) wb_session_resume(bundle, workspace) else
  wb_session_create(bundle, workspace_root = session_root,
    session_id = session_id)
if (!is.null(session$durable$promotion_receipt)) stop("Review session was promoted")
score_runs <- names(session$runs)[vapply(session$runs, function(x)
  identical(x$op_name, "ucell"), logical(1))]
matching <- score_runs[vapply(score_runs, function(id) {
  tab <- session$runs[[id]]$outputs$signature_summary$table
  setequal(as.character(tab$signature), names(signatures))
}, logical(1))]
if (!length(matching)) {
  if (length(score_runs)) stop("Unexpected UCell run in review session")
  log_event("Scoring", length(signatures),
    "passing sets with native scWorkbench UCell")
  scored <- wb_session_run_op(session, bundle, "ucell", inputs = list(
    matrix = "RNA_counts", cell_set = "all", cell_ids = cells,
    min_cells = 10L, features = signatures, name = session_id,
    max_rank = 1500L, ncores = ncores, parallel_backend = "auto",
    chunk_size = 400L, missing_features_policy = "error"),
    notes = "Unpromoted native Bdev UCell scores for passing all-cytokine refined sets")
  session <- scored$session
  run <- scored$run
} else {
  stopifnot(length(matching) == 1L)
  run <- session$runs[[matching[[1L]]]]
  log_event("Reusing completed native UCell score run")
}
native_summary <- run$outputs$signature_summary$table
native_scores <- run$outputs$signature_scores
scores <- data.frame(key = as.character(native_scores$key),
  cell_id = as.character(native_scores$cell_id),
  value = as.numeric(native_scores$value), stringsAsFactors = FALSE)
stopifnot(nrow(scores) == length(cells) * length(signatures),
  setequal(unique(scores$key), names(signatures)), all(is.finite(scores$value)),
  all(scores$value >= 0 & scores$value <= 1),
  all(native_summary$n_features_missing == 0L))
for (key in names(signatures)) {
  local <- scores[scores$key == key, , drop = FALSE]
  stopifnot(identical(local$cell_id, cells), !anyDuplicated(local$cell_id))
}
saveRDS(native_scores, file.path(stage_dir, "ucell_scores.rds"))
write_tsv_gz(scores, file.path(stage_dir, "ucell_cell_scores.tsv.gz"))
write_tsv(native_summary, file.path(stage_dir, "native_ucell_summary.tsv"))

label_info <- passing[, c("signature_id", "cytokine", "direction", "final_size",
  "final_mean_pass", "final_responder_pass"), drop = FALSE]
label_info$panel_label <- paste0(label_info$cytokine, " ",
  toupper(label_info$direction), " — ", label_info$final_size, " genes")
label_map <- setNames(label_info$panel_label, label_info$signature_id)
direction_map <- setNames(label_info$direction, label_info$signature_id)
cytokine_map <- setNames(label_info$cytokine, label_info$signature_id)

# Existing pseudotime is used exactly; missing values remain missing.
pt <- data.frame(cell_id = cells[is.finite(metadata$pt)],
  value = as.numeric(metadata$pt[is.finite(metadata$pt)]),
  stringsAsFactors = FALSE)
change_labels <- sort(unique(na.omit(as.character(metadata$bb_up))))
stopifnot(identical(change_labels, c("up1", "up2")))
change_points <- do.call(rbind, lapply(change_labels, function(label) {
  selected <- !is.na(metadata$bb_up) & metadata$bb_up == label &
    is.finite(metadata$pt)
  values <- as.numeric(metadata$pt[selected])
  data.frame(metadata_field = "bb_up", label = label,
    n_cells = length(values), pt_min = min(values), pt_median = median(values),
    pt_max = max(values), stringsAsFactors = FALSE)
}))
log_event("Calculating 5% sliding-window summaries")
trajectory <- sliding_window_signature_summary(scores, pt,
  width_fractions = 0.05, step_fraction = 0.01, min_cells = 20L)
summaries <- trajectory$summaries
windows <- trajectory$windows
summaries$window_midpoint <- (summaries$window_left + summaries$window_right) / 2
windows$window_midpoint <- (windows$window_left + windows$window_right) / 2
summaries$direction <- factor(direction_map[summaries$feature_set],
  levels = c("up", "down"))
summaries$cytokine <- cytokine_map[summaries$feature_set]
summaries$panel_label <- factor(label_map[summaries$feature_set],
  levels = label_info$panel_label)
summaries <- summaries[order(summaries$feature_set,
  summaries$window_index), , drop = FALSE]
summaries$supported_segment <- ave(!summaries$eligible_window,
  summaries$feature_set, FUN = cumsum)
limits <- range(pt$value)
unsupported <- windows[!windows$eligible_window, , drop = FALSE]
unsupported$region_left <- unsupported$window_left
unsupported$region_right <- unsupported$window_right
palette <- c(up = "#D73027", down = "#2166AC")

trajectory_plot <- function(data, title, facet_formula, ncol = NULL) {
  p <- ggplot(data[data$eligible_window, ], aes(window_midpoint,
      score_median, group = interaction(feature_set, supported_segment),
      color = direction, fill = direction)) +
    geom_rect(data = unsupported,
      aes(xmin = region_left, xmax = region_right, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "#64748B", alpha = 0.13) +
    geom_vline(data = change_points, aes(xintercept = pt_median),
      inherit.aes = FALSE, linewidth = 0.55, linetype = "dotted",
      color = "#111827") +
    geom_ribbon(aes(ymin = score_q25, ymax = score_q75),
      alpha = 0.22, colour = NA) +
    geom_line(linewidth = 0.85) +
    scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
    labs(title = title, x = "Existing Bdev pseudotime (pt)",
      y = "UCell score", color = "Direction", fill = "Direction",
      caption = paste("5% pt-range windows, 1% steps, minimum 20 cells;",
        "median and type-7 IQR. Gray regions are unsupported; no smoothing.",
        "Dotted lines are median original pt for bb_up groups.")) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
      legend.position = "bottom", strip.text = element_text(face = "bold"))
  p + facet_wrap(facet_formula, scales = "free_y", ncol = ncol)
}

trajectory_dir <- file.path(stage_dir, "trajectory_figures")
dir.create(trajectory_dir, showWarnings = FALSE)
overview <- trajectory_plot(summaries,
  "Passing competitor-refined cytokine signals along Bdev pseudotime",
  ~panel_label, ncol = 4L)
trajectory_overview <- "all_passing_signatures_5pct_iqr.png"
ggsave(file.path(trajectory_dir, trajectory_overview), overview,
  width = 16, height = max(9, 2.6 * ceiling(nrow(label_info) / 4)),
  dpi = 220, units = "in", bg = "white", limitsize = FALSE)
trajectory_rows <- list(data.frame(cytokine = "all", direction = "all",
  file = trajectory_overview, stringsAsFactors = FALSE))
for (cytokine in unique(label_info$cytokine)) {
  local <- summaries[summaries$cytokine == cytokine, , drop = FALSE]
  p <- trajectory_plot(local, paste(cytokine,
    "competitor-refined UCell signal"), ~direction, ncol = 2L)
  file <- paste0(slug(cytokine), "_5pct_iqr.png")
  ggsave(file.path(trajectory_dir, file), p, width = 11, height = 6.5,
    dpi = 240, units = "in", bg = "white")
  trajectory_rows[[length(trajectory_rows) + 1L]] <- data.frame(
    cytokine = cytokine, direction = "up/down", file = file,
    stringsAsFactors = FALSE)
}
trajectory_inventory <- do.call(rbind, trajectory_rows)
write_tsv(summaries, file.path(stage_dir, "sliding_window_summary.tsv"))
write_tsv(windows, file.path(stage_dir, "window_support.tsv"))
write_tsv(change_points, file.path(stage_dir, "bb_up_change_points.tsv"))
write_tsv(data.frame(cell_id = pt$cell_id, pt = pt$value),
  file.path(stage_dir, "source_pseudotime.tsv"))
write_tsv(trajectory_inventory,
  file.path(stage_dir, "trajectory_figure_inventory.tsv"))
saveRDS(trajectory$membership, file.path(stage_dir, "window_membership.rds"))
write_json(trajectory$parameters, file.path(stage_dir,
  "trajectory_parameters.json"))

# Raw scores are preserved; only the UMAP display color is panel-normalized.
space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells))
  coordinates <- data.frame(cell_id = cells, x = as.numeric(space$x),
    y = as.numeric(space$y), stringsAsFactors = FALSE)
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  coordinates <- data.frame(cell_id = cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]), stringsAsFactors = FALSE)
}
score_summary <- do.call(rbind, lapply(names(signatures), function(key) {
  values <- scores$value[scores$key == key]
  positive <- values > 0
  cap <- if (any(positive)) as.numeric(quantile(values[positive], 0.99,
    names = FALSE, type = 7)) else NA_real_
  data.frame(signature_id = key, cytokine = cytokine_map[[key]],
    direction = direction_map[[key]], genes = length(signatures[[key]]),
    cells = length(values), positive_cells = sum(positive),
    positive_percent = 100 * mean(positive), score_min = min(values),
    score_median = median(values), score_mean = mean(values),
    score_max = max(values), positive_p99_display_cap = cap,
    stringsAsFactors = FALSE)
}))
scores$display_score <- NA_real_
for (key in names(signatures)) {
  ix <- scores$key == key & scores$value > 0
  cap <- score_summary$positive_p99_display_cap[
    score_summary$signature_id == key]
  if (is.finite(cap) && cap > 0) scores$display_score[ix] <-
    pmin(scores$value[ix], cap) / cap
}
umap_label_map <- setNames(vapply(names(signatures), function(key) {
  stat <- score_summary[score_summary$signature_id == key, ]
  paste0(label_map[[key]], "\n",
    sprintf("%.1f%% cells with UCell > 0", stat$positive_percent))
}, character(1)), names(signatures))
scores$panel_label <- factor(umap_label_map[scores$key],
  levels = umap_label_map[names(signatures)])
plot_data <- merge(scores, coordinates, by = "cell_id", sort = FALSE)
plot_data <- plot_data[order(plot_data$panel_label,
  ifelse(is.na(plot_data$display_score), -Inf, plot_data$display_score)), ]
stopifnot(nrow(plot_data) == length(cells) * length(signatures),
  all(is.finite(plot_data$x)), all(is.finite(plot_data$y)))

add_umap <- function(data, title, facet = TRUE, ncol = 4L) {
  p <- ggplot(data, aes(x, y)) +
    geom_point(data = data[is.na(data$display_score), ], color = "#CDD2DA",
      size = 0.12, alpha = 0.7, stroke = 0) +
    geom_point(data = data[!is.na(data$display_score), ],
      aes(color = display_score), size = 0.28, alpha = 1, stroke = 0) +
    scale_color_gradientn(colours = c("#0066FF", "#00D9FF", "#00F060",
      "#FFF000", "#FF5A00", "#FF00A8"), limits = c(0, 1),
      oob = scales::squish, na.value = "#CDD2DA", name = "Panel-scaled\nUCell") +
    coord_equal() + theme_void(base_size = 9) +
    theme(strip.text = element_text(face = "bold"), legend.position = "bottom") +
    labs(title = title, caption = paste("Native UCell on Bdev RNA_counts;",
      "each panel is scaled to its positive-score 99th percentile.",
      "Colors cannot compare absolute activity between signatures."))
  if (facet) p + facet_wrap(~panel_label, ncol = ncol) else p
}
umap_dir <- file.path(stage_dir, "umap_figures")
dir.create(umap_dir, showWarnings = FALSE)
umap_overview <- "all_passing_signatures_atac_umap.png"
ggsave(file.path(umap_dir, umap_overview), add_umap(plot_data,
  "Passing competitor-refined cytokine signals", TRUE, 4L), width = 16,
  height = max(9, 3.1 * ceiling(nrow(label_info) / 4)), dpi = 220,
  units = "in", bg = "white", limitsize = FALSE)
umap_rows <- list(data.frame(signature_id = "all", cytokine = "all",
  direction = "all", file = umap_overview, stringsAsFactors = FALSE))
for (cytokine in unique(label_info$cytokine)) {
  keys <- label_info$signature_id[label_info$cytokine == cytokine]
  local <- plot_data[plot_data$key %in% keys, , drop = FALSE]
  file <- paste0(slug(cytokine), "_up_down_atac_umap.png")
  ggsave(file.path(umap_dir, file), add_umap(local,
    paste(cytokine, "competitor-refined UCell signals"), TRUE, 2L),
    width = 12, height = 6.8, dpi = 240, units = "in", bg = "white")
  umap_rows[[length(umap_rows) + 1L]] <- data.frame(
    signature_id = "cytokine_group", cytokine = cytokine,
    direction = "up/down", file = file, stringsAsFactors = FALSE)
}
for (key in names(signatures)) {
  local <- plot_data[plot_data$key == key, , drop = FALSE]
  stat <- score_summary[score_summary$signature_id == key, ]
  title <- umap_label_map[[key]]
  file <- paste0(slug(paste(stat$cytokine, stat$direction)),
    "_atac_umap.png")
  ggsave(file.path(umap_dir, file), add_umap(local, title, FALSE),
    width = 8.5, height = 7.5, dpi = 260, units = "in", bg = "white")
  umap_rows[[length(umap_rows) + 1L]] <- data.frame(signature_id = key,
    cytokine = stat$cytokine, direction = stat$direction, file = file,
    stringsAsFactors = FALSE)
}
umap_inventory <- do.call(rbind, umap_rows)
write_tsv(score_summary, file.path(stage_dir, "ucell_score_summary.tsv"))
write_tsv(umap_inventory, file.path(stage_dir, "umap_figure_inventory.tsv"))
write_tsv_gz(scores[, c("cell_id", "key", "value", "display_score")],
  file.path(stage_dir, "ucell_cell_scores_with_display.tsv.gz"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after),
  all(file.exists(file.path(trajectory_dir, trajectory_inventory$file))),
  all(file.exists(file.path(umap_dir, umap_inventory$file))))
input_paths <- c(file.path(bundle_dir, "manifest.json"), source_completion_path,
  file.path(parse_report, "output_manifest.tsv"),
  file.path(parse_report, "METHODS.md"), summary_path, membership_path,
  script_path,
  "figure_scripts/bdev_reference/run_08_all_cytokine_refined_ucell_report.sh",
  "src/sliding_window_signature_summary.R", "config/bdev_reference.yml")
input_manifest <- data.frame(path = normalizePath(input_paths, mustWork = TRUE),
  sha256 = vapply(input_paths, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))

report <- c("# All-cytokine competitor-refined Bdev signals", "",
  paste0("Native scWorkbench UCell scores were calculated for **",
    length(signatures), "** Parse-passing signatures across **",
    length(unique(label_info$cytokine)), "** cytokines."), "",
  "- Input: exact final Parse memberships passing mean and/or responder discrimination.",
  "- Score: Bdev raw `RNA_counts`, `max_rank = 1500`, no graph smoothing.",
  "- Display: unchanged `umap.atac`; per-panel p99 color normalization.",
  "- Trajectory: unchanged source `pt`, 5% windows, 1% steps, minimum 20 cells, median and IQR.",
  "- Missing pseudotime values remain missing; no refitting, rescaling, imputation, smoothing, or promotion.", "",
  "The patterns indicate resemblance to Parse-derived response programs, not cytokine exposure or causality.")
writeLines(report, file.path(stage_dir, "REPORT.md"))
completion <- list(status = "complete", promoted = FALSE,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  source_parse_report = parse_report,
  source_parse_fingerprint = source_completion$parse_bundle_fingerprint,
  bdev_bundle = bundle_dir, bdev_bundle_fingerprint = fingerprint_before,
  bdev_bundle_unchanged = TRUE, passing_signatures = length(signatures),
  cytokines = length(unique(label_info$cytokine)), cells = length(cells),
  score = "native scWorkbench UCell on RNA_counts", max_rank = 1500L,
  ncores = ncores,
  missing_features_policy = "error", graph_smoothing = FALSE,
  embedding = "umap.atac", defined_pt_cells = nrow(pt),
  missing_pt_cells = sum(is.na(metadata$pt)), width_fraction = 0.05,
  step_fraction = 0.01, min_cells = 20L, result_promotions = 0L,
  scoring_session_workspace = workspace,
  scworkbench_version = as.character(packageVersion("scworkbench")))
write_json(completion, file.path(stage_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
output_files <- list.files(stage_dir, recursive = TRUE, full.names = TRUE)
output_files <- output_files[basename(output_files) != "output_manifest.tsv" &
  !grepl("(^|/)sessions(/|$)", substring(output_files,
    nchar(stage_dir) + 2L))]
output_manifest <- data.frame(file = substring(output_files,
  nchar(stage_dir) + 2L), bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE)
write_tsv(output_manifest[order(output_manifest$file), ],
  file.path(stage_dir, "output_manifest.tsv"))
if (!file.rename(stage_dir, output_dir)) stop("Failed to finalize output directory")
log_event("Completed all-cytokine Bdev UCell report at", output_dir)
