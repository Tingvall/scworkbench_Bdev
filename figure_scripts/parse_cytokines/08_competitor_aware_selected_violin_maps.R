#!/usr/bin/env Rscript

# Render the same cell-level, replicate-separated UCell violin maps used for
# the previous rest-only candidates, now for all 52 competitor-aware selected
# gene sets. Completed native UCell scores are reused; nothing is promoted.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 08_competitor_aware_selected_violin_maps.R <new_output_directory>")
}

output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
if (file.exists(output_dir)) stop("Output path already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir)) stop("Incomplete output path already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
  library(yaml)
})
source("src/bundle_config.R")
source("src/bundle_validation.R")

registry <- load_bundle_registry()
cfg <- read_project_yaml("config/parse_cytokines_rest_only.yml")
runtime <- registry$bundles$parse_cytokines$runtime
stopifnot(identical(as.character(packageVersion("scworkbench")), runtime$scworkbench_version))

review_dir <- file.path(cfg$storage_root, "reports", cfg$run_id,
  "competitor_aware_signature_review_v3")
input_paths <- c(
  signature_summary = file.path(review_dir, "signature_summary.tsv"),
  final_membership = file.path(review_dir, "final_gene_membership.tsv"),
  core_membership = file.path(review_dir, "core_gene_membership.tsv"),
  analysis_completion = file.path(review_dir, "completion.json"),
  analysis_methods = file.path(review_dir, "METHODS.md"),
  project_config = "config/parse_cytokines_rest_only.yml",
  bundle_registry = "config/bundles.yml",
  renderer = "figure_scripts/parse_cytokines/08_competitor_aware_selected_violin_maps.R",
  bundle_manifest = file.path(cfg$bundle_dir, "manifest.json")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) stop("Missing inputs: ", paste(missing_inputs, collapse = ", "))

read_tsv <- function(path) read.delim(path, check.names = FALSE, na.strings = c("NA", ""))
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
fmt <- function(x) ifelse(is.finite(x), formatC(x, digits = 3L, format = "f"), "NA")

signatures <- read_tsv(input_paths[["signature_summary"]])
members <- read_tsv(input_paths[["final_membership"]])
cores <- read_tsv(input_paths[["core_membership"]])
analysis_completion <- jsonlite::read_json(input_paths[["analysis_completion"]],
  simplifyVector = TRUE)
stopifnot(nrow(signatures) == 52L, !anyDuplicated(signatures$signature_id),
  all(signatures$full_status %in% c("PASS", "FAIL")),
  all(signatures$direction %in% c("up", "down")),
  nrow(members) == sum(signatures$signature_size),
  all(table(members$signature_id) == signatures$signature_size[
    match(names(table(members$signature_id)), signatures$signature_id)]),
  isTRUE(analysis_completion$source_bundle_unchanged),
  identical(analysis_completion$source_bundle, cfg$bundle_dir))

candidate_sets <- setNames(lapply(signatures$signature_id, function(id) {
  local <- members[members$signature_id == id, , drop = FALSE]
  local <- local[order(local$rank, local$feature_id), , drop = FALSE]
  genes <- as.character(local$feature_id)
  expected <- signatures$signature_size[signatures$signature_id == id]
  if (length(genes) != expected || anyDuplicated(genes)) {
    stop("Invalid final membership for ", id)
  }
  genes
}), signatures$signature_id)

parent_check <- validate_registered_bundle(registry, "parse_cytokines",
  run_recipient_validation = FALSE)
stopifnot(isTRUE(jsonlite::read_json(registry$bundles$parse_cytokines$validation_report,
  simplifyVector = TRUE)$ok))
bundle <- wb_open_bundle(cfg$bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
bundle_fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(bundle_fingerprint_before,
  analysis_completion$source_bundle_fingerprint))

cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "sample_id", "condition", "replicate",
  "treatment_category"
))
stopifnot(length(cells) == 154209L, identical(as.character(metadata$cell_id), cells),
  !anyNA(metadata))
sample_metadata <- unique(metadata[, c("sample_id", "condition", "replicate",
  "treatment_category")])
sample_metadata <- sample_metadata[match(unique(metadata$sample_id), sample_metadata$sample_id), ]
stopifnot(nrow(sample_metadata) == 81L, all(table(sample_metadata$condition) == 3L),
  all(table(sample_metadata$replicate) == 27L))

# Reassemble each full-data selected score from the three replicate-specific
# sessions used by the competitor-aware analysis. UCell is cell-wise, so this
# is the same score vector as scoring the set over all cells in one operation.
score_matrix <- matrix(NA_real_, nrow = length(cells), ncol = nrow(signatures),
  dimnames = list(NULL, signatures$signature_id))
native_rows <- list()
session_workspaces <- analysis_completion$full_session_workspaces
stopifnot(length(session_workspaces) == 6L)

for (replicate in sort(unique(as.character(metadata$replicate)))) {
  for (direction in c("up", "down")) {
    batch_id <- paste0("rep", replicate, "_", direction)
    workspace <- as.character(session_workspaces[[batch_id]])
    if (!dir.exists(workspace)) stop("Missing completed session: ", workspace)
    session <- wb_session_resume(bundle, workspace)
    runs <- session$runs
    score_runs <- names(runs)[vapply(runs, function(x) {
      identical(x$op_name, "ucell") &&
        all(c("signature_scores", "signature_summary") %in% names(x$outputs))
    }, logical(1))]
    if (length(score_runs) != 1L) stop("Expected one UCell run in ", workspace)
    run <- runs[[score_runs]]
    score <- run$outputs$signature_scores
    native <- run$outputs$signature_summary$table
    batch <- signatures[signatures$direction == direction, , drop = FALSE]
    expected_cells <- cells[as.character(metadata$replicate) == replicate]

    groups <- rle(as.character(score$key))
    ends <- cumsum(groups$lengths)
    starts <- ends - groups$lengths + 1L
    index_by_key <- setNames(Map(seq.int, starts, ends), groups$values)
    stopifnot(all(batch$score_key %in% names(index_by_key)))
    for (i in seq_len(nrow(batch))) {
      key <- batch$score_key[[i]]
      id <- batch$signature_id[[i]]
      index <- index_by_key[[key]]
      score_cells <- as.character(score$cell_id[index])
      values <- as.numeric(score$value[index])
      stopifnot(length(index) == length(expected_cells),
        identical(score_cells, expected_cells), all(is.finite(values)),
        all(values >= 0 & values <= 1))
      cell_index <- match(score_cells, cells)
      score_matrix[cell_index, id] <- values

      row <- native[native$signature == key, , drop = FALSE]
      stopifnot(nrow(row) == 1L, row$n_cells == length(expected_cells),
        row$n_features_missing == 0L,
        row$n_features_matched == length(candidate_sets[[id]]),
        setequal(strsplit(row$matched_features, ";", fixed = TRUE)[[1L]],
          candidate_sets[[id]]))
      row$replicate <- replicate
      row$direction <- direction
      row$signature_id <- id
      row$selected_size <- batch$signature_size[[i]]
      native_rows[[length(native_rows) + 1L]] <- row
    }
    try(wb_close_bundle(wb_session_view(session)), silent = TRUE)
    rm(session, run, score, native, index_by_key)
    invisible(gc())
  }
}
stopifnot(!anyNA(score_matrix), all(is.finite(score_matrix)),
  all(score_matrix >= 0 & score_matrix <= 1))
native_selected <- do.call(rbind, native_rows)
write_tsv(native_selected, file.path(stage_dir, "native_selected_score_summary.tsv"))

sample_rows <- vector("list", nrow(signatures))
for (i in seq_len(nrow(signatures))) {
  values <- score_matrix[, signatures$signature_id[[i]]]
  medians <- tapply(values, metadata$sample_id, median)
  means <- tapply(values, metadata$sample_id, mean)
  local <- sample_metadata[match(names(medians), sample_metadata$sample_id), , drop = FALSE]
  stopifnot(nrow(local) == 81L, identical(names(medians), names(means)), !anyNA(local))
  local$signature_id <- signatures$signature_id[[i]]
  local$cytokine <- signatures$cytokine[[i]]
  local$direction <- signatures$direction[[i]]
  local$sample_median_ucell <- as.numeric(medians)
  local$sample_mean_ucell <- as.numeric(means)
  sign <- if (signatures$direction[[i]] == "up") 1 else -1
  local$oriented_sample_mean <- sign * local$sample_mean_ucell
  sample_rows[[i]] <- local
}
sample_summaries <- do.call(rbind, sample_rows)

category_order <- c("full_pass_up", "full_pass_down", "full_fail_up", "full_fail_down")
signatures$figure_category <- paste0("full_", tolower(signatures$full_status), "_",
  signatures$direction)
stopifnot(all(signatures$figure_category %in% category_order))
order_index <- order(match(signatures$figure_category, category_order),
  signatures$cytokine, signatures$direction)
signatures <- signatures[order_index, , drop = FALSE]
signatures$page <- seq_len(nrow(signatures))
sample_summaries$figure_category <- signatures$figure_category[
  match(sample_summaries$signature_id, signatures$signature_id)]
write_tsv(signatures, file.path(stage_dir, "selected_signature_summary.tsv"))
write_tsv(sample_summaries, file.path(stage_dir, "selected_sample_summaries.tsv"))
write_tsv(members, file.path(stage_dir, "selected_gene_sets.tsv"))
write_tsv(cores, file.path(stage_dir, "selected_core_gene_sets.tsv"))

plot_candidate <- function(info) {
  id <- info$signature_id[[1L]]
  target <- info$cytokine[[1L]]
  direction <- info$direction[[1L]]
  values <- score_matrix[, id]
  condition_order <- c("control", target,
    sort(setdiff(unique(metadata$condition), c("control", target))))
  local_samples <- sample_metadata[order(match(sample_metadata$condition, condition_order),
    suppressWarnings(as.numeric(sample_metadata$replicate)), sample_metadata$sample_id), , drop = FALSE]
  ymax <- max(0.02, max(values) * 1.06)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(10.5, 4.8, 5.3, 1.3), mgp = c(2.8, 0.8, 0), xpd = FALSE)
  stable <- info$holdout_mean_pass[[1L]] || info$holdout_responder_pass[[1L]]
  core_label <- if (is.na(info$core_size[[1L]])) "-" else info$core_size[[1L]]
  plot(NA_real_, NA_real_, xlim = c(0.35, length(condition_order) + 0.65),
    ylim = c(0, ymax), xaxt = "n", xlab = "", ylab = "Per-cell UCell score",
    las = 1, bty = "l",
    main = paste0(target, " ", toupper(direction), " — full-data ", info$full_status[[1L]]),
    sub = paste0(id, " | ", info$signature_size[[1L]], " genes | selected by ",
      info$selection_route[[1L]], " | mean gap ", fmt(info$full_mean_median_standardized_gap[[1L]]),
      " | responder gap ", fmt(info$full_responder_median_standardized_gap[[1L]]),
      " | holdout ", if (stable) "stable" else "not stable", " | core ", core_label),
    cex.main = 1.05, cex.sub = 0.77)
  grid(nx = NA, ny = NULL, col = "#E8EDF0", lty = 1)

  sample_offsets <- c(-0.22, 0, 0.22)
  for (condition_index in seq_along(condition_order)) {
    condition <- condition_order[[condition_index]]
    condition_samples <- local_samples[local_samples$condition == condition, , drop = FALSE]
    stopifnot(nrow(condition_samples) == 3L)
    fill <- if (condition == "control") "#B8B8B8" else if (condition == target) {
      if (info$full_status[[1L]] == "PASS") "#D95F02" else "#E6AB02"
    } else "#80B1D3"
    edge <- if (condition == "control") "#4D4D4D" else if (condition == target) {
      "#8C2D04"
    } else "#2166AC"
    for (sample_index in seq_len(nrow(condition_samples))) {
      sample_id <- condition_samples$sample_id[[sample_index]]
      sample_values <- values[metadata$sample_id == sample_id]
      x <- condition_index + sample_offsets[[sample_index]]
      if (diff(range(sample_values)) > .Machine$double.eps^0.5) {
        den <- density(sample_values, n = 256L, bw = "nrd0", from = 0, to = ymax)
        keep <- which(is.finite(den$y) & den$y >= max(den$y, na.rm = TRUE) * 0.002)
        if (length(keep) >= 2L) {
          keep <- seq.int(min(keep), max(keep))
          half_width <- 0.09 * den$y[keep] / max(den$y[keep])
          polygon(c(x - half_width, rev(x + half_width)),
            c(den$x[keep], rev(den$x[keep])),
            col = adjustcolor(fill, alpha.f = 0.62),
            border = adjustcolor(edge, alpha.f = 0.88), lwd = 0.65)
        }
      } else {
        segments(x - 0.07, sample_values[[1L]], x + 0.07, sample_values[[1L]],
          col = edge, lwd = 3)
      }
      points(x, median(sample_values), pch = 21, bg = "white", col = "#111111",
        cex = 0.62, lwd = 0.75)
    }
  }
  axis(1, at = seq_along(condition_order), labels = FALSE, tck = -0.012)
  usr <- par("usr")
  label_y <- usr[[3L]] - 0.047 * diff(usr[3:4])
  label_col <- ifelse(condition_order == target, "#8C2D04",
    ifelse(condition_order == "control", "#222222", "#333333"))
  text(seq_along(condition_order), label_y, condition_order, srt = 55, adj = 1,
    cex = 0.70, xpd = NA, col = label_col)
  mtext("Condition (three biological samples per condition)", side = 1,
    line = 8.9, cex = 0.88)
  expected <- if (direction == "up") "higher target score expected" else
    "lower target score expected"
  legend("topright", inset = 0.005, bty = "n", cex = 0.76,
    legend = c("control", "target cytokine", "other cytokines", "sample median", expected),
    pch = c(15, 15, 15, 21, NA),
    pt.bg = c("#B8B8B8", if (info$full_status[[1L]] == "PASS") "#D95F02" else
      "#E6AB02", "#80B1D3", "white", NA),
    col = c("#4D4D4D", "#8C2D04", "#2166AC", "#111111", "#333333"),
    lty = c(NA, NA, NA, NA, 1), lwd = c(NA, NA, NA, NA, 1))
}

render_pdf <- function(path, rows) {
  pdf(path, width = 14, height = 8.5, onefile = TRUE, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  if (!length(rows)) {
    plot.new()
    text(0.5, 0.55, "No signatures in this category", cex = 1.4, font = 2)
  } else {
    for (i in rows) plot_candidate(signatures[i, , drop = FALSE])
  }
}

combined_pdf <- file.path(stage_dir, "all_52_competitor_aware_ucell_violin_maps.pdf")
render_pdf(combined_pdf, seq_len(nrow(signatures)))
category_files <- setNames(character(length(category_order)), category_order)
for (category in category_order) {
  rows <- which(signatures$figure_category == category)
  category_files[[category]] <- paste0(category, "_violin_maps.pdf")
  render_pdf(file.path(stage_dir, category_files[[category]]), rows)
}

category_summary <- do.call(rbind, lapply(category_order, function(category) {
  rows <- signatures$figure_category == category
  data.frame(category = category, signatures = sum(rows),
    pdf = category_files[[category]], stringsAsFactors = FALSE)
}))
write_tsv(category_summary, file.path(stage_dir, "category_summary.tsv"))

input_manifest <- data.frame(path = normalizePath(input_paths, mustWork = TRUE),
  sha256 = vapply(input_paths, wb_source_sha256, character(1)), stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))

report_lines <- c(
  "# Competitor-aware selected-set UCell violin maps", "",
  paste0("All 52 full-data-selected cytokine-direction gene sets are shown with the same layout as the prior candidate plots. ",
    "Each page contains raw cell-level UCell distributions for control and all 26 cytokines. Three narrow violins represent the ",
    "biological replicates, and white points mark their cell-level medians."), "",
  paste0("The four categories use the competitor-aware full-data result: ",
    sum(signatures$full_status == "PASS"), " pass and ", sum(signatures$full_status == "FAIL"),
    " fail. The plotted medians describe the distributions; the formal pass uses the prespecified sample-mean or responder-fraction ",
    "gap and requires the target to rank strictly first in at least two of three replicates."), "",
  paste0("For down sets, lower raw UCell scores indicate the requested response. Control is displayed and was included as a competitor. ",
    "All scores are reused from the six completed native scWorkbench sessions, and all memberships are the selected full-data memberships."),
  "", "## Outputs", "",
  "- `all_52_competitor_aware_ucell_violin_maps.pdf`: all sets ordered pass/up, pass/down, fail/up, fail/down.",
  "- `full_pass_up_violin_maps.pdf`: passing up sets.",
  "- `full_pass_down_violin_maps.pdf`: passing down sets.",
  "- `full_fail_up_violin_maps.pdf`: failed up sets.",
  "- `full_fail_down_violin_maps.pdf`: failed down sets.",
  "- `selected_signature_summary.tsv`: selected size, full/holdout results, category, and page.",
  "- `selected_sample_summaries.tsv`: sample means and medians for every plotted set.",
  "- `selected_gene_sets.tsv` and `selected_core_gene_sets.tsv`: exact plotted and core memberships.", "",
  "No score or gene set was promoted, and the source bundle is unchanged."
)
writeLines(report_lines, file.path(stage_dir, "REPORT.md"))

bundle_fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(bundle_fingerprint_before, bundle_fingerprint_after),
  nrow(signatures) == 52L, nrow(sample_summaries) == 52L * 81L,
  all(table(sample_summaries$signature_id) == 81L),
  all(table(interaction(sample_summaries$signature_id, sample_summaries$condition,
    drop = TRUE)) == 3L),
  nrow(native_selected) == 52L * 3L,
  file.exists(combined_pdf), all(file.exists(file.path(stage_dir, category_files))))

completion <- list(
  status = "complete", created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  dataset_id = cfg$bundle_id, source_run = cfg$run_id,
  source_analysis = review_dir, signatures = nrow(signatures),
  full_pass = sum(signatures$full_status == "PASS"),
  full_fail = sum(signatures$full_status == "FAIL"),
  category_counts = setNames(as.list(category_summary$signatures), category_summary$category),
  cells_per_signature = length(cells), conditions = length(unique(metadata$condition)),
  biological_samples = nrow(sample_metadata), native_sessions_reused = length(session_workspaces),
  result_promotions = 0L, input_bundle_state_unchanged = TRUE,
  bundle_fingerprint = bundle_fingerprint_before, accepted_parent = parent_check,
  runtime = runtime
)
write_json(completion, file.path(stage_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

output_files <- list.files(stage_dir, full.names = TRUE, recursive = TRUE)
output_files <- output_files[basename(output_files) != "output_manifest.tsv"]
output_manifest <- data.frame(
  path = substring(output_files, nchar(stage_dir) + 2L),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(output_manifest[order(output_manifest$path), ],
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not atomically finalize output directory")
cat("Completed competitor-aware selected-set violin maps:", output_dir, "\n")
