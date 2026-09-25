#!/usr/bin/env Rscript

# Render cell-level UCell violin maps for all 51 rest-only directions with at
# least five eligible genes. The 47 discovery-passing memberships reuse their
# completed UCell scores. The best-sized memberships for the four discovery
# failures are scored in a separate, unpromoted scWorkbench session.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 07_rest_only_all_candidate_violin_maps.R <new_output_directory>")

output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
if (file.exists(output_dir)) stop("Output path already exists; refusing to overwrite: ", output_dir)
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

base_report <- file.path(cfg$storage_root, "reports", cfg$run_id)
base_work <- file.path(cfg$storage_root, "work", cfg$run_id)
input_paths <- c(
  discovery_summary = file.path(base_report, "discovery_summary.tsv"),
  cross_validation = file.path(base_report, "cross_validation_performance.tsv"),
  eligible_rows = file.path(base_report, "eligible_rows.tsv.gz"),
  selected_genes = file.path(base_report, "selected_gene_evidence.tsv"),
  scored_status = file.path(base_report, "all_signature_status.tsv"),
  existing_scores = file.path(base_work, "ucell_scores.rds"),
  project_config = "config/parse_cytokines_rest_only.yml",
  bundle_registry = "config/bundles.yml",
  renderer = "figure_scripts/parse_cytokines/07_rest_only_all_candidate_violin_maps.R"
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) stop("Missing inputs: ", paste(missing_inputs, collapse = ", "))

read_tsv <- function(path) read.delim(path, check.names = FALSE, na.strings = c("NA", ""))
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE,
  digits = NA, null = "null", na = "null")
auc <- function(values, positive) {
  stopifnot(length(values) == length(positive), all(is.finite(values)), !anyNA(positive))
  n1 <- sum(positive); n0 <- sum(!positive)
  stopifnot(n1 > 0L, n0 > 0L)
  (sum(rank(values, ties.method = "average")[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

summary <- read_tsv(input_paths[["discovery_summary"]])
performance <- read_tsv(input_paths[["cross_validation"]])
eligible <- read_tsv(input_paths[["eligible_rows"]])
selected <- read_tsv(input_paths[["selected_genes"]])
prior_status <- read_tsv(input_paths[["scored_status"]])
existing_scores <- readRDS(input_paths[["existing_scores"]])

stopifnot(nrow(summary) == 52L, sum(summary$eligible_genes >= cfg$analysis$minimum_genes) == 51L,
  sum(summary$discovery_status == "discovery_pass") == 47L,
  length(unique(existing_scores$key)) == 47L,
  setequal(unique(existing_scores$key), summary$signature_id[summary$discovery_status == "discovery_pass"]))

# Recover the prespecified best size for every sufficient direction. Discovery
# failures previously retained no membership; these four sets are exploratory.
candidates <- summary[summary$eligible_genes >= cfg$analysis$minimum_genes, , drop = FALSE]
candidates$candidate_size <- candidates$selected_genes
candidates$membership_scope <- ifelse(candidates$discovery_status == "discovery_pass",
  "pipeline_selected", "exploratory_best_size_after_discovery_failure")
for (i in which(candidates$discovery_status != "discovery_pass")) {
  id <- candidates$signature_id[[i]]
  local <- performance[
    performance$signature_id == id & performance$signature_size <= candidates$eligible_genes[[i]],
    , drop = FALSE
  ]
  if (!nrow(local)) stop("No candidate-size performance rows for ", id)
  near <- local[local$cv_rest_auc >= max(local$cv_rest_auc) - 0.02, , drop = FALSE]
  chosen <- near[order(-near$fold_gene_jaccard, near$signature_size), , drop = FALSE][1L, ]
  candidates$candidate_size[[i]] <- chosen$signature_size
  stopifnot(isTRUE(all.equal(candidates$cv_rest_auc[[i]], chosen$cv_rest_auc, tolerance = 1e-12)),
    isTRUE(all.equal(candidates$fold_gene_jaccard[[i]], chosen$fold_gene_jaccard, tolerance = 1e-12)))
}
stopifnot(nrow(candidates) == 51L, all(candidates$candidate_size >= 5L),
  sum(candidates$membership_scope == "exploratory_best_size_after_discovery_failure") == 4L)

candidate_sets <- setNames(lapply(candidates$signature_id, function(id) {
  info <- candidates[candidates$signature_id == id, , drop = FALSE]
  if (info$discovery_status == "discovery_pass") {
    genes <- selected$feature_id[selected$signature_id == id]
  } else {
    rows <- eligible[eligible$signature_id == id, , drop = FALSE]
    rows <- rows[order(rows$rank, rows$feature_id), , drop = FALSE]
    genes <- head(rows$feature_id, info$candidate_size)
  }
  genes <- as.character(genes)
  if (length(genes) != info$candidate_size || anyDuplicated(genes)) stop("Invalid candidate membership for ", id)
  genes
}), candidates$signature_id)

candidate_members <- do.call(rbind, lapply(seq_len(nrow(candidates)), function(i) {
  id <- candidates$signature_id[[i]]
  data.frame(
    signature_id = id, cytokine = candidates$cytokine[[i]], direction = candidates$direction[[i]],
    discovery_status = candidates$discovery_status[[i]], membership_scope = candidates$membership_scope[[i]],
    rank_in_candidate = seq_along(candidate_sets[[id]]), feature_id = candidate_sets[[id]],
    stringsAsFactors = FALSE
  )
}))
write_tsv(candidate_members, file.path(stage_dir, "candidate_gene_sets.tsv"))

parent_check <- validate_registered_bundle(registry, "parse_cytokines", run_recipient_validation = FALSE)
stopifnot(isTRUE(jsonlite::read_json(registry$bundles$parse_cytokines$validation_report,
  simplifyVector = TRUE)$ok))

bundle <- wb_open_bundle(cfg$bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
bundle_fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "sample_id", "condition", "replicate", "treatment_category"
))
stopifnot(length(cells) == 154209L, identical(as.character(metadata$cell_id), cells), !anyNA(metadata),
  sum(metadata$treatment_category == "single") == 149393L)
sample_metadata <- unique(metadata[, c("sample_id", "condition", "replicate", "treatment_category")])
stopifnot(nrow(sample_metadata) == 81L, all(table(sample_metadata$condition) == 3L),
  sum(sample_metadata$treatment_category == "single") == 78L,
  sum(sample_metadata$treatment_category == "control") == 3L)

# Score only the four previously unscored exploratory sets in an isolated
# session. No output is promoted to the review bundle or accepted parent.
extra_ids <- candidates$signature_id[candidates$discovery_status != "discovery_pass"]
extra_sets <- candidate_sets[extra_ids]
session_root <- file.path(cfg$storage_root, "sessions", cfg$run_id)
session_id <- "all_51_ucell_violin_extra"
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) wb_session_resume(bundle, workspace) else
  wb_session_create(bundle, workspace_root = session_root, session_id = session_id)
on.exit(try(wb_close_bundle(wb_session_view(session)), silent = TRUE), add = TRUE)
runs <- session$runs
score_runs <- names(runs)[vapply(runs, function(x) identical(x$op_name, "ucell"), logical(1))]
if (!length(score_runs) && length(runs)) {
  score_runs <- names(runs)[vapply(runs, function(x) "signature_scores" %in% names(x$outputs), logical(1))]
}
if (!length(score_runs)) {
  scored <- wb_session_run_op(session, bundle, "ucell", inputs = list(
    matrix = "RNA_counts", cell_set = "all", min_cells = 10L, features = extra_sets,
    name = paste0(cfg$namespace, ".figure_only"), max_rank = cfg$analysis$max_rank,
    ncores = 4L, parallel_backend = "auto", chunk_size = cfg$analysis$chunk_size,
    missing_features_policy = "error"
  ), notes = paste(
    "Figure-only UCell scoring for four best-sized rest-only candidates that failed the discovery gate.",
    "No outputs may be promoted; results support the requested complete 51-candidate violin inventory."
  ))
  session <- scored$session
  extra_score <- scored$run$outputs$signature_scores
  extra_native_summary <- scored$run$outputs$signature_summary$table
} else {
  stopifnot(length(score_runs) == 1L)
  score_run <- session$runs[[score_runs]]
  extra_score <- score_run$outputs$signature_scores
  extra_native_summary <- score_run$outputs$signature_summary$table
}
stopifnot(setequal(unique(extra_score$key), extra_ids),
  length(extra_score$value) == length(cells) * length(extra_ids),
  all(is.finite(extra_score$value)), all(extra_score$value >= 0 & extra_score$value <= 1),
  setequal(extra_native_summary$signature, extra_ids), all(extra_native_summary$n_features_missing == 0L),
  all(extra_native_summary$n_cells == length(cells)))
for (id in extra_ids) {
  row <- extra_native_summary[extra_native_summary$signature == id, , drop = FALSE]
  stopifnot(row$n_features_matched == length(candidate_sets[[id]]),
    setequal(strsplit(row$matched_features, ";", fixed = TRUE)[[1L]], candidate_sets[[id]]))
}
saveRDS(extra_score, file.path(stage_dir, "exploratory_discovery_failure_ucell_scores.rds"), compress = "gzip")
write_tsv(extra_native_summary, file.path(stage_dir, "exploratory_discovery_failure_ucell_summary.tsv"))

primitive_index <- function(x) {
  stopifnot(length(x$key) == length(x$cell_id), length(x$key) == length(x$value))
  groups <- rle(as.character(x$key)); ends <- cumsum(groups$lengths); starts <- ends - groups$lengths + 1L
  setNames(Map(seq.int, starts, ends), groups$values)
}
existing_index <- primitive_index(existing_scores)
extra_index <- primitive_index(extra_score)
score_values <- function(id) {
  source <- if (id %in% names(existing_index)) existing_scores else extra_score
  index <- if (id %in% names(existing_index)) existing_index[[id]] else extra_index[[id]]
  score_cells <- as.character(source$cell_id[index])
  if (length(score_cells) != length(cells) || anyDuplicated(score_cells) || !setequal(score_cells, cells)) {
    stop("Invalid score coverage for ", id)
  }
  values <- as.numeric(source$value[index])[match(cells, score_cells)]
  if (anyNA(values) || any(!is.finite(values))) stop("Invalid score values for ", id)
  values
}

sample_rows <- vector("list", nrow(candidates))
ucell_auc <- target_group_median <- rest_group_median <- control_group_median <- numeric(nrow(candidates))
for (i in seq_len(nrow(candidates))) {
  id <- candidates$signature_id[[i]]
  values <- score_values(id)
  medians <- tapply(values, metadata$sample_id, median)
  local <- sample_metadata[match(names(medians), sample_metadata$sample_id), , drop = FALSE]
  stopifnot(nrow(local) == 81L, !anyNA(local))
  local$signature_id <- id
  local$cytokine <- candidates$cytokine[[i]]
  local$direction <- candidates$direction[[i]]
  local$sample_median_ucell <- as.numeric(medians)
  direction_sign <- if (candidates$direction[[i]] == "up") 1 else -1
  local$oriented_sample_median <- direction_sign * local$sample_median_ucell
  single <- local$treatment_category == "single"
  target <- single & local$condition == candidates$cytokine[[i]]
  rest <- single & !target
  stopifnot(sum(target) == 3L, sum(rest) == 75L)
  ucell_auc[[i]] <- auc(local$oriented_sample_median[single], target[single])
  target_group_median[[i]] <- median(local$sample_median_ucell[target])
  rest_group_median[[i]] <- median(local$sample_median_ucell[rest])
  control_group_median[[i]] <- median(local$sample_median_ucell[local$treatment_category == "control"])
  sample_rows[[i]] <- local
}
sample_medians <- do.call(rbind, sample_rows)
candidates$ucell_rest_auc <- ucell_auc
candidates$ucell_check <- ifelse(candidates$ucell_rest_auc >= cfg$analysis$ucell_rest_auc_minimum, "PASS", "FAIL")
candidates$target_sample_median <- target_group_median
candidates$other_single_sample_median <- rest_group_median
candidates$control_sample_median <- control_group_median
candidates$figure_category <- paste("UCell", tolower(candidates$ucell_check), candidates$direction, sep = "_")
category_order <- c("UCell_pass_up", "UCell_pass_down", "UCell_fail_up", "UCell_fail_down")
candidates$figure_category <- factor(candidates$figure_category, levels = category_order)
stopifnot(!anyNA(candidates$figure_category))

# The 47 existing candidates must reproduce their recorded AUC and decisions.
prior <- candidates$discovery_status == "discovery_pass"
matching <- match(candidates$signature_id[prior], prior_status$signature_id)
stopifnot(!anyNA(matching),
  max(abs(candidates$ucell_rest_auc[prior] - prior_status$ucell_rest_auc[matching])) < 1e-12,
  identical(candidates$ucell_check[prior], prior_status$ucell_check[matching]))

order_index <- order(match(candidates$figure_category, category_order), candidates$cytokine, candidates$direction)
candidates <- candidates[order_index, , drop = FALSE]
candidates$page <- seq_len(nrow(candidates))
candidates$figure_category <- as.character(candidates$figure_category)
sample_medians$figure_category <- candidates$figure_category[match(sample_medians$signature_id, candidates$signature_id)]
write_tsv(candidates, file.path(stage_dir, "candidate_ucell_summary.tsv"))
write_tsv(sample_medians, file.path(stage_dir, "candidate_sample_medians.tsv"))

plot_candidate <- function(info) {
  id <- info$signature_id[[1L]]; target <- info$cytokine[[1L]]; direction <- info$direction[[1L]]
  values <- score_values(id)
  conditions <- unique(metadata$condition)
  condition_order <- c("control", target, sort(setdiff(conditions, c("control", target))))
  local_samples <- sample_metadata[order(match(sample_metadata$condition, condition_order),
    suppressWarnings(as.numeric(sample_metadata$replicate)), sample_metadata$sample_id), , drop = FALSE]
  ymax <- max(0.02, max(values) * 1.06)
  old_par <- par(no.readonly = TRUE); on.exit(par(old_par), add = TRUE)
  par(mar = c(10.5, 4.8, 5.3, 1.3), mgp = c(2.8, 0.8, 0), xpd = FALSE)
  plot(NA_real_, NA_real_, xlim = c(0.35, length(condition_order) + 0.65), ylim = c(0, ymax),
    xaxt = "n", xlab = "", ylab = "Per-cell UCell score", las = 1, bty = "l",
    main = paste0(target, " ", toupper(direction), " — ", gsub("_", " / ", info$figure_category[[1L]])),
    sub = paste0(id, " | ", info$candidate_size[[1L]], " genes | rest AUC ",
      formatC(info$ucell_rest_auc[[1L]], digits = 3L, format = "f"), " | discovery: ",
      info$discovery_status[[1L]], if (info$membership_scope[[1L]] != "pipeline_selected")
      " (exploratory membership)" else "",
      " | target n=3, other-single n=75; control displayed but excluded from AUC"),
    cex.main = 1.05, cex.sub = 0.77)
  grid(nx = NA, ny = NULL, col = "#E8EDF0", lty = 1)

  sample_offsets <- c(-0.22, 0, 0.22)
  for (condition_index in seq_along(condition_order)) {
    condition <- condition_order[[condition_index]]
    condition_samples <- local_samples[local_samples$condition == condition, , drop = FALSE]
    stopifnot(nrow(condition_samples) == 3L)
    fill <- if (condition == "control") "#B8B8B8" else if (condition == target) {
      if (info$ucell_check[[1L]] == "PASS") "#D95F02" else "#E6AB02"
    } else "#80B1D3"
    edge <- if (condition == "control") "#4D4D4D" else if (condition == target) "#8C2D04" else "#2166AC"
    for (sample_index in seq_len(nrow(condition_samples))) {
      sample_id <- condition_samples$sample_id[[sample_index]]
      sample_values <- values[metadata$sample_id == sample_id]
      x <- condition_index + sample_offsets[[sample_index]]
      if (diff(range(sample_values)) > .Machine$double.eps^0.5) {
        den <- density(sample_values, n = 256L, bw = "nrd0", from = 0, to = ymax)
        keep <- which(is.finite(den$y) & den$y >= max(den$y, na.rm = TRUE) * 0.002)
        if (length(keep) >= 2L) {
          keep <- seq.int(min(keep), max(keep)); half_width <- 0.09 * den$y[keep] / max(den$y[keep])
          polygon(c(x - half_width, rev(x + half_width)), c(den$x[keep], rev(den$x[keep])),
            col = adjustcolor(fill, alpha.f = 0.62), border = adjustcolor(edge, alpha.f = 0.88), lwd = 0.65)
        }
      } else {
        segments(x - 0.07, sample_values[[1L]], x + 0.07, sample_values[[1L]], col = edge, lwd = 3)
      }
      points(x, median(sample_values), pch = 21, bg = "white", col = "#111111", cex = 0.62, lwd = 0.75)
    }
  }
  axis(1, at = seq_along(condition_order), labels = FALSE, tck = -0.012)
  usr <- par("usr"); label_y <- usr[[3L]] - 0.047 * diff(usr[3:4])
  label_col <- ifelse(condition_order == target, "#8C2D04",
    ifelse(condition_order == "control", "#222222", "#333333"))
  text(seq_along(condition_order), label_y, condition_order, srt = 55, adj = 1,
    cex = 0.70, xpd = NA, col = label_col)
  mtext("Condition (three biological samples per condition)", side = 1, line = 8.9, cex = 0.88)
  expected <- if (direction == "up") "higher target score expected" else "lower target score expected"
  legend("topright", inset = 0.005, bty = "n", cex = 0.76,
    legend = c("control", "target cytokine", "other cytokines", "sample median", expected),
    pch = c(15, 15, 15, 21, NA),
    pt.bg = c("#B8B8B8", if (info$ucell_check[[1L]] == "PASS") "#D95F02" else "#E6AB02",
      "#80B1D3", "white", NA),
    col = c("#4D4D4D", "#8C2D04", "#2166AC", "#111111", "#333333"),
    lty = c(NA, NA, NA, NA, 1), lwd = c(NA, NA, NA, NA, 1))
}

render_pdf <- function(path, rows) {
  pdf(path, width = 14, height = 8.5, onefile = TRUE, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  if (!length(rows)) {
    plot.new()
    text(0.5, 0.55, "No candidates in this category", cex = 1.4, font = 2)
  } else {
    for (i in rows) plot_candidate(candidates[i, , drop = FALSE])
  }
}

combined_pdf <- file.path(stage_dir, "all_51_candidate_ucell_violin_maps.pdf")
render_pdf(combined_pdf, seq_len(nrow(candidates)))
category_files <- setNames(character(length(category_order)), category_order)
for (category in category_order) {
  rows <- which(candidates$figure_category == category)
  category_files[[category]] <- paste0(tolower(category), "_violin_maps.pdf")
  render_pdf(file.path(stage_dir, category_files[[category]]), rows)
}

category_summary <- do.call(rbind, lapply(category_order, function(category) {
  rows <- candidates$figure_category == category
  data.frame(category = category, candidates = sum(rows),
    discovery_pass = sum(rows & candidates$discovery_status == "discovery_pass"),
    discovery_fail_exploratory = sum(rows & candidates$discovery_status != "discovery_pass"),
    pdf = category_files[[category]], stringsAsFactors = FALSE)
}))
write_tsv(category_summary, file.path(stage_dir, "category_summary.tsv"))

input_manifest <- data.frame(path = normalizePath(input_paths, mustWork = TRUE),
  sha256 = vapply(input_paths, wb_source_sha256, character(1)), stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))

report_lines <- c(
  "# All 51 one-versus-all candidate UCell violin maps", "",
  paste0("The 51 directions with at least five eligible genes are grouped by UCell confirmation and direction. ",
    "Each PDF page shows cell-level raw UCell distributions for control and all 26 individual cytokine conditions. ",
    "The three narrow violins within each condition are its biological samples; white points mark sample medians."), "",
  paste0("UCell confirmation uses the direction-oriented AUC of three target sample medians versus 75 other-single sample medians. ",
    "PASS requires AUC >= ", format(cfg$analysis$ucell_rest_auc_minimum), ". Control is displayed but excluded from this test."), "",
  paste0("The completed pipeline selected and scored 47 discovery-passing sets. For this requested complete figure, ",
    "the same prespecified size-selection rule recovered best-sized memberships for the four discovery failures; those four were ",
    "scored with native scWorkbench UCell in an isolated session and are labeled exploratory. Their UCell result does not override their discovery failure."),
  "", "## Outputs", "",
  "- `all_51_candidate_ucell_violin_maps.pdf`: all candidates ordered as pass/up, pass/down, fail/up, fail/down.",
  paste0("- `", category_files[["UCell_pass_up"]], "`: UCell pass, up direction."),
  paste0("- `", category_files[["UCell_pass_down"]], "`: UCell pass, down direction."),
  paste0("- `", category_files[["UCell_fail_up"]], "`: UCell fail, up direction."),
  paste0("- `", category_files[["UCell_fail_down"]], "`: UCell fail, down direction."),
  "- `candidate_ucell_summary.tsv`: candidate sizes, discovery status, UCell AUC, category, and page number.",
  "- `candidate_sample_medians.tsv`: exact sample medians used for confirmation.",
  "- `candidate_gene_sets.tsv`: exact membership of all 51 plotted candidates.", "",
  "No score, gene set, or figure-only artifact was promoted. The accepted Parse bundle and staged 37-set rest-only collection are unchanged."
)
writeLines(report_lines, file.path(stage_dir, "REPORT.md"))

stopifnot(identical(bundle_fingerprint_before, wb_manifest(bundle)$current_state$fingerprint),
  nrow(candidates) == 51L, nrow(sample_medians) == 51L * 81L,
  all(table(sample_medians$signature_id) == 81L),
  all(table(interaction(sample_medians$signature_id, sample_medians$condition, drop = TRUE)) == 3L),
  file.exists(combined_pdf), all(file.exists(file.path(stage_dir, category_files))))

completion <- list(
  status = "complete", created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  dataset_id = cfg$bundle_id, source_run = cfg$run_id, candidates = nrow(candidates),
  discovery_passing_reused = sum(candidates$discovery_status == "discovery_pass"),
  discovery_failures_scored_exploratorily = sum(candidates$discovery_status != "discovery_pass"),
  category_counts = setNames(as.list(category_summary$candidates), category_summary$category),
  cells_per_candidate = length(cells), conditions = length(unique(metadata$condition)),
  biological_samples = nrow(sample_metadata), result_promotions = 0L,
  input_bundle_state_unchanged = TRUE, bundle_fingerprint = bundle_fingerprint_before,
  accepted_parent = parent_check, runtime = runtime, figure_only_session = workspace
)
write_json(completion, file.path(stage_dir, "completion.json"))

output_files <- list.files(stage_dir, full.names = TRUE, recursive = TRUE)
output_files <- output_files[basename(output_files) != "output_manifest.tsv"]
output_manifest <- data.frame(path = substring(output_files, nchar(stage_dir) + 2L),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, wb_source_sha256, character(1)), stringsAsFactors = FALSE)
write_tsv(output_manifest, file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle); bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not atomically finalize output directory")
cat("Completed all-candidate UCell violin maps:", output_dir, "\n")
