#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config"))
  stop("Usage: 01_export_margin_decisions.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(digest)
})

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "06_manifests")), c("tables", "figures", "reports", "logs",
  "manifests"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_export_margin_decisions.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
for (entry in cfg$inputs) stopifnot(identical(sha256(entry$path), entry$sha256))
stopifnot(identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256))

scores <- readRDS(cfg$inputs$raw_signature_scores$path)
stopifnot(inherits(scores, "wb_primitive"), identical(scores$kind,
  "cell_vector_family"), length(scores$key) == length(scores$cell_id),
  length(scores$key) == length(scores$value), all(is.finite(scores$value)))
score_table <- data.frame(cell_id = scores$cell_id, program = scores$key,
  value = scores$value, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(score_table[c("cell_id", "program")]))
cell_groups <- split(score_table, score_table$cell_id)
expected_cells <- as.integer(cfg$acceptance_boundary$expected_cells)
expected_programs <- as.integer(cfg$acceptance_boundary$expected_programs)
stopifnot(length(cell_groups) == expected_cells,
  all(vapply(cell_groups, nrow, integer(1L)) == expected_programs))

tie_tolerance <- as.numeric(cfg$comparison$tie_tolerance)
winner_rows <- lapply(cell_groups, function(local) {
  ordering <- order(-local$value, local$program)
  top <- local$value[ordering[[1L]]]
  runner <- local$value[ordering[[2L]]]
  data.frame(
    cell_id = local$cell_id[[1L]],
    winning_program = local$program[ordering[[1L]]],
    top_score = top,
    runner_up_score = runner,
    score_margin = top - runner,
    top_tie = sum(abs(local$value - top) <= tie_tolerance) > 1L,
    stringsAsFactors = FALSE)
})
decisions <- do.call(rbind, winner_rows)
rownames(decisions) <- NULL

baseline <- read.delim(cfg$inputs$decisions_margin001$path,
  check.names = FALSE, quote = "\"")
coordinates <- read.delim(gzfile(cfg$inputs$coordinates$path),
  check.names = FALSE, quote = "\"")
stopifnot(nrow(baseline) == expected_cells, !anyDuplicated(baseline$cell_id),
  nrow(coordinates) == expected_cells, !anyDuplicated(coordinates$cell_id),
  setequal(decisions$cell_id, baseline$cell_id),
  setequal(decisions$cell_id, coordinates$cell_id))
decisions <- decisions[match(baseline$cell_id, decisions$cell_id), , drop = FALSE]
coordinates <- coordinates[match(baseline$cell_id, coordinates$cell_id),
  , drop = FALSE]
stopifnot(
  all(abs(decisions$top_score - baseline$top_score) < 1e-12),
  all(abs(decisions$runner_up_score - baseline$runner_up_score) < 1e-12),
  all(abs(decisions$score_margin - baseline$score_margin) < 1e-12))

unassigned <- cfg$comparison$unassigned_label
rules <- cfg$comparison$rules
for (rule in rules) {
  field <- paste0("label_", rule$rule_id)
  threshold <- as.numeric(rule$minimum_raw_margin)
  decisions[[field]] <- ifelse(!decisions$top_tie &
    decisions$score_margin >= threshold, decisions$winning_program,
    unassigned)
}
stopifnot(identical(decisions$label_margin001, as.character(baseline$label)))

sample_field <- cfg$comparison$sample_field
decisions[[sample_field]] <- baseline[[sample_field]]
decisions$pt <- baseline$pt
decisions$UMAP_1 <- coordinates$UMAP_1
decisions$UMAP_2 <- coordinates$UMAP_2

cell_path <- file.path(dirs[["tables"]], "01_cell_margin_comparison.tsv.gz")
connection <- gzfile(cell_path, open = "wt")
write.table(decisions, connection, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
close(connection)

threshold_rows <- lapply(rules, function(rule) {
  field <- paste0("label_", rule$rule_id)
  assigned <- decisions[[field]] != unassigned
  data.frame(
    rule_id = rule$rule_id,
    display_label = rule$display_label,
    minimum_raw_margin = as.numeric(rule$minimum_raw_margin),
    assigned = sum(assigned),
    unassigned = sum(!assigned),
    assigned_labels = length(unique(decisions[[field]][assigned])),
    stringsAsFactors = FALSE)
})
threshold_summary <- do.call(rbind, threshold_rows)
write.table(threshold_summary, file.path(dirs[["tables"]],
  "02_threshold_summary.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

programs <- sort(unique(decisions$winning_program))
label_rows <- lapply(programs, function(program) {
  values <- vapply(rules, function(rule) {
    field <- paste0("label_", rule$rule_id)
    sum(decisions[[field]] == program)
  }, integer(1L))
  data.frame(program = program, t(values), stringsAsFactors = FALSE)
})
label_counts <- do.call(rbind, label_rows)
names(label_counts)[-1L] <- vapply(rules, function(rule)
  paste0("assigned_", rule$rule_id), character(1L))
write.table(label_counts, file.path(dirs[["tables"]],
  "03_label_counts_by_rule.tsv"), sep = "\t", quote = TRUE,
  row.names = FALSE)

sample_ids <- sort(unique(as.character(decisions[[sample_field]])))
sample_rows <- lapply(sample_ids, function(sample_id) {
  local <- decisions[as.character(decisions[[sample_field]]) == sample_id,
    , drop = FALSE]
  values <- vapply(rules, function(rule) {
    field <- paste0("label_", rule$rule_id)
    sum(local[[field]] != unassigned)
  }, integer(1L))
  output <- data.frame(sample_id = sample_id, cells = nrow(local),
    t(values), stringsAsFactors = FALSE)
  names(output)[-(1:2)] <- vapply(rules, function(rule)
    paste0("assigned_", rule$rule_id), character(1L))
  output
})
sample_counts <- do.call(rbind, sample_rows)
write.table(sample_counts, file.path(dirs[["tables"]],
  "04_sample_counts_by_rule.tsv"), sep = "\t", quote = TRUE,
  row.names = FALSE)

expected_assigned <- unlist(cfg$acceptance_boundary$expected_assigned)
expected_labels <- unlist(cfg$acceptance_boundary$expected_assigned_labels)
stopifnot(
  identical(as.integer(threshold_summary$assigned),
    as.integer(expected_assigned[threshold_summary$rule_id])),
  identical(as.integer(threshold_summary$assigned_labels),
    as.integer(expected_labels[threshold_summary$rule_id])),
  sum(decisions$top_tie) == as.integer(cfg$acceptance_boundary$expected_ties))

receipt <- list(
  schema_version = "bdev.simple_mean_contrast_margin_export_receipt.v1",
  run_identity = cfg$run_identity,
  config_path = cfg_path,
  config_sha256 = sha256(cfg_path),
  source_score_run_id = cfg$scope$source_score_run_id,
  raw_signature_scores_sha256 = sha256(cfg$inputs$raw_signature_scores$path),
  cells = nrow(decisions),
  programs = expected_programs,
  exact_top_ties = sum(decisions$top_tie),
  baseline_decisions_reproduced = TRUE)
writeLines(jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE),
  file.path(dirs[["manifests"]], "01_export_receipt.json"))
print(threshold_summary)
