#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) {
  stop("Usage: 03_score_hierarchical_fixture.R --config <config.yaml>")
}

suppressPackageStartupMessages({
  library(scworkbench)
  library(yaml)
  library(jsonlite)
  library(digest)
})

config_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(config_path)
input_cfg <- yaml::read_yaml(cfg$inputs_manifest)
root <- cfg$storage$results_root
bundle_dir <- cfg$storage$session_bundle
session_root <- cfg$storage$sessions
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports", "04_logs", "06_manifests")),
  c("tables", "figures", "reports", "logs", "manifests"))
invisible(lapply(c(dirs, sessions = session_root), dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_hierarchical_scoring.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
write_tsv <- function(x, path) {
  connection <- if (grepl("[.]gz$", path)) gzfile(path, "wt") else path
  on.exit(if (inherits(connection, "connection")) close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  output <- lapply(fields, function(x) unique(x[seq.int(3L, length(x))]))
  names(output) <- vapply(fields, `[[`, character(1), 1L)
  stopifnot(!anyDuplicated(names(output)))
  output
}
run_ref <- function(run, output) wb_ref_run(run$run_id, output)
primitive_long <- function(primitive, variant, level, value_name) {
  output <- data.frame(
    variant = variant,
    level = level,
    cell_id = as.character(primitive$cell_id),
    key = as.character(primitive$key),
    value = as.numeric(primitive$value),
    stringsAsFactors = FALSE)
  names(output)[names(output) == "value"] <- value_name
  output
}
pt_binary_sha256 <- function(cell_id, pt) {
  ordering <- order(cell_id, method = "radix")
  cell_id <- as.character(cell_id[ordering]); pt <- as.numeric(pt[ordering])
  connection <- rawConnection(raw(), open = "wb")
  on.exit(close(connection))
  for (i in seq_along(cell_id)) {
    encoded <- charToRaw(enc2utf8(cell_id[[i]]))
    writeBin(as.integer(length(encoded)), connection, size = 4L, endian = "big")
    writeBin(encoded, connection)
    if (is.na(pt[[i]])) {
      writeBin(as.raw(0L), connection)
    } else {
      writeBin(as.raw(1L), connection)
      writeBin(pt[[i]], connection, size = 8L, endian = "big")
    }
  }
  digest::digest(rawConnectionValue(connection), algo = "sha256", serialize = FALSE)
}

fixture <- read.delim(file.path(dirs[["tables"]], "01_fixture_cells.tsv"), check.names = FALSE)
stopifnot(nrow(fixture) == as.integer(cfg$fixture$expected_cells), !anyDuplicated(fixture$cell_id))
broad_map <- read.delim(file.path(dirs[["tables"]], "03_broad_lineage_set_map.tsv"), check.names = FALSE)
population_map <- read.delim(input_cfg$population_collection$target_summary, check.names = FALSE)
lineage_by_population <- setNames(population_map$lineage, population_map$geneset_id)
lineage_by_broad <- setNames(broad_map$lineage, broad_map$broad_set_id)

broad_gmt <- file.path(dirs[["tables"]], "02_broad_lineage_genesets.gmt")
filtered_gmt <- file.path(dirs[["tables"]], "04_shared_gene_filtered_population_genesets.gmt")
original_gmt <- input_cfg$population_collection$source_gmt
sets <- list(broad = read_gmt(broad_gmt), original = read_gmt(original_gmt),
  shared_gene_filtered = read_gmt(filtered_gmt))
stopifnot(length(sets$broad) == 9L, length(sets$original) == 64L,
  length(sets$shared_gene_filtered) == 64L,
  identical(sort(names(sets$original)), sort(names(sets$shared_gene_filtered))))

log_event("Loading scwbGenePrograms and opening the isolated session bundle")
wb_operation_pack_load("scwbGenePrograms")
pack_record <- wb_operation_packs()
stopifnot(nrow(pack_record) == 1L,
  identical(as.character(pack_record$manifest_sha256[[1L]]),
    input_cfg$runtime$scwb_gene_programs_manifest_sha256))
bundle <- wb_open_bundle(bundle_dir, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_cells <- as.character(wb_cells(bundle, "all"))
stopifnot(length(all_cells) == input_cfg$bdev_reference$expected_cells)
full_pt <- wb_meta(bundle, cells = sort(all_cells, method = "radix"),
  cols = list(wb_ref_metadata_field("core:cell_id"), cfg$fixture$pseudotime_field))
observed_pt_sha <- pt_binary_sha256(full_pt$cell_id, full_pt[[cfg$fixture$pseudotime_field]])
stopifnot(identical(observed_pt_sha,
  input_cfg$bdev_reference$authoritative_pseudotime$binary_sha256),
  sum(is.finite(full_pt[[cfg$fixture$pseudotime_field]])) ==
    input_cfg$bdev_reference$authoritative_pseudotime$defined_cells,
  sum(is.na(full_pt[[cfg$fixture$pseudotime_field]])) ==
    input_cfg$bdev_reference$authoritative_pseudotime$missing_cells)
fixture_md <- wb_meta(bundle, cells = fixture$cell_id)
stopifnot(identical(as.character(fixture_md$cell_id), fixture$cell_id))

analysis_session <- wb_session_create(bundle, workspace_root = session_root,
  session_id = "hierarchical_scoring_fixture_v1")
session <- analysis_session
score_results <- list(); calibration_results <- list(); score_audits <- list()
for (variant in names(sets)) {
  log_event("Session-local UCell scoring without import or promotion:", variant)
  score_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::ucell",
    inputs = list(
      matrix = cfg$scoring$matrix_id,
      cell_ids = fixture$cell_id,
      min_cells = 10L,
      features = sets[[variant]],
      name = paste0("IMMGEN_HIER_", toupper(variant)),
      max_rank = as.integer(cfg$scoring$ucell$max_rank),
      ncores = as.integer(cfg$scoring$ucell$ncores),
      parallel_backend = cfg$scoring$ucell$parallel_backend,
      chunk_size = as.integer(cfg$scoring$ucell$chunk_size),
      missing_features_policy = cfg$scoring$ucell$missing_features_policy),
    notes = paste("Native UCell scores for hierarchical variant:", variant))
  session <- score_result$session
  calibration_result <- wb_session_run_op(session, bundle,
    "scwbgeneprograms::cell_vector_family_calibrate", inputs = list(
      signal = run_ref(score_result$run, "signature_scores"),
      cell_ids = fixture$cell_id,
      min_cells = 10L,
      method = cfg$scoring$calibration$method,
      ties_method = cfg$scoring$calibration$ties_method,
      na_rm = TRUE),
    notes = paste("Percentile calibration for hierarchical variant:", variant))
  session <- calibration_result$session
  score_results[[variant]] <- score_result
  calibration_results[[variant]] <- calibration_result
  score_audit <- score_result$run$outputs$signature_summary$table
  score_audit$detected_fraction <- score_audit$n_features_matched / score_audit$n_features_input
  score_audit$screening_status <- ifelse(
    score_audit$n_features_matched >= as.integer(cfg$scoring$feature_set_qc$minimum_mapped_genes) &
      score_audit$detected_fraction >= as.numeric(cfg$scoring$feature_set_qc$minimum_detected_fraction),
    "include", "exclude")
  stopifnot(all(score_audit$screening_status == "include"))
  score_audits[[variant]] <- score_audit
  write_tsv(score_audit, file.path(dirs[["tables"]], paste0("07_score_overlap_audit_", variant, ".tsv")))
}
log_event("Native scoring complete; retaining ordinary exploratory tables without promotion")

raw_scores <- do.call(rbind, lapply(names(sets), function(variant)
  primitive_long(score_results[[variant]]$run$outputs$signature_scores,
    variant, if (variant == "broad") "broad_state" else "population", "raw_score")))
calibrated_scores <- do.call(rbind, lapply(names(sets), function(variant)
  primitive_long(calibration_results[[variant]]$run$outputs$calibrated_signal,
    variant, if (variant == "broad") "broad_state" else "population", "calibrated_score")))
write_tsv(raw_scores, file.path(dirs[["tables"]], "10_raw_scores.tsv.gz"))
write_tsv(calibrated_scores, file.path(dirs[["tables"]], "11_calibrated_scores.tsv.gz"))

broad_long <- calibrated_scores[calibrated_scores$variant == "broad", , drop = FALSE]
broad_split <- split(broad_long, broad_long$cell_id)
broad_min <- as.numeric(cfg$hierarchy$broad_state$minimum_calibrated_score)
broad_delta <- as.numeric(cfg$hierarchy$broad_state$maximum_delta_from_top)
broad_decisions <- do.call(rbind, lapply(fixture$cell_id, function(cell) {
  rows <- broad_split[[cell]]
  rows <- rows[order(-rows$calibrated_score, rows$key), , drop = FALSE]
  top <- rows$calibrated_score[[1L]]
  eligible <- rows$calibrated_score >= broad_min & (top - rows$calibrated_score) <= broad_delta + 1e-12
  states <- sort(unname(lineage_by_broad[rows$key[eligible]]))
  data.frame(cell_id = cell,
    broad_top_state = unname(lineage_by_broad[rows$key[[1L]]]),
    broad_top_score = top,
    broad_runner_up_score = rows$calibrated_score[[2L]],
    broad_margin = top - rows$calibrated_score[[2L]],
    eligible_broad_states = if (length(states)) paste(states, collapse = ";") else cfg$hierarchy$broad_state$unknown_label,
    eligible_broad_state_count = length(states),
    broad_abstain_reason = if (length(states)) "none" else "min_score",
    stringsAsFactors = FALSE)
}))
write_tsv(broad_decisions, file.path(dirs[["tables"]], "12_broad_state_decisions.tsv"))

hierarchical_reduce <- function(variant) {
  score_long <- calibrated_scores[calibrated_scores$variant == variant, , drop = FALSE]
  score_split <- split(score_long, score_long$cell_id)
  min_score <- as.numeric(cfg$hierarchy$population$minimum_calibrated_score)
  min_margin <- as.numeric(cfg$hierarchy$population$minimum_calibrated_margin)
  tie_tol <- as.numeric(cfg$hierarchy$population$tie_tolerance)
  output <- do.call(rbind, lapply(seq_len(nrow(broad_decisions)), function(index) {
    broad_row <- broad_decisions[index, ]
    states <- if (broad_row$eligible_broad_state_count == 0L) character() else
      strsplit(broad_row$eligible_broad_states, ";", fixed = TRUE)[[1L]]
    rows <- score_split[[broad_row$cell_id]]
    rows <- rows[lineage_by_population[rows$key] %in% states, , drop = FALSE]
    if (!nrow(rows)) return(data.frame(
      variant = variant, cell_id = broad_row$cell_id, label = cfg$hierarchy$population$unknown_label,
      label_lineage = NA_character_, eligible_population_count = 0L,
      top_score = NA_real_, runner_up_score = NA_real_, score_margin = NA_real_,
      abstain_reason = "broad_unassigned", stringsAsFactors = FALSE))
    rows <- rows[order(-rows$calibrated_score, rows$key), , drop = FALSE]
    top <- rows$calibrated_score[[1L]]
    runner <- if (nrow(rows) > 1L) rows$calibrated_score[[2L]] else NA_real_
    margin <- if (is.finite(runner)) top - runner else Inf
    reason <- if (top < min_score) "min_score" else if (is.finite(margin) && margin <= tie_tol) "tie" else if (margin < min_margin) "min_margin" else "none"
    label <- if (reason == "none") rows$key[[1L]] else cfg$hierarchy$population$unknown_label
    data.frame(variant = variant, cell_id = broad_row$cell_id, label = label,
      label_lineage = if (reason == "none") unname(lineage_by_population[label]) else NA_character_,
      eligible_population_count = nrow(rows), top_score = top,
      runner_up_score = runner, score_margin = margin, abstain_reason = reason,
      stringsAsFactors = FALSE)
  }))
  output
}

fine_decisions <- rbind(hierarchical_reduce("original"), hierarchical_reduce("shared_gene_filtered"))
fine_decisions <- merge(fine_decisions, broad_decisions, by = "cell_id", all.x = TRUE, sort = FALSE)
fine_decisions <- merge(fine_decisions,
  fixture_md[, c("cell_id", cfg$fixture$cluster_field, cfg$fixture$pseudotime_field, "Phase", "orig.ident"), drop = FALSE],
  by = "cell_id", all.x = TRUE, sort = FALSE)
fine_decisions$pt_defined <- ifelse(is.finite(fine_decisions[[cfg$fixture$pseudotime_field]]), "defined", "missing")
write_tsv(fine_decisions, file.path(dirs[["tables"]], "13_hierarchical_population_decisions.tsv"))

count_rows <- function(data, fields) {
  output <- as.data.frame(table(data[fields], useNA = "ifany"), stringsAsFactors = FALSE)
  names(output)[ncol(output)] <- "n_cells"
  output[output$n_cells > 0L, , drop = FALSE]
}
write_tsv(count_rows(broad_decisions, c("eligible_broad_states", "eligible_broad_state_count")),
  file.path(dirs[["tables"]], "14_broad_state_summary.tsv"))
write_tsv(count_rows(fine_decisions, c("variant", "label", "abstain_reason")),
  file.path(dirs[["tables"]], "15_population_label_summary.tsv"))
write_tsv(count_rows(fine_decisions, c("variant", cfg$fixture$cluster_field, "label")),
  file.path(dirs[["tables"]], "16_population_labels_by_cluster.tsv"))
write_tsv(count_rows(fine_decisions, c("variant", "pt_defined", "label")),
  file.path(dirs[["tables"]], "17_population_labels_by_pt_status.tsv"))

original_decision <- fine_decisions[fine_decisions$variant == "original", c("cell_id", "label", "abstain_reason")]
filtered_decision <- fine_decisions[fine_decisions$variant == "shared_gene_filtered", c("cell_id", "label", "abstain_reason")]
comparison <- merge(original_decision, filtered_decision, by = "cell_id", suffixes = c("_original", "_filtered"))
comparison$label_agreement <- comparison$label_original == comparison$label_filtered
comparison$assignment_change <- ifelse(comparison$label_original == "unassigned" & comparison$label_filtered != "unassigned", "gained",
  ifelse(comparison$label_original != "unassigned" & comparison$label_filtered == "unassigned", "lost",
    ifelse(comparison$label_original != comparison$label_filtered, "changed", "unchanged")))
write_tsv(comparison, file.path(dirs[["tables"]], "18_variant_cell_comparison.tsv"))
write_tsv(count_rows(comparison, c("assignment_change", "label_agreement")),
  file.path(dirs[["tables"]], "19_variant_comparison_summary.tsv"))

broad_counts <- sort(table(broad_decisions$eligible_broad_states), decreasing = TRUE)
png(file.path(dirs[["figures"]], "01_eligible_broad_states.png"), width = 1900, height = 1100, res = 180)
par(mar = c(11, 5, 3, 1)); barplot(broad_counts, las = 2, col = "#3A6EA5", ylab = "Fixture cells", main = "Eligible broad-state combinations")
dev.off()
label_counts <- table(fine_decisions$variant, fine_decisions$label)
keep <- names(head(sort(colSums(label_counts), decreasing = TRUE), 25L))
png(file.path(dirs[["figures"]], "02_hierarchical_label_counts.png"), width = 2100, height = 1200, res = 180)
par(mar = c(12, 5, 3, 1)); barplot(label_counts[, keep, drop = FALSE], beside = TRUE, las = 2,
  col = c("#386CB0", "#F0027F"), ylab = "Fixture cells", main = "Original versus shared-gene-filtered hierarchical calls")
legend("topright", legend = rownames(label_counts), fill = c("#386CB0", "#F0027F"), bty = "n")
dev.off()

manifest_final <- wb_manifest(bundle)
stopifnot(identical(manifest_final$current_state$fingerprint,
  input_cfg$parent_execution$bundle_fingerprint))
archive_manifest <- file.path(input_cfg$bdev_reference$archived_bundle_path, "manifest.json")
stopifnot(digest::digest(file = archive_manifest, algo = "sha256", serialize = FALSE) ==
  input_cfg$bdev_reference$archived_manifest_sha256)
assigned_by_variant <- tapply(fine_decisions$label != cfg$hierarchy$population$unknown_label,
  fine_decisions$variant, sum)
summary <- list(
  schema_version = "bdev.immgen_bdev_hierarchical_labeling_summary.v1",
  status = "diagnostic_hierarchical_comparison_complete",
  run_identity = cfg$run_identity,
  fixture_cells = nrow(fixture),
  exact_pt_preserved = identical(observed_pt_sha, input_cfg$bdev_reference$authoritative_pseudotime$binary_sha256),
  broad_sets_scored = length(sets$broad),
  population_sets_scored = length(sets$original),
  shared_genes_removed = nrow(read.delim(file.path(dirs[["tables"]], "06_removed_shared_genes.tsv"))),
  broad_single_state_cells = sum(broad_decisions$eligible_broad_state_count == 1L),
  broad_multi_state_cells = sum(broad_decisions$eligible_broad_state_count > 1L),
  broad_unassigned_cells = sum(broad_decisions$eligible_broad_state_count == 0L),
  assigned_original = unname(assigned_by_variant[["original"]]),
  assigned_shared_gene_filtered = unname(assigned_by_variant[["shared_gene_filtered"]]),
  label_agreement_cells = sum(comparison$label_agreement),
  accepted_labels_emitted = FALSE,
  winner_cell_sets_registered = FALSE,
  implementation_stage_unlocked = FALSE,
  parent_bundle_fingerprint = input_cfg$parent_execution$bundle_fingerprint,
  session_bundle_fingerprint = manifest_final$current_state$fingerprint,
  session_bundle_registry_unchanged = TRUE,
  scworkbench_artifacts_promoted = FALSE,
  operation_pack = as.list(pack_record[1L, , drop = FALSE]))
write_json(summary, file.path(dirs[["reports"]], "02_hierarchical_summary.json"))
report <- c(
  "# Bounded hierarchical ImmGen labeling qualification",
  "",
  paste0("- Status: **", summary$status, "**"),
  paste0("- Fixture cells: ", summary$fixture_cells),
  paste0("- Broad state: ", summary$broad_single_state_cells, " single-state; ", summary$broad_multi_state_cells, " multi-state; ", summary$broad_unassigned_cells, " unassigned"),
  paste0("- Original population panels assigned: ", summary$assigned_original),
  paste0("- Shared-gene-filtered population panels assigned: ", summary$assigned_shared_gene_filtered),
  paste0("- Cross-lineage shared genes removed: ", summary$shared_genes_removed),
  paste0("- Variant label agreement: ", summary$label_agreement_cells, "/", summary$fixture_cells),
  "",
  "The broad gate uses globally contrasted lineage panels and retains multiple states within the declared score window. Population competition is restricted to those eligible lineages. Both variants remain diagnostic; no accepted annotations or winner cell sets were emitted.")
writeLines(report, file.path(dirs[["reports"]], "03_hierarchical_qualification_report.md"))
write_json(list(
  schema_version = "bdev.immgen_hierarchy_session_receipt.v1",
  status = "complete_without_promotion",
  session_id = "hierarchical_scoring_fixture_v1",
  score_run_ids = lapply(score_results, function(x) x$run$run_id),
  calibration_run_ids = lapply(calibration_results, function(x) x$run$run_id),
  scworkbench_artifacts_promoted = FALSE),
  file.path(dirs[["manifests"]], "03_session_receipt.json"))
log_event("Completed bounded hierarchical scoring comparison")
