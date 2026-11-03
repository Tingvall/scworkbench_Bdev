#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) {
  stop("Usage: 03_score_and_label_fixture.R --config <config.yaml>")
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
results_root <- cfg$storage$results_root
clone_dir <- cfg$storage$execution_clone
session_root <- cfg$storage$sessions
dirs <- setNames(file.path(results_root, c(
  "01_tables", "02_figures", "03_reports", "04_logs", "06_manifests")),
  c("tables", "figures", "reports", "logs", "manifests"))
invisible(lapply(c(dirs, sessions = session_root), dir.create, recursive = TRUE,
  showWarnings = FALSE))
log_path <- file.path(dirs[["logs"]], "01_score_and_label.log")
sink(log_path, split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  out <- lapply(fields, function(x) unique(x[seq.int(3L, length(x))]))
  names(out) <- vapply(fields, `[[`, character(1), 1L)
  stopifnot(length(out) == 64L, !anyDuplicated(names(out)))
  out
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
  value <- rawConnectionValue(connection)
  digest::digest(value, algo = "sha256", serialize = FALSE)
}
run_ref <- function(run, output) wb_ref_run(run$run_id, output)
artifact_id_for <- function(promotion, run_id, output) {
  rows <- promotion$promoted[
    promotion$promoted$run_id == run_id &
      promotion$promoted$output == output &
      promotion$promoted$status == "ok", , drop = FALSE]
  stopifnot(nrow(rows) == 1L, nzchar(rows$artifact_id[[1L]]))
  as.character(rows$artifact_id[[1L]])
}
safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 10L || length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

stopifnot(dir.exists(clone_dir))
preflight <- jsonlite::read_json(file.path(dirs[["reports"]],
  "01_fixture_preflight.json"), simplifyVector = TRUE)
fixture <- read.delim(file.path(dirs[["tables"]], "01_fixture_cells.tsv"),
  check.names = FALSE)
stopifnot(identical(nrow(fixture), as.integer(cfg$fixture$expected_cells)),
  !anyDuplicated(fixture$cell_id))

log_event("Activating scwbGenePrograms")
wb_operation_pack_load("scwbGenePrograms")
pack_record <- wb_operation_packs()
stopifnot(nrow(pack_record) == 1L,
  identical(as.character(pack_record$manifest_sha256[[1L]]),
    input_cfg$runtime$scwb_gene_programs_manifest_sha256))

for (op in c("feature_set_import", "scwbgeneprograms::feature_set_overlap_qc",
    "scwbgeneprograms::ucell", "scwbgeneprograms::cell_vector_family_qc",
    "scwbgeneprograms::cell_vector_family_calibrate",
    "scwbgeneprograms::cell_vector_family_reduce")) {
  write_json(wb_op_interface(op), file.path(dirs[["manifests"]],
    paste0("interface_", gsub("[^A-Za-z0-9]+", "_", op), ".json")))
}

log_event("Opening writable execution clone; schema migration is clone-only")
bundle <- wb_open_bundle(clone_dir, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest_after_migration <- wb_manifest(bundle)
all_cells <- as.character(wb_cells(bundle, "all"))
stopifnot(length(all_cells) == input_cfg$bdev_reference$expected_cells,
  all(fixture$cell_id %in% all_cells))
full_pt <- wb_meta(bundle, cells = sort(all_cells, method = "radix"), cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$fixture$pseudotime_field))
stopifnot(sum(is.finite(full_pt[[cfg$fixture$pseudotime_field]])) ==
    input_cfg$bdev_reference$authoritative_pseudotime$defined_cells,
  sum(is.na(full_pt[[cfg$fixture$pseudotime_field]])) ==
    input_cfg$bdev_reference$authoritative_pseudotime$missing_cells)
observed_pt_sha <- pt_binary_sha256(full_pt$cell_id,
  full_pt[[cfg$fixture$pseudotime_field]])
stopifnot(identical(observed_pt_sha, preflight$source_pt_binary_sha256))
log_event("Exact authoritative pt values and missing entries verified after migration")

fixture_md <- wb_meta(bundle, cells = fixture$cell_id)
stopifnot(identical(as.character(fixture_md$cell_id), fixture$cell_id))
grouping <- wb_prim_cell_vector(stats::setNames(
  as.character(fixture_md[[cfg$fixture$cluster_field]]), fixture_md$cell_id),
  "categorical")

gmt_path <- input_cfg$immgen_collection$candidate_gmt
gene_sets <- read_gmt(gmt_path)
namespace <- cfg$scoring$namespace
import_workspace <- file.path(session_root, "import_64_population_collection")
registered_sets <- wb_feature_sets(bundle, matrix_id = cfg$scoring$matrix_id,
  namespace = namespace)
import_reused <- nrow(registered_sets) == 64L
if (import_reused) {
  log_event("Reusing the already promoted 64-population feature-set namespace")
  import_manifest <- read.delim(file.path(dirs[["tables"]],
    "03_feature_set_import_manifest.tsv"), check.names = FALSE)
  stopifnot(nrow(import_manifest) == 64L,
    sum(import_manifest$status == "imported") == 64L)
  import_result <- NULL
  import_promotion <- NULL
} else {
  stopifnot(nrow(registered_sets) == 0L, !dir.exists(import_workspace))
  import_session <- wb_session_create(bundle, workspace_root = session_root,
    session_id = "import_64_population_collection")
  log_event("Importing pinned 64-population GMT")
  import_result <- wb_session_run_op(import_session, bundle, "feature_set_import",
    inputs = list(
      matrix = cfg$scoring$matrix_id,
      cell_ids = fixture$cell_id,
      min_cells = 10L,
      source_path = gmt_path,
      source_format = "gmt",
      namespace = namespace,
      collection = "ImmGen GSE109125 vM25 all-population Parse candidates",
      species = "Mus musculus",
      gene_id_type = "gene_symbol",
      overwrite = TRUE,
      strict = FALSE,
      created_by = "scworkbench_Bdev bounded development acceptance",
      source_id = input_cfg$immgen_collection$run_identity,
      source_version = "20261102",
      origin = "diagnostic_candidate",
      is_readonly = TRUE),
    notes = "Pinned 64-population diagnostic ImmGen collection; development acceptance only.")
  import_manifest <- import_result$run$outputs$import_manifest$table
  stopifnot(nrow(import_manifest) == 64L,
    sum(import_manifest$status == "imported") == 64L)
  write_tsv(import_manifest, file.path(dirs[["tables"]],
    "03_feature_set_import_manifest.tsv"))
  import_promotion <- wb_session_promote_current(bundle, import_result$session,
    outputs = c("import_summary", "import_manifest", "rejected_sets"),
    name_prefix = "bdev_immgen64_import",
    description = "Diagnostic 64-population ImmGen collection for bounded Bdev labeling acceptance.",
    created_by = "scworkbench_Bdev development acceptance", fail_fast = TRUE)
}

analysis_workspace <- file.path(session_root, "score_and_reduce_fixture_v1")
stopifnot(!dir.exists(analysis_workspace))
session <- wb_session_create(bundle, workspace_root = session_root,
  session_id = "score_and_reduce_fixture_v1")
log_event("Running feature-set overlap QC")
overlap_result <- wb_session_run_op(session, bundle,
  "scwbgeneprograms::feature_set_overlap_qc", inputs = list(
    matrix = cfg$scoring$matrix_id,
    cell_ids = fixture$cell_id,
    min_cells = 10L,
    namespace = namespace,
    min_set_size = as.integer(cfg$scoring$feature_set_qc$minimum_mapped_genes),
    min_detected_fraction = as.numeric(
      cfg$scoring$feature_set_qc$minimum_detected_fraction),
    deduplicate_sets = isTRUE(cfg$scoring$feature_set_qc$deduplicate_sets)),
  notes = "Screen the pinned 64-population collection against Bdev RNA_counts.")
session <- overlap_result$session
overlap <- overlap_result$run$outputs$overlap_qc$table
write_tsv(overlap, file.path(dirs[["tables"]], "04_feature_set_overlap_qc.tsv"))
write_tsv(overlap_result$run$outputs$namespace_summary$table,
  file.path(dirs[["tables"]], "05_feature_set_overlap_summary.tsv"))
included <- as.character(overlap$feature_set[overlap$screening_status == "include"])
stopifnot(length(included) > 1L, all(included %in% names(gene_sets)))
included_sets <- gene_sets[included]
log_event("Scoring", length(included_sets), "of 64 sets after declared overlap QC")

ucell_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::ucell",
  inputs = list(
    matrix = cfg$scoring$matrix_id,
    cell_ids = fixture$cell_id,
    min_cells = 10L,
    features = included_sets,
    name = "IMMGEN_POP",
    max_rank = as.integer(cfg$scoring$ucell$max_rank),
    ncores = as.integer(cfg$scoring$ucell$ncores),
    parallel_backend = cfg$scoring$ucell$parallel_backend,
    chunk_size = as.integer(cfg$scoring$ucell$chunk_size),
    missing_features_policy = cfg$scoring$ucell$missing_features_policy),
  notes = "Raw UCell scores for included ImmGen population sets on the bounded Bdev fixture.")
session <- ucell_result$session

log_event("Running raw-family QC with cluster grouping")
qc_result <- wb_session_run_op(session, bundle,
  "scwbgeneprograms::cell_vector_family_qc", inputs = list(
    signal = run_ref(ucell_result$run, "signature_scores"),
    grouping = grouping,
    cell_ids = fixture$cell_id,
    min_cells = 10L,
    correlation_methods = unlist(cfg$scoring$raw_family_qc$correlation_methods),
    tie_tol = as.numeric(cfg$scoring$raw_family_qc$tie_tolerance)),
  notes = "Raw ImmGen score-family QC grouped by frozen Bdev cluster metadata.")
session <- qc_result$session

log_event("Percentile-calibrating each population score")
calibration_result <- wb_session_run_op(session, bundle,
  "scwbgeneprograms::cell_vector_family_calibrate", inputs = list(
    signal = run_ref(ucell_result$run, "signature_scores"),
    cell_ids = fixture$cell_id,
    min_cells = 10L,
    method = cfg$scoring$calibration$method,
    ties_method = cfg$scoring$calibration$ties_method,
    na_rm = TRUE),
  notes = "Endorsed percentile calibration before cross-key label reduction.")
session <- calibration_result$session

log_event("Reducing calibrated scores with conservative abstention")
reduction_result <- wb_session_run_op(session, bundle,
  "scwbgeneprograms::cell_vector_family_reduce", inputs = list(
    signal = run_ref(calibration_result$run, "calibrated_signal"),
    cell_ids = fixture$cell_id,
    min_cells = 10L,
    rule = cfg$scoring$reduction$rule,
    min_score = as.numeric(cfg$scoring$reduction$minimum_calibrated_score),
    min_margin = as.numeric(cfg$scoring$reduction$minimum_calibrated_margin),
    unknown_label = cfg$scoring$reduction$unknown_label,
    ties = cfg$scoring$reduction$ties,
    tie_tol = as.numeric(cfg$scoring$reduction$tie_tolerance)),
  notes = "Diagnostic provisional labels with explicit unassigned state; no winner cell-set registration.")
session <- reduction_result$session

promotion_outputs <- c(
  lapply(c("overlap_qc", "namespace_summary"), function(x)
    run_ref(overlap_result$run, x)),
  lapply(c("signature_scores", "signature_summary"), function(x)
    run_ref(ucell_result$run, x)),
  lapply(c("key_summary", "pairwise_summary", "group_summary"), function(x)
    run_ref(qc_result$run, x)),
  lapply(c("calibrated_signal", "calibration_summary"), function(x)
    run_ref(calibration_result$run, x)),
  lapply(c("label", "top_score", "runner_up_score", "score_margin",
    "abstain_reason", "decision_table", "label_summary"), function(x)
    run_ref(reduction_result$run, x)))
log_event("Promoting the complete diagnostic evidence chain into the execution clone")
promotion <- wb_session_promote_current(bundle, session,
  outputs = promotion_outputs,
  name_prefix = "bdev_immgen64_labeling_acceptance",
  description = paste("Bounded development acceptance for the diagnostic",
    "64-population ImmGen collection; raw scores, QC, calibration, labels,",
    "margins and abstention evidence; no accepted annotations."),
  created_by = "scworkbench_Bdev development acceptance", fail_fast = TRUE)

signature_summary <- ucell_result$run$outputs$signature_summary$table
key_summary <- qc_result$run$outputs$key_summary$table
pairwise_summary <- qc_result$run$outputs$pairwise_summary$table
group_summary <- qc_result$run$outputs$group_summary$table
calibration_summary <- calibration_result$run$outputs$calibration_summary$table
decision <- reduction_result$run$outputs$decision_table$table
label_summary <- reduction_result$run$outputs$label_summary$table
write_tsv(signature_summary, file.path(dirs[["tables"]], "06_ucell_signature_summary.tsv"))
write_tsv(key_summary, file.path(dirs[["tables"]], "07_raw_score_key_qc.tsv"))
write_tsv(pairwise_summary, file.path(dirs[["tables"]], "08_raw_score_pairwise_qc.tsv"))
write_tsv(group_summary, file.path(dirs[["tables"]], "09_raw_score_cluster_qc.tsv"))
write_tsv(calibration_summary, file.path(dirs[["tables"]], "10_calibration_summary.tsv"))
write_tsv(decision, file.path(dirs[["tables"]], "11_cell_label_decisions.tsv"))
write_tsv(label_summary, file.path(dirs[["tables"]], "12_label_summary.tsv"))

audit <- merge(decision, fixture_md, by = "cell_id", all.x = TRUE, sort = FALSE)
audit$pt_defined <- ifelse(is.finite(audit[[cfg$fixture$pseudotime_field]]),
  "defined", "missing")
target_summary <- read.delim(input_cfg$immgen_collection$target_summary,
  check.names = FALSE)
evidence <- target_summary[, c("geneset_id", "target", "source_role",
  "lineage", "evidence_tier", "status"), drop = FALSE]
audit <- merge(audit, evidence, by.x = "label", by.y = "geneset_id",
  all.x = TRUE, sort = FALSE)
audit$evidence_tier[is.na(audit$evidence_tier)] <- "unassigned"
write_tsv(audit, file.path(dirs[["tables"]], "13_label_decisions_with_metadata.tsv"))

count_table <- function(data, fields) {
  out <- as.data.frame(table(data[fields], useNA = "ifany"), stringsAsFactors = FALSE)
  names(out)[ncol(out)] <- "n_cells"
  out[out$n_cells > 0L, , drop = FALSE]
}
write_tsv(count_table(audit, c(cfg$fixture$cluster_field, "label")),
  file.path(dirs[["tables"]], "14_labels_by_cluster.tsv"))
write_tsv(count_table(audit, c("pt_defined", "label")),
  file.path(dirs[["tables"]], "15_labels_by_pt_status.tsv"))
write_tsv(count_table(audit, c("Phase", "label")),
  file.path(dirs[["tables"]], "16_labels_by_cell_cycle_phase.tsv"))
write_tsv(count_table(audit, c("orig.ident", "label")),
  file.path(dirs[["tables"]], "17_labels_by_source_identity.tsv"))
write_tsv(count_table(audit, c("evidence_tier", "label")),
  file.path(dirs[["tables"]], "18_labels_by_evidence_tier.tsv"))

score_primitive <- ucell_result$run$outputs$signature_scores
score_long <- data.frame(key = as.character(score_primitive$key),
  cell_id = as.character(score_primitive$cell_id),
  score = as.numeric(score_primitive$value), stringsAsFactors = FALSE)
score_md <- merge(score_long, fixture_md, by = "cell_id", all.x = TRUE,
  sort = FALSE)
confounders <- unlist(cfg$evaluation$confounder_fields)
correlation_rows <- do.call(rbind, lapply(split(score_md, score_md$key), function(x)
  do.call(rbind, lapply(confounders, function(field) data.frame(
    key = x$key[[1L]], field = field, n_finite = sum(is.finite(x$score) &
      is.finite(as.numeric(x[[field]]))), spearman_rho = safe_cor(x$score,
      as.numeric(x[[field]])), stringsAsFactors = FALSE)))))
write_tsv(correlation_rows, file.path(dirs[["tables"]],
  "19_score_confounder_correlations.tsv"))

label_counts <- sort(table(decision$label), decreasing = TRUE)
png(file.path(dirs[["figures"]], "01_label_counts.png"), width = 1800,
  height = 1100, res = 180)
par(mar = c(10, 5, 3, 1))
barplot(head(label_counts, 25L), las = 2, col = "#32688E",
  ylab = "Fixture cells", main = "Top provisional ImmGen labels (including abstention)")
dev.off()

abstain_counts <- sort(table(decision$abstain_reason), decreasing = TRUE)
png(file.path(dirs[["figures"]], "02_abstention_reasons.png"), width = 1400,
  height = 900, res = 180)
par(mar = c(6, 5, 3, 1))
barplot(abstain_counts, las = 2, col = "#D1495B", ylab = "Fixture cells",
  main = "Label reduction outcomes")
dev.off()

png(file.path(dirs[["figures"]], "03_score_margin_diagnostic.png"), width = 1400,
  height = 1000, res = 180)
cols <- ifelse(decision$label == cfg$scoring$reduction$unknown_label,
  grDevices::adjustcolor("grey50", alpha.f = 0.35),
  grDevices::adjustcolor("#2A9D8F", alpha.f = 0.45))
plot(decision$top_score, decision$score_margin, pch = 16, cex = 0.55,
  col = cols, xlab = "Top percentile-calibrated score", ylab = "Winner margin",
  main = "Conservative score and margin gates")
abline(v = cfg$scoring$reduction$minimum_calibrated_score,
  h = cfg$scoring$reduction$minimum_calibrated_margin, lty = 2, col = "#D1495B")
dev.off()

cluster_label <- table(audit[[cfg$fixture$cluster_field]], audit$label)
keep_labels <- names(head(sort(colSums(cluster_label), decreasing = TRUE), 20L))
heat <- cluster_label[, keep_labels, drop = FALSE]
heat <- heat / pmax(1, rowSums(heat))
png(file.path(dirs[["figures"]], "04_cluster_label_composition.png"), width = 2200,
  height = 1500, res = 180)
par(mar = c(12, 7, 4, 2))
image(seq_len(ncol(heat)), seq_len(nrow(heat)), t(heat), axes = FALSE,
  col = hcl.colors(64, "YlOrRd", rev = FALSE), xlab = "", ylab = "",
  main = "Cluster composition of the 20 most frequent provisional labels")
axis(1, at = seq_len(ncol(heat)), labels = colnames(heat), las = 2, cex.axis = 0.65)
axis(2, at = seq_len(nrow(heat)), labels = rownames(heat), las = 2, cex.axis = 0.8)
box()
dev.off()

wb_finalize_recipient_bundle(bundle)
manifest_final <- wb_manifest(bundle)
promotion_record <- list(
  import = if (import_reused) list(reused_existing_promotion = TRUE) else
    import_promotion$promoted,
  analysis = promotion$promoted,
  import_snapshot_id = if (import_reused) NULL else import_promotion$snapshot_id,
  analysis_snapshot_id = promotion$snapshot_id,
  import_manifest_artifact_id = if (import_reused) NULL else
    artifact_id_for(import_promotion, import_result$run$run_id, "import_manifest"),
  label_artifact_id = artifact_id_for(promotion,
    reduction_result$run$run_id, "label"),
  decision_table_artifact_id = artifact_id_for(promotion,
    reduction_result$run$run_id, "decision_table"))
write_json(promotion_record, file.path(dirs[["manifests"]],
  "03_promotion_receipt.json"))

assigned <- sum(decision$label != cfg$scoring$reduction$unknown_label)
pack_record <- wb_operation_packs()
summary <- list(
  schema_version = "bdev.immgen_bdev_labeling_acceptance_summary.v1",
  status = "technical_chain_complete_biological_review_required",
  run_identity = cfg$run_identity,
  source_bundle_mutated = FALSE,
  execution_clone = clone_dir,
  source_cells = length(all_cells),
  fixture_cells = nrow(fixture),
  exact_pt_binary_sha256 = observed_pt_sha,
  exact_pt_preserved = identical(observed_pt_sha,
    preflight$source_pt_binary_sha256),
  collection_sets = length(gene_sets),
  imported_sets = sum(import_manifest$status == "imported"),
  scored_sets = length(included_sets),
  excluded_sets = sum(overlap$screening_status == "exclude"),
  assigned_cells = assigned,
  unassigned_cells = nrow(decision) - assigned,
  assigned_fraction = assigned / nrow(decision),
  reduction = cfg$scoring$reduction,
  accepted_labels_emitted = FALSE,
  winner_cell_sets_registered = FALSE,
  implementation_stage_unlocked = FALSE,
  bundle_fingerprint_after_migration =
    manifest_after_migration$current_state$fingerprint,
  bundle_fingerprint_after_promotion = manifest_final$current_state$fingerprint,
  operation_pack = as.list(pack_record[1L, , drop = FALSE]))
write_json(summary, file.path(dirs[["reports"]], "02_labeling_summary.json"))

report_lines <- c(
  "# Bounded 64-population ImmGen labeling acceptance",
  "",
  paste0("- Status: **", summary$status, "**"),
  paste0("- Fixture: ", summary$fixture_cells, " of ", summary$source_cells,
    " Bdev cells"),
  paste0("- Collection: ", summary$imported_sets, " imported; ",
    summary$scored_sets, " passed overlap QC; ", summary$excluded_sets, " excluded"),
  paste0("- Provisional assignments: ", summary$assigned_cells, " (",
    sprintf("%.1f%%", 100 * summary$assigned_fraction), ")"),
  paste0("- Unassigned: ", summary$unassigned_cells),
  paste0("- Exact authoritative `pt` preservation: ", summary$exact_pt_preserved),
  "",
  "These are diagnostic development labels, not accepted Bdev annotations.",
  "Review overlap exclusions, evidence tiers, cluster localization, abstention",
  "reasons, and confounder correlations before changing thresholds or promoting",
  "the candidate collection.")
writeLines(report_lines, file.path(dirs[["reports"]],
  "03_labeling_acceptance_report.md"))
log_event("Completed bounded labeling acceptance; biological review remains required")
