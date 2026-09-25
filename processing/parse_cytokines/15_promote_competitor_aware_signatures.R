#!/usr/bin/env Rscript

# Promote the reviewed competitor-aware passing full sets and their three-fold
# cores into distinct read-only feature-set namespaces. Per-cell scores remain
# in durable sessions; failed candidates remain report-only.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 15_promote_competitor_aware_signatures.R <new_report_directory>")
}

report_dir <- normalizePath(args[[1L]], mustWork = FALSE)
if (file.exists(report_dir)) stop("Report path already exists: ", report_dir)
stage_dir <- paste0(report_dir, ".incomplete")
if (file.exists(stage_dir)) stop("Incomplete report path already exists: ", stage_dir)
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
stopifnot(identical(as.character(packageVersion("scworkbench")), runtime$scworkbench_version),
  isTRUE(cfg$publication$promote_results), isTRUE(cfg$publication$report_before_promotion))

full_namespace <- cfg$publication$full_feature_set_namespace
core_namespace <- cfg$publication$core_feature_set_namespace
review_dir <- cfg$publication$source_review
work_dir <- file.path(cfg$storage_root, "work", cfg$run_id,
  "competitor_aware_promotion_inputs")
session_root <- file.path(cfg$storage_root, "sessions", cfg$run_id)
session_id <- "competitor_aware_feature_set_promotion_20260925"
session_workspace <- file.path(session_root, session_id)
if (dir.exists(session_workspace)) stop("Promotion session already exists: ", session_workspace)
if (dir.exists(work_dir)) stop("Promotion input directory already exists: ", work_dir)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- c(
  signature_summary = file.path(review_dir, "signature_summary.tsv"),
  final_membership = file.path(review_dir, "final_gene_membership.tsv"),
  core_membership = file.path(review_dir, "core_gene_membership.tsv"),
  analysis_completion = file.path(review_dir, "completion.json"),
  analysis_methods = file.path(review_dir, "METHODS.md"),
  figure_completion = file.path(review_dir, "selected_ucell_violins", "completion.json"),
  project_config = "config/parse_cytokines_rest_only.yml",
  parent_config = "config/parse_cytokines.yml",
  bundle_registry = "config/bundles.yml",
  script = "processing/parse_cytokines/15_promote_competitor_aware_signatures.R",
  bundle_manifest_before = file.path(cfg$bundle_dir, "manifest.json")
)
missing <- input_paths[!file.exists(input_paths)]
if (length(missing)) stop("Missing inputs: ", paste(missing, collapse = ", "))

read_tsv <- function(path) read.delim(path, check.names = FALSE, na.strings = c("NA", ""))
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}
artifact_id_for <- function(promotion, run_id, output_name) {
  rows <- promotion$promoted[
    promotion$promoted$run_id == run_id & promotion$promoted$output == output_name &
      promotion$promoted$status == "ok", , drop = FALSE]
  if (nrow(rows) != 1L) stop("Expected one promoted artifact for ", run_id, "/", output_name)
  as.character(rows$artifact_id[[1L]])
}

# Capture input identities before the bundle manifest changes.
input_manifest <- data.frame(
  path = normalizePath(input_paths, mustWork = TRUE),
  sha256 = vapply(input_paths, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE
)

summary <- read_tsv(input_paths[["signature_summary"]])
final_members <- read_tsv(input_paths[["final_membership"]])
core_members <- read_tsv(input_paths[["core_membership"]])
analysis_completion <- jsonlite::read_json(input_paths[["analysis_completion"]],
  simplifyVector = TRUE)
figure_completion <- jsonlite::read_json(input_paths[["figure_completion"]],
  simplifyVector = TRUE)
passing <- summary[summary$full_overall_pass, , drop = FALSE]
passing_ids <- as.character(passing$signature_id)
stopifnot(nrow(summary) == 52L, nrow(passing) == 33L,
  identical(analysis_completion$full_passes, 33L),
  identical(figure_completion$full_pass, 33L),
  all(passing$full_status == "PASS"), all(passing$core_size > 0L),
  setequal(unique(core_members$signature_id), passing_ids))

full_rows <- final_members[final_members$signature_id %in% passing_ids, , drop = FALSE]
full_rows <- full_rows[order(match(full_rows$signature_id, passing_ids), full_rows$rank), ]
full_import <- data.frame(
  feature_set = full_rows$signature_id,
  set_name = paste0(full_rows$cytokine, " ", toupper(full_rows$direction),
    " competitor-aware signature"),
  feature_id = full_rows$feature_id,
  description = paste0(
    "Reviewed Parse target-versus-other-cytokines ", full_rows$direction,
    " signature for ", full_rows$cytokine,
    "; full-data mean-or-responder competitor rule passed; selected size ",
    full_rows$selected_size, "."),
  stringsAsFactors = FALSE
)

core_rows <- core_members[order(match(core_members$signature_id, passing_ids),
  core_members$core_rank), ]
core_id_map <- setNames(paste0(passing_ids, ".core"), passing_ids)
core_import <- data.frame(
  feature_set = unname(core_id_map[core_rows$signature_id]),
  set_name = paste0(core_rows$cytokine, " ", toupper(core_rows$direction),
    " competitor-aware stable core"),
  feature_id = core_rows$feature_id,
  description = paste0(
    "Stable core of reviewed Parse target-versus-other-cytokines ",
    core_rows$direction, " signature for ", core_rows$cytokine,
    "; gene selected in the full set and all three leave-one-replicate-out rankings."),
  stringsAsFactors = FALSE
)
stopifnot(length(unique(full_import$feature_set)) == 33L,
  length(unique(core_import$feature_set)) == 33L,
  nrow(full_import) == sum(passing$signature_size),
  nrow(core_import) == sum(passing$core_size),
  all(table(full_import$feature_set) == passing$signature_size[
    match(names(table(full_import$feature_set)), passing$signature_id)]),
  all(table(core_import$feature_set) == passing$core_size[
    match(sub("\\.core$", "", names(table(core_import$feature_set))), passing$signature_id)]))

full_csv <- file.path(work_dir, "passing_full_feature_sets.csv")
core_csv <- file.path(work_dir, "passing_core_feature_sets.csv")
write.csv(full_import, full_csv, row.names = FALSE)
write.csv(core_import, core_csv, row.names = FALSE)
stopifnot(identical(read.csv(full_csv, check.names = FALSE), full_import),
  identical(read.csv(core_csv, check.names = FALSE), core_import))
write_tsv(full_import, file.path(stage_dir, "promoted_full_feature_sets.tsv"))
write_tsv(core_import, file.path(stage_dir, "promoted_core_feature_sets.tsv"))
writeLines(vapply(unique(full_import$feature_set), function(id) {
  x <- full_import[full_import$feature_set == id, ]
  paste(c(id, x$description[[1L]], x$feature_id), collapse = "\t")
}, character(1)), file.path(stage_dir, "promoted_full_feature_sets.gmt"))
writeLines(vapply(unique(core_import$feature_set), function(id) {
  x <- core_import[core_import$feature_set == id, ]
  paste(c(id, x$description[[1L]], x$feature_id), collapse = "\t")
}, character(1)), file.path(stage_dir, "promoted_core_feature_sets.gmt"))

parent_check <- validate_registered_bundle(registry, "parse_cytokines",
  run_recipient_validation = FALSE)
stopifnot(isTRUE(jsonlite::read_json(registry$bundles$parse_cytokines$validation_report,
  simplifyVector = TRUE)$ok))
bundle <- wb_open_bundle(cfg$bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells_before <- as.character(wb_cells(bundle, "all"))
cell_sets_before <- sort(as.character(wb_cell_sets(bundle)))
stopifnot(identical(fingerprint_before, analysis_completion$source_bundle_fingerprint),
  length(cells_before) == 154209L,
  nrow(wb_feature_sets(bundle, matrix_id = cfg$matrix_id,
    namespace = full_namespace)) == 0L,
  nrow(wb_feature_sets(bundle, matrix_id = cfg$matrix_id,
    namespace = core_namespace)) == 0L)

log_event("Staging", length(unique(full_import$feature_set)), "full feature sets")
session <- wb_session_create(bundle, workspace_root = session_root, session_id = session_id)
full_run <- wb_session_run_op(session, bundle, "feature_set_import", inputs = list(
  matrix = cfg$matrix_id, cell_set = "all", min_cells = 10L,
  source_path = full_csv, source_format = "table", namespace = full_namespace,
  collection = "competitor_aware_full", species = "Mus musculus",
  organism = "Mus musculus", gene_id_type = "bundle_feature_id",
  set_id_col = "feature_set", set_name_col = "set_name", gene_col = "feature_id",
  description_col = "description", split_by_collection = FALSE,
  overwrite = FALSE, strict = TRUE,
  created_by = "scworkbench_Bdev competitor-aware workflow",
  source_id = "parse_cytokines_rest_only_competitor_aware",
  source_version = "2026-09-25", origin = "derived", is_readonly = TRUE
), notes = paste(
  "User-authorized promotion after full-data competitor analysis, holdout stability review,",
  "core-set review, and inspection of all 52 violin maps. Only 33 full-data passing sets are imported."
))
full_manifest <- full_run$run$outputs$import_manifest$table
stopifnot(nrow(full_manifest) == 33L,
  setequal(as.character(full_manifest$feature_set), unique(full_import$feature_set)),
  all(full_manifest$status == "imported"),
  all(full_manifest$n_genes_raw == full_manifest$n_genes_mapped_to_bundle))

log_event("Staging", length(unique(core_import$feature_set)), "core feature sets")
core_run <- wb_session_run_op(full_run$session, bundle, "feature_set_import", inputs = list(
  matrix = cfg$matrix_id, cell_set = "all", min_cells = 10L,
  source_path = core_csv, source_format = "table", namespace = core_namespace,
  collection = "competitor_aware_core", species = "Mus musculus",
  organism = "Mus musculus", gene_id_type = "bundle_feature_id",
  set_id_col = "feature_set", set_name_col = "set_name", gene_col = "feature_id",
  description_col = "description", split_by_collection = FALSE,
  overwrite = FALSE, strict = TRUE,
  created_by = "scworkbench_Bdev competitor-aware workflow",
  source_id = "parse_cytokines_rest_only_competitor_aware_core",
  source_version = "2026-09-25", origin = "derived", is_readonly = TRUE
), notes = paste(
  "Stable core sets for the 33 promoted full signatures. Each core gene belongs to the",
  "full-data membership and all three selected-size leave-one-replicate-out memberships."
))
core_manifest <- core_run$run$outputs$import_manifest$table
stopifnot(nrow(core_manifest) == 33L,
  setequal(as.character(core_manifest$feature_set), unique(core_import$feature_set)),
  all(core_manifest$status == "imported"),
  all(core_manifest$n_genes_raw == core_manifest$n_genes_mapped_to_bundle))

dag <- wb_session_dag(core_run$session, bundle = bundle,
  include_artifacts = TRUE, include_active = TRUE)
write_json(dag, file.path(stage_dir, "session_dag_before_promotion.json"))
write_tsv(full_manifest, file.path(stage_dir, "full_import_manifest.tsv"))
write_tsv(core_manifest, file.path(stage_dir, "core_import_manifest.tsv"))

log_event("Promoting both reviewed collections in one closure-aware receipt")
promotion <- wb_session_promote_current(bundle, core_run$session, outputs = list(
  wb_ref_run(full_run$run$run_id, "import_summary"),
  wb_ref_run(full_run$run$run_id, "import_manifest"),
  wb_ref_run(full_run$run$run_id, "rejected_sets"),
  wb_ref_run(core_run$run$run_id, "import_summary"),
  wb_ref_run(core_run$run$run_id, "import_manifest"),
  wb_ref_run(core_run$run$run_id, "rejected_sets")
), name_prefix = "parse_rest_only_competitor_aware",
  description = "Reviewed competitor-aware Parse cytokine signatures and stable three-fold cores.",
  created_by = "scworkbench_Bdev competitor-aware workflow", fail_fast = TRUE)

promotion_record <- list(
  promoted = promotion$promoted,
  snapshot_id = if (is.null(promotion$snapshot_id)) NA_character_ else promotion$snapshot_id,
  commit_id = if (is.null(promotion$commit_id)) NA_character_ else promotion$commit_id,
  durable_session_state = if (is.null(promotion$session$durable$state)) "promoted" else
    promotion$session$durable$state,
  session_workspace = session_workspace,
  full_run_id = full_run$run$run_id,
  core_run_id = core_run$run$run_id,
  full_manifest_artifact_id = artifact_id_for(promotion, full_run$run$run_id,
    "import_manifest"),
  core_manifest_artifact_id = artifact_id_for(promotion, core_run$run$run_id,
    "import_manifest")
)
saveRDS(promotion_record, file.path(stage_dir, "promotion_receipt.rds"))
write_json(promotion_record, file.path(stage_dir, "promotion_receipt.json"))

registered_full <- wb_feature_sets(bundle, matrix_id = cfg$matrix_id,
  namespace = full_namespace)
registered_core <- wb_feature_sets(bundle, matrix_id = cfg$matrix_id,
  namespace = core_namespace)
stopifnot(nrow(registered_full) == 33L, nrow(registered_core) == 33L,
  setequal(registered_full$feature_set, unique(full_import$feature_set)),
  setequal(registered_core$feature_set, unique(core_import$feature_set)),
  all(registered_full$is_readonly == 1L), all(registered_core$is_readonly == 1L),
  identical(as.character(wb_cells(bundle, "all")), cells_before),
  identical(sort(as.character(wb_cell_sets(bundle))), cell_sets_before))

membership_checks <- do.call(rbind, c(
  lapply(unique(full_import$feature_set), function(id) {
    expected <- full_import$feature_id[full_import$feature_set == id]
    observed <- as.character(wb_feature_set_get(bundle, id, cfg$matrix_id))
    data.frame(namespace = full_namespace, feature_set = id,
      expected_genes = length(expected), observed_genes = length(observed),
      ordered_membership_exact = identical(observed, expected), stringsAsFactors = FALSE)
  }),
  lapply(unique(core_import$feature_set), function(id) {
    expected <- core_import$feature_id[core_import$feature_set == id]
    observed <- as.character(wb_feature_set_get(bundle, id, cfg$matrix_id))
    data.frame(namespace = core_namespace, feature_set = id,
      expected_genes = length(expected), observed_genes = length(observed),
      ordered_membership_exact = identical(observed, expected), stringsAsFactors = FALSE)
  })
))
stopifnot(nrow(membership_checks) == 66L,
  all(membership_checks$ordered_membership_exact))
write_tsv(membership_checks, file.path(stage_dir, "feature_set_membership_checks.tsv"))

log_event("Finalizing recipient readiness")
readiness <- wb_finalize_recipient_bundle(bundle)
stopifnot(isTRUE(readiness$ready))
fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(!identical(fingerprint_before, fingerprint_after))
wb_close_bundle(bundle)
bundle <- NULL

log_event("Running full recipient validation")
validation <- wb_bundle_validation_report(cfg$bundle_dir, level = "recipient", progress = TRUE)
stopifnot(isTRUE(validation$ok),
  identical(validation$pre_validation_fingerprint, validation$post_validation_fingerprint))
write_json(validation, file.path(stage_dir, "validation.json"))

recipient <- wb_open_bundle(cfg$bundle_dir, read_only = TRUE, validate = FALSE)
recipient_manifest <- wb_manifest(recipient)
recipient_readiness <- wb_recipient_readiness(recipient)
stopifnot(identical(recipient_manifest$current_state$fingerprint, fingerprint_after),
  isTRUE(recipient_readiness$ready),
  nrow(wb_feature_sets(recipient, matrix_id = cfg$matrix_id,
    namespace = full_namespace)) == 33L,
  nrow(wb_feature_sets(recipient, matrix_id = cfg$matrix_id,
    namespace = core_namespace)) == 33L)
wb_close_bundle(recipient)
recipient <- NULL

write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))
write_json(recipient_manifest, file.path(stage_dir, "bundle_manifest_after.json"))
write_json(recipient_readiness, file.path(stage_dir, "recipient_readiness.json"))

completion <- list(
  status = "complete", created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  dataset_id = cfg$bundle_id, bundle_dir = cfg$bundle_dir, source_review = review_dir,
  full_namespace = full_namespace, core_namespace = core_namespace,
  promoted_full_sets = nrow(registered_full), promoted_core_sets = nrow(registered_core),
  promoted_full_genes = nrow(full_import), promoted_core_genes = nrow(core_import),
  failed_candidates_promoted = 0L, per_cell_scores_promoted = FALSE,
  bundle_fingerprint_before = fingerprint_before,
  bundle_state_fingerprint = fingerprint_after,
  recipient_validation_fingerprint = validation$post_validation_fingerprint,
  validation_ok = validation$ok, validation_status = validation$status,
  ordered_memberships_exact = all(membership_checks$ordered_membership_exact),
  cells_unchanged = TRUE, cell_sets_unchanged = TRUE,
  parent = parent_check, runtime = runtime, promotion = promotion_record
)
write_json(completion, file.path(stage_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

report_lines <- c(
  "# Competitor-aware signature promotion", "",
  paste0("Promotion and recipient validation completed for `", cfg$bundle_id, "`."), "",
  paste0("- Full passing signatures: **", nrow(registered_full), "** in `", full_namespace, "`."),
  paste0("- Stable core signatures: **", nrow(registered_core), "** in `", core_namespace, "`."),
  paste0("- Exact ordered membership checks: **", sum(membership_checks$ordered_membership_exact),
    "/", nrow(membership_checks), "**."),
  "- Failed candidates promoted: **0**.",
  "- Per-cell score families promoted: **no**; they remain in the completed durable sessions.", "",
  paste0("Bundle state fingerprint: `", fingerprint_after, "`"),
  paste0("Recipient validation fingerprint: `", validation$post_validation_fingerprint, "`"),
  paste0("Validation status: `", validation$status, "`"), "",
  "The original 154,209-cell universe and all cell-set memberships are unchanged. Full and core sets are read-only and use distinct namespaces and identifiers."
)
writeLines(report_lines, file.path(stage_dir, "REPORT.md"))

output_files <- list.files(stage_dir, recursive = TRUE, full.names = TRUE)
output_files <- output_files[basename(output_files) != "output_manifest.tsv"]
output_manifest <- data.frame(
  path = substring(output_files, nchar(stage_dir) + 2L),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(output_manifest[order(output_manifest$path), ],
  file.path(stage_dir, "output_manifest.tsv"))

if (!file.rename(stage_dir, report_dir)) stop("Failed to finalize promotion report")
log_event("Competitor-aware feature-set promotion complete at", report_dir)
