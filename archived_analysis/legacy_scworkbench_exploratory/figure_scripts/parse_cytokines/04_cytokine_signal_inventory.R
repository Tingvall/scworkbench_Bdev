#!/usr/bin/env Rscript
# Read-only inventory of saved Parse discovery evidence and native score QC.
# No DE fits, cell scoring, new hypotheses, or bundle promotions are performed.
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(scworkbench); library(jsonlite); library(yaml)})
source("src/bundle_config.R")
source("src/bundle_validation.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- normalizePath(args[[1]], mustWork = TRUE)
stopifnot(length(list.files(out, all.files = TRUE, no.. = TRUE)) == 0L)
registry <- load_bundle_registry()
cfg <- read_project_yaml("config/parse_cytokines.yml")
sig_cfg <- cfg$signature_analysis
discovery <- sig_cfg$discovery$report
confirmation <- sig_cfg$final$report
inputs <- character()
read_source <- function(path) {
  inputs <<- unique(c(inputs, normalizePath(path, mustWork = TRUE)))
  read.delim(path, check.names = FALSE, na.strings = c("NA", ""))
}
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = "\t",
  row.names = FALSE, quote = TRUE, na = "NA")
write_json <- function(x, name) jsonlite::write_json(x, file.path(out, name),
  pretty = TRUE, auto_unbox = TRUE, digits = NA, na = "null", null = "null")
log_event <- function(...) {cat(format(Sys.time(), tz = "UTC"), ..., "\n"); flush.console()}
near <- function(x, y) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-10))

main <- function() {
  log_event("Checking accepted Parse bundle and source tables")
  checked <- validate_registered_bundle(registry, "parse_cytokines", run_recipient_validation = FALSE)
  b <- wb_open_bundle(checked$bundle_dir, read_only = TRUE, validate = FALSE)
  on.exit(wb_close_bundle(b), add = TRUE)
  before <- wb_manifest(b)$current_state$fingerprint
  validation_path <- registry$bundles$parse_cytokines$validation_report
  stopifnot(isTRUE(read_json(validation_path, simplifyVector = TRUE)$ok))
  status <- read_source(file.path(confirmation, "all_signature_status.tsv"))
  reasons <- read_source(file.path(confirmation, "signature_exclusion_reasons.tsv"))
  members <- read_source(file.path(discovery, "selected_signature_members.tsv"))
  eligible <- read_source(file.path(discovery, "eligible_signature_candidates.tsv.gz"))
  final_sets <- read_source(sig_cfg$final$feature_set_tsv)
  scores <- read_source(file.path(confirmation, "ucell_sample_medians.tsv"))
  profile <- read_source(file.path(confirmation, "ucell_condition_specificity_profile.tsv"))
  expression_profile <- read_source(file.path(discovery, "selected_signature_offtarget_profile.tsv"))
  quantiles <- read_source(file.path(confirmation, "ucell_cell_distribution_quantiles.tsv"))
  conditions <- read_source(cfg$treatment_classification$condition_table)
  sample_metadata <- read_source(file.path(discovery, "sample_metadata.tsv"))
  parameters_path <- file.path(discovery, "analysis_parameters.json")
  parameters <- read_json(parameters_path, simplifyVector = TRUE)
  single <- conditions$condition[conditions$treatment_category == "single"]
  accepted <- status$signature_id[status$final_status == "accepted"]
  scored <- status$signature_id[status$selected_genes > 0L]
  stopifnot(nrow(status) == 52L, length(single) == 26L, length(accepted) == 23L,
    length(scored) == 31L, !anyDuplicated(status$signature_id),
    !anyDuplicated(status[c("cytokine", "direction")]),
    all(table(status$cytokine) == 2L), setequal(status$cytokine, single),
    all(status$direction %in% c("up", "down")),
    setequal(reasons$signature_id, status$signature_id),
    setequal(unique(members$signature_id), scored),
    !anyDuplicated(members[c("signature_id", "feature_id")]),
    all(members$selected), all(eligible$fold_direction_stable),
    setequal(unique(scores$signature_id), scored),
    !anyDuplicated(scores[c("signature_id", "sample_id")]),
    all(table(scores$signature_id) == 81L),
    !anyDuplicated(profile[c("signature_id", "evaluated_condition")]),
    nrow(profile) == 31L * 26L, nrow(expression_profile) == 31L * 26L)
  registered <- wb_feature_sets(b, matrix_id = "RNA_counts", namespace = sig_cfg$final$bundle_namespace)
  stopifnot(setequal(registered$feature_set, accepted), setequal(final_sets$feature_set, accepted))
  for (id in accepted) {
    genes <- final_sets$feature_id[final_sets$feature_set == id]
    stopifnot(identical(as.character(wb_feature_set_get(b, id, "RNA_counts")), genes),
      identical(genes, members$feature_id[members$signature_id == id]))
  }

  cells <- wb_cells(b, "set:parse_qc_control_single")
  md <- wb_meta(b, cells = cells, cols = list(wb_ref_metadata_field("core:cell_id"),
    "sample_id", "condition", "replicate", "treatment_category"))
  stopifnot(length(cells) == 154209L, !anyDuplicated(cells), identical(md$cell_id, cells),
    !anyNA(md), setequal(md$sample_id, sample_metadata$sample_id))
  observed_samples <- unique(md[, names(sample_metadata)])
  observed_samples <- observed_samples[match(sample_metadata$sample_id, observed_samples$sample_id), ]
  stopifnot(nrow(observed_samples) == 81L)
  for (field in names(sample_metadata)) stopifnot(identical(
    as.character(observed_samples[[field]]), as.character(sample_metadata[[field]])))
  score_md <- sample_metadata[match(scores$sample_id, sample_metadata$sample_id), ]
  for (field in c("condition", "replicate", "treatment_category")) stopifnot(
    identical(as.character(scores[[field]]), as.character(score_md[[field]])))
  for (condition in c(single, "control")) stopifnot(
    sum(md$condition == condition) == conditions$qc_accepted[conditions$condition == condition])
  stopifnot(near(scores$oriented_sample_median,
    scores$sample_median_ucell * ifelse(scores$direction == "up", 1, -1)))
  sample_counts <- as.data.frame(table(md$sample_id), stringsAsFactors = FALSE)
  names(sample_counts) <- c("sample_id", "analysis_cells")
  samples <- merge(sample_metadata, sample_counts, by = "sample_id", sort = FALSE)
  write_tsv(samples, "analysis_samples.tsv")
  write_tsv(conditions, "all_conditions_scope.tsv")

  log_event("Loading saved scWorkbench QC artifacts")
  artifact_ids <- unlist(sig_cfg$ucell_transfer[c("score_summary_artifact_id",
    "qc_key_artifact_id", "qc_pairwise_artifact_id")])
  native <- lapply(artifact_ids, function(id) wb_art_load(b, id, expected_kind = "obs_table")$table)
  names(native) <- c("score_summary", "key_qc", "pairwise")
  for (name in names(native)) write_tsv(native[[name]], paste0("native_", name, ".tsv"))
  pairwise <- native$pairwise
  key_qc <- native$key_qc
  stopifnot(nrow(pairwise) == choose(length(scored), 2L),
    setequal(c(pairwise$key_a, pairwise$key_b), scored), all(pairwise$n_overlap == length(cells)),
    setequal(key_qc$key, scored), all(key_qc$n_finite == length(cells)))
  genesets <- split(members$feature_id, members$signature_id)
  pairwise$signature_a_status <- status$final_status[match(pairwise$key_a, status$signature_id)]
  pairwise$signature_b_status <- status$final_status[match(pairwise$key_b, status$signature_id)]
  pairwise$genes_a <- lengths(genesets[pairwise$key_a])
  pairwise$genes_b <- lengths(genesets[pairwise$key_b])
  pairwise$shared_genes <- mapply(function(a, z) length(intersect(genesets[[a]], genesets[[z]])),
    pairwise$key_a, pairwise$key_b)
  pairwise$gene_jaccard <- with(pairwise, shared_genes / (genes_a + genes_b - shared_genes))
  pairwise$overlap_fraction_smaller_set <- with(pairwise, shared_genes / pmin(genes_a, genes_b))
  pairwise$shared_gene_ids <- mapply(function(a, z) paste(intersect(genesets[[a]], genesets[[z]]), collapse = ";"),
    pairwise$key_a, pairwise$key_b)
  write_tsv(pairwise, "signal_pairwise_comparisons.tsv")

  inventory <- reasons[match(status$signature_id, reasons$signature_id), ]
  inventory$program_scope <- status$program_scope
  inventory$registered_in_parse <- inventory$final_status == "accepted"
  inventory$selected_gene_ids <- vapply(inventory$signature_id, function(id) {
    if (is.null(genesets[[id]])) NA_character_ else paste(genesets[[id]], collapse = ";")
  }, character(1))
  inventory$top_genes <- status$top_genes
  inventory$discovery_status <- inventory$status
  inventory$status <- NULL
  inventory$gene_set_available <- NULL
  inventory$violin_available <- NULL
  # Saved gene-level tests: counts retain technical exclusions; no DE is rerun.
  log_event("Summarizing 26 saved differential-expression tables")
  de_by_condition <- lapply(single, function(condition) {
    id <- status$signature_id[status$cytokine == condition & status$direction == "up"]
    slug <- sub("\\.up$", "", sub("^parse_cytokine\\.", "", id))
    read_source(file.path(discovery, "de_full", paste0(slug, ".tsv.gz")))
  })
  names(de_by_condition) <- single
  for (i in seq_len(nrow(inventory))) {
    id <- inventory$signature_id[i]; cytokine <- inventory$cytokine[i]
    sign <- if (inventory$direction[i] == "up") 1 else -1
    de <- de_by_condition[[cytokine]]
    usable <- de$signature_eligible_feature
    response <- sign * de$response_logFC >= parameters$modelling$response_absolute_log2fc_minimum &
      de$response_FDR <= parameters$modelling$full_data_fdr_maximum_for_both_contrasts
    specificity <- sign * de$specificity_logFC >= parameters$modelling$specificity_absolute_log2fc_minimum &
      de$specificity_FDR <= parameters$modelling$full_data_fdr_maximum_for_both_contrasts
    stopifnot(!anyNA(c(response, specificity, usable)), !anyDuplicated(de$feature_id))
    inventory$modelled_genes[i] <- nrow(de)
    inventory$membership_eligible_modelled_genes[i] <- sum(usable)
    inventory$de_genes_vs_control[i] <- sum(usable & response)
    inventory$de_genes_vs_rest[i] <- sum(usable & specificity)
    inventory$de_genes_both_contrasts[i] <- sum(usable & response & specificity)
    e <- eligible[eligible$signature_id == id, ]
    m <- members[members$signature_id == id, ]
    stopifnot(nrow(e) == inventory$full_eligible_genes[i], nrow(m) == inventory$selected_genes[i],
      all(e$feature_id %in% de$feature_id[usable & response & specificity]))
    if (nrow(m)) {
      matched <- de[match(m$feature_id, de$feature_id), ]
      for (field in c("response_logFC", "response_p_value", "response_FDR",
                      "specificity_logFC", "specificity_p_value", "specificity_FDR"))
        stopifnot(near(m[[field]], matched[[field]]))
    }
    for (field in c("response_logFC", "specificity_logFC", "state_concordance", "target_dominance_ratio")) {
      inventory[i, paste0("selected_gene_median_", field)] <- if (nrow(m)) median(m[[field]], na.rm = TRUE) else NA_real_
    }
    for (field in c("response_FDR", "specificity_FDR")) {
      inventory[i, paste0("selected_gene_max_", field)] <- if (nrow(m)) max(m[[field]]) else NA_real_
    }
    inventory$target_samples[i] <- sum(samples$condition == cytokine)
    inventory$control_samples[i] <- sum(samples$condition == "control")
    inventory$rest_samples[i] <- sum(samples$treatment_category == "single" & samples$condition != cytokine)
    inventory$target_cells[i] <- sum(md$condition == cytokine)
    inventory$control_cells[i] <- sum(md$condition == "control")
    inventory$rest_cells[i] <- sum(md$treatment_category == "single" & md$condition != cytokine)
    for (group in c("target", "control", "rest")) {
      s <- scores[scores$signature_id == id, ]
      use <- switch(group, target = s$condition == cytokine, control = s$condition == "control",
        rest = s$treatment_category == "single" & s$condition != cytokine)
      vals <- s$sample_median_ucell[use]
      inventory[i, paste0(group, "_sample_median_iqr")] <- if (length(vals)) IQR(vals) else NA_real_
      old <- switch(group, target = "target_group_median", control = "control_group_median", rest = "other_single_group_median")
      if (length(vals)) stopifnot(near(median(vals), inventory[[old]][i]))
    }
    p <- profile[profile$signature_id == id, ]
    if (nrow(p)) {
      stopifnot(setequal(p$evaluated_condition, single), sum(p$is_target) == 1L)
      target <- p$oriented_shift_from_control[p$is_target]
      off <- p[!p$is_target, ]; off <- off[order(-off$oriented_shift_from_control, off$evaluated_condition), ]
      eps <- 1e-12
      equal <- abs(off$oriented_shift_from_control - target) <= eps
      inventory$target_oriented_shift[i] <- target
      inventory$top_offtarget_condition[i] <- off$evaluated_condition[1]
      inventory$top_offtarget_oriented_shift[i] <- off$oriented_shift_from_control[1]
      inventory$target_minus_top_offtarget_shift[i] <- target - off$oriented_shift_from_control[1]
      inventory$other_conditions_stronger[i] <- sum(off$oriented_shift_from_control > target + eps)
      inventory$other_conditions_tied[i] <- sum(equal)
      inventory$target_response_rank[i] <- 1L + inventory$other_conditions_stronger[i]
      ties <- off$evaluated_condition[abs(off$oriented_shift_from_control - off$oriented_shift_from_control[1]) <= eps]
      inventory$top_offtarget_ties[i] <- paste(ties, collapse = ";")
      q <- key_qc[key_qc$key == id, ]
      inventory$pooled_cell_score_iqr[i] <- q$score_q75 - q$score_q25
      inventory$pooled_cell_score_sd[i] <- q$score_sd
      pp <- pairwise[pairwise$key_a == id | pairwise$key_b == id, ]
      other <- ifelse(pp$key_a == id, pp$key_b, pp$key_a)
      # Include failures here to expose redundancies across all 31 candidate sets.
      best <- order(-pp$gene_jaccard, other)[1]
      inventory$most_overlapping_signature[i] <- other[best]
      inventory$maximum_gene_jaccard[i] <- pp$gene_jaccard[best]
      best <- order(-abs(pp$spearman_cor), other, na.last = TRUE)[1]
      inventory$most_correlated_signature[i] <- other[best]
      inventory$strongest_absolute_spearman_signed[i] <- pp$spearman_cor[best]
    } else {
      for (field in c("target_oriented_shift", "top_offtarget_oriented_shift", "target_minus_top_offtarget_shift",
        "other_conditions_stronger", "other_conditions_tied", "target_response_rank", "pooled_cell_score_iqr",
        "pooled_cell_score_sd", "maximum_gene_jaccard", "strongest_absolute_spearman_signed")) inventory[i, field] <- NA_real_
      for (field in c("top_offtarget_condition", "top_offtarget_ties", "most_overlapping_signature", "most_correlated_signature"))
        inventory[i, field] <- NA_character_
    }
  }
  # Preserve original statuses while adding report context to downloadable evidence.
  for (x in c("members", "eligible", "profile", "expression_profile", "quantiles", "scores")) {
    value <- get(x)
    value$final_status <- status$final_status[match(value$signature_id, status$signature_id)]
    stopifnot(!anyNA(value$final_status))
    if (x == "eligible") {
      con <- gzfile(file.path(out, "eligible_gene_evidence.tsv.gz"), "wt")
      write.table(value, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA"); close(con)
    } else write_tsv(value, switch(x, members = "selected_gene_evidence.tsv", profile = "condition_score_profiles.tsv",
      expression_profile = "condition_expression_profiles.tsv", quantiles = "cell_score_quantiles.tsv",
      scores = "sample_median_scores.tsv"))
  }
  write_tsv(inventory, "cytokine_signal_inventory.tsv")
  write.csv(inventory, file.path(out, "cytokine_signal_inventory.csv"), row.names = FALSE, na = "NA")
  write_tsv(inventory[inventory$registered_in_parse, ], "accepted_signals.tsv")
  log_event("Rendering inventory and comparison figures")
  source("figure_scripts/parse_cytokines/04_render_signal_inventory.R", local = TRUE)
  render_inventory(inventory, pairwise, conditions, out, checked)
  source("figure_scripts/parse_cytokines/04_render_signal_checks.R", local = TRUE)
  render_signal_checks(inventory, members, eligible, out, parameters, sig_cfg$ucell_transfer$acceptance)
  stopifnot(identical(before, wb_manifest(b)$current_state$fingerprint))
  inputs <<- unique(c(inputs, normalizePath(c("config/bundles.yml", "config/parse_cytokines.yml",
    "figure_scripts/parse_cytokines/04_cytokine_signal_inventory.R",
    "figure_scripts/parse_cytokines/04_render_signal_inventory.R",
    "figure_scripts/parse_cytokines/04_render_signal_checks.R",
    "figure_scripts/parse_cytokines/run_04_cytokine_signal_inventory.sh", "src/bundle_config.R",
    "src/bundle_validation.R", "src/parse_runtime.sh")), parameters_path, validation_path,
    registry$bundles$parse_cytokines$runtime$build_provenance))
  input_manifest <- data.frame(path = inputs, sha256 = vapply(inputs, wb_source_sha256, character(1)))
  write_tsv(input_manifest, "input_manifest.tsv")
  write_json(list(status = "complete", created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    bundle = checked, bundle_access = "read_only", bundle_fingerprint_unchanged = TRUE,
    runtime = registry$bundles$parse_cytokines$runtime, artifacts = as.list(artifact_ids),
    report_method = "descriptive inventory of saved evidence; no new fits, scores, or hypothesis tests",
    seed = "not applicable; deterministic report", original_discovery_parameters = parameters,
    counts = list(individual_conditions = length(single), directions = nrow(inventory), scored = length(scored),
      accepted = length(accepted), pairwise_comparisons = nrow(pairwise), analysis_cells = length(cells),
      samples = nrow(samples), selected_gene_memberships = nrow(members), registered_gene_memberships = nrow(final_sets)),
    checks = c("accepted bundle identity and fingerprint", "retained recipient validation passed",
      "all 23 registered gene sets match source membership and order", "all 52 direction statuses and eligible counts",
      "selected gene effects, raw p-values and FDR match full DE tables",
      "all eight check outcomes reconstruct the saved discovery and final statuses",
      "cell IDs, sample IDs, condition and replicate alignment",
      "all 31 score keys cover 81 biological samples", "all 465 native pairwise comparisons cover 154209 cells",
      "raw and oriented score sign consistency", "condition cell counts and saved group medians")), "validation.json")
  capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
  output_files <- list.files(out, full.names = TRUE)
  write_tsv(data.frame(file = basename(output_files), sha256 = vapply(output_files, wb_source_sha256, character(1))),
    "output_manifest.tsv")
  log_event("Complete:", file.path(out, "CHECKS.html"))
}
main()
