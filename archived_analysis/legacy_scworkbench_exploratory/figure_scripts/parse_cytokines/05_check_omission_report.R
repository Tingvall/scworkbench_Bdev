#!/usr/bin/env Rscript
# Counterfactual counts from saved Parse evidence, with memberships held fixed
# for score checks and explicitly separate gene-eligibility sensitivity counts.
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(scworkbench); library(jsonlite); library(yaml)})
source("src/bundle_config.R")
source("src/bundle_validation.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- normalizePath(args[[1]], mustWork = TRUE)
stopifnot(!length(list.files(out, all.files = TRUE, no.. = TRUE)))
registry <- load_bundle_registry()
cfg <- read_project_yaml("config/parse_cytokines.yml")
sa <- cfg$signature_analysis
discovery <- sa$discovery$report
input_paths <- c("config/bundles.yml", "config/parse_cytokines.yml", "src/bundle_config.R",
  "src/bundle_validation.R", "src/parse_runtime.sh", "figure_scripts/parse_cytokines/05_check_omission_report.R",
  "figure_scripts/parse_cytokines/run_05_check_omission_report.sh")
read_source <- function(path, type = "tsv") {
  input_paths <<- unique(c(input_paths, normalizePath(path, mustWork = TRUE)))
  switch(type, tsv = read.delim(path, check.names = FALSE, na.strings = c("NA", "")),
    rds = readRDS(path), json = jsonlite::read_json(path, simplifyVector = TRUE))
}
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
main <- function() {
  checked <- validate_registered_bundle(registry, "parse_cytokines", run_recipient_validation = FALSE)
  validation_path <- registry$bundles$parse_cytokines$validation_report
  stopifnot(isTRUE(read_source(validation_path, "json")$ok))
  status <- read_source(file.path(sa$final$report, "all_signature_status.tsv"))
  pars <- read_source(file.path(discovery, "analysis_parameters.json"), "json")
  selected <- read_source(file.path(discovery, "selected_signature_members.tsv"))
  fits <- read_source(file.path(discovery, "limma_voom_fits.rds"), "rds")
  filters <- read_source(file.path(discovery, "feature_filter.tsv"))
  stopifnot(nrow(status) == 52L, !anyDuplicated(status$signature_id), length(fits$folds) == 3L)

  score_fields <- c("discovery_cv_response_auc", "discovery_cv_specificity_auc", "discovery_cv_balanced_auc",
    "discovery_fold_gene_jaccard", "ucell_response_auc", "ucell_specificity_auc", "ucell_balanced_auc")
  score_names <- c("Discovery: treatment vs control", "Discovery: treatment vs rest", "Discovery: mean AUC",
    "Discovery: fold-gene Jaccard", "UCell: treatment vs control", "UCell: treatment vs rest", "UCell: mean AUC")
  thresholds <- c(.60, .60, .65, .10, sa$ucell_transfer$acceptance$response_auc_minimum,
    sa$ucell_transfer$acceptance$specificity_auc_minimum, sa$ucell_transfer$acceptance$balanced_auc_minimum)
  passed <- cbind(minimum_eligible_genes = status$full_eligible_genes >= pars$validation$minimum_signature_genes,
    sweep(as.matrix(status[, score_fields]), 2, thresholds, ">="))
  colnames(passed) <- c("minimum_eligible_genes", score_fields)
  rownames(passed) <- status$signature_id
  baseline_ids <- status$signature_id[status$final_status == "accepted"]
  baseline <- apply(passed, 1, function(x) all(x %in% TRUE))
  stopifnot(setequal(rownames(passed)[baseline], baseline_ids), sum(baseline) == 23L)
  omissions <- setNames(lapply(seq_len(ncol(passed)), function(i) colnames(passed)[i]),
    c("Minimum five eligible genes", score_names))
  omissions <- c(list("None (baseline)" = character()), omissions,
    list("All UCell confirmation checks" = score_fields[5:7],
      "Both control AUC checks (mean checks retained)" = score_fields[c(1, 5)],
      "Control AUC and both mean-AUC checks" = score_fields[c(1, 3, 5, 7)],
      "Both rest AUC checks (mean checks retained)" = score_fields[c(2, 6)],
      "Rest AUC and both mean-AUC checks" = score_fields[c(2, 3, 6, 7)]))
  scenario_states <- list()
  counts <- do.call(rbind, lapply(names(omissions), function(name) {
    remaining <- passed[, setdiff(colnames(passed), omissions[[name]]), drop = FALSE]
    state <- apply(remaining, 1, function(x) {
      if (any(x == FALSE, na.rm = TRUE)) "excluded" else if (anyNA(x)) "not_evaluated" else "included"
    })
    included <- names(state)[state == "included"]
    stopifnot(all(baseline_ids %in% included))
    scenario_states[[name]] <<- data.frame(scenario = name, signature_id = status$signature_id,
      cytokine = status$cytokine, direction = status$direction, selected_genes = status$selected_genes,
      result = unname(state), original_final_status = status$final_status)
    data.frame(omitted_check = name, included_existing_sets = sum(state == "included"),
      added_existing_sets = length(setdiff(included, baseline_ids)),
      excluded_by_remaining_checks = sum(state == "excluded"),
      not_evaluated_on_remaining_checks = sum(state == "not_evaluated"),
      omitted_fields = paste(omissions[[name]], collapse = ";"),
      added_signature_ids = paste(setdiff(included, baseline_ids), collapse = ";"))
  }))
  write_tsv(counts, "existing_set_check_omissions.tsv")
  write_tsv(do.call(rbind, scenario_states), "existing_set_omission_details.tsv")

  # Re-evaluate gene eligibility against the ORIGINAL stored fits. This neither
  # reranks signatures nor reuses the old AUCs for changed gene memberships.
  features <- rownames(fits$full$coefficients)
  stopifnot(identical(features, fits$filtered_genes), !anyDuplicated(features))
  filter_rows <- filters[match(features, filters$feature_id), ]
  stopifnot(!anyNA(filter_rows$feature_id), identical(filter_rows$feature_id, features))
  allowed <- filter_rows$eligible_for_signature_membership
  gene_omissions <- list(
    "None (baseline)" = character(),
    "Entire treatment-vs-control gene filter" = c("control_direction", "control_effect", "control_FDR", "control_fold_sign"),
    "Entire treatment-vs-rest gene filter" = c("rest_direction", "rest_effect", "rest_FDR", "rest_fold_sign"),
    "Control gene FDR cutoff" = "control_FDR",
    "Rest gene FDR cutoff" = "rest_FDR",
    "Control gene effect-size cutoff (sign retained)" = "control_effect",
    "Rest gene effect-size cutoff (sign retained)" = "rest_effect",
    "Three-fold gene-direction stability" = c("control_fold_sign", "rest_fold_sign"),
    "Minimum five genes (nonempty sets only)" = character())
  gene_rows <- list()
  for (i in seq_len(nrow(status))) {
    s <- if (status$direction[i] == "up") 1 else -1
    r <- paste0("response::", status$cytokine[i]); z <- paste0("specificity::", status$cytokine[i])
    response <- s * fits$full$coefficients[, r]; specificity <- s * fits$full$coefficients[, z]
    fold_ok <- function(contrast) Reduce(`&`, lapply(fits$folds, function(f) {
      stopifnot(identical(rownames(f$coefficients), features))
      s * f$coefficients[, contrast] > 0
    }))
    gates <- cbind(control_direction = response > 0, rest_direction = specificity > 0,
      control_effect = response >= pars$modelling$response_absolute_log2fc_minimum,
      rest_effect = specificity >= pars$modelling$specificity_absolute_log2fc_minimum,
      control_FDR = fits$full$adjusted_p_value[, r] <= pars$modelling$full_data_fdr_maximum_for_both_contrasts,
      rest_FDR = fits$full$adjusted_p_value[, z] <= pars$modelling$full_data_fdr_maximum_for_both_contrasts,
      control_fold_sign = fold_ok(r), rest_fold_sign = fold_ok(z))
    stopifnot(!anyNA(gates), all(is.finite(fits$full$t[, c(r, z)])))
    current_genes <- features[allowed & rowSums(!gates) == 0L]
    stopifnot(length(current_genes) == status$full_eligible_genes[i],
      all(selected$feature_id[selected$signature_id == status$signature_id[i]] %in% current_genes))
    for (name in names(gene_omissions)) {
      keep <- allowed & rowSums(!gates[, setdiff(colnames(gates), gene_omissions[[name]]), drop = FALSE]) == 0L
      minimum <- if (name == "Minimum five genes (nonempty sets only)") 1L else pars$validation$minimum_signature_genes
      stopifnot(all(current_genes %in% features[keep]))
      gene_rows[[length(gene_rows) + 1L]] <- data.frame(omitted_check = name,
        signature_id = status$signature_id[i], cytokine = status$cytokine[i], direction = status$direction[i],
        eligible_genes = sum(keep), minimum_genes_required = minimum, candidate_eligible = sum(keep) >= minimum,
        baseline_eligible_genes = length(current_genes), original_final_status = status$final_status[i],
        final_acceptance_after_reselection = NA_character_)
    }
  }
  gene_details <- do.call(rbind, gene_rows)
  baseline_candidates <- status$signature_id[status$full_eligible_genes >= pars$validation$minimum_signature_genes]
  gene_counts <- do.call(rbind, lapply(names(gene_omissions), function(name) {
    d <- gene_details[gene_details$omitted_check == name & gene_details$candidate_eligible, ]
    added <- setdiff(d$signature_id, baseline_candidates)
    data.frame(omitted_gene_filter = name, eligible_candidate_sets = nrow(d),
      additional_candidates = length(added), up = sum(d$direction == "up"), down = sum(d$direction == "down"),
      cytokines = length(unique(d$cytokine)), total_eligible_gene_memberships = sum(d$eligible_genes),
      final_accepted_sets_after_reselection = if (name == "None (baseline)") 23L else NA_integer_,
      added_signature_ids = paste(added, collapse = ";"))
  }))
  stopifnot(gene_counts$eligible_candidate_sets[1] == 31L,
    all(gene_counts$eligible_candidate_sets >= 31L))
  write_tsv(gene_counts, "gene_filter_omissions.tsv")
  write_tsv(gene_details, "gene_filter_omission_details.tsv")

  md_table <- function(d) c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |")))
  text <- c("# Parse: how many sets remain when checks are omitted?", "",
    "Baseline: 52 cytokine/direction possibilities, 31 selected and scored sets, 23 accepted sets.", "",
    "## Existing gene sets: omit checks without changing membership", "",
    "Each row removes only the named requirement(s); all remaining checks still apply. NOT EVALUATED is not treated as PASS.", "",
    md_table(counts[, 1:5]), "",
    "The eight UCell failures fail the control, rest and mean-AUC cutoffs simultaneously. Removing any single one therefore rescues none. Removing all three includes the eight previously excluded scored candidates, increasing the total to 31.", "",
    "Removing the minimum-five-genes rule alone leaves 23 confirmed existing sets and 21 directions with missing downstream evidence. Eighteen of those directions have zero eligible genes and cannot form a nonempty signature; only three have 1–4 eligible genes. These are not automatically accepted.", "",
    "Mean AUC contains both control and rest comparisons. Grouped scenarios show both retaining the mean checks and removing them with the omitted comparison; neither choice rescues the eight failures while the other UCell contrast remains required.", "",
    "## Gene discovery filters: candidate counts before reselection and scoring", "",
    md_table(gene_counts[, 1:6]), "",
    "These are eligible candidate directions, NOT final accepted signatures. Counts use the saved full-data and three-fold model coefficients. Omitting an entire contrast removes its FDR, effect-size, direction and fold-sign requirements. The other contrast and its fold-sign requirements remain. Individual effect-size omissions retain the up/down sign; minimum-size omission still requires at least one gene.", "",
    "Changing these filters changes gene membership and ranking. The old AUCs cannot validate the revised sets. Final accepted totals for rebuilt sets require reselection, cross-validation and native UCell rescoring and are deliberately left NA. No comparison was refitted, no revised sets were scored, and no bundle was promoted.", "",
    "[Existing-set counts](existing_set_check_omissions.tsv) · [Which existing sets are included](existing_set_omission_details.tsv) · [Candidate counts](gene_filter_omissions.tsv) · [Candidate genes/counts by direction](gene_filter_omission_details.tsv)", "",
    paste("Accepted Parse fingerprint:", checked$bundle_state_fingerprint))
  writeLines(text, file.path(out, "REPORT.md"))
  escaped <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)
  html_table <- function(d) paste0("<div class='scroll'><table><thead><tr>",
    paste0("<th>", escaped(names(d)), "</th>", collapse = ""), "</tr></thead><tbody>",
    paste(apply(d, 1, function(x) paste0("<tr>", paste0("<td>", escaped(x), "</td>", collapse = ""), "</tr>")), collapse = ""),
    "</tbody></table></div>")
  writeLines(c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Parse check omissions</title>',
    '<style>body{font:15px/1.5 system-ui,sans-serif;max-width:1250px;margin:30px auto;padding:0 20px;color:#193249}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:9px;text-align:left;border-bottom:1px solid #d8e0e8}th{background:#eaf0f6}.scroll{overflow:auto}a{color:#076690}</style></head><body>',
    '<h1>Parse: counts with checks omitted</h1><p><b>Baseline: 31 scored gene sets; 23 accepted.</b></p>',
    '<h2>Existing sets, unchanged memberships</h2>', html_table(counts[, 1:5]),
    '<p>All eight UCell failures fail all three UCell cutoffs. Omitting one cutoff therefore changes no acceptance decisions. Removing the minimum gene-count rule leaves 21 directions with unassessed downstream checks, not 21 additional accepted sets.</p>',
    '<h2>Gene-selection changes: eligible candidates only</h2>', html_table(gene_counts[, 1:6]),
    '<p>These counts are candidates before reselection and scoring. The entire-contrast scenarios drop that contrast’s FDR, effect-size, direction and fold-sign rules. The other contrast remains. Existing AUCs do not apply to revised memberships; final acceptance is not yet determined. This report uses saved fits and scores and makes no bundle changes.</p>',
    '<p><a href="REPORT.md">Full definitions and interpretation</a> · <a href="existing_set_check_omissions.tsv">Existing-set counts</a> · <a href="existing_set_omission_details.tsv">Existing-set details</a> · <a href="gene_filter_omissions.tsv">Candidate counts</a> · <a href="gene_filter_omission_details.tsv">Candidate details</a> · <a href="validation.json">Validation</a></p></body></html>'), file.path(out, "REPORT.html"))
  after <- validate_registered_bundle(registry, "parse_cytokines", run_recipient_validation = FALSE)
  stopifnot(identical(after$bundle_state_fingerprint, checked$bundle_state_fingerprint))
  write_tsv(data.frame(path = normalizePath(input_paths), sha256 = vapply(input_paths, wb_source_sha256, character(1))), "input_manifest.tsv")
  jsonlite::write_json(list(status = "complete", bundle = checked, bundle_fingerprint_unchanged = TRUE,
    method = "saved-check logical omissions and saved-fit gene eligibility masks; no fitting, selection or scoring",
    scope = "pooled Parse control/single dataset only", source_analysis_cells = pars$analysis_scope$analysis_cells,
    existing_set_scenarios = nrow(counts), gene_filter_scenarios = nrow(gene_counts),
    checks = c("original 23 acceptance decisions reproduced", "all 52 baseline eligible-gene counts reproduced",
      "all original selected genes pass reconstructed masks", "full and fold gene alignment verified",
      "omissions never remove baseline eligible genes or accepted sets", "missing checks never count as passes"),
    seed = "not applicable; deterministic saved-result summaries", runtime = registry$bundles$parse_cytokines$runtime),
    file.path(out, "validation.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
  capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
  files <- list.files(out, full.names = TRUE)
  write_tsv(data.frame(file = basename(files), sha256 = vapply(files, wb_source_sha256, character(1))), "output_manifest.tsv")
  print(counts[, 1:5]); print(gene_counts[, 1:6])
  cat("Report:", file.path(out, "REPORT.html"), "\n")
}
main()
