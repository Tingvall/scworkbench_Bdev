# Shared configuration and checks for the competitor-aware Parse-to-Bdev review.
source("src/integration_workflow.R")

rest_cfg <- read_project_yaml("config/parse_cytokines_rest_only.yml")
competitor_cfg <- cfg$competitor_aware_cytokine_analysis
source_review_dir <- rest_cfg$publication$source_review
source_promotion_dir <- rest_cfg$promotion$report
source_membership_tsv <- file.path(source_promotion_dir,
  "promoted_full_feature_sets.tsv")
source_signature_summary_tsv <- file.path(source_review_dir,
  "signature_summary.tsv")

check_competitor_aware_parents <- function() {
  ids <- c("bdev_reference", "parse_cytokines_rest_only")
  out <- lapply(ids, function(id) {
    checked <- validate_registered_bundle(registry, id,
      run_recipient_validation = FALSE)
    path <- if (id == "bdev_reference") ref_cfg$validation$report else
      registry$bundles[[id]]$validation_report
    validation <- jsonlite::read_json(path, simplifyVector = FALSE)
    stopifnot(isTRUE(validation$ok))
    checked$accepted_validation_report <- path
    checked$accepted_validation_sha256 <- wb_source_sha256(path)
    checked
  })
  names(out) <- ids
  out
}

read_source_memberships <- function() {
  sets <- read.delim(source_membership_tsv, check.names = FALSE,
    stringsAsFactors = FALSE)
  summary <- read.delim(source_signature_summary_tsv, check.names = FALSE,
    stringsAsFactors = FALSE)
  passing <- summary[summary$full_overall_pass, , drop = FALSE]
  stopifnot(nrow(passing) == 33L, nrow(sets) == 490L,
    setequal(unique(sets$feature_set), passing$signature_id))
  sets <- sets[order(match(sets$feature_set, passing$signature_id)), , drop = FALSE]
  list(sets = sets, summary = passing)
}

competitor_label <- function(ids, coverage) {
  paste(coverage$cytokine[match(ids, coverage$feature_set)],
    toupper(coverage$direction[match(ids, coverage$feature_set)]))
}
