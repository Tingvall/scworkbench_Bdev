#!/usr/bin/env Rscript

# Bdev-feature-compatible competitor-aware signatures with two retained
# variants per cytokine direction: a flexible 25-50 gene full set and a fixed
# top-25 set. This is a private review branch; no outputs are promoted.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 16_competitor_aware_bdev25_dual.R <new_output_directory>")
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
parent_cfg <- read_project_yaml("config/parse_cytokines.yml")
cfg <- read_project_yaml("config/parse_cytokines_rest_only.yml")
ref_cfg <- read_project_yaml("config/bdev_reference.yml")
runtime <- registry$bundles$parse_cytokines$runtime
stopifnot(identical(as.character(packageVersion("scworkbench")),
  runtime$scworkbench_version))

analysis_id <- "competitor_aware_bdev25_dual_v1"
namespace <- "parse_cytokine_rest_only_bdev25_v1"
log2fc_minimum <- 0.5
fdr_maximum <- 0.05
minimum_size <- 25L
maximum_size <- 50L
max_rank <- 1500L
chunk_size <- 400L

source_dir <- parent_cfg$signature_analysis$discovery$report
input_paths <- c(
  fits = file.path(source_dir, "limma_voom_fits.rds"),
  feature_filter = file.path(source_dir, "feature_filter.tsv"),
  parent_config = "config/parse_cytokines.yml",
  review_config = "config/parse_cytokines_rest_only.yml",
  reference_config = "config/bdev_reference.yml",
  bundle_registry = "config/bundles.yml",
  methods = "processing/parse_cytokines/COMPETITOR_AWARE_BDEV25_DUAL_WORKFLOW.md",
  script = "processing/parse_cytokines/16_competitor_aware_bdev25_dual.R",
  source_bundle_manifest = file.path(cfg$bundle_dir, "manifest.json"),
  bdev_bundle_manifest = file.path(registry$bundles$bdev_reference$bundle_dir,
    "manifest.json")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) stop("Missing inputs: ", paste(missing_inputs, collapse = ", "))

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
slug <- function(x) gsub("(^_+|_+$)", "", tolower(gsub("[^A-Za-z0-9]+", "_", x)))
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}
fmt_num <- function(x, digits = 4L) ifelse(is.finite(x),
  formatC(x, digits = digits, format = "f"), "NA")
bind_rows_fill <- function(rows, label) {
  stopifnot(length(rows) > 0L)
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    for (column in setdiff(columns, names(x))) x[[column]] <- NA
    x[, columns, drop = FALSE]
  })
  out <- do.call(rbind, unname(rows))
  rownames(out) <- NULL
  log_event(label, nrow(out), "rows")
  out
}

# Validate all accepted inputs before reading their public state.
for (id in c("parse_cytokines", "parse_cytokines_rest_only", "bdev_reference")) {
  checked <- validate_registered_bundle(registry, id,
    run_recipient_validation = FALSE)
  stopifnot(isTRUE(checked$identity_ok), isTRUE(checked$fingerprint_ok))
  validation_path <- if (id == "bdev_reference") ref_cfg$validation$report else
    registry$bundles[[id]]$validation_report
  stopifnot(isTRUE(jsonlite::read_json(validation_path,
    simplifyVector = TRUE)$ok))
}

bdev <- wb_open_bundle(registry$bundles$bdev_reference$bundle_dir,
  read_only = TRUE, validate = FALSE)
bdev_manifest <- wb_manifest(bdev)
bdev_features <- as.character(wb_feature_ids(bdev, "RNA_counts"))
stopifnot(length(bdev_features) == 32285L, !anyDuplicated(bdev_features))
wb_close_bundle(bdev)

fits <- readRDS(input_paths[["fits"]])
filters <- read_tsv(input_paths[["feature_filter"]])
features <- as.character(fits$filtered_genes)
allowed <- as.character(filters$feature_id[
  as.logical(filters$eligible_for_signature_membership)])
allowed_flag <- features %in% allowed
bdev_flag <- features %in% bdev_features
folds <- sort(names(fits$folds))
specificity_columns <- grep("^specificity::", colnames(fits$full$coefficients),
  value = TRUE)
conditions <- sub("^specificity::", "", specificity_columns)
stopifnot(length(features) == nrow(fits$full$coefficients),
  length(conditions) == 26L, identical(folds, c("1", "2", "3")),
  identical(rownames(fits$full$coefficients), features),
  all(vapply(fits$folds, function(x) identical(rownames(x$coefficients),
    features), logical(1))))

rank_genes <- function(stats, target, direction) {
  sign <- if (direction == "up") 1 else -1
  column <- paste0("specificity::", target)
  coefficient <- as.numeric(stats$coefficients[, column])
  statistic <- as.numeric(stats$t[, column])
  p_value <- as.numeric(stats$p_value[, column])
  fdr <- as.numeric(stats$adjusted_p_value[, column])
  oriented_log2fc <- sign * coefficient
  oriented_t <- sign * statistic
  keep <- allowed_flag & bdev_flag & is.finite(oriented_log2fc) &
    is.finite(oriented_t) & oriented_log2fc > log2fc_minimum
  index <- which(keep)
  strict <- is.finite(fdr[index]) & fdr[index] < fdr_maximum
  index <- index[order(-as.integer(strict), -oriented_t[index],
    -oriented_log2fc[index], fdr[index], features[index])]
  strict <- is.finite(fdr[index]) & fdr[index] < fdr_maximum
  data.frame(
    feature_id = features[index],
    specificity_log2fc = coefficient[index],
    oriented_log2fc = oriented_log2fc[index],
    specificity_t = statistic[index],
    oriented_t = oriented_t[index],
    specificity_p_value = p_value[index],
    specificity_fdr = fdr[index],
    strict_significant = strict,
    rank = seq_along(index),
    stringsAsFactors = FALSE
  )
}

construct_full <- function(ranking) {
  n_eligible <- nrow(ranking)
  n_significant <- sum(ranking$strict_significant)
  if (n_eligible < minimum_size) {
    return(list(ok = FALSE, size = NA_integer_, mode = "insufficient_effect_genes",
      genes = character()))
  }
  if (n_significant >= maximum_size) {
    size <- maximum_size
    mode <- "significant_capped_at_50"
  } else if (n_significant >= minimum_size) {
    size <- n_significant
    mode <- "all_significant"
  } else {
    size <- minimum_size
    mode <- "effect_size_fallback_to_25"
  }
  list(ok = TRUE, size = as.integer(size), mode = mode,
    genes = ranking$feature_id[seq_len(size)])
}

construct_top25 <- function(ranking) {
  if (nrow(ranking) < minimum_size) {
    return(list(ok = FALSE, size = NA_integer_, mode = "insufficient_effect_genes",
      genes = character()))
  }
  list(ok = TRUE, size = minimum_size, mode = if (
    sum(ranking$strict_significant) >= minimum_size) "significant_top25" else
      "effect_size_fallback_to_25",
    genes = ranking$feature_id[seq_len(minimum_size)])
}

signature_rows <- list()
full_rankings <- list()
fold_rankings <- list()
counter <- 0L
for (target in conditions) for (direction in c("up", "down")) {
  counter <- counter + 1L
  base_id <- paste(namespace, slug(target), direction, sep = ".")
  full_rankings[[base_id]] <- rank_genes(fits$full, target, direction)
  fold_rankings[[base_id]] <- setNames(lapply(folds, function(fold) {
    rank_genes(fits$folds[[fold]], target, direction)
  }), folds)
  ranking <- full_rankings[[base_id]]
  signature_rows[[counter]] <- data.frame(
    base_signature_id = base_id, cytokine = target, direction = direction,
    full_effect_eligible_genes = nrow(ranking),
    full_significant_genes = sum(ranking$strict_significant),
    fold1_effect_eligible_genes = nrow(fold_rankings[[base_id]][["1"]]),
    fold2_effect_eligible_genes = nrow(fold_rankings[[base_id]][["2"]]),
    fold3_effect_eligible_genes = nrow(fold_rankings[[base_id]][["3"]]),
    fold1_significant_genes = sum(fold_rankings[[base_id]][["1"]]$strict_significant),
    fold2_significant_genes = sum(fold_rankings[[base_id]][["2"]]$strict_significant),
    fold3_significant_genes = sum(fold_rankings[[base_id]][["3"]]$strict_significant),
    stringsAsFactors = FALSE)
}
signature_info <- do.call(rbind, signature_rows)
stopifnot(nrow(signature_info) == 52L, !anyDuplicated(signature_info$base_signature_id))

construction_rows <- list()
member_rows <- list()
full_sets <- list()
row_counter <- 0L
for (i in seq_len(nrow(signature_info))) {
  info <- signature_info[i, ]
  ranking <- full_rankings[[info$base_signature_id]]
  constructions <- list(full = construct_full(ranking),
    top25 = construct_top25(ranking))
  for (variant in names(constructions)) {
    row_counter <- row_counter + 1L
    built <- constructions[[variant]]
    set_id <- paste(info$base_signature_id, variant, sep = ".")
    score_key <- paste(set_id, "all", sep = ".")
    construction_rows[[row_counter]] <- data.frame(
      signature_id = set_id, base_signature_id = info$base_signature_id,
      cytokine = info$cytokine, direction = info$direction, variant = variant,
      signature_size = built$size, construction_ok = built$ok,
      construction_mode = built$mode,
      full_effect_eligible_genes = info$full_effect_eligible_genes,
      full_significant_genes = info$full_significant_genes,
      score_key = if (built$ok) score_key else NA_character_,
      stringsAsFactors = FALSE)
    if (built$ok) {
      full_sets[[score_key]] <- built$genes
      local <- ranking[match(built$genes, ranking$feature_id), , drop = FALSE]
      local$score_key <- score_key
      local$signature_id <- set_id
      local$base_signature_id <- info$base_signature_id
      local$cytokine <- info$cytokine
      local$direction <- info$direction
      local$variant <- variant
      local$signature_size <- built$size
      local$construction_mode <- built$mode
      member_rows[[length(member_rows) + 1L]] <- local[, c("score_key",
        "signature_id", "base_signature_id", "cytokine", "direction", "variant",
        "signature_size", "construction_mode", "rank", "feature_id",
        "strict_significant", "specificity_log2fc", "oriented_log2fc",
        "specificity_t", "oriented_t", "specificity_p_value", "specificity_fdr")]
    }
  }
}
construction_inventory <- do.call(rbind, construction_rows)
full_membership <- do.call(rbind, member_rows)
full_map <- construction_inventory[construction_inventory$construction_ok,
  c("score_key", "signature_id", "base_signature_id", "cytokine", "direction",
    "variant", "signature_size", "construction_mode"), drop = FALSE]
stopifnot(nrow(construction_inventory) == 104L,
  length(full_sets) == nrow(full_map),
  identical(names(full_sets), full_map$score_key),
  nrow(full_membership) == sum(full_map$signature_size))

full_ranking_evidence <- do.call(rbind, lapply(seq_len(nrow(signature_info)),
  function(i) {
    info <- signature_info[i, ]
    x <- full_rankings[[info$base_signature_id]]
    x$base_signature_id <- info$base_signature_id
    x$cytokine <- info$cytokine
    x$direction <- info$direction
    x[, c("base_signature_id", "cytokine", "direction", "rank", "feature_id",
      "strict_significant", "specificity_log2fc", "oriented_log2fc",
      "specificity_t", "oriented_t", "specificity_p_value", "specificity_fdr")]
  }))

b <- wb_open_bundle(cfg$bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(b), silent = TRUE), add = TRUE)
bundle_fingerprint_before <- wb_manifest(b)$current_state$fingerprint
cells <- as.character(wb_cells(b, "all"))
md <- wb_meta(b, cells = cells, cols = list(wb_ref_metadata_field("core:cell_id"),
  "sample_id", "condition", "replicate", "treatment_category"))
stopifnot(length(cells) == 154209L, identical(as.character(md$cell_id), cells),
  !anyNA(md))
sample_md <- unique(md[, c("sample_id", "condition", "replicate",
  "treatment_category")])
sample_md <- sample_md[match(unique(md$sample_id), sample_md$sample_id), ]
stopifnot(nrow(sample_md) == 81L, all(table(sample_md$condition) == 3L),
  setequal(unique(sample_md$condition), c("control", conditions)),
  all(table(sample_md$replicate) == 27L))
replicate_levels <- sort(unique(as.character(md$replicate)))
write_json(wb_op_interface("ucell", bundle = b),
  file.path(stage_dir, "ucell_interface.json"))

session_root <- file.path(cfg$storage_root, "sessions", cfg$run_id)
open_or_run_ucell <- function(session_id, signatures, note, score_cells) {
  workspace <- file.path(session_root, session_id)
  session <- if (dir.exists(workspace)) wb_session_resume(b, workspace) else
    wb_session_create(b, workspace_root = session_root, session_id = session_id)
  if (!is.null(session$durable$promotion_receipt)) {
    stop("Review session was unexpectedly promoted: ", workspace)
  }
  runs <- session$runs
  score_runs <- names(runs)[vapply(runs,
    function(x) identical(x$op_name, "ucell"), logical(1))]
  matching <- score_runs[vapply(score_runs, function(id) {
    tab <- runs[[id]]$outputs$signature_summary$table
    setequal(as.character(tab$signature), names(signatures))
  }, logical(1))]
  if (!length(matching)) {
    if (length(score_runs)) stop("Unexpected UCell run in ", workspace)
    log_event("Native UCell scoring", length(signatures), "sets in", session_id)
    result <- wb_session_run_op(session, b, "ucell", inputs = list(
      matrix = "RNA_counts", cell_set = "all", cell_ids = score_cells,
      min_cells = 10L, features = signatures, name = session_id,
      max_rank = max_rank, ncores = 4L, parallel_backend = "auto",
      chunk_size = chunk_size, missing_features_policy = "error"), notes = note)
    session <- result$session
    run <- result$run
  } else {
    stopifnot(length(matching) == 1L)
    log_event("Reusing completed UCell run in", session_id)
    run <- runs[[matching[[1L]]]]
  }
  native <- run$outputs$signature_summary$table
  score <- run$outputs$signature_scores
  stopifnot(setequal(as.character(native$signature), names(signatures)),
    all(native$n_features_missing == 0L),
    all(native$n_cells == length(score_cells)),
    length(score$value) == length(score_cells) * length(signatures),
    setequal(unique(score$key), names(signatures)), all(is.finite(score$value)),
    all(score$value >= 0 & score$value <= 1))
  for (signature in names(signatures)) {
    row <- native[native$signature == signature, , drop = FALSE]
    stopifnot(nrow(row) == 1L,
      row$n_features_matched == length(signatures[[signature]]),
      setequal(strsplit(row$matched_features, ";", fixed = TRUE)[[1L]],
        signatures[[signature]]))
  }
  list(session = session, score = score, native = native,
    workspace = workspace)
}

candidate_replicate_summaries <- function(score, candidate_map, score_cells,
                                          heldout_only = FALSE) {
  cell_index <- match(score_cells, cells)
  stopifnot(!anyNA(cell_index), identical(cells[cell_index], score_cells))
  local_md <- md[cell_index, , drop = FALSE]
  local_sample_md <- unique(local_md[, c("sample_id", "condition", "replicate",
    "treatment_category")])
  local_sample_md <- local_sample_md[match(unique(local_md$sample_id),
    local_sample_md$sample_id), ]
  sample_factor <- factor(local_md$sample_id, levels = local_sample_md$sample_id)
  sample_counts <- tabulate(as.integer(sample_factor), nbins = nrow(local_sample_md))
  local_replicates <- sort(unique(as.character(local_md$replicate)))
  key_rle <- rle(as.character(score$key))
  stopifnot(length(key_rle$values) == nrow(candidate_map),
    all(key_rle$lengths == length(score_cells)),
    setequal(key_rle$values, candidate_map$score_key))
  ends <- cumsum(key_rle$lengths)
  starts <- ends - key_rle$lengths + 1L
  result <- vector("list", length(key_rle$values) * 2L * length(local_replicates))
  out_counter <- 0L
  for (j in seq_along(key_rle$values)) {
    key <- key_rle$values[[j]]
    entry <- candidate_map[candidate_map$score_key == key, , drop = FALSE]
    stopifnot(nrow(entry) == 1L)
    index <- starts[[j]]:ends[[j]]
    stopifnot(identical(as.character(score$cell_id[index]), score_cells))
    values <- as.numeric(score$value[index])
    sample_mean <- as.numeric(rowsum(values, sample_factor, reorder = FALSE)) /
      sample_counts
    thresholds <- setNames(vapply(local_replicates, function(replicate) {
      use <- local_md$condition == "control" &
        as.character(local_md$replicate) == replicate
      probability <- if (entry$direction == "up") 0.95 else 0.05
      as.numeric(quantile(values[use], probability, names = FALSE, type = 7))
    }, numeric(1)), local_replicates)
    responding <- if (entry$direction == "up") {
      values > thresholds[as.character(local_md$replicate)]
    } else values < thresholds[as.character(local_md$replicate)]
    responder_fraction <- as.numeric(rowsum(as.numeric(responding), sample_factor,
      reorder = FALSE)) / sample_counts
    local <- local_sample_md
    local$sample_mean <- sample_mean
    local$responder_fraction <- responder_fraction
    reps <- if (heldout_only) as.character(entry$held_out_replicate) else
      local_replicates
    for (replicate in reps) {
      block <- local[as.character(local$replicate) == replicate, , drop = FALSE]
      stopifnot(nrow(block) == 27L, !anyDuplicated(block$condition),
        sum(block$condition == entry$cytokine) == 1L)
      for (route in c("mean", "responder")) {
        response <- if (route == "mean") {
          if (entry$direction == "up") block$sample_mean else -block$sample_mean
        } else block$responder_fraction
        target_index <- which(block$condition == entry$cytokine)
        competitor_index <- which(block$condition != entry$cytokine)
        competitor_max <- max(response[competitor_index])
        tied <- block$condition[competitor_index][
          response[competitor_index] == competitor_max]
        spread <- sd(response)
        gap <- if (is.finite(spread) && spread > 0) {
          (response[target_index] - competitor_max) / spread
        } else 0
        out_counter <- out_counter + 1L
        result[[out_counter]] <- data.frame(
          score_key = key, signature_id = entry$signature_id,
          base_signature_id = entry$base_signature_id,
          cytokine = entry$cytokine, direction = entry$direction,
          variant = entry$variant, signature_size = entry$signature_size,
          held_out_replicate = if (heldout_only)
            as.character(entry$held_out_replicate) else NA_character_,
          replicate = replicate, route = route,
          target_response = response[target_index],
          strongest_competitor_response = competitor_max,
          response_sd_across_27_conditions = spread,
          standardized_gap = gap,
          target_rank_min = 1L + sum(response[competitor_index] >
            response[target_index]),
          strict_rank1 = isTRUE(response[target_index] > competitor_max),
          strongest_competitor = paste(sort(tied), collapse = ";"),
          responder_threshold = if (route == "responder")
            thresholds[[replicate]] else NA_real_, stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, result[seq_len(out_counter)])
}

collapse_competitors <- function(x) paste(paste0("rep", x$replicate, ":",
  x$strongest_competitor), collapse = " | ")
summarize_routes <- function(replicate_rows, group_column) {
  groups <- unique(as.character(replicate_rows[[group_column]]))
  rows <- lapply(groups, function(group) {
    local <- replicate_rows[as.character(replicate_rows[[group_column]]) == group,
      , drop = FALSE]
    first <- local[1L, ]
    mean_rows <- local[local$route == "mean", , drop = FALSE]
    responder_rows <- local[local$route == "responder", , drop = FALSE]
    stopifnot(nrow(mean_rows) == 3L, nrow(responder_rows) == 3L,
      setequal(mean_rows$replicate, replicate_levels),
      setequal(responder_rows$replicate, replicate_levels))
    mean_gap <- median(mean_rows$standardized_gap)
    responder_gap <- median(responder_rows$standardized_gap)
    data.frame(group_id = group, signature_id = first$signature_id,
      base_signature_id = first$base_signature_id, cytokine = first$cytokine,
      direction = first$direction, variant = first$variant,
      signature_size = first$signature_size,
      mean_median_standardized_gap = mean_gap,
      mean_rank1_replicates = sum(mean_rows$strict_rank1),
      mean_pass = isTRUE(mean_gap > 0),
      mean_competitors = collapse_competitors(mean_rows[order(mean_rows$replicate), ]),
      responder_median_standardized_gap = responder_gap,
      responder_rank1_replicates = sum(responder_rows$strict_rank1),
      responder_pass = isTRUE(responder_gap > 0),
      responder_competitors = collapse_competitors(
        responder_rows[order(responder_rows$replicate), ]),
      stringsAsFactors = FALSE)
  })
  bind_rows_fill(rows, paste0("Route summaries by ", group_column))
}

# Full-data memberships are evaluated independently for all three replicates.
full_native_batches <- list()
full_replicate_batches <- list()
full_session_workspaces <- character()
for (replicate in replicate_levels) for (direction in c("up", "down")) {
  use <- full_map$direction == direction
  batch_map <- full_map[use, , drop = FALSE]
  batch_sets <- full_sets[batch_map$score_key]
  score_cells <- cells[as.character(md$replicate) == replicate]
  batch_id <- paste0("rep", replicate, "_", direction)
  session_id <- paste0("bdev25_v1_full_", batch_id)
  scored <- open_or_run_ucell(session_id, batch_sets,
    paste("Bdev-feature-compatible full and top25 memberships; replicate",
      replicate, "direction", direction, "; private review only."), score_cells)
  full_native_batches[[batch_id]] <- scored$native
  full_replicate_batches[[batch_id]] <- candidate_replicate_summaries(
    scored$score, batch_map, score_cells, heldout_only = FALSE)
  full_session_workspaces[[batch_id]] <- scored$workspace
  try(wb_close_bundle(wb_session_view(scored$session)), silent = TRUE)
  rm(scored)
  invisible(gc())
}
full_native <- bind_rows_fill(full_native_batches, "Full native summaries")
full_replicate_results <- bind_rows_fill(full_replicate_batches,
  "Full replicate results")
full_results <- summarize_routes(full_replicate_results, "score_key")
names(full_results)[names(full_results) == "group_id"] <- "score_key"
full_results$full_overall_pass <- full_results$mean_pass |
  full_results$responder_pass
full_results$full_status <- ifelse(full_results$full_overall_pass, "PASS", "FAIL")

# Rebuild each variant at its fixed full-data size in every training fold.
holdout_inventory_rows <- list()
holdout_member_rows <- list()
holdout_sets <- list()
hold_counter <- 0L
for (i in seq_len(nrow(full_map))) {
  entry <- full_map[i, ]
  for (fold in folds) {
    hold_counter <- hold_counter + 1L
    ranking <- fold_rankings[[entry$base_signature_id]][[fold]]
    ok <- nrow(ranking) >= entry$signature_size
    score_key <- paste(entry$signature_id, paste0("fold_", fold), sep = ".")
    mode <- if (!ok) "insufficient_effect_genes" else if (
      sum(ranking$strict_significant) >= entry$signature_size) {
      "significant_fixed_size"
    } else "effect_size_fallback_fixed_size"
    holdout_inventory_rows[[hold_counter]] <- data.frame(
      score_key = if (ok) score_key else NA_character_,
      signature_id = entry$signature_id,
      base_signature_id = entry$base_signature_id,
      cytokine = entry$cytokine, direction = entry$direction,
      variant = entry$variant, signature_size = entry$signature_size,
      held_out_replicate = fold, construction_ok = ok,
      construction_mode = mode, fold_effect_eligible_genes = nrow(ranking),
      fold_significant_genes = sum(ranking$strict_significant),
      stringsAsFactors = FALSE)
    if (ok) {
      genes <- ranking$feature_id[seq_len(entry$signature_size)]
      holdout_sets[[score_key]] <- genes
      local <- ranking[seq_len(entry$signature_size), , drop = FALSE]
      local$score_key <- score_key
      local$signature_id <- entry$signature_id
      local$base_signature_id <- entry$base_signature_id
      local$cytokine <- entry$cytokine
      local$direction <- entry$direction
      local$variant <- entry$variant
      local$signature_size <- entry$signature_size
      local$held_out_replicate <- fold
      local$construction_mode <- mode
      holdout_member_rows[[length(holdout_member_rows) + 1L]] <- local[, c(
        "score_key", "signature_id", "base_signature_id", "cytokine", "direction",
        "variant", "signature_size", "held_out_replicate", "construction_mode",
        "rank", "feature_id", "strict_significant", "specificity_log2fc",
        "oriented_log2fc", "specificity_t", "oriented_t",
        "specificity_p_value", "specificity_fdr")]
    }
  }
}
holdout_inventory <- do.call(rbind, holdout_inventory_rows)
holdout_membership <- do.call(rbind, holdout_member_rows)
holdout_map <- holdout_inventory[holdout_inventory$construction_ok,
  c("score_key", "signature_id", "base_signature_id", "cytokine", "direction",
    "variant", "signature_size", "held_out_replicate"), drop = FALSE]
stopifnot(length(holdout_sets) == nrow(holdout_map),
  nrow(holdout_inventory) == nrow(full_map) * 3L)

valid_cv_ids <- names(which(table(holdout_map$signature_id) == 3L))
holdout_native_batches <- list()
holdout_replicate_batches <- list()
holdout_session_workspaces <- character()
for (fold in folds) {
  use <- holdout_map$held_out_replicate == fold &
    holdout_map$signature_id %in% valid_cv_ids
  batch_map <- holdout_map[use, , drop = FALSE]
  if (!nrow(batch_map)) next
  batch_sets <- holdout_sets[batch_map$score_key]
  score_cells <- cells[as.character(md$replicate) == fold]
  session_id <- paste0("bdev25_v1_holdout_fold_", fold)
  scored <- open_or_run_ucell(session_id, batch_sets,
    paste("Fixed-size full and top25 training memberships scored in held-out replicate",
      fold, "; private stability review only."), score_cells)
  holdout_native_batches[[fold]] <- scored$native
  holdout_replicate_batches[[fold]] <- candidate_replicate_summaries(
    scored$score, batch_map, score_cells, heldout_only = TRUE)
  holdout_session_workspaces[[fold]] <- scored$workspace
  try(wb_close_bundle(wb_session_view(scored$session)), silent = TRUE)
  rm(scored)
  invisible(gc())
}
holdout_native <- bind_rows_fill(holdout_native_batches,
  "Holdout native summaries")
holdout_replicate_results <- bind_rows_fill(holdout_replicate_batches,
  "Holdout replicate results")
holdout_results <- summarize_routes(holdout_replicate_results, "signature_id")
names(holdout_results)[names(holdout_results) == "group_id"] <- "holdout_group_id"

jaccard <- function(a, b) {
  union_size <- length(union(a, b))
  if (!union_size) 1 else length(intersect(a, b)) / union_size
}
pairwise_fold_jaccard <- function(sets) {
  pairs <- combn(seq_along(sets), 2L)
  mean(apply(pairs, 2L, function(pair) jaccard(sets[[pair[[1L]]]],
    sets[[pair[[2L]]]])))
}

final_member_rows <- list()
core_member_rows <- list()
stability_rows <- list()
for (i in seq_len(nrow(full_map))) {
  entry <- full_map[i, ]
  final_genes <- full_sets[[entry$score_key]]
  fold_rows <- holdout_map[holdout_map$signature_id == entry$signature_id,
    , drop = FALSE]
  cv_evaluable <- nrow(fold_rows) == 3L
  fold_genes <- if (cv_evaluable) setNames(lapply(folds, function(fold) {
    key <- fold_rows$score_key[fold_rows$held_out_replicate == fold]
    stopifnot(length(key) == 1L)
    holdout_sets[[key]]
  }), folds) else setNames(vector("list", length(folds)), folds)
  core_genes <- if (cv_evaluable) {
    core <- Reduce(intersect, c(list(final_genes), fold_genes))
    final_genes[final_genes %in% core]
  } else character()
  local <- full_membership[full_membership$score_key == entry$score_key,
    , drop = FALSE]
  stopifnot(identical(local$feature_id, final_genes))
  for (fold in folds) {
    local[[paste0("selected_in_fold", fold)]] <- if (cv_evaluable)
      local$feature_id %in% fold_genes[[fold]] else NA
  }
  local$fold_selection_count <- if (cv_evaluable)
    rowSums(local[, paste0("selected_in_fold", folds), drop = FALSE]) else NA_integer_
  local$is_core <- if (cv_evaluable) local$fold_selection_count == 3L else FALSE
  final_member_rows[[i]] <- local
  if (length(core_genes)) {
    core_local <- local[match(core_genes, local$feature_id), , drop = FALSE]
    core_local$core_rank <- seq_along(core_genes)
    core_local$core_feature_set_id <- paste0(entry$signature_id, ".core")
    core_member_rows[[length(core_member_rows) + 1L]] <- core_local[, c(
      "core_feature_set_id", "signature_id", "base_signature_id", "cytokine",
      "direction", "variant", "signature_size", "core_rank", "feature_id",
      "strict_significant", "specificity_log2fc", "specificity_t",
      "specificity_p_value", "specificity_fdr")]
  }
  stability_rows[[i]] <- data.frame(signature_id = entry$signature_id,
    cv_evaluable = cv_evaluable,
    fold_pairwise_jaccard = if (cv_evaluable)
      pairwise_fold_jaccard(fold_genes) else NA_real_,
    final_fold1_jaccard = if (cv_evaluable)
      jaccard(final_genes, fold_genes[["1"]]) else NA_real_,
    final_fold2_jaccard = if (cv_evaluable)
      jaccard(final_genes, fold_genes[["2"]]) else NA_real_,
    final_fold3_jaccard = if (cv_evaluable)
      jaccard(final_genes, fold_genes[["3"]]) else NA_real_,
    core_size = if (cv_evaluable) length(core_genes) else NA_integer_,
    stringsAsFactors = FALSE)
}
final_membership <- do.call(rbind, final_member_rows)
core_membership <- if (length(core_member_rows)) do.call(rbind, core_member_rows) else
  data.frame(core_feature_set_id = character(), signature_id = character(),
    base_signature_id = character(), cytokine = character(), direction = character(),
    variant = character(), signature_size = integer(), core_rank = integer(),
    feature_id = character(), strict_significant = logical(),
    specificity_log2fc = numeric(), specificity_t = numeric(),
    specificity_p_value = numeric(), specificity_fdr = numeric())
stability <- do.call(rbind, stability_rows)

# Assemble one row per separately retained full/top25 variant.
full_eval <- full_results[, c("signature_id", "mean_median_standardized_gap",
  "mean_rank1_replicates", "mean_pass", "mean_competitors",
  "responder_median_standardized_gap", "responder_rank1_replicates",
  "responder_pass", "responder_competitors")]
names(full_eval)[-1L] <- paste0("full_", names(full_eval)[-1L])
cv_eval <- holdout_results[, c("signature_id", "mean_median_standardized_gap",
  "mean_rank1_replicates", "mean_pass", "mean_competitors",
  "responder_median_standardized_gap", "responder_rank1_replicates",
  "responder_pass", "responder_competitors")]
names(cv_eval)[-1L] <- paste0("holdout_", names(cv_eval)[-1L])
summary_table <- merge(construction_inventory, full_eval, by = "signature_id",
  all.x = TRUE, sort = FALSE)
summary_table <- merge(summary_table, cv_eval, by = "signature_id",
  all.x = TRUE, sort = FALSE)
summary_table <- merge(summary_table, stability, by = "signature_id",
  all.x = TRUE, sort = FALSE)
summary_table <- summary_table[match(construction_inventory$signature_id,
  summary_table$signature_id), ]
summary_table$full_overall_pass <- with(summary_table,
  ifelse(construction_ok, full_mean_pass | full_responder_pass, FALSE))
summary_table$full_status <- ifelse(!summary_table$construction_ok,
  "CONSTRUCTION_FAIL", ifelse(summary_table$full_overall_pass, "PASS", "FAIL"))
summary_table$holdout_stable <- with(summary_table,
  ifelse(!cv_evaluable, NA, holdout_mean_pass | holdout_responder_pass))
stopifnot(nrow(summary_table) == 104L, !anyDuplicated(summary_table$signature_id),
  all(summary_table$signature_size[summary_table$construction_ok] >= minimum_size),
  all(summary_table$signature_size[summary_table$construction_ok] <= maximum_size))

write_sets <- function(members, path_tsv, path_gmt, core = FALSE) {
  if (!nrow(members)) {
    table <- data.frame(feature_set = character(), set_name = character(),
      feature_id = character(), description = character())
    write_tsv(table, path_tsv)
    writeLines(character(), path_gmt)
    return(invisible(NULL))
  }
  ids <- if (core) members$core_feature_set_id else members$signature_id
  table <- data.frame(feature_set = ids, set_name = ids,
    feature_id = members$feature_id,
    description = if (core) paste0("Stable three-fold core of ", members$variant,
      " Bdev-compatible cytokine signature; review only") else paste0(
      members$variant, " Bdev-compatible cytokine signature; review only"),
    stringsAsFactors = FALSE)
  write_tsv(table, path_tsv)
  lines <- vapply(unique(table$feature_set), function(id) {
    x <- table[table$feature_set == id, , drop = FALSE]
    paste(c(id, x$description[[1L]], x$feature_id), collapse = "\t")
  }, character(1))
  writeLines(lines, path_gmt)
}

# Evidence tables and separate exports for both retained variants.
write_tsv(signature_info, file.path(stage_dir, "signature_inventory.tsv"))
write_tsv(construction_inventory, file.path(stage_dir, "construction_inventory.tsv"))
write_tsv(summary_table, file.path(stage_dir, "signature_summary.tsv"))
write_tsv(summary_table[summary_table$variant == "full", ],
  file.path(stage_dir, "full_signature_summary.tsv"))
write_tsv(summary_table[summary_table$variant == "top25", ],
  file.path(stage_dir, "top25_signature_summary.tsv"))
write_tsv(final_membership, file.path(stage_dir, "final_gene_membership.tsv"))
write_tsv(final_membership[final_membership$variant == "full", ],
  file.path(stage_dir, "full_gene_membership.tsv"))
write_tsv(final_membership[final_membership$variant == "top25", ],
  file.path(stage_dir, "top25_gene_membership.tsv"))
write_tsv(core_membership, file.path(stage_dir, "core_gene_membership.tsv"))
write_tsv(core_membership[core_membership$variant == "full", ],
  file.path(stage_dir, "full_core_gene_membership.tsv"))
write_tsv(core_membership[core_membership$variant == "top25", ],
  file.path(stage_dir, "top25_core_gene_membership.tsv"))
write_tsv(holdout_inventory, file.path(stage_dir, "holdout_construction_inventory.tsv"))
write_tsv(holdout_membership, file.path(stage_dir, "holdout_gene_membership.tsv"))
write_tsv(full_replicate_results, file.path(stage_dir, "full_replicate_results.tsv"))
write_tsv(holdout_replicate_results,
  file.path(stage_dir, "holdout_replicate_results.tsv"))
write_tsv(full_native, file.path(stage_dir, "native_ucell_full_summary.tsv"))
write_tsv(holdout_native, file.path(stage_dir, "native_ucell_holdout_summary.tsv"))
write_tsv_gz(full_ranking_evidence,
  file.path(stage_dir, "full_gene_rankings.tsv.gz"))

for (variant in c("full", "top25")) {
  members <- final_membership[final_membership$variant == variant, ]
  cores <- core_membership[core_membership$variant == variant, ]
  write_sets(members,
    file.path(stage_dir, paste0("constructed_", variant, "_gene_sets.tsv")),
    file.path(stage_dir, paste0("constructed_", variant, "_gene_sets.gmt")))
  write_sets(cores,
    file.path(stage_dir, paste0(variant, "_core_gene_sets.tsv")),
    file.path(stage_dir, paste0(variant, "_core_gene_sets.gmt")), core = TRUE)
}

file.copy(input_paths[["methods"]], file.path(stage_dir, "METHODS.md"),
  overwrite = FALSE)
input_manifest <- data.frame(path = normalizePath(input_paths, mustWork = TRUE),
  sha256 = vapply(input_paths, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))

ordered <- summary_table[order(match(summary_table$cytokine, conditions),
  match(summary_table$direction, c("up", "down")),
  match(summary_table$variant, c("full", "top25"))), ]
report <- c("# Bdev-compatible full and top25 cytokine signatures", "",
  "This is a separate, unpromoted review branch. The earlier competitor-aware analysis and accepted feature-set namespaces are unchanged.", "",
  paste0("Bdev contributes only its ", length(bdev_features),
    " RNA feature identifiers. Gene eligibility requires direction-oriented log2FC > ",
    log2fc_minimum, "; strict significance requires FDR < ", fdr_maximum, "."), "",
  paste0("Constructed full sets: **", sum(summary_table$variant == "full" &
    summary_table$construction_ok), "/52**; constructed top25 sets: **",
    sum(summary_table$variant == "top25" & summary_table$construction_ok),
    "/52**. No winner is selected between variants."), "",
  paste0("Full-data passes: full **", sum(summary_table$variant == "full" &
    summary_table$full_overall_pass, na.rm = TRUE), "**, top25 **",
    sum(summary_table$variant == "top25" & summary_table$full_overall_pass,
      na.rm = TRUE), "**. Holdout-stable: full **",
    sum(summary_table$variant == "full" & summary_table$holdout_stable,
      na.rm = TRUE), "**, top25 **",
    sum(summary_table$variant == "top25" & summary_table$holdout_stable,
      na.rm = TRUE), "**."), "",
  "| Cytokine | Direction | Variant | Size | Construction | Significant genes | Mean | Responder | Full status | Holdout | Core |",
  "|---|---|---|---:|---|---:|---|---|---|---|---:|")
for (i in seq_len(nrow(ordered))) {
  x <- ordered[i, ]
  report <- c(report, paste0("| ", x$cytokine, " | ", x$direction,
    " | ", x$variant, " | ", ifelse(is.na(x$signature_size), "-", x$signature_size),
    " | ", x$construction_mode, " | ", x$full_significant_genes,
    " | ", ifelse(isTRUE(x$full_mean_pass), "PASS", "FAIL"),
    " | ", ifelse(isTRUE(x$full_responder_pass), "PASS", "FAIL"),
    " | **", x$full_status, "** | ",
    ifelse(is.na(x$holdout_stable), "not evaluable",
      ifelse(x$holdout_stable, "stable", "not stable")),
    " | ", ifelse(is.na(x$core_size), "-", x$core_size), " |"))
}
report <- c(report, "", "## Main files", "",
  "- `signature_summary.tsv`: all 104 separately retained variant results.",
  "- `full_signature_summary.tsv` and `top25_signature_summary.tsv`: variant-specific results.",
  "- `full_gene_membership.tsv` and `top25_gene_membership.tsv`: exact full-data memberships.",
  "- `full_core_gene_membership.tsv` and `top25_core_gene_membership.tsv`: separate stable cores.",
  "- `full_replicate_results.tsv` and `holdout_replicate_results.tsv`: competitor-aware evidence.",
  "- `METHODS.md`: exact construction and interpretation rules.")
writeLines(report, file.path(stage_dir, "REPORT.md"))

bundle_fingerprint_after <- wb_manifest(b)$current_state$fingerprint
stopifnot(identical(bundle_fingerprint_before, bundle_fingerprint_after),
  all(final_membership$feature_id %in% bdev_features),
  all(core_membership$feature_id %in% final_membership$feature_id),
  nrow(final_membership) == sum(full_map$signature_size),
  all(table(final_membership$signature_id) == full_map$signature_size[
    match(names(table(final_membership$signature_id)), full_map$signature_id)]))

completion <- list(status = "complete", promoted = FALSE,
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  analysis_id = analysis_id, namespace = namespace,
  source_bundle = cfg$bundle_dir,
  source_bundle_fingerprint = bundle_fingerprint_before,
  source_bundle_unchanged = TRUE,
  bdev_bundle = registry$bundles$bdev_reference$bundle_dir,
  bdev_bundle_fingerprint = bdev_manifest$current_state$fingerprint,
  bdev_use = "RNA_counts feature identifiers only",
  candidate_directions = 52L, variants = c("full", "top25"),
  constructed_full_sets = sum(summary_table$variant == "full" &
    summary_table$construction_ok),
  constructed_top25_sets = sum(summary_table$variant == "top25" &
    summary_table$construction_ok),
  full_variant_passes = sum(summary_table$variant == "full" &
    summary_table$full_overall_pass, na.rm = TRUE),
  top25_variant_passes = sum(summary_table$variant == "top25" &
    summary_table$full_overall_pass, na.rm = TRUE),
  parameters = list(log2fc_minimum_strictly_greater = log2fc_minimum,
    fdr_maximum_strictly_less = fdr_maximum, minimum_size = minimum_size,
    maximum_size = maximum_size,
    ranking = "strict FDR class, directional moderated t, directional log2FC, FDR, feature ID",
    full_construction = "all significant genes from 25 through 49; cap at 50; effect-size fallback to 25",
    top25_construction = "first 25 genes in the same ranking",
    holdout_size = "fixed at the corresponding full-data variant size",
    holdout_role = "stability annotation only",
    core = "full-data membership intersected with all three fixed-size training-fold memberships"),
  runtime = runtime,
  full_session_workspaces = as.list(full_session_workspaces),
  holdout_session_workspaces = as.list(holdout_session_workspaces),
  result_promotions = 0L)
write_json(completion, file.path(stage_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

output_files <- list.files(stage_dir, recursive = TRUE, full.names = TRUE)
output_files <- output_files[basename(output_files) != "output_manifest.tsv"]
output_manifest <- data.frame(
  file = substring(output_files, nchar(stage_dir) + 2L),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE)
write_tsv(output_manifest[order(output_manifest$file), ],
  file.path(stage_dir, "output_manifest.tsv"))

if (!file.rename(stage_dir, output_dir)) stop("Failed to finalize output directory")
log_event("Bdev-compatible dual-signature review complete at", output_dir)
