#!/usr/bin/env Rscript

# Full-data competitor-aware size selection followed by selected-size
# leave-one-replicate-out stability analysis. All scoring uses native
# scWorkbench UCell in private sessions; no result is promoted.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 14_competitor_aware_signature_review.R <new_output_directory>")
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
runtime <- registry$bundles$parse_cytokines$runtime
stopifnot(identical(as.character(packageVersion("scworkbench")), runtime$scworkbench_version))

candidate_sizes <- c(5L, 10L, 15L, 20L, 30L, 40L, 50L)
log2fc_minimum <- 0.25
max_rank <- 1500L
chunk_size <- 400L
random_seed <- 1729L
set.seed(random_seed)

source_dir <- parent_cfg$signature_analysis$discovery$report
input_paths <- c(
  fits = file.path(source_dir, "limma_voom_fits.rds"),
  feature_filter = file.path(source_dir, "feature_filter.tsv"),
  parent_config = "config/parse_cytokines.yml",
  review_config = "config/parse_cytokines_rest_only.yml",
  bundle_registry = "config/bundles.yml",
  methods = "processing/parse_cytokines/COMPETITOR_AWARE_SIGNATURE_WORKFLOW.md",
  script = "processing/parse_cytokines/14_competitor_aware_signature_review.R",
  bundle_manifest = file.path(cfg$bundle_dir, "manifest.json")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) stop("Missing inputs: ", paste(missing_inputs, collapse = ", "))

read_tsv <- function(path) read.delim(path, check.names = FALSE, na.strings = c("NA", ""))
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
fmt_num <- function(x, digits = 4L) ifelse(is.finite(x), formatC(x, digits = digits,
  format = "f"), "NA")
bind_rows_fill <- function(rows, label) {
  stopifnot(length(rows) > 0L)
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  missing <- lapply(rows, function(x) setdiff(columns, names(x)))
  if (any(lengths(missing) > 0L)) {
    detail <- paste(vapply(seq_along(rows), function(i) {
      if (!length(missing[[i]])) return(NA_character_)
      paste0(names(rows)[[i]], "=[", paste(missing[[i]], collapse = ","), "]")
    }, character(1)), collapse = "; ")
    log_event(label, "tables have optional columns; filling missing values:", detail)
  }
  rows <- lapply(rows, function(x) {
    for (column in setdiff(columns, names(x))) x[[column]] <- NA
    x[, columns, drop = FALSE]
  })
  do.call(rbind, unname(rows))
}

fits <- readRDS(input_paths[["fits"]])
filters <- read_tsv(input_paths[["feature_filter"]])
features <- as.character(fits$filtered_genes)
allowed <- as.character(filters$feature_id[as.logical(filters$eligible_for_signature_membership)])
allowed_flag <- features %in% allowed
folds <- sort(names(fits$folds))
specificity_columns <- grep("^specificity::", colnames(fits$full$coefficients), value = TRUE)
conditions <- sub("^specificity::", "", specificity_columns)
stopifnot(length(features) == nrow(fits$full$coefficients), length(conditions) == 26L,
  length(folds) == 3L, all(folds %in% c("1", "2", "3")),
  identical(rownames(fits$full$coefficients), features),
  all(vapply(fits$folds, function(x) identical(rownames(x$coefficients), features), logical(1))))

rank_genes <- function(stats, target, direction) {
  sign <- if (direction == "up") 1 else -1
  column <- paste0("specificity::", target)
  coefficient <- as.numeric(stats$coefficients[, column])
  statistic <- as.numeric(stats$t[, column])
  p_value <- as.numeric(stats$p_value[, column])
  fdr <- as.numeric(stats$adjusted_p_value[, column])
  oriented_log2fc <- sign * coefficient
  oriented_t <- sign * statistic
  keep <- allowed_flag & is.finite(oriented_log2fc) & is.finite(oriented_t) &
    oriented_log2fc >= log2fc_minimum
  index <- which(keep)
  index <- index[order(-oriented_t[index], -oriented_log2fc[index], fdr[index], features[index])]
  data.frame(
    feature_id = features[index],
    specificity_log2fc = coefficient[index],
    oriented_log2fc = oriented_log2fc[index],
    specificity_t = statistic[index],
    oriented_t = oriented_t[index],
    specificity_p_value = p_value[index],
    specificity_fdr = fdr[index],
    rank = seq_along(index),
    stringsAsFactors = FALSE
  )
}

signature_rows <- list()
full_rankings <- list()
fold_rankings <- list()
counter <- 0L
for (target in conditions) for (direction in c("up", "down")) {
  counter <- counter + 1L
  signature_id <- paste(cfg$namespace, "competitor_aware", slug(target), direction, sep = ".")
  full_rankings[[signature_id]] <- rank_genes(fits$full, target, direction)
  fold_rankings[[signature_id]] <- setNames(lapply(folds, function(fold) {
    rank_genes(fits$folds[[fold]], target, direction)
  }), folds)
  signature_rows[[counter]] <- data.frame(
    signature_id = signature_id,
    cytokine = target,
    direction = direction,
    full_eligible_genes = nrow(full_rankings[[signature_id]]),
    fold1_eligible_genes = nrow(fold_rankings[[signature_id]][["1"]]),
    fold2_eligible_genes = nrow(fold_rankings[[signature_id]][["2"]]),
    fold3_eligible_genes = nrow(fold_rankings[[signature_id]][["3"]]),
    stringsAsFactors = FALSE
  )
}
signature_info <- do.call(rbind, signature_rows)
stopifnot(nrow(signature_info) == 52L, all(signature_info$full_eligible_genes >= max(candidate_sizes)),
  all(signature_info[, grep("^fold[123]_eligible", names(signature_info))] >= max(candidate_sizes)))

full_map_rows <- list()
full_member_rows <- list()
full_sets <- list()
map_counter <- member_counter <- 0L
for (i in seq_len(nrow(signature_info))) {
  info <- signature_info[i, ]
  ranking <- full_rankings[[info$signature_id]]
  for (size in candidate_sizes) {
    map_counter <- map_counter + 1L
    score_key <- paste(info$signature_id, paste0("size_", size), "full", sep = ".")
    genes <- ranking$feature_id[seq_len(size)]
    full_sets[[score_key]] <- genes
    full_map_rows[[map_counter]] <- data.frame(
      score_key = score_key, signature_id = info$signature_id,
      cytokine = info$cytokine, direction = info$direction,
      signature_size = size, stringsAsFactors = FALSE
    )
    local <- ranking[seq_len(size), ]
    local$score_key <- score_key
    local$signature_id <- info$signature_id
    local$cytokine <- info$cytokine
    local$direction <- info$direction
    local$signature_size <- size
    full_member_rows[[map_counter]] <- local[, c("score_key", "signature_id", "cytokine",
      "direction", "signature_size", "rank", "feature_id", "specificity_log2fc",
      "specificity_t", "specificity_p_value", "specificity_fdr")]
  }
}
full_map <- do.call(rbind, full_map_rows)
full_membership <- do.call(rbind, full_member_rows)
stopifnot(nrow(full_map) == 52L * length(candidate_sizes), length(full_sets) == nrow(full_map),
  identical(names(full_sets), full_map$score_key), !anyDuplicated(full_map$score_key))

full_ranking_evidence <- do.call(rbind, lapply(seq_len(nrow(signature_info)), function(i) {
  info <- signature_info[i, ]
  x <- full_rankings[[info$signature_id]]
  x$signature_id <- info$signature_id
  x$cytokine <- info$cytokine
  x$direction <- info$direction
  x[, c("signature_id", "cytokine", "direction", "rank", "feature_id",
    "specificity_log2fc", "oriented_log2fc", "specificity_t", "oriented_t",
    "specificity_p_value", "specificity_fdr")]
}))

parent_check <- validate_registered_bundle(registry, "parse_cytokines", run_recipient_validation = FALSE)
stopifnot(isTRUE(jsonlite::read_json(registry$bundles$parse_cytokines$validation_report,
  simplifyVector = TRUE)$ok))
b <- wb_open_bundle(cfg$bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(b), silent = TRUE), add = TRUE)
bundle_fingerprint_before <- wb_manifest(b)$current_state$fingerprint
cells <- as.character(wb_cells(b, "all"))
md <- wb_meta(b, cells = cells, cols = list(wb_ref_metadata_field("core:cell_id"),
  "sample_id", "condition", "replicate", "treatment_category"))
stopifnot(length(cells) == 154209L, identical(as.character(md$cell_id), cells), !anyNA(md))
sample_md <- unique(md[, c("sample_id", "condition", "replicate", "treatment_category")])
sample_md <- sample_md[match(unique(md$sample_id), sample_md$sample_id), ]
stopifnot(nrow(sample_md) == 81L, all(table(sample_md$condition) == 3L),
  setequal(unique(sample_md$condition), c("control", conditions)),
  all(table(sample_md$replicate) == 27L))

replicate_levels <- sort(unique(as.character(md$replicate)))

open_or_run_ucell <- function(session_id, signatures, note, score_cells = cells) {
  workspace <- file.path(cfg$storage_root, "sessions", cfg$run_id, session_id)
  session <- if (dir.exists(workspace)) wb_session_resume(b, workspace) else
    wb_session_create(b, workspace_root = file.path(cfg$storage_root, "sessions", cfg$run_id),
      session_id = session_id)
  runs <- session$runs
  score_runs <- names(runs)[vapply(runs, function(x) identical(x$op_name, "ucell"), logical(1))]
  matching <- score_runs[vapply(score_runs, function(id) {
    tab <- runs[[id]]$outputs$signature_summary$table
    setequal(as.character(tab$signature), names(signatures))
  }, logical(1))]
  if (!length(matching)) {
    if (length(score_runs)) stop("Session contains an unexpected UCell run: ", workspace)
    log_event("Native UCell scoring", length(signatures), "sets in", session_id)
    result <- wb_session_run_op(session, b, "ucell", inputs = list(
      matrix = "RNA_counts", cell_set = "all", cell_ids = score_cells, min_cells = 10L,
      features = signatures, name = session_id, max_rank = max_rank,
      ncores = 4L, parallel_backend = "auto", chunk_size = chunk_size,
      missing_features_policy = "error"), notes = note)
    session <- result$session
    run <- result$run
  } else {
    stopifnot(length(matching) == 1L)
    log_event("Reusing completed UCell run in", session_id)
    run <- runs[[matching]]
  }
  native <- run$outputs$signature_summary$table
  score <- run$outputs$signature_scores
  stopifnot(setequal(as.character(native$signature), names(signatures)),
    all(native$n_features_missing == 0L), all(native$n_cells == length(score_cells)),
    length(score$value) == length(score_cells) * length(signatures),
    setequal(unique(score$key), names(signatures)), all(is.finite(score$value)),
    all(score$value >= 0 & score$value <= 1))
  for (signature in names(signatures)) {
    row <- native[native$signature == signature, , drop = FALSE]
    stopifnot(nrow(row) == 1L, row$n_features_matched == length(signatures[[signature]]),
      setequal(strsplit(row$matched_features, ";", fixed = TRUE)[[1L]], signatures[[signature]]))
  }
  list(session = session, run = run, score = score, native = native, workspace = workspace)
}

candidate_replicate_summaries <- function(score, candidate_map, score_cells,
                                          heldout_only = FALSE) {
  cell_index <- match(score_cells, cells)
  stopifnot(!anyNA(cell_index), identical(cells[cell_index], score_cells))
  local_md <- md[cell_index, , drop = FALSE]
  local_sample_md <- unique(local_md[, c("sample_id", "condition", "replicate",
    "treatment_category")])
  local_sample_md <- local_sample_md[match(unique(local_md$sample_id), local_sample_md$sample_id), ]
  local_sample_factor <- factor(local_md$sample_id, levels = local_sample_md$sample_id)
  local_sample_counts <- tabulate(as.integer(local_sample_factor), nbins = nrow(local_sample_md))
  local_replicates <- sort(unique(as.character(local_md$replicate)))
  key_rle <- rle(as.character(score$key))
  stopifnot(length(key_rle$values) == nrow(candidate_map),
    all(key_rle$lengths == length(score_cells)), setequal(key_rle$values, candidate_map$score_key))
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
    sample_mean <- as.numeric(rowsum(values, local_sample_factor, reorder = FALSE)) /
      local_sample_counts
    thresholds <- setNames(vapply(local_replicates, function(replicate) {
      use <- local_md$condition == "control" & as.character(local_md$replicate) == replicate
      probability <- if (entry$direction == "up") 0.95 else 0.05
      as.numeric(quantile(values[use], probability, names = FALSE, type = 7))
    }, numeric(1)), local_replicates)
    cell_threshold <- thresholds[as.character(local_md$replicate)]
    responding <- if (entry$direction == "up") values > cell_threshold else values < cell_threshold
    responder_fraction <- as.numeric(rowsum(as.numeric(responding), local_sample_factor,
      reorder = FALSE)) / local_sample_counts
    local <- local_sample_md
    local$sample_mean <- sample_mean
    local$responder_fraction <- responder_fraction
    reps <- if (heldout_only) as.character(entry$held_out_replicate) else local_replicates
    stopifnot(all(reps %in% local_replicates))
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
        tied <- block$condition[competitor_index][response[competitor_index] == competitor_max]
        spread <- sd(response)
        gap <- if (is.finite(spread) && spread > 0) {
          (response[target_index] - competitor_max) / spread
        } else 0
        out_counter <- out_counter + 1L
        result[[out_counter]] <- data.frame(
          score_key = key,
          signature_id = entry$signature_id,
          cytokine = entry$cytokine,
          direction = entry$direction,
          signature_size = entry$signature_size,
          held_out_replicate = if (heldout_only) as.character(entry$held_out_replicate) else NA_character_,
          replicate = replicate,
          route = route,
          target_response = response[target_index],
          strongest_competitor_response = competitor_max,
          response_sd_across_27_conditions = spread,
          standardized_gap = gap,
          target_rank_min = 1L + sum(response[competitor_index] > response[target_index]),
          strict_rank1 = isTRUE(response[target_index] > competitor_max),
          strongest_competitor = paste(sort(tied), collapse = ";"),
          responder_threshold = if (route == "responder") thresholds[[replicate]] else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
    if (j %% 25L == 0L || j == length(key_rle$values)) {
      log_event("Summarized", j, "of", length(key_rle$values), "UCell sets")
    }
  }
  do.call(rbind, result[seq_len(out_counter)])
}

collapse_competitors <- function(x) paste(paste0("rep", x$replicate, ":", x$strongest_competitor),
  collapse = " | ")

summarize_routes <- function(replicate_rows, group_column) {
  groups <- unique(as.character(replicate_rows[[group_column]]))
  rows <- lapply(groups, function(group) {
    local <- replicate_rows[as.character(replicate_rows[[group_column]]) == group, , drop = FALSE]
    first <- local[1L, ]
    mean_rows <- local[local$route == "mean", , drop = FALSE]
    responder_rows <- local[local$route == "responder", , drop = FALSE]
    stopifnot(nrow(mean_rows) == 3L, nrow(responder_rows) == 3L,
      setequal(mean_rows$replicate, replicate_levels), setequal(responder_rows$replicate, replicate_levels))
    mean_gap <- median(mean_rows$standardized_gap)
    responder_gap <- median(responder_rows$standardized_gap)
    data.frame(
      group_id = group,
      signature_id = first$signature_id,
      cytokine = first$cytokine,
      direction = first$direction,
      signature_size = first$signature_size,
      mean_median_standardized_gap = mean_gap,
      mean_rank1_replicates = sum(mean_rows$strict_rank1),
      mean_pass = isTRUE(mean_gap > 0),
      mean_competitors = collapse_competitors(mean_rows[order(mean_rows$replicate), ]),
      responder_median_standardized_gap = responder_gap,
      responder_rank1_replicates = sum(responder_rows$strict_rank1),
      responder_pass = isTRUE(responder_gap > 0),
      responder_competitors = collapse_competitors(responder_rows[order(responder_rows$replicate), ]),
      stringsAsFactors = FALSE
    )
  })
  bind_rows_fill(rows, paste0("Route summaries grouped by ", group_column))
}

log_event("Scoring full-data size grid", length(full_sets), "sets in replicate-direction batches")
full_native_batches <- list()
full_replicate_batches <- list()
full_session_workspaces <- character()
for (replicate in replicate_levels) for (direction in c("up", "down")) {
  use <- full_map$direction == direction
  batch_map <- full_map[use, , drop = FALSE]
  batch_sets <- full_sets[batch_map$score_key]
  score_cells <- cells[as.character(md$replicate) == replicate]
  batch_id <- paste0("rep", replicate, "_", direction)
  session_id <- paste0("competitor_aware_full_", batch_id)
  full_scored <- open_or_run_ucell(
    session_id,
    batch_sets,
    paste0("Full-data competitor-aware candidate grid for replicate ", replicate,
      ", direction ", direction, ". Gene FDR is not an eligibility gate; private review only."),
    score_cells = score_cells
  )
  full_native_batches[[batch_id]] <- full_scored$native
  full_replicate_batches[[batch_id]] <- candidate_replicate_summaries(
    full_scored$score, batch_map, score_cells = score_cells, heldout_only = FALSE)
  full_session_workspaces[[batch_id]] <- full_scored$workspace
  try(wb_close_bundle(wb_session_view(full_scored$session)), silent = TRUE)
  rm(full_scored)
  invisible(gc())
}
log_event("Combining full native UCell summaries")
full_native <- bind_rows_fill(full_native_batches, "Full native UCell summary")
log_event("Combining full replicate results")
full_replicate_results <- bind_rows_fill(full_replicate_batches,
  "Full candidate replicate results")
log_event("Summarizing full candidate routes")
full_size_results <- summarize_routes(full_replicate_results, "score_key")
names(full_size_results)[names(full_size_results) == "group_id"] <- "score_key"
full_size_results <- merge(full_map, full_size_results[, setdiff(names(full_size_results),
  c("signature_id", "cytokine", "direction", "signature_size"))], by = "score_key",
  all.x = TRUE, sort = FALSE)
full_size_results <- full_size_results[match(full_map$score_key, full_size_results$score_key), ]
stopifnot(nrow(full_size_results) == nrow(full_map), !anyNA(full_size_results$mean_pass),
  !anyNA(full_size_results$responder_pass))

select_size <- function(local) {
  mean_candidates <- local[local$mean_pass, , drop = FALSE]
  if (nrow(mean_candidates)) {
    chosen <- mean_candidates[order(-mean_candidates$mean_median_standardized_gap,
      mean_candidates$signature_size), , drop = FALSE][1L, ]
    route <- "mean"
    gap <- chosen$mean_median_standardized_gap
  } else {
    responder_candidates <- local[local$responder_pass, , drop = FALSE]
    if (nrow(responder_candidates)) {
      chosen <- responder_candidates[order(-responder_candidates$responder_median_standardized_gap,
        responder_candidates$signature_size), , drop = FALSE][1L, ]
      route <- "responder"
      gap <- chosen$responder_median_standardized_gap
    } else {
      best_failed_gap <- pmax(local$mean_median_standardized_gap,
        local$responder_median_standardized_gap)
      chosen <- local[order(-best_failed_gap, local$signature_size), , drop = FALSE][1L, ]
      route <- if (chosen$mean_median_standardized_gap >=
        chosen$responder_median_standardized_gap) "mean" else "responder"
      gap <- max(chosen$mean_median_standardized_gap,
        chosen$responder_median_standardized_gap)
    }
  }
  chosen$selection_route <- route
  chosen$selection_gap <- gap
  chosen
}

selected_rows <- lapply(signature_info$signature_id, function(id) {
  select_size(full_size_results[full_size_results$signature_id == id, , drop = FALSE])
})
selected <- do.call(rbind, selected_rows)
selected$full_overall_pass <- selected$mean_pass | selected$responder_pass
selected$full_status <- ifelse(selected$full_overall_pass, "PASS", "FAIL")
stopifnot(nrow(selected) == 52L, !anyDuplicated(selected$signature_id),
  all(selected$signature_size %in% candidate_sizes))
full_size_results$selected_size <- full_size_results$score_key %in% selected$score_key

holdout_map_rows <- list()
holdout_member_rows <- list()
holdout_sets <- list()
holdout_counter <- 0L
for (i in seq_len(nrow(selected))) {
  entry <- selected[i, ]
  size <- entry$signature_size
  for (fold in folds) {
    holdout_counter <- holdout_counter + 1L
    ranking <- fold_rankings[[entry$signature_id]][[fold]]
    genes <- ranking$feature_id[seq_len(size)]
    score_key <- paste(entry$signature_id, paste0("size_", size), paste0("fold_", fold), sep = ".")
    holdout_sets[[score_key]] <- genes
    holdout_map_rows[[holdout_counter]] <- data.frame(
      score_key = score_key, signature_id = entry$signature_id,
      cytokine = entry$cytokine, direction = entry$direction,
      signature_size = size, held_out_replicate = fold,
      stringsAsFactors = FALSE
    )
    local <- ranking[seq_len(size), ]
    local$score_key <- score_key
    local$signature_id <- entry$signature_id
    local$cytokine <- entry$cytokine
    local$direction <- entry$direction
    local$signature_size <- size
    local$held_out_replicate <- fold
    holdout_member_rows[[holdout_counter]] <- local[, c("score_key", "signature_id",
      "cytokine", "direction", "signature_size", "held_out_replicate", "rank",
      "feature_id", "specificity_log2fc", "specificity_t", "specificity_p_value",
      "specificity_fdr")]
  }
}
holdout_map <- do.call(rbind, holdout_map_rows)
holdout_membership <- do.call(rbind, holdout_member_rows)
stopifnot(nrow(holdout_map) == 52L * 3L, length(holdout_sets) == nrow(holdout_map),
  identical(names(holdout_sets), holdout_map$score_key))

log_event("Scoring selected-size holdout grid", length(holdout_sets), "sets in fold-specific batches")
holdout_native_batches <- list()
holdout_replicate_batches <- list()
holdout_session_workspaces <- character()
for (fold in folds) {
  use <- holdout_map$held_out_replicate == fold
  batch_map <- holdout_map[use, , drop = FALSE]
  batch_sets <- holdout_sets[batch_map$score_key]
  session_id <- paste0("competitor_aware_selected_holdout_fold_", fold)
  holdout_scored <- open_or_run_ucell(
    session_id,
    batch_sets,
    paste0("Selected full-data size with fold-specific training membership, held-out replicate ",
      fold, ". Stability only; results remain private and are not promoted."),
    score_cells = cells[as.character(md$replicate) == fold]
  )
  holdout_native_batches[[fold]] <- holdout_scored$native
  holdout_replicate_batches[[fold]] <- candidate_replicate_summaries(
    holdout_scored$score, batch_map,
    score_cells = cells[as.character(md$replicate) == fold], heldout_only = TRUE)
  holdout_session_workspaces[[fold]] <- holdout_scored$workspace
  try(wb_close_bundle(wb_session_view(holdout_scored$session)), silent = TRUE)
  rm(holdout_scored)
  invisible(gc())
}
log_event("Combining holdout native UCell summaries")
holdout_native <- bind_rows_fill(holdout_native_batches, "Holdout native UCell summary")
log_event("Combining holdout replicate results")
holdout_replicate_results <- bind_rows_fill(holdout_replicate_batches,
  "Holdout replicate results")
log_event("Summarizing holdout routes")
holdout_summary <- summarize_routes(holdout_replicate_results, "signature_id")
names(holdout_summary)[names(holdout_summary) == "group_id"] <- "holdout_group_id"

jaccard <- function(a, b) {
  union_size <- length(union(a, b))
  if (!union_size) 1 else length(intersect(a, b)) / union_size
}
pairwise_fold_jaccard <- function(sets) {
  pairs <- combn(seq_along(sets), 2L)
  mean(apply(pairs, 2L, function(pair) jaccard(sets[[pair[[1L]]]], sets[[pair[[2L]]]])))
}

final_member_rows <- list()
core_member_rows <- list()
stability_rows <- list()
for (i in seq_len(nrow(selected))) {
  entry <- selected[i, ]
  final_genes <- full_sets[[entry$score_key]]
  fold_genes <- setNames(lapply(folds, function(fold) {
    key <- holdout_map$score_key[holdout_map$signature_id == entry$signature_id &
      holdout_map$held_out_replicate == fold]
    stopifnot(length(key) == 1L)
    holdout_sets[[key]]
  }), folds)
  core_genes <- Reduce(intersect, c(list(final_genes), fold_genes))
  core_genes <- final_genes[final_genes %in% core_genes]
  ranking <- full_rankings[[entry$signature_id]]
  local <- ranking[match(final_genes, ranking$feature_id), ]
  local$signature_id <- entry$signature_id
  local$cytokine <- entry$cytokine
  local$direction <- entry$direction
  local$selected_size <- entry$signature_size
  local$selected_route <- entry$selection_route
  local$full_status <- entry$full_status
  local$selected_in_fold1 <- local$feature_id %in% fold_genes[["1"]]
  local$selected_in_fold2 <- local$feature_id %in% fold_genes[["2"]]
  local$selected_in_fold3 <- local$feature_id %in% fold_genes[["3"]]
  local$fold_selection_count <- local$selected_in_fold1 + local$selected_in_fold2 +
    local$selected_in_fold3
  local$is_core <- entry$full_overall_pass & local$fold_selection_count == 3L
  final_member_rows[[i]] <- local[, c("signature_id", "cytokine", "direction",
    "selected_size", "selected_route", "full_status", "rank", "feature_id",
    "specificity_log2fc", "specificity_t", "specificity_p_value", "specificity_fdr",
    "selected_in_fold1", "selected_in_fold2", "selected_in_fold3",
    "fold_selection_count", "is_core")]
  if (entry$full_overall_pass && length(core_genes)) {
    core_local <- local[local$feature_id %in% core_genes, , drop = FALSE]
    core_local$core_rank <- match(core_local$feature_id, core_genes)
    core_member_rows[[length(core_member_rows) + 1L]] <- core_local[, c("signature_id",
      "cytokine", "direction", "selected_size", "core_rank", "feature_id",
      "specificity_log2fc", "specificity_t", "specificity_p_value", "specificity_fdr")]
  }
  stability_rows[[i]] <- data.frame(
    signature_id = entry$signature_id,
    fold_pairwise_jaccard = pairwise_fold_jaccard(fold_genes),
    final_fold1_jaccard = jaccard(final_genes, fold_genes[["1"]]),
    final_fold2_jaccard = jaccard(final_genes, fold_genes[["2"]]),
    final_fold3_jaccard = jaccard(final_genes, fold_genes[["3"]]),
    core_size = if (entry$full_overall_pass) length(core_genes) else NA_integer_,
    stringsAsFactors = FALSE
  )
}
final_membership <- do.call(rbind, final_member_rows)
core_membership <- if (length(core_member_rows)) do.call(rbind, core_member_rows) else
  data.frame(signature_id = character(), cytokine = character(), direction = character(),
    selected_size = integer(), core_rank = integer(), feature_id = character(),
    specificity_log2fc = numeric(), specificity_t = numeric(),
    specificity_p_value = numeric(), specificity_fdr = numeric())
stability <- do.call(rbind, stability_rows)

holdout_for_merge <- holdout_summary[, c("signature_id", "mean_median_standardized_gap",
  "mean_rank1_replicates", "mean_pass", "mean_competitors",
  "responder_median_standardized_gap", "responder_rank1_replicates",
  "responder_pass", "responder_competitors")]
names(holdout_for_merge)[-1L] <- paste0("holdout_", names(holdout_for_merge)[-1L])
summary_table <- merge(selected, holdout_for_merge, by = "signature_id", all.x = TRUE, sort = FALSE)
summary_table <- merge(summary_table, stability, by = "signature_id", all.x = TRUE, sort = FALSE)
summary_table <- summary_table[match(signature_info$signature_id, summary_table$signature_id), ]
summary_table <- merge(signature_info, summary_table[, setdiff(names(summary_table),
  c("cytokine", "direction"))], by = "signature_id", all.x = TRUE, sort = FALSE)
summary_table <- summary_table[match(signature_info$signature_id, summary_table$signature_id), ]
names(summary_table)[names(summary_table) == "mean_median_standardized_gap"] <-
  "full_mean_median_standardized_gap"
names(summary_table)[names(summary_table) == "mean_rank1_replicates"] <-
  "full_mean_rank1_replicates"
names(summary_table)[names(summary_table) == "mean_pass"] <- "full_mean_pass"
names(summary_table)[names(summary_table) == "mean_competitors"] <- "full_mean_competitors"
names(summary_table)[names(summary_table) == "responder_median_standardized_gap"] <-
  "full_responder_median_standardized_gap"
names(summary_table)[names(summary_table) == "responder_rank1_replicates"] <-
  "full_responder_rank1_replicates"
names(summary_table)[names(summary_table) == "responder_pass"] <- "full_responder_pass"
names(summary_table)[names(summary_table) == "responder_competitors"] <-
  "full_responder_competitors"
stopifnot(nrow(summary_table) == 52L, !anyDuplicated(summary_table$signature_id),
  all(summary_table$full_overall_pass ==
    (summary_table$full_mean_pass | summary_table$full_responder_pass)),
  all(summary_table$core_size[!summary_table$full_overall_pass] %in% NA_integer_))

write_sets <- function(members, path_tsv, path_gmt, description) {
  if (!nrow(members)) {
    write_tsv(data.frame(feature_set = character(), set_name = character(),
      feature_id = character(), description = character()), path_tsv)
    writeLines(character(), path_gmt)
    return(invisible(NULL))
  }
  table <- data.frame(feature_set = members$signature_id,
    set_name = members$signature_id, feature_id = members$feature_id,
    description = description, stringsAsFactors = FALSE)
  write_tsv(table, path_tsv)
  ids <- unique(table$feature_set)
  lines <- vapply(ids, function(id) paste(c(id, description,
    table$feature_id[table$feature_set == id]), collapse = "\t"), character(1))
  writeLines(lines, path_gmt)
}

write_tsv(signature_info, file.path(stage_dir, "signature_inventory.tsv"))
write_tsv(full_size_results, file.path(stage_dir, "all_candidate_size_results.tsv"))
write_tsv(full_replicate_results, file.path(stage_dir, "all_candidate_replicate_results.tsv"))
write_tsv(full_membership, file.path(stage_dir, "all_candidate_membership.tsv"))
write_tsv_gz(full_ranking_evidence, file.path(stage_dir, "full_gene_rankings.tsv.gz"))
write_tsv(summary_table, file.path(stage_dir, "signature_summary.tsv"))
write_tsv(final_membership, file.path(stage_dir, "final_gene_membership.tsv"))
write_tsv(holdout_replicate_results, file.path(stage_dir, "holdout_replicate_results.tsv"))
write_tsv(holdout_membership, file.path(stage_dir, "holdout_gene_membership.tsv"))
write_tsv(core_membership, file.path(stage_dir, "core_gene_membership.tsv"))
write_tsv(full_native, file.path(stage_dir, "native_ucell_full_candidate_summary.tsv"))
write_tsv(holdout_native, file.path(stage_dir, "native_ucell_holdout_summary.tsv"))

passing_members <- final_membership[final_membership$full_status == "PASS", , drop = FALSE]
write_sets(passing_members,
  file.path(stage_dir, "passing_gene_sets.tsv"),
  file.path(stage_dir, "passing_gene_sets.gmt"),
  "Full-data competitor-aware Parse cytokine signature; pending review")
write_sets(core_membership,
  file.path(stage_dir, "passing_core_gene_sets.tsv"),
  file.path(stage_dir, "passing_core_gene_sets.gmt"),
  "Core genes selected in the final and all three holdout-training rankings; pending review")
file.copy(input_paths[["methods"]], file.path(stage_dir, "METHODS.md"), overwrite = FALSE)

input_manifest <- data.frame(
  path = normalizePath(input_paths, mustWork = TRUE),
  sha256 = vapply(input_paths, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))

ordered_summary <- summary_table[order(match(summary_table$cytokine, conditions),
  match(summary_table$direction, c("up", "down"))), ]
report_lines <- c(
  "# Competitor-aware cytokine signatures", "",
  paste0("This private review analyzes all 52 cytokine-direction combinations. Full-data ",
    "performance selects the size and controls the final pass; selected-size leave-one-replicate-out ",
    "performance is reported as stability evidence and never changes the final status."), "",
  paste0("Candidate sizes: `", paste(candidate_sizes, collapse = ", "), "`. Gene eligibility requires ",
    "directional specificity log2FC >= 0.25; no gene-level FDR cutoff is used."), "",
  paste0("Full-data passes: **", sum(ordered_summary$full_overall_pass), "/52**; mean passes: **",
    sum(ordered_summary$full_mean_pass), "**; responder passes: **",
    sum(ordered_summary$full_responder_pass), "**. Holdout-stable by either route: **",
    sum(ordered_summary$holdout_mean_pass | ordered_summary$holdout_responder_pass), "/52**."), "",
  "## Signature summary", "",
  "| Cytokine | Direction | Size | Selected by | Full mean gap | Full responder gap | Mean | Responder | Final | Holdout mean | Holdout responder | Core genes |",
  "|---|---:|---:|---|---:|---:|---|---|---|---|---|---:|"
)
for (i in seq_len(nrow(ordered_summary))) {
  x <- ordered_summary[i, ]
  report_lines <- c(report_lines, paste0(
    "| ", x$cytokine, " | ", x$direction, " | ", x$signature_size, " | ",
    x$selection_route, " | ", fmt_num(x$full_mean_median_standardized_gap), " | ",
    fmt_num(x$full_responder_median_standardized_gap), " | ",
    ifelse(x$full_mean_pass, "PASS", "FAIL"), " | ",
    ifelse(x$full_responder_pass, "PASS", "FAIL"), " | **", x$full_status, "** | ",
    ifelse(x$holdout_mean_pass, "stable", "not stable"), " | ",
    ifelse(x$holdout_responder_pass, "stable", "not stable"), " | ",
    ifelse(is.na(x$core_size), "-", x$core_size), " |"))
}
report_lines <- c(report_lines, "", "## Files", "",
  "- `all_candidate_size_results.tsv`: every tested size, route result, and strongest competitor.",
  "- `signature_summary.tsv`: one row per cytokine-direction combination.",
  "- `final_gene_membership.tsv`: final membership and fold-selection frequency.",
  "- `core_gene_membership.tsv`: stable core membership for full-data passing sets.",
  "- `holdout_replicate_results.tsv`: selected-size cross-fitted stability details.",
  "- `METHODS.md`: exact prespecified method and interpretation limits.", "",
  "## Interpretation", "",
  paste0("The final status is internally optimized on this dataset. The holdout result tests whether ",
    "training-pair gene rankings at the chosen size reproduce in the omitted replicate, but the size itself ",
    "was selected on the full dataset. Exact final memberships therefore require independent validation."))
writeLines(report_lines, file.path(stage_dir, "REPORT.md"))

bundle_fingerprint_after <- wb_manifest(b)$current_state$fingerprint
stopifnot(identical(bundle_fingerprint_before, bundle_fingerprint_after),
  nrow(wb_feature_sets(b, namespace = cfg$namespace)) == 0L,
  nrow(final_membership) == sum(summary_table$signature_size),
  all(core_membership$signature_id %in% summary_table$signature_id[summary_table$full_overall_pass]),
  all(core_membership$feature_id %in% final_membership$feature_id))

completion <- list(
  status = "complete",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  dataset_id = cfg$bundle_id,
  source_bundle = cfg$bundle_dir,
  source_bundle_fingerprint = bundle_fingerprint_before,
  source_bundle_unchanged = identical(bundle_fingerprint_before, bundle_fingerprint_after),
  candidate_directions = nrow(summary_table),
  candidate_sizes = candidate_sizes,
  full_candidate_sets_scored = length(full_sets),
  holdout_sets_scored = length(holdout_sets),
  cells_per_set = length(cells),
  full_passes = sum(summary_table$full_overall_pass),
  full_mean_passes = sum(summary_table$full_mean_pass),
  full_responder_passes = sum(summary_table$full_responder_pass),
  holdout_stable_any_route = sum(summary_table$holdout_mean_pass |
    summary_table$holdout_responder_pass),
  passing_core_sets = length(unique(core_membership$signature_id)),
  passing_sets_with_empty_core = sum(summary_table$full_overall_pass & summary_table$core_size == 0L),
  parameters = list(
    gene_log2fc_minimum = log2fc_minimum,
    gene_fdr_cutoff = NULL,
    ranking = "directional moderated t, directional log2FC, FDR, feature ID",
    max_rank = max_rank,
    responder_up = "strictly above matched-control 95th percentile",
    responder_down = "strictly below matched-control 5th percentile",
    size_selection = "largest full-data median standardized mean gap among mean passes; responder fallback; exact ties use smaller size",
    full_pass = "median standardized gap > 0 for mean or responder",
    holdout_role = "stability annotation only",
    core = "final membership intersected with all three selected-size fold memberships"
  ),
  runtime = runtime,
  full_session_workspaces = as.list(full_session_workspaces),
  holdout_session_workspaces = as.list(holdout_session_workspaces),
  result_promotions = 0L
)
write_json(completion, file.path(stage_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

output_files <- list.files(stage_dir, recursive = TRUE, full.names = TRUE)
output_files <- output_files[basename(output_files) != "output_manifest.tsv"]
output_manifest <- data.frame(
  file = substring(output_files, nchar(stage_dir) + 2L),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(output_manifest[order(output_manifest$file), ], file.path(stage_dir, "output_manifest.tsv"))

if (!file.rename(stage_dir, output_dir)) stop("Failed to finalize output directory")
log_event("Competitor-aware review complete at", output_dir)
