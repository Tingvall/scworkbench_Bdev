#!/usr/bin/env Rscript

# Replicate-aware Parse-only cytokine-signature discovery. This is the
# explicitly approved outside-scWorkbench statistical step: BPCells performs
# exact sample/state count aggregation and limma-voom performs sample-level
# differential modelling. No Bdev or combined-stimulation cells are used.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 08_discover_specific_signatures.R BUNDLE_DIR STATE_REPORT_DIR OUTPUT_DIR", call. = FALSE)
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
state_report_dir <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
if (file.exists(output_dir) || dir.exists(output_dir)) {
  stop("Output destination already exists; refusing to overwrite: ", output_dir, call. = FALSE)
}

options(width = 180, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(BPCells)
  library(limma)
  library(jsonlite)
})

analysis_set <- "set:parse_qc_control_single"
counts_matrix_id <- "RNA_counts"
threads <- as.integer(Sys.getenv("PARSE_THREADS", unset = "8"))
seed <- 1729L
set.seed(seed)

response_logfc_min <- 0.5
specificity_logfc_min <- 0.25
full_fdr_max <- 0.05
fold_p_max <- 0.05
min_signature_genes <- 5L
candidate_sizes <- c(5L, 10L, 15L, 20L, 25L, 30L, 40L, 50L)
min_state_cells_per_sample <- 20L

write_json_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  jsonlite::write_json(value, temporary_path, pretty = TRUE, auto_unbox = TRUE,
                       null = "null", na = "null", digits = NA)
  if (!file.rename(temporary_path, path)) stop("Could not publish JSON: ", path, call. = FALSE)
  invisible(path)
}

write_tsv_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  utils::write.table(value, temporary_path, sep = "\t", quote = TRUE,
                     row.names = FALSE, col.names = TRUE, na = "NA")
  if (!file.rename(temporary_path, path)) stop("Could not publish TSV: ", path, call. = FALSE)
  invisible(path)
}

write_tsv_gz_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  connection <- gzfile(temporary_path, open = "wt", compression = 6L)
  connection_open <- TRUE
  on.exit({
    if (isTRUE(connection_open)) try(close(connection), silent = TRUE)
    if (file.exists(temporary_path)) unlink(temporary_path)
  }, add = TRUE)
  utils::write.table(value, connection, sep = "\t", quote = TRUE,
                     row.names = FALSE, col.names = TRUE, na = "NA")
  close(connection)
  connection_open <- FALSE
  if (!file.rename(temporary_path, path)) stop("Could not publish compressed TSV: ", path, call. = FALSE)
  invisible(path)
}

write_lines_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  writeLines(value, temporary_path, useBytes = TRUE)
  if (!file.rename(temporary_path, path)) stop("Could not publish text: ", path, call. = FALSE)
  invisible(path)
}

slugify <- function(x) {
  value <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("(^_+|_+$)", "", value)
}

markdown_table <- function(value, digits = 4L) {
  if (!nrow(value)) return("_No rows._")
  formatted <- lapply(value, function(x) {
    if (is.numeric(x)) x <- ifelse(is.na(x), "NA", format(signif(x, digits), trim = TRUE))
    x <- ifelse(is.na(x), "NA", as.character(x))
    gsub("|", "\\\\|", x, fixed = TRUE)
  })
  rows <- vapply(seq_len(nrow(value)), function(i) {
    paste0("| ", paste(vapply(formatted, `[[`, character(1L), i), collapse = " | "), " |")
  }, character(1L))
  paste(c(
    paste0("| ", paste(names(value), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(value)), collapse = " | "), " |"),
    rows
  ), collapse = "\n")
}

condition_design <- function(sample_metadata, conditions) {
  design <- matrix(0, nrow = nrow(sample_metadata), ncol = length(conditions),
                   dimnames = list(sample_metadata$sample_id, conditions))
  design[cbind(seq_len(nrow(sample_metadata)), match(sample_metadata$condition, conditions))] <- 1
  design
}

contrast_matrix <- function(conditions, control_condition, single_conditions) {
  out <- matrix(0, nrow = length(conditions), ncol = 2L * length(single_conditions),
                dimnames = list(conditions, character(2L * length(single_conditions))))
  column <- 0L
  for (cytokine in single_conditions) {
    column <- column + 1L
    out[cytokine, column] <- 1
    out[control_condition, column] <- -1
    colnames(out)[column] <- paste0("response::", cytokine)
    column <- column + 1L
    out[cytokine, column] <- 1
    others <- setdiff(single_conditions, cytokine)
    out[others, column] <- -1 / length(others)
    colnames(out)[column] <- paste0("specificity::", cytokine)
  }
  out
}

fit_voom_contrasts <- function(counts, sample_metadata, conditions, contrasts, library_sizes) {
  stopifnot(identical(colnames(counts), sample_metadata$sample_id))
  design <- condition_design(sample_metadata, conditions)
  voomed <- limma::voom(counts, design = design, lib.size = library_sizes,
                        normalize.method = "none", plot = FALSE)
  fitted <- limma::lmFit(voomed, design)
  fitted <- limma::contrasts.fit(fitted, contrasts)
  fitted <- limma::eBayes(fitted, robust = TRUE)
  adjusted <- apply(fitted$p.value, 2L, p.adjust, method = "BH")
  dimnames(adjusted) <- dimnames(fitted$p.value)
  response_columns <- grep("^response::", colnames(fitted$coefficients), value = TRUE)
  absolute_response <- abs(fitted$coefficients[, response_columns, drop = FALSE])
  maximum_response <- apply(absolute_response, 1L, max)
  maximum_index <- max.col(absolute_response, ties.method = "first")
  without_maximum <- absolute_response
  without_maximum[cbind(seq_len(nrow(without_maximum)), maximum_index)] <- -Inf
  second_maximum <- apply(without_maximum, 1L, max)
  off_target_maximum <- matrix(
    maximum_response, nrow = nrow(absolute_response), ncol = ncol(absolute_response),
    dimnames = dimnames(absolute_response)
  )
  off_target_maximum[cbind(seq_len(nrow(off_target_maximum)), maximum_index)] <- second_maximum
  list(
    coefficients = fitted$coefficients,
    t = fitted$t,
    p_value = fitted$p.value,
    adjusted_p_value = adjusted,
    off_target_maximum = off_target_maximum,
    df_residual = fitted$df.residual,
    voom_design = design
  )
}

rank_candidates <- function(stats, target, direction,
                            use_adjusted = FALSE, state_concordance = NULL) {
  response_name <- paste0("response::", target)
  specificity_name <- paste0("specificity::", target)
  response <- stats$coefficients[, response_name]
  specificity <- stats$coefficients[, specificity_name]
  response_t <- stats$t[, response_name]
  specificity_t <- stats$t[, specificity_name]
  probability_response <- if (use_adjusted) stats$adjusted_p_value[, response_name] else stats$p_value[, response_name]
  probability_specificity <- if (use_adjusted) stats$adjusted_p_value[, specificity_name] else stats$p_value[, specificity_name]
  probability_cutoff <- if (use_adjusted) full_fdr_max else fold_p_max
  off_target <- stats$off_target_maximum[, paste0("response::", target)]
  dominance <- abs(response) / (off_target + 0.25)
  eligible <- direction * response >= response_logfc_min &
    direction * specificity >= specificity_logfc_min &
    probability_response <= probability_cutoff &
    probability_specificity <= probability_cutoff
  fallback <- FALSE
  if (!use_adjusted && sum(eligible, na.rm = TRUE) < max(candidate_sizes)) {
    eligible <- direction * response >= response_logfc_min &
      direction * specificity >= specificity_logfc_min &
      probability_response <= 0.20 & probability_specificity <= 0.20
    fallback <- TRUE
  }
  indices <- which(eligible & is.finite(response_t) & is.finite(specificity_t))
  if (!length(indices)) {
    return(list(table = data.frame(), genes = character(), fallback = fallback))
  }
  concordance <- if (is.null(state_concordance)) rep(NA_real_, length(response)) else state_concordance
  concordance_for_score <- ifelse(is.finite(concordance), concordance, 0)
  score <- direction * response_t + direction * specificity_t +
    log2(1 + dominance) + concordance_for_score
  order_index <- indices[order(-score[indices], -dominance[indices], rownames(stats$coefficients)[indices])]
  table <- data.frame(
    feature_id = rownames(stats$coefficients)[order_index],
    response_logFC = response[order_index],
    response_t = response_t[order_index],
    response_p_value = stats$p_value[order_index, response_name],
    response_FDR = stats$adjusted_p_value[order_index, response_name],
    specificity_logFC = specificity[order_index],
    specificity_t = specificity_t[order_index],
    specificity_p_value = stats$p_value[order_index, specificity_name],
    specificity_FDR = stats$adjusted_p_value[order_index, specificity_name],
    maximum_other_response_abs_logFC = off_target[order_index],
    target_dominance_ratio = dominance[order_index],
    state_concordance = concordance[order_index],
    ranking_score = score[order_index],
    stringsAsFactors = FALSE
  )
  list(table = table, genes = table$feature_id, fallback = fallback)
}

auc_binary <- function(score, positive) {
  keep <- is.finite(score) & !is.na(positive)
  score <- score[keep]
  positive <- as.logical(positive[keep])
  n_positive <- sum(positive)
  n_negative <- sum(!positive)
  if (!n_positive || !n_negative) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[positive]) - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)
}

jaccard_three <- function(sets) {
  pairs <- combn(seq_along(sets), 2L)
  mean(apply(pairs, 2L, function(pair) {
    union_size <- length(union(sets[[pair[[1L]]]], sets[[pair[[2L]]]]))
    if (!union_size) return(NA_real_)
    length(intersect(sets[[pair[[1L]]]], sets[[pair[[2L]]]])) / union_size
  }), na.rm = TRUE)
}

main <- function() {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "de_full"), showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Could not create output directory", call. = FALSE)

  state_completion_path <- file.path(state_report_dir, "completion.json")
  if (!file.exists(state_completion_path)) stop("State report is incomplete", call. = FALSE)
  state_completion <- jsonlite::read_json(state_completion_path, simplifyVector = TRUE)
  stopifnot(identical(state_completion$status, "complete"), isFALSE(state_completion$bdev_used))
  cluster_artifact_id <- as.character(state_completion$cluster_artifact_id)

  bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
  on.exit(if (!is.null(bundle)) try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  manifest <- wb_manifest(bundle)
  analysis_cells <- as.character(wb_cells(bundle, cell_set = analysis_set))
  imported_cells <- as.character(wb_cells(bundle, cell_set = "set:parse_cytokines"))
  stopifnot(length(analysis_cells) == 154209L, length(imported_cells) == 251599L)

  metadata <- wb_meta(
    bundle, cells = analysis_cells,
    cols = list(wb_ref_metadata_field("core:cell_id"), "sample_id", "condition",
                "replicate", "treatment_category")
  )
  stopifnot(
    identical(as.character(metadata$cell_id), analysis_cells),
    !anyNA(metadata),
    all(metadata$treatment_category %in% c("control", "single"))
  )
  sample_metadata <- unique(metadata[, c("sample_id", "condition", "replicate", "treatment_category")])
  sample_metadata <- sample_metadata[order(sample_metadata$condition, sample_metadata$replicate), , drop = FALSE]
  rownames(sample_metadata) <- NULL
  stopifnot(
    nrow(sample_metadata) == 81L,
    !anyDuplicated(sample_metadata$sample_id),
    all(table(sample_metadata$condition) == 3L),
    length(unique(sample_metadata$replicate)) == 3L
  )
  control_condition <- unique(sample_metadata$condition[sample_metadata$treatment_category == "control"])
  single_conditions <- sort(unique(sample_metadata$condition[sample_metadata$treatment_category == "single"]))
  conditions <- c(control_condition, single_conditions)
  stopifnot(length(control_condition) == 1L, length(single_conditions) == 26L)
  sample_metadata <- sample_metadata[order(match(sample_metadata$condition, conditions), sample_metadata$replicate), , drop = FALSE]

  clusters <- wb_art_load(bundle, cluster_artifact_id, expected_kind = "cell_vector")
  cluster_index <- match(analysis_cells, as.character(clusters$cell_id))
  stopifnot(!anyNA(cluster_index))
  metadata$state <- as.character(clusters$value[cluster_index])
  states <- sort(unique(metadata$state))
  stopifnot(length(states) >= 2L)

  counts <- wb_matrix(bundle, counts_matrix_id)
  matrix_cell_index <- match(analysis_cells, colnames(counts))
  stopifnot(!anyNA(matrix_cell_index))
  analysis_counts <- counts[, matrix_cell_index, drop = FALSE]
  stopifnot(identical(colnames(analysis_counts), analysis_cells))

  sample_factor <- factor(metadata$sample_id, levels = sample_metadata$sample_id)
  sample_pseudobulk <- BPCells::pseudobulk_matrix(
    analysis_counts, cell_groups = sample_factor, method = "sum", threads = threads
  )
  sample_library_sizes <- as.numeric(colSums(sample_pseudobulk))
  names(sample_library_sizes) <- colnames(sample_pseudobulk)
  stopifnot(identical(colnames(sample_pseudobulk), sample_metadata$sample_id), all(sample_library_sizes > 0))

  state_group <- paste(metadata$sample_id, metadata$state, sep = "||")
  state_group_metadata <- unique(data.frame(
    group = state_group,
    sample_id = as.character(metadata$sample_id),
    state = as.character(metadata$state),
    stringsAsFactors = FALSE
  ))
  state_group_metadata <- merge(
    state_group_metadata, sample_metadata, by = "sample_id", all.x = TRUE, sort = FALSE
  )
  state_group_metadata <- state_group_metadata[
    order(match(state_group_metadata$state, states), match(state_group_metadata$sample_id, sample_metadata$sample_id)),
    , drop = FALSE
  ]
  state_factor <- factor(state_group, levels = state_group_metadata$group)
  state_pseudobulk <- BPCells::pseudobulk_matrix(
    analysis_counts, cell_groups = state_factor, method = "sum", threads = threads
  )
  stopifnot(identical(colnames(state_pseudobulk), state_group_metadata$group))
  state_group_metadata$cells <- as.integer(table(state_factor)[state_group_metadata$group])
  state_library_sizes <- as.numeric(colSums(state_pseudobulk))
  names(state_library_sizes) <- colnames(state_pseudobulk)

  sample_dense_all <- as.matrix(sample_pseudobulk)
  rownames(sample_dense_all) <- rownames(sample_pseudobulk)
  cpm <- 1e6 * sweep(sample_dense_all, 2L, sample_library_sizes, "/")
  expression_keep <- rowSums(cpm >= 1) >= 3L
  technical_exclusion <- grepl("^(mt-|Rpl|Rps)", rownames(sample_dense_all), perl = TRUE)
  analysis_keep <- expression_keep
  signature_eligible_feature <- expression_keep & !technical_exclusion
  stopifnot(sum(analysis_keep) > 1000L, sum(signature_eligible_feature) > 1000L)
  sample_dense <- sample_dense_all[analysis_keep, , drop = FALSE]
  filtered_genes <- rownames(sample_dense)
  signature_eligible_genes <- rownames(sample_dense_all)[signature_eligible_feature]
  sample_logcpm <- log2(1e6 * sweep(sample_dense + 0.5, 2L, sample_library_sizes + 0.5 * nrow(sample_dense), "/"))
  sample_percentile_rank <- apply(sample_logcpm, 2L, function(x) {
    (rank(x, ties.method = "average") - 1) / (length(x) - 1)
  })
  rownames(sample_percentile_rank) <- filtered_genes
  colnames(sample_percentile_rank) <- colnames(sample_dense)
  rm(cpm)

  contrasts <- contrast_matrix(conditions, control_condition, single_conditions)
  full_stats <- fit_voom_contrasts(
    sample_dense,
    sample_metadata[match(colnames(sample_dense), sample_metadata$sample_id), , drop = FALSE],
    conditions, contrasts, sample_library_sizes[colnames(sample_dense)]
  )

  replicate_levels <- sort(unique(as.character(sample_metadata$replicate)))
  fold_stats <- vector("list", length(replicate_levels))
  names(fold_stats) <- replicate_levels
  for (held_out in replicate_levels) {
    training_samples <- sample_metadata$sample_id[sample_metadata$replicate != held_out]
    training_metadata <- sample_metadata[match(training_samples, sample_metadata$sample_id), , drop = FALSE]
    fold_stats[[held_out]] <- fit_voom_contrasts(
      sample_dense[, training_samples, drop = FALSE], training_metadata,
      conditions, contrasts, sample_library_sizes[training_samples]
    )
  }

  state_dense <- as.matrix(state_pseudobulk[analysis_keep, , drop = FALSE])
  rownames(state_dense) <- filtered_genes
  state_logcpm <- log2(1e6 * sweep(state_dense + 0.5, 2L,
                                   state_library_sizes + 0.5 * nrow(state_dense), "/"))
  state_response <- vector("list", length(states))
  state_specificity <- vector("list", length(states))
  names(state_response) <- names(state_specificity) <- states
  for (state in states) {
    columns <- which(state_group_metadata$state == state &
                       state_group_metadata$cells >= min_state_cells_per_sample)
    local_metadata <- state_group_metadata[columns, , drop = FALSE]
    condition_means <- matrix(NA_real_, nrow = nrow(state_logcpm), ncol = length(conditions),
                              dimnames = list(filtered_genes, conditions))
    for (condition in conditions) {
      condition_columns <- columns[local_metadata$condition == condition]
      if (length(condition_columns) >= 2L) {
        condition_means[, condition] <- rowMeans(state_logcpm[, condition_columns, drop = FALSE])
      }
    }
    response <- matrix(NA_real_, nrow = nrow(condition_means), ncol = length(single_conditions),
                       dimnames = list(filtered_genes, single_conditions))
    specificity <- response
    for (cytokine in single_conditions) {
      response[, cytokine] <- condition_means[, cytokine] - condition_means[, control_condition]
      other_means <- condition_means[, setdiff(single_conditions, cytokine), drop = FALSE]
      specificity[, cytokine] <- condition_means[, cytokine] - rowMeans(other_means, na.rm = TRUE)
      specificity[rowSums(is.finite(other_means)) < 10L, cytokine] <- NA_real_
    }
    state_response[[state]] <- response
    state_specificity[[state]] <- specificity
  }

  fold_rankings <- list()
  fold_diagnostic_rows <- list()
  fold_membership_rows <- list()
  oof_rows <- list()
  row_counter <- 0L
  eligible_gene_flag <- filtered_genes %in% signature_eligible_genes
  for (cytokine in single_conditions) {
    for (direction_label in c("up", "down")) {
      direction <- if (direction_label == "up") 1 else -1
      signature_key <- paste(cytokine, direction_label, sep = "::")
      fold_rankings[[signature_key]] <- list()
      for (held_out in replicate_levels) {
        ranked <- rank_candidates(fold_stats[[held_out]], cytokine, direction,
                                  use_adjusted = FALSE)
        ranked$genes <- ranked$genes[ranked$genes %in% signature_eligible_genes]
        fold_rankings[[signature_key]][[held_out]] <- ranked$genes
        fold_diagnostic_rows[[length(fold_diagnostic_rows) + 1L]] <- data.frame(
          cytokine = cytokine,
          direction = direction_label,
          held_out_replicate = held_out,
          ranked_training_genes = length(ranked$genes),
          relaxed_training_p_value_to_0_20 = ranked$fallback,
          stringsAsFactors = FALSE
        )
        held_samples <- sample_metadata$sample_id[sample_metadata$replicate == held_out]
        available_sizes <- candidate_sizes[candidate_sizes <= length(ranked$genes)]
        for (size in available_sizes) {
          genes <- utils::head(ranked$genes, size)
          gene_indices <- match(genes, filtered_genes)
          raw_score <- colMeans(sample_percentile_rank[gene_indices, held_samples, drop = FALSE])
          oriented_score <- direction * raw_score
          held_metadata <- sample_metadata[match(held_samples, sample_metadata$sample_id), , drop = FALSE]
          row_counter <- row_counter + 1L
          oof_rows[[row_counter]] <- data.frame(
            cytokine = cytokine,
            direction = direction_label,
            signature_size = size,
            held_out_replicate = held_out,
            sample_id = held_samples,
            condition = held_metadata$condition,
            treatment_category = held_metadata$treatment_category,
            raw_rank_score = as.numeric(raw_score),
            oriented_rank_score = as.numeric(oriented_score),
            stringsAsFactors = FALSE
          )
          fold_membership_rows[[length(fold_membership_rows) + 1L]] <- data.frame(
            cytokine = cytokine, direction = direction_label,
            signature_size = size, held_out_replicate = held_out,
            feature_id = genes, stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  oof_scores <- do.call(rbind, oof_rows)
  fold_membership <- do.call(rbind, fold_membership_rows)
  fold_diagnostics <- do.call(rbind, fold_diagnostic_rows)

  performance_rows <- list()
  performance_counter <- 0L
  for (cytokine in single_conditions) {
    for (direction_label in c("up", "down")) {
      signature_key <- paste(cytokine, direction_label, sep = "::")
      for (size in candidate_sizes) {
        rows <- oof_scores$cytokine == cytokine & oof_scores$direction == direction_label &
          oof_scores$signature_size == size
        values <- oof_scores[rows, , drop = FALSE]
        if (nrow(values) != 81L) next
        response_rows <- values$condition %in% c(control_condition, cytokine)
        auc_response <- auc_binary(values$oriented_rank_score[response_rows],
                                   values$condition[response_rows] == cytokine)
        specificity_rows <- values$treatment_category == "single"
        auc_specificity <- auc_binary(values$oriented_rank_score[specificity_rows],
                                      values$condition[specificity_rows] == cytokine)
        sets <- lapply(replicate_levels, function(fold) {
          utils::head(fold_rankings[[signature_key]][[fold]], size)
        })
        performance_counter <- performance_counter + 1L
        performance_rows[[performance_counter]] <- data.frame(
          cytokine = cytokine,
          direction = direction_label,
          signature_size = size,
          auc_response_vs_control = auc_response,
          auc_specificity_vs_other_singles = auc_specificity,
          balanced_auc = mean(c(auc_response, auc_specificity)),
          fold_gene_jaccard = jaccard_three(sets),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  cv_performance <- do.call(rbind, performance_rows)

  final_candidate_rows <- list()
  signature_rows <- list()
  member_rows <- list()
  off_target_rows <- list()
  candidate_counter <- signature_counter <- member_counter <- off_counter <- 0L
  for (cytokine in single_conditions) {
    for (direction_label in c("up", "down")) {
      direction <- if (direction_label == "up") 1 else -1
      state_concordance <- rep(NA_real_, length(filtered_genes))
      best_state <- rep(NA_character_, length(filtered_genes))
      state_metric <- matrix(NA_real_, nrow = length(filtered_genes), ncol = length(states),
                             dimnames = list(filtered_genes, states))
      for (state in states) {
        response <- state_response[[state]][, cytokine]
        specificity <- state_specificity[[state]][, cytokine]
        state_metric[, state] <- pmin(direction * response, direction * specificity)
      }
      state_concordance <- rowMeans(state_metric > 0, na.rm = TRUE)
      state_concordance[rowSums(is.finite(state_metric)) == 0L] <- NA_real_
      finite_state_count <- rowSums(is.finite(state_metric))
      state_metric_for_max <- replace(state_metric, !is.finite(state_metric), -Inf)
      best_state <- states[max.col(state_metric_for_max, ties.method = "first")]
      best_state[finite_state_count == 0L] <- NA_character_

      ranked <- rank_candidates(full_stats, cytokine, direction,
                                use_adjusted = TRUE, state_concordance = state_concordance)
      if (nrow(ranked$table)) {
        ranked$table <- ranked$table[ranked$table$feature_id %in% signature_eligible_genes, , drop = FALSE]
      }
      fold_direction_stable <- if (nrow(ranked$table)) {
        vapply(ranked$table$feature_id, function(gene) {
          all(vapply(fold_stats, function(fold) {
            direction * fold$coefficients[gene, paste0("response::", cytokine)] > 0 &&
              direction * fold$coefficients[gene, paste0("specificity::", cytokine)] > 0
          }, logical(1L)))
        }, logical(1L))
      } else logical()
      if (nrow(ranked$table)) {
        ranked$table$fold_direction_stable <- fold_direction_stable
        ranked$table$best_state <- best_state[match(ranked$table$feature_id, filtered_genes)]
        ranked$table <- ranked$table[ranked$table$fold_direction_stable, , drop = FALSE]
        ranked$table$rank <- seq_len(nrow(ranked$table))
      }

      performance <- cv_performance[
        cv_performance$cytokine == cytokine & cv_performance$direction == direction_label &
          cv_performance$signature_size <= nrow(ranked$table), , drop = FALSE
      ]
      status <- "insufficient_genes"
      selected_size <- NA_integer_
      selected_performance <- NULL
      if (nrow(ranked$table) >= min_signature_genes && nrow(performance)) {
        best_auc <- max(performance$balanced_auc, na.rm = TRUE)
        near_best <- performance[performance$balanced_auc >= best_auc - 0.02, , drop = FALSE]
        near_best <- near_best[order(-near_best$fold_gene_jaccard, near_best$signature_size), , drop = FALSE]
        selected_performance <- near_best[1L, , drop = FALSE]
        selected_size <- as.integer(selected_performance$signature_size)
        validated <- selected_performance$auc_response_vs_control >= 0.60 &&
          selected_performance$auc_specificity_vs_other_singles >= 0.60 &&
          selected_performance$balanced_auc >= 0.65 &&
          is.finite(selected_performance$fold_gene_jaccard) &&
          selected_performance$fold_gene_jaccard >= 0.10
        status <- if (validated) "validated" else "failed_held_out_validation"
      }

      signature_id <- paste0("parse_cytokine.", slugify(cytokine), ".", direction_label)
      selected_genes <- if (identical(status, "validated")) {
        utils::head(ranked$table$feature_id, selected_size)
      } else character()
      if (nrow(ranked$table)) {
        ranked$table$cytokine <- cytokine
        ranked$table$direction <- direction_label
        ranked$table$signature_id <- signature_id
        ranked$table$selected <- ranked$table$feature_id %in% selected_genes
        final_candidate_rows[[length(final_candidate_rows) + 1L]] <- ranked$table
      }

      selected_state_concordance <- if (!length(selected_genes)) numeric() else
        ranked$table$state_concordance[match(selected_genes, ranked$table$feature_id)]
      median_state_concordance <- if (!length(selected_state_concordance) ||
                                      !any(is.finite(selected_state_concordance))) {
        NA_real_
      } else {
        median(selected_state_concordance, na.rm = TRUE)
      }
      signature_counter <- signature_counter + 1L
      signature_rows[[signature_counter]] <- data.frame(
        cytokine = cytokine,
        direction = direction_label,
        signature_id = signature_id,
        status = status,
        full_eligible_genes = nrow(ranked$table),
        selected_genes = length(selected_genes),
        selected_size = selected_size,
        auc_response_vs_control = if (is.null(selected_performance)) NA_real_ else selected_performance$auc_response_vs_control,
        auc_specificity_vs_other_singles = if (is.null(selected_performance)) NA_real_ else selected_performance$auc_specificity_vs_other_singles,
        balanced_auc = if (is.null(selected_performance)) NA_real_ else selected_performance$balanced_auc,
        fold_gene_jaccard = if (is.null(selected_performance)) NA_real_ else selected_performance$fold_gene_jaccard,
        program_scope = if (!length(selected_genes) || !is.finite(median_state_concordance)) {
          NA_character_
        } else if (median_state_concordance >= 0.5) {
          "broad_across_states"
        } else {
          "state_enriched"
        },
        stringsAsFactors = FALSE
      )

      if (length(selected_genes)) {
        selected_table <- ranked$table[match(selected_genes, ranked$table$feature_id), , drop = FALSE]
        selected_table$member_rank <- seq_len(nrow(selected_table))
        member_rows[[length(member_rows) + 1L]] <- selected_table
        for (condition in single_conditions) {
          off_counter <- off_counter + 1L
          effect <- full_stats$coefficients[selected_genes, paste0("response::", condition)]
          off_target_rows[[off_counter]] <- data.frame(
            signature_id = signature_id,
            cytokine = cytokine,
            direction = direction_label,
            evaluated_condition = condition,
            is_target = condition == cytokine,
            mean_oriented_response_logFC = mean(direction * effect),
            median_oriented_response_logFC = median(direction * effect),
            fraction_genes_expected_direction = mean(direction * effect > 0),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    response_name <- paste0("response::", cytokine)
    specificity_name <- paste0("specificity::", cytokine)
    de_table <- data.frame(
      feature_id = filtered_genes,
      response_logFC = full_stats$coefficients[, response_name],
      response_t = full_stats$t[, response_name],
      response_p_value = full_stats$p_value[, response_name],
      response_FDR = full_stats$adjusted_p_value[, response_name],
      specificity_logFC = full_stats$coefficients[, specificity_name],
      specificity_t = full_stats$t[, specificity_name],
      specificity_p_value = full_stats$p_value[, specificity_name],
      specificity_FDR = full_stats$adjusted_p_value[, specificity_name],
      signature_eligible_feature = eligible_gene_flag,
      stringsAsFactors = FALSE
    )
    write_tsv_gz_atomic(de_table, file.path(output_dir, "de_full", paste0(slugify(cytokine), ".tsv.gz")))
  }

  signature_summary <- do.call(rbind, signature_rows)
  final_candidates <- if (length(final_candidate_rows)) do.call(rbind, final_candidate_rows) else data.frame()
  signature_members <- if (length(member_rows)) do.call(rbind, member_rows) else data.frame()
  off_target_profile <- if (length(off_target_rows)) do.call(rbind, off_target_rows) else data.frame()
  validated_summary <- signature_summary[signature_summary$status == "validated", , drop = FALSE]
  if (!nrow(validated_summary)) stop("No signatures passed held-out validation", call. = FALSE)

  feature_set_table <- signature_members[, c("signature_id", "feature_id"), drop = FALSE]
  names(feature_set_table)[[1L]] <- "feature_set"
  feature_set_table$set_name <- feature_set_table$feature_set
  description_map <- setNames(
    paste0("Parse individual-cytokine ", validated_summary$direction,
           " signature for ", validated_summary$cytokine,
           "; replicate-aware response and specificity contrasts with held-out validation."),
    validated_summary$signature_id
  )
  feature_set_table$description <- unname(description_map[feature_set_table$feature_set])
  feature_set_table <- feature_set_table[, c("feature_set", "set_name", "feature_id", "description")]

  gmt_lines <- vapply(unique(feature_set_table$feature_set), function(set_id) {
    genes <- feature_set_table$feature_id[feature_set_table$feature_set == set_id]
    paste(c(set_id, description_map[[set_id]], genes), collapse = "\t")
  }, character(1L))

  write_tsv_atomic(sample_metadata, file.path(output_dir, "sample_metadata.tsv"))
  write_tsv_atomic(state_group_metadata, file.path(output_dir, "pseudobulk_state_cell_counts.tsv"))
  write_tsv_atomic(data.frame(
    feature_id = rownames(sample_dense_all),
    cpm_at_least_1_in_three_samples = expression_keep,
    excluded_from_signatures_as_mitochondrial_or_ribosomal = technical_exclusion,
    used_for_modelling = analysis_keep,
    eligible_for_signature_membership = signature_eligible_feature,
    stringsAsFactors = FALSE
  ), file.path(output_dir, "feature_filter.tsv"))
  write_tsv_atomic(cv_performance, file.path(output_dir, "cross_validation_performance.tsv"))
  write_tsv_gz_atomic(oof_scores, file.path(output_dir, "cross_validation_oof_scores.tsv.gz"))
  write_tsv_gz_atomic(fold_membership, file.path(output_dir, "cross_validation_fold_membership.tsv.gz"))
  write_tsv_atomic(fold_diagnostics, file.path(output_dir, "cross_validation_fold_diagnostics.tsv"))
  write_tsv_atomic(signature_summary, file.path(output_dir, "signature_summary.tsv"))
  write_tsv_gz_atomic(final_candidates, file.path(output_dir, "eligible_signature_candidates.tsv.gz"))
  write_tsv_atomic(signature_members, file.path(output_dir, "selected_signature_members.tsv"))
  write_tsv_atomic(off_target_profile, file.path(output_dir, "selected_signature_offtarget_profile.tsv"))
  write_tsv_atomic(feature_set_table, file.path(output_dir, "validated_feature_sets.tsv"))
  write_lines_atomic(gmt_lines, file.path(output_dir, "validated_feature_sets.gmt"))

  saveRDS(list(
    counts = sample_dense_all,
    sample_metadata = sample_metadata,
    state_counts = state_dense,
    state_group_metadata = state_group_metadata
  ), file.path(output_dir, "pseudobulk_counts.rds"), compress = "gzip")
  saveRDS(list(full = full_stats, folds = fold_stats, contrasts = contrasts,
               filtered_genes = filtered_genes),
          file.path(output_dir, "limma_voom_fits.rds"), compress = "gzip")

  write_json_atomic(list(
    analysis_scope = list(
      bundle_id = "parse_cytokines",
      bundle_state_fingerprint = manifest$current_state$fingerprint,
      input_cell_set = analysis_set,
      analysis_cells = length(analysis_cells),
      imported_cells_retained = length(imported_cells),
      samples = nrow(sample_metadata),
      control_samples = sum(sample_metadata$treatment_category == "control"),
      single_cytokine_samples = sum(sample_metadata$treatment_category == "single"),
      single_cytokines = length(single_conditions),
      combined_cells_used = FALSE,
      bdev_used = FALSE,
      cluster_artifact_id = cluster_artifact_id,
      state_report_dir = state_report_dir
    ),
    aggregation = list(
      method = "exact raw-count sum with BPCells::pseudobulk_matrix",
      primary_unit = "biological sample",
      state_unit = "biological sample by Parse-only Leiden state",
      minimum_cells_per_sample_state_for_state_effect = min_state_cells_per_sample
    ),
    modelling = list(
      method = "limma-voom with robust empirical Bayes",
      count_normalization = "library-size CPM; voom normalize.method=none",
      gene_filter = "CPM >= 1 in at least 3 of 81 samples",
      signature_technical_exclusions = "feature IDs beginning mt-, Rpl, or Rps",
      response_contrast = "target single cytokine minus unstimulated control",
      specificity_contrast = "target single cytokine minus the equal-weight mean of the other 25 single cytokines",
      response_absolute_log2fc_minimum = response_logfc_min,
      specificity_absolute_log2fc_minimum = specificity_logfc_min,
      full_data_fdr_maximum_for_both_contrasts = full_fdr_max
    ),
    validation = list(
      folds = "leave one replicate label out; 3 folds",
      candidate_sizes = candidate_sizes,
      score = "mean within-sample expression percentile rank on held-out sample pseudobulk",
      primary_unit = "held-out biological sample",
      selection_rule = "highest fold-membership Jaccard among sizes within 0.02 of maximum balanced AUC; smaller size breaks ties",
      acceptance = "response AUC >= 0.60, specificity AUC >= 0.60, balanced AUC >= 0.65, fold-gene Jaccard >= 0.10",
      minimum_signature_genes = min_signature_genes
    ),
    seed = seed
  ), file.path(output_dir, "analysis_parameters.json"))

  status_counts <- as.data.frame(table(status = signature_summary$status), stringsAsFactors = FALSE)
  direction_counts <- aggregate(
    rep(1L, nrow(validated_summary)),
    by = list(direction = validated_summary$direction), sum
  )
  names(direction_counts)[[2L]] <- "validated_signatures"
  report <- paste0(
    "# Parse individual-cytokine signature discovery\n\n",
    "This report covers discovery and held-out validation before final scWorkbench UCell scoring. ",
    "The analysis used only `", analysis_set, "` (", format(length(analysis_cells), big.mark = ","),
    " QC-accepted control/single cells; 81 samples, 26 individual cytokines plus control). ",
    "Combined-stimulation cells and Bdev were not used. All ", format(length(imported_cells), big.mark = ","),
    " imported cells remain retained in the source bundle.\n\n",
    "## Statistical design\n\n",
    "Raw counts were summed by biological sample. limma-voom tested each cytokine against control and, separately, ",
    "against the equal-weight mean of the other 25 single-cytokine conditions. Candidate genes had to pass both ",
    "contrasts, effect-size cutoffs, FDR 0.05, and direction stability in all three leave-one-replicate-out fits. ",
    "Mitochondrial and ribosomal-prefix genes were excluded from signature membership.\n\n",
    "Compact size was chosen from 5-50 genes by sample-level held-out rank enrichment and fold-membership stability. A direction was accepted only ",
    "with response AUC >=0.60, specificity AUC >=0.60, balanced AUC >=0.65, and fold-gene Jaccard >=0.10; unsupported directions were retained ",
    "in `signature_summary.tsv` without a forced gene set.\n\n",
    "## Outcomes\n\n",
    markdown_table(status_counts), "\n\n",
    markdown_table(direction_counts), "\n\n",
    "The exact validated memberships are in `validated_feature_sets.tsv` and `validated_feature_sets.gmt`. ",
    "Complete response/specificity statistics are in `de_full/`; cross-validation scores, fold memberships, ",
    "off-target profiles, state robustness annotations, parameters, and serialized pseudobulk/fit objects are retained beside this report.\n"
  )
  write_lines_atomic(report, file.path(output_dir, "signature_discovery_report.md"))
  capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
  write_json_atomic(list(
    status = "complete",
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    validated_signatures = nrow(validated_summary),
    validated_up_signatures = sum(validated_summary$direction == "up"),
    validated_down_signatures = sum(validated_summary$direction == "down"),
    cytokines_with_validated_up_signature = length(unique(validated_summary$cytokine[validated_summary$direction == "up"])),
    cytokines_with_any_validated_signature = length(unique(validated_summary$cytokine)),
    combined_cells_used = FALSE,
    bdev_used = FALSE,
    feature_set_table = file.path(output_dir, "validated_feature_sets.tsv")
  ), file.path(output_dir, "completion.json"))
}

main()
