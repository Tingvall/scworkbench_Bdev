#!/usr/bin/env Rscript
# Re-run membership selection and held-out scoring using only saved specificity
# statistics from the previously authorized limma-voom discovery workflow.
source("src/parse_rest_only_workflow.R")
main <- function() {
  stopifnot(!file.exists(file.path(report, "discovery.json")))
  parent <- parent_check()
  dir.create(discovery_dir, recursive = TRUE, showWarnings = FALSE)
  source_dir <- parent_cfg$signature_analysis$discovery$report
  paths <- file.path(source_dir, c("limma_voom_fits.rds", "pseudobulk_counts.rds", "feature_filter.tsv"))
  log_event("Loading saved specificity statistics and exact pseudobulk counts")
  fits <- readRDS(paths[1]); pb <- readRDS(paths[2]); filters <- read_tsv(paths[3])
  md <- pb$sample_metadata
  stopifnot(identical(colnames(pb$counts), md$sample_id), !anyDuplicated(md$sample_id), nrow(md) == 81L)
  singles <- sort(unique(md$condition[md$treatment_category == "single"]))
  folds <- sort(unique(as.character(md$replicate)))
  features <- fits$filtered_genes
  allowed <- filters$feature_id[filters$eligible_for_signature_membership]
  stopifnot(length(singles) == 26L, length(folds) == 3L,
    all(table(md$condition) == 3L), identical(features, rownames(fits$full$coefficients)))
  # Delete every control contrast and dominance quantity from the working fit
  # objects, so subsequent selection code can only access specificity evidence.
  columns <- paste0("specificity::", singles)
  reduce_fit <- function(x) lapply(x[c("coefficients", "t", "p_value", "adjusted_p_value")], function(m) m[, columns, drop = FALSE])
  full <- reduce_fit(fits$full); training <- lapply(fits$folds, reduce_fit); rm(fits)
  libraries <- colSums(pb$counts)
  counts <- pb$counts[match(features, rownames(pb$counts)), , drop = FALSE]
  logcpm <- log2(1e6 * sweep(counts + .5, 2, libraries + .5 * nrow(counts), "/"))
  ranks <- apply(logcpm, 2, function(x) (rank(x, ties.method = "average") - 1) / (length(x) - 1))
  rownames(ranks) <- features; colnames(ranks) <- md$sample_id
  stopifnot(all(is.finite(ranks)), all(libraries > 0))

  # Specificity-only state concordance replaces the old two-contrast term.
  # Use total state library sizes from recorded cell/sample pseudobulk evidence:
  # the saved state count matrix contains the modelled genes only, as in stage 08.
  state_md <- pb$state_group_metadata
  stopifnot(identical(colnames(pb$state_counts), state_md$group))
  # Reconstruct the original per-state library sizes from the parent raw counts;
  # do not normalize state evidence using only modelled genes.
  b <- wb_open_bundle(parent$bundle_dir, read_only = TRUE, validate = FALSE)
  on.exit(wb_close_bundle(b), add = TRUE)
  cells <- wb_cells(b, cfg$source_cell_set)
  cell_md <- wb_meta(b, cells = cells, cols = list(wb_ref_metadata_field("core:cell_id"), "sample_id"))
  state_artifact <- parent_cfg$signature_analysis$state_construction$cluster_artifact_id
  state <- wb_art_load(b, state_artifact, expected_kind = "cell_vector")
  state_values <- state$value[match(cells, state$cell_id)]
  stopifnot(length(cells) == 154209L, identical(cell_md$cell_id, cells), !anyNA(state_values))
  raw <- wb_matrix(b, "RNA_counts")
  raw <- raw[, match(cells, colnames(raw)), drop = FALSE]
  totals <- as.numeric(colSums(raw))
  group <- paste(cell_md$sample_id, state_values, sep = "||")
  state_libraries <- tapply(totals, group, sum)[state_md$group]
  stopifnot(!anyNA(state_libraries), all(state_libraries > 0), sum(state_md$cells) == length(cells))
  state_log <- log2(1e6 * sweep(pb$state_counts + .5, 2, state_libraries + .5 * nrow(pb$state_counts), "/"))
  state_effects <- lapply(sort(unique(state_md$state)), function(st) {
    means <- sapply(singles, function(condition) {
      use <- state_md$state == st & state_md$condition == condition & state_md$cells >= 20L
      if (sum(use) >= 2L) rowMeans(state_log[, use, drop = FALSE]) else rep(NA_real_, nrow(state_log))
    })
    out <- means
    for (condition in singles) {
      other <- means[, setdiff(singles, condition), drop = FALSE]
      out[, condition] <- means[, condition] - rowMeans(other, na.rm = TRUE)
      out[rowSums(is.finite(other)) < 10L, condition] <- NA_real_
    }
    out
  })
  stopifnot(identical(rownames(pb$state_counts), features))
  rank_genes <- function(stats, target, sign, adjusted, concordance = NULL) {
    col <- paste0("specificity::", target)
    prob <- stats[[if (adjusted) "adjusted_p_value" else "p_value"]][, col]
    cutoff <- if (adjusted) cfg$analysis$gene_fdr_maximum else cfg$analysis$training_gene_p_maximum
    use <- sign * stats$coefficients[, col] >= cfg$analysis$gene_log2fc_minimum & prob <= cutoff
    fallback <- FALSE
    if (!adjusted && sum(use) < max(cfg$analysis$candidate_sizes)) {
      use <- sign * stats$coefficients[, col] >= cfg$analysis$gene_log2fc_minimum & prob <= cfg$analysis$training_gene_p_fallback
      fallback <- TRUE
    }
    use <- use & features %in% allowed & is.finite(stats$t[, col])
    if (is.null(concordance)) concordance <- rep(0, length(features))
    score <- sign * stats$t[, col] + ifelse(is.finite(concordance), concordance, 0)
    ids <- which(use); ids <- ids[order(-score[ids], features[ids])]
    list(genes = features[ids], score = score, fallback = fallback)
  }
  summary_rows <- member_rows <- eligible_rows <- cv_rows <- oof_rows <- fold_member_rows <- diagnostics <- list()
  for (target in singles) for (direction in c("up", "down")) {
    sign <- if (direction == "up") 1 else -1
    id <- paste(cfg$namespace, slug(target), direction, sep = ".")
    col <- paste0("specificity::", target)
    effects <- sapply(state_effects, function(m) m[, target])
    concordance <- rowMeans(sign * effects > 0, na.rm = TRUE)
    concordance[!is.finite(concordance)] <- NA_real_
    ranked <- rank_genes(full, target, sign, TRUE, concordance)
    stable <- Reduce(`&`, lapply(training, function(f) sign * f$coefficients[, col] > 0))
    genes <- ranked$genes[stable[match(ranked$genes, features)]]
    fold_sets <- setNames(lapply(folds, function(fold) {
      x <- rank_genes(training[[fold]], target, sign, FALSE)
      diagnostics[[length(diagnostics) + 1L]] <<- data.frame(signature_id = id, held_out_replicate = fold,
        available_genes = length(x$genes), relaxed_p_to_0_20 = x$fallback)
      x$genes
    }), folds)
    perf <- list()
    for (size in cfg$analysis$candidate_sizes) {
      if (any(lengths(fold_sets) < size)) next
      fold_scores <- list()
      for (fold in folds) {
        local_genes <- head(fold_sets[[fold]], size)
        use <- md$replicate == fold & md$treatment_category == "single"
        values <- sign * colMeans(ranks[local_genes, md$sample_id[use], drop = FALSE])
        fold_scores[[fold]] <- data.frame(signature_id = id, signature_size = size, held_out_replicate = fold,
          sample_id = md$sample_id[use], condition = md$condition[use], oriented_rank_score = values)
        fold_member_rows[[length(fold_member_rows) + 1L]] <- data.frame(signature_id = id,
          signature_size = size, held_out_replicate = fold, feature_id = local_genes)
      }
      local <- do.call(rbind, fold_scores)
      stopifnot(nrow(local) == 78L, !anyDuplicated(local$sample_id), all(local$condition != "control"))
      oof_rows[[length(oof_rows) + 1L]] <- local
      row <- data.frame(signature_id = id, signature_size = size,
        cv_rest_auc = auc(local$oriented_rank_score, local$condition == target),
        fold_gene_jaccard = jaccard(lapply(fold_sets, head, size)))
      perf[[length(perf) + 1L]] <- row
    }
    performance <- if (length(perf)) do.call(rbind, perf) else data.frame()
    if (nrow(performance)) cv_rows[[length(cv_rows) + 1L]] <- performance
    valid <- if (nrow(performance)) performance[performance$signature_size <= length(genes), ] else performance
    decision <- "insufficient_genes"; chosen <- NULL; selected_genes <- character()
    if (length(genes) >= cfg$analysis$minimum_genes) {
      decision <- "insufficient_training_genes"
      if (nrow(valid)) {
        near <- valid[valid$cv_rest_auc >= max(valid$cv_rest_auc) - .02, ]
        chosen <- near[order(-near$fold_gene_jaccard, near$signature_size), ][1, ]
        pass <- chosen$cv_rest_auc >= cfg$analysis$discovery_rest_auc_minimum && chosen$fold_gene_jaccard >= cfg$analysis$fold_gene_jaccard_minimum
        decision <- if (pass) "discovery_pass" else "discovery_fail"
        if (pass) selected_genes <- head(genes, chosen$signature_size)
      }
    }
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(signature_id = id, cytokine = target, direction = direction,
      eligible_genes = length(genes), selected_genes = length(selected_genes),
      cv_rest_auc = if (is.null(chosen)) NA_real_ else chosen$cv_rest_auc,
      fold_gene_jaccard = if (is.null(chosen)) NA_real_ else chosen$fold_gene_jaccard,
      discovery_status = decision)
    if (length(genes)) {
      j <- match(genes, features)
      table <- data.frame(signature_id = id, cytokine = target, direction = direction, feature_id = genes,
        specificity_logFC = full$coefficients[j, col], specificity_t = full$t[j, col],
        specificity_p_value = full$p_value[j, col], specificity_FDR = full$adjusted_p_value[j, col],
        specificity_state_concordance = concordance[j], ranking_score = ranked$score[j],
        all_three_folds_same_direction = stable[j], rank = seq_along(genes), selected = genes %in% selected_genes)
      eligible_rows[[length(eligible_rows) + 1L]] <- table
      if (length(selected_genes)) member_rows[[length(member_rows) + 1L]] <- table[table$selected, ]
    }
    log_event(target, direction, "eligible", length(genes), "selected", length(selected_genes), decision)
  }
  summary <- do.call(rbind, summary_rows); members <- do.call(rbind, member_rows)
  stopifnot(nrow(summary) == 52L, sum(summary$eligible_genes >= 5L) == 51L, nrow(members) > 0L,
    all(members$specificity_FDR <= .05), all(members$all_three_folds_same_direction))
  write_tsv(summary, "discovery_summary.tsv")
  write_tsv(members, "selected_gene_evidence.tsv")
  write_tsv(do.call(rbind, cv_rows), "cross_validation_performance.tsv")
  write_tsv(do.call(rbind, diagnostics), "fold_diagnostics.tsv")
  for (name in c("eligible_rows", "oof_rows", "fold_member_rows")) {
    con <- gzfile(file.path(report, paste0(name, ".tsv.gz")), "wt")
    write.table(do.call(rbind, get(name)), con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA"); close(con)
  }
  write_sets(members, summary$signature_id[summary$discovery_status == "discovery_pass"], "discovery_gene_sets")
  saveRDS(list(summary = summary, members = members), file.path(discovery_dir, "discovery.rds"))
  save_provenance("discovery", c(paths, "processing/parse_cytokines/11_discover_rest_only_signatures.R"),
    list(candidate_directions = 52L, sufficient_genes = 51L, discovery_pass = sum(summary$discovery_status == "discovery_pass"),
      discovery_up = sum(summary$direction == "up" & summary$discovery_status == "discovery_pass"),
      discovery_down = sum(summary$direction == "down" & summary$discovery_status == "discovery_pass"),
      response_contrast_statistics_accessed = FALSE, source_state_artifact = state_artifact,
      source_states_note = "Existing Parse states reused; their original HVG selection used control cells, but no treatment-control effect enters this workflow."))
  log_event("Discovery complete", sum(summary$discovery_status == "discovery_pass"), "sets")
}
main()
