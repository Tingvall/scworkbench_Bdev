#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(jsonlite)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop("Missing argument: ", flag)
  args[[index + 1L]]
}

config_path <- arg_value("--config")
output_root <- arg_value("--output-root")
config <- yaml::read_yaml(config_path)
if (!identical(config$schema_version, "bdev.immgen_marker_workflow_qualification.v1")) {
  stop("Unexpected qualification config schema")
}

resolve_output <- function(relative) file.path(output_root, relative)
fixture <- config$outputs$fixture
model_outputs <- config$outputs$model

counts_frame <- read.delim(resolve_output(fixture$counts), check.names = FALSE, stringsAsFactors = FALSE)
gene_symbols <- counts_frame[[1L]]
counts <- as.matrix(counts_frame[, -1L, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- gene_symbols

metadata <- read.delim(resolve_output(fixture$metadata), check.names = FALSE, stringsAsFactors = FALSE)
annotation <- read.delim(resolve_output(fixture$annotation), check.names = FALSE, stringsAsFactors = FALSE)
plan <- read.delim(resolve_output(fixture$contrast_plan), check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(identical(colnames(counts), metadata$library_id))
stopifnot(identical(rownames(counts), annotation$gene_symbol))
stopifnot(ncol(counts) == 126L, sum(metadata$role == "discovery") == 121L)
stopifnot(all(counts >= 0), all(counts == round(counts)))

aggregate_by_population <- function(input_counts, input_metadata) {
  populations <- unique(input_metadata$population)
  aggregated <- vapply(
    populations,
    function(population) rowSums(input_counts[, input_metadata$population == population, drop = FALSE]),
    numeric(nrow(input_counts))
  )
  rownames(aggregated) <- rownames(input_counts)
  population_metadata <- do.call(
    rbind,
    lapply(populations, function(population) {
      rows <- input_metadata[input_metadata$population == population, , drop = FALSE]
      data.frame(
        library_id = population,
        population = population,
        table_s1_lineage = unique(rows$table_s1_lineage),
        superlineage = unique(rows$superlineage),
        b_branch = unique(rows$b_branch),
        role = unique(rows$role),
        stringsAsFactors = FALSE
      )
    })
  )
  stopifnot(ncol(aggregated) == nrow(population_metadata))
  list(counts = aggregated, metadata = population_metadata)
}

auc_rank <- function(values, positive) {
  positive <- as.logical(positive)
  n_pos <- sum(positive)
  n_neg <- sum(!positive)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  ranks <- rank(values, ties.method = "average")
  (sum(ranks[positive]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

weights_for <- function(text, levels) {
  parsed <- jsonlite::fromJSON(text)
  result <- setNames(rep(0, length(levels)), levels)
  unknown <- setdiff(names(parsed), levels)
  if (length(unknown)) stop("Contrast contains absent groups: ", paste(unknown, collapse = ","))
  result[names(parsed)] <- as.numeric(parsed)
  result
}

run_family <- function(family) {
  family_plan <- plan[plan$family == family, , drop = FALSE]
  stopifnot(nrow(family_plan) > 0L)
  scopes <- unique(family_plan$sample_scope)
  group_columns <- unique(family_plan$group_column)
  stopifnot(length(scopes) == 1L, length(group_columns) == 1L)
  scope <- scopes[[1L]]
  group_column <- group_columns[[1L]]

  selected <- metadata$role == "discovery"
  if (scope %in% c("b_discovery", "b_local_discovery")) selected <- selected & metadata$table_s1_lineage == "B"
  if (scope == "b_local_discovery") {
    local_populations <- unique(c(family_plan$target, unlist(lapply(family_plan$weights_json, function(value) names(jsonlite::fromJSON(value))))))
    selected <- selected & metadata$population %in% local_populations
  }
  family_counts <- counts[, selected, drop = FALSE]
  family_metadata <- metadata[selected, , drop = FALSE]
  if (family != "b_local_stage") {
    aggregated <- aggregate_by_population(family_counts, family_metadata)
    family_counts <- aggregated$counts
    family_metadata <- aggregated$metadata
  }

  groups <- factor(family_metadata[[group_column]])
  design <- model.matrix(~ 0 + groups)
  group_levels <- levels(groups)
  colnames(design) <- group_levels
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient design for ", family)

  y_all <- DGEList(counts = family_counts)
  keep <- filterByExpr(y_all, design = design, min.count = 10, min.total.count = 15)
  if (sum(keep) < 2L) stop("Too few genes pass filtering for ", family)
  y <- y_all[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y_all$samples$norm.factors <- y$samples$norm.factors
  y <- estimateDisp(y, design = design)
  fit <- glmQLFit(y, design = design, robust = TRUE)
  normalized_cpm <- cpm(y_all, normalized.lib.sizes = TRUE)

  design_audit <- data.frame(
    family = family,
    unit = if (family == "b_local_stage") "rna_library" else "population_profile",
    samples = ncol(family_counts),
    groups = nlevels(groups),
    design_columns = ncol(design),
    design_rank = qr(design)$rank,
    residual_df = nrow(design) - qr(design)$rank,
    filtered_in = sum(keep),
    filtered_out = sum(!keep),
    stringsAsFactors = FALSE
  )

  statistics <- lapply(seq_len(nrow(family_plan)), function(index) {
    plan_row <- family_plan[index, , drop = FALSE]
    contrast <- weights_for(plan_row$weights_json, group_levels)
    test <- glmTreat(fit, contrast = contrast, lfc = 0.5)
    table <- topTags(test, n = Inf, sort.by = "none")$table
    tested <- data.frame(
      gene_symbol = rownames(table),
      logFC = table$logFC,
      logCPM = table$logCPM,
      PValue = table$PValue,
      FDR = table$FDR,
      stringsAsFactors = FALSE
    )
    complete <- merge(data.frame(gene_symbol = rownames(y_all), filter_pass = keep, stringsAsFactors = FALSE), tested, by = "gene_symbol", all.x = TRUE, sort = FALSE)
    complete <- complete[match(rownames(y_all), complete$gene_symbol), , drop = FALSE]

    parsed_weights <- jsonlite::fromJSON(plan_row$weights_json)
    competitor_groups <- names(parsed_weights)[as.numeric(parsed_weights) < 0]
    target_samples <- groups == plan_row$target
    competitor_samples <- groups %in% competitor_groups
    target_expression_fraction <- rowMeans(normalized_cpm[, target_samples, drop = FALSE] >= 1)
    competitor_median <- apply(normalized_cpm[, competitor_samples, drop = FALSE], 1L, median)
    target_direction_fraction <- rowMeans(normalized_cpm[, target_samples, drop = FALSE] > competitor_median)
    auc <- vapply(
      seq_len(nrow(normalized_cpm)),
      function(gene_index) auc_rank(normalized_cpm[gene_index, target_samples | competitor_samples], target_samples[target_samples | competitor_samples]),
      numeric(1L)
    )

    complete$family <- family
    complete$target <- plan_row$target
    complete$test_id <- plan_row$test_id
    complete$test_type <- plan_row$test_type
    complete$competitor_groups <- paste(competitor_groups, collapse = ";")
    complete$target_expression_fraction <- target_expression_fraction
    complete$target_direction_fraction <- target_direction_fraction
    complete$target_support_fraction <- pmin(target_expression_fraction, target_direction_fraction)
    complete$target_versus_competitor_auroc <- auc
    complete$stability_fraction <- 1.0
    complete$stability_status <- "qualification_placeholder_not_production"
    complete <- merge(complete, annotation, by = "gene_symbol", all.x = TRUE, sort = FALSE)
    complete <- complete[match(rownames(y_all), complete$gene_symbol), , drop = FALSE]
    complete[, c(
      "family", "target", "test_id", "test_type", "competitor_groups", "gene_symbol",
      "filter_pass", "logFC", "logCPM", "PValue", "FDR", "target_expression_fraction",
      "target_direction_fraction", "target_support_fraction", "target_versus_competitor_auroc",
      "stability_fraction", "stability_status", "gene_class", "annotation_eligible",
      "state_flag", "mapping_status"
    )]
  })
  list(statistics = do.call(rbind, statistics), audit = design_audit)
}

families <- c("superlineage", "broad_lineage", "b_branch", "b_local_stage")
results <- lapply(families, run_family)
complete_statistics <- do.call(rbind, lapply(results, `[[`, "statistics"))
design_audit <- do.call(rbind, lapply(results, `[[`, "audit"))

statistics_path <- resolve_output(model_outputs$complete_statistics)
audit_path <- resolve_output(model_outputs$design_audit)
dir.create(dirname(statistics_path), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(statistics_path, open = "wt", compression = 9)
write.table(complete_statistics, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
close(connection)
write.table(design_audit, audit_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat(jsonlite::toJSON(list(
  genes = nrow(counts),
  libraries = ncol(counts),
  contrast_tests = nrow(plan),
  complete_statistic_rows = nrow(complete_statistics),
  designs = design_audit
), auto_unbox = TRUE, pretty = TRUE), "\n")
