#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L, 5L)) stop(paste(
  "Usage: 32_significant_region_residual_ucell_kmeans.R",
  "<residual_ucell_report> <bb_up_change_points.tsv> <output_dir>",
  "[bp=K,mf=K,cc=K] [region_count_stratified]"))
report_root <- normalizePath(args[[1L]], mustWork = TRUE)
bb_path <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(cluster)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

source_completion <- read_json(file.path(report_root, "completion.json"),
  simplifyVector = TRUE)
stopifnot(identical(source_completion$status, "complete"),
  isTRUE(source_completion$complexity_residualization_enabled),
  source_completion$permutations == 1000L,
  isTRUE(source_completion$groupwise_multiple_testing))
bb <- fread(bb_path)
stopifnot(nrow(bb) == 2L,
  identical(as.character(bb$label), c("up1", "up2")),
  all(is.finite(bb$pt_median)))
bb[, boundary_label := c("BB1", "BB2")]

groups <- c("bp", "mf", "cc")
group_labels <- c(bp = "GO Biological Process", mf = "GO Molecular Function",
  cc = "GO Cellular Component")
candidate_k <- list(bp = 4:20, mf = 3:15, cc = 3:12)
fixed_k <- NULL
if (length(args) >= 4L) {
  pieces <- strsplit(args[[4L]], ",", fixed = TRUE)[[1L]]
  parsed <- strsplit(pieces, "=", fixed = TRUE)
  stopifnot(length(parsed) == 3L,
    all(vapply(parsed, length, integer(1)) == 2L))
  fixed_k <- setNames(as.integer(vapply(parsed, `[[`, character(1), 2L)),
    vapply(parsed, `[[`, character(1), 1L))
  stopifnot(setequal(names(fixed_k), groups),
    all(is.finite(fixed_k)), all(fixed_k >= 2L))
  fixed_k <- fixed_k[groups]
}
topology_aware <- length(args) == 5L
if (topology_aware)
  stopifnot(identical(args[[5L]], "region_count_stratified"),
    !is.null(fixed_k))
seed <- 314159L

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
pretty_term <- function(x) {
  x <- sub("^GO(BP|MF|CC)_", "", x)
  tools::toTitleCase(tolower(gsub("_", " ", x, fixed = TRUE)))
}
safe_cor <- function(x, y) {
  if (sd(x) <= .Machine$double.eps || sd(y) <= .Machine$double.eps)
    return(NA_real_)
  cor(x, y)
}
correlation_order <- function(x) {
  if (nrow(x) <= 1L) return(seq_len(nrow(x)))
  similarity <- cor(t(x), use = "pairwise.complete.obs")
  similarity[!is.finite(similarity)] <- 0
  diag(similarity) <- 1
  distance <- pmax(1 - similarity, 0)
  hclust(as.dist(distance), method = "average")$order
}
common_keywords <- function(x, maximum = 8L) {
  stop_words <- c("GO", "GOBP", "GOMF", "GOCC", "OF", "IN", "TO", "BY",
    "AND", "OR", "THE", "A", "AN", "PROCESS", "CELLULAR", "REGULATION",
    "POSITIVE", "NEGATIVE", "ACTIVITY", "PROTEIN", "DEPENDENT", "INVOLVED",
    "COMPLEX", "BINDING")
  words <- lapply(x, function(z) {
    y <- unlist(strsplit(gsub("[^A-Z0-9]+", " ", toupper(z)), " +"))
    unique(y[nchar(y) >= 3L & !y %in% stop_words])
  })
  counts <- sort(table(unlist(words)), decreasing = TRUE)
  if (!length(counts)) return("")
  counts <- head(counts, maximum)
  paste0(names(counts), " (", as.integer(counts), ")", collapse = "; ")
}

# Broad, predeclared labels make the cluster summaries reproducible. They are
# descriptive lexical annotations of GO names, not a second enrichment test.
theme_patterns <- c(
  "cell cycle, chromosome segregation, and DNA replication" =
    "CELL_CYCLE|MITOT|MEIOT|CHROMOSOM|CHROMATID|CENTROMER|KINETOCHORE|SPINDLE|DNA_REPLICATION|DNA_REPAIR|CYTOKINESIS|TELOMERE",
  "RNA processing, translation, and ribosome biology" =
    "RNA_|MRNA|TRNA|RRNA|SPLIC|RIBOSOM|TRANSLAT|CODON|RIBONUCLEOPROTEIN",
  "chromatin, transcription, and epigenetic control" =
    "CHROMATIN|HISTONE|NUCLEOSOM|EPIGEN|TRANSCRIPTION|POLYMERASE|DNA_METHYL|HETEROCHROMATIN",
  "B-cell, lymphocyte, antigen, and immune programs" =
    "B_CELL|LYMPH|IMMUN|ANTIGEN|MHC_|LEUKOCYTE|INTERFERON|T_CELL|V_D_J|RECOMBINATION|RECEPTOR_SIGNALING",
  "cytokine, growth-factor, and intracellular signaling" =
    "CYTOKINE|INTERLEUKIN|CHEMOKINE|GROWTH_FACTOR|SIGNAL|KINASE|PHOSPHATASE|GTPASE|SECOND_MESSENGER",
  "development, differentiation, and morphogenesis" =
    "DEVELOP|DIFFERENTIATION|MORPHOGEN|ORGANOGEN|STEM_CELL|LINEAGE|FATE|MATURATION",
  "metabolism, mitochondria, and bioenergetics" =
    "METABOL|BIOSYNTH|CATABOL|MITOCHON|RESPIRATORY|OXIDATIVE|GLYCOL|LIPID|FATTY_ACID|NUCLEOTIDE|ATP",
  "transport, vesicles, secretion, and localization" =
    "TRANSPORT|VESICLE|ENDOSOM|LYSOSOM|GOLGI|SECRETION|LOCALIZATION|MEMBRANE|IMPORT|EXPORT",
  "adhesion, migration, and cytoskeletal organization" =
    "ADHES|MIGRAT|CYTOSKELE|ACTIN|MICROTUBULE|MOTILITY|EXTRACELLULAR_MATRIX|CELL_JUNCTION",
  "stress responses, apoptosis, and proteostasis" =
    "STRESS|APOPT|CELL_DEATH|UNFOLDED|PROTEAS|AUTOPHAG|DAMAGE|HYPOX|OXIDATIVE_STRESS",
  "enzyme, receptor, and molecular binding functions" =
    "BINDING|CATALYTIC|ENZYME|RECEPTOR|TRANSFERASE|HYDROLASE|LIGASE|OXIDOREDUCTASE|CHANNEL|TRANSPORTER",
  "organelle, membrane, and structural compartments" =
    "ORGANELLE|MEMBRANE|NUCLEUS|NUCLEAR|CYTOPLAS|CYTOSOL|MITOCHON|RIBOSOM|CHROMOSOM|VESICLE|EXTRACELLULAR|CORTEX"
)
annotate_themes <- function(terms) {
  assigned <- rep("other or mixed GO processes", length(terms))
  for (theme in names(theme_patterns)) {
    take <- assigned == "other or mixed GO processes" &
      grepl(theme_patterns[[theme]], terms, perl = TRUE)
    assigned[take] <- theme
  }
  counts <- sort(table(assigned), decreasing = TRUE)
  total <- length(terms)
  display <- paste0(names(head(counts, 3L)), " (",
    as.integer(head(counts, 3L)), ", ",
    sprintf("%.1f%%", 100 * as.numeric(head(counts, 3L)) / total), ")",
    collapse = "; ")
  named <- counts[names(counts) != "other or mixed GO processes"]
  label <- if (length(named)) names(named)[[1L]] else "mixed GO processes"
  if (length(named) >= 2L && as.numeric(named[[2L]]) >= 0.6 * as.numeric(named[[1L]]))
    label <- paste(label, names(named)[[2L]], sep = "; ")
  list(primary = label, display = display, assignments = assigned)
}

all_membership <- list(); all_centroids <- list(); all_summaries <- list()
all_silhouette <- list(); completion_groups <- list()
report <- c(if (topology_aware)
    "# Region-count-stratified significant-region residual-UCell k-means clusters" else if (is.null(fixed_k))
    "# Significant-region residual-UCell k-means clusters" else
    "# Fine-resolution significant-region residual-UCell k-means clusters", "",
  paste0("Only GO sets with at least one permutation-supported positive region were clustered. ",
    "For each set, its 97-window RNA-feature-residualized UCell median was min-max scaled to 0–1, ",
    "then set to zero outside supported regions. Thus k-means uses only localized supported signal; ",
    "unsupported trajectory values and the original gene-expression z-score do not enter clustering."), "",
  if (topology_aware) paste0(
    "Visual review showed that temporal k-means could mix sets with one versus two supported regions. ",
    "Sets are therefore stratified before clustering as exactly 1, exactly 2, or 3+ supported regions. ",
    "The requested fine total resolution is retained (BP=", fixed_k[["bp"]], ", MF=",
    fixed_k[["mf"]], ", CC=", fixed_k[["cc"]], "), allocated across strata in proportion to the ",
    "square root of stratum size with at least one cluster per nonempty stratum. No final cluster can ",
    "mix region-count strata. Within each stratum, k-means uses all eligible sets with 100 starts.") else if (is.null(fixed_k)) paste0(
    "The number of clusters was selected independently by the highest mean silhouette over a ",
    "predeclared range (BP 4–20, MF 3–15, CC 3–12), using a deterministic sample of at most 1,000 ",
    "sets. Final k-means used all eligible sets with 100 starts. C01, C02, … are ordered by centroid ",
    "peak pseudotime from early to late; identifiers are descriptive, not cross-ontology matches.") else paste0(
    "Fine resolution was explicitly requested after visual review of the coarse silhouette solution. ",
    "Cluster counts are fixed at BP=", fixed_k[["bp"]], ", MF=", fixed_k[["mf"]],
    ", and CC=", fixed_k[["cc"]], ". Silhouette values are retained as diagnostics rather than used ",
    "for selection. Final k-means used all eligible sets with 100 starts. C01, C02, … are ordered by ",
    "centroid peak pseudotime from early to late; identifiers are descriptive, not cross-ontology matches."), "",
  paste0("Biological labels below are reproducible lexical summaries of member GO names. They are ",
    "interpretive aids rather than additional enrichment tests. Representative terms are the members ",
    "most correlated with the cluster centroid."), "")

for (group_index in seq_along(groups)) {
  g <- groups[[group_index]]
  message("Processing ", g)
  group_dir <- file.path(report_root, g)
  curves <- read_tsv_gz(file.path(group_dir, "window_curves.tsv.gz"))
  regions <- fread(file.path(group_dir, "significant_positive_regions.tsv"))
  passing <- fread(file.path(group_dir, "significant_gene_sets.tsv"))
  group_completion <- read_json(file.path(group_dir, "completion.json"),
    simplifyVector = TRUE)
  selected <- sort(unique(regions$feature_set))
  stopifnot(length(selected) > max(candidate_k[[g]]),
    all(selected %chin% passing$feature_set),
    group_completion$significant_positive_regions == nrow(regions))
  curves <- curves[feature_set %chin% selected]
  setorder(curves, feature_set, window_index)
  window_map <- unique(curves[, .(window_index, window_center)])
  setorder(window_map, window_index)
  wide <- dcast(curves, feature_set ~ window_index, value.var = "ucell_median")
  feature_sets <- wide$feature_set
  raw <- as.matrix(wide[, -"feature_set"]); rownames(raw) <- feature_sets
  stopifnot(ncol(raw) == 97L, !anyNA(raw), all(is.finite(raw)))
  minima <- apply(raw, 1L, min); spans <- apply(raw, 1L, function(x) diff(range(x)))
  stopifnot(all(spans > .Machine$double.eps))
  scaled <- (raw - minima) / spans
  mask <- matrix(FALSE, nrow(scaled), ncol(scaled), dimnames = dimnames(scaled))
  for (i in seq_along(feature_sets)) {
    rr <- regions[feature_set == feature_sets[[i]]]
    for (j in seq_len(nrow(rr))) mask[i, ] <- mask[i, ] |
      (window_map$window_center >= rr$region_left[[j]] &
        window_map$window_center <= rr$region_right[[j]])
  }
  masked <- scaled; masked[!mask] <- 0
  stopifnot(all(rowSums(mask) > 0L), all(rowSums(masked) > 0))
  region_count <- regions[, .N, by = feature_set]$N[
    match(feature_sets, regions[, .N, by = feature_set]$feature_set)]
  stopifnot(!anyNA(region_count), all(region_count >= 1L))
  region_count_stratum <- pmin(region_count, 3L)

  set.seed(seed + group_index)
  sample_index <- if (nrow(masked) <= 1000L) seq_len(nrow(masked)) else
    sort(sample(seq_len(nrow(masked)), 1000L))
  sample_matrix <- masked[sample_index, , drop = FALSE]
  sample_distance <- dist(sample_matrix)
  silhouette_table <- rbindlist(lapply(candidate_k[[g]], function(k) {
    set.seed(seed + group_index * 100L + k)
    fit <- kmeans(sample_matrix, centers = k, nstart = 30L, iter.max = 200L)
    sil <- silhouette(fit$cluster, sample_distance)
    data.table(analysis_group = g, k = k,
      mean_silhouette = mean(sil[, "sil_width"]),
      median_silhouette = median(sil[, "sil_width"]),
      minimum_cluster_size = min(table(fit$cluster)),
      sampled_gene_sets = nrow(sample_matrix))
  }))
  chosen_k <- if (is.null(fixed_k))
    silhouette_table$k[[which.max(silhouette_table$mean_silhouette)]] else
    fixed_k[[g]]
  if (!chosen_k %in% silhouette_table$k) {
    set.seed(seed + group_index * 100L + chosen_k)
    diagnostic_fit <- kmeans(sample_matrix, centers = chosen_k, nstart = 30L,
      iter.max = 200L)
    diagnostic_sil <- silhouette(diagnostic_fit$cluster, sample_distance)
    silhouette_table <- rbind(silhouette_table, data.table(
      analysis_group = g, k = chosen_k,
      mean_silhouette = mean(diagnostic_sil[, "sil_width"]),
      median_silhouette = median(diagnostic_sil[, "sil_width"]),
      minimum_cluster_size = min(table(diagnostic_fit$cluster)),
      sampled_gene_sets = nrow(sample_matrix)))
    setorder(silhouette_table, k)
  }
  silhouette_table[, selected := k == chosen_k]
  all_silhouette[[g]] <- silhouette_table

  set.seed(seed + group_index)
  stratum_allocation <- NULL
  if (topology_aware) {
    stratum_sizes <- table(factor(region_count_stratum, levels = 1:3))
    active <- which(stratum_sizes > 0L)
    weights <- sqrt(as.numeric(stratum_sizes[active]))
    allocation <- pmax(1L, floor(chosen_k * weights / sum(weights)))
    allocation <- pmin(allocation, as.integer(stratum_sizes[active]))
    while (sum(allocation) < chosen_k) {
      capacity <- as.integer(stratum_sizes[active]) - allocation
      add <- which.max(ifelse(capacity > 0L,
        chosen_k * weights / sum(weights) - allocation, -Inf))
      allocation[[add]] <- allocation[[add]] + 1L
    }
    while (sum(allocation) > chosen_k) {
      remove <- which.max(ifelse(allocation > 1L,
        allocation - chosen_k * weights / sum(weights), -Inf))
      allocation[[remove]] <- allocation[[remove]] - 1L
    }
    stratum_allocation <- setNames(allocation, active)
    raw_cluster <- integer(nrow(masked)); center_parts <- list(); offset <- 0L
    for (stratum in active) {
      idx <- which(region_count_stratum == stratum)
      stratum_k <- stratum_allocation[[as.character(stratum)]]
      set.seed(seed + group_index * 1000L + stratum)
      stratum_fit <- kmeans(masked[idx, , drop = FALSE], centers = stratum_k,
        nstart = 100L, iter.max = 300L)
      raw_cluster[idx] <- offset + stratum_fit$cluster
      centers <- stratum_fit$centers
      if (is.null(dim(centers))) centers <- matrix(centers, nrow = 1L)
      center_parts[[as.character(stratum)]] <- centers
      offset <- offset + stratum_k
    }
    raw_centers <- do.call(rbind, center_parts)
    stopifnot(offset == chosen_k, all(raw_cluster > 0L),
      nrow(raw_centers) == chosen_k)
  } else {
    fit <- kmeans(masked, centers = chosen_k, nstart = 100L, iter.max = 300L)
    raw_cluster <- fit$cluster
    raw_centers <- fit$centers
  }
  peak_col <- max.col(raw_centers, ties.method = "first")
  weighted_pt <- apply(raw_centers, 1L, function(x)
    weighted.mean(window_map$window_center, pmax(x, 0)))
  raw_order <- order(window_map$window_center[peak_col], weighted_pt)
  raw_to_display <- setNames(sprintf("C%02d", seq_len(chosen_k)), raw_order)
  cluster_id <- unname(raw_to_display[as.character(raw_cluster)])
  display_ids <- sprintf("C%02d", seq_len(chosen_k))

  ordered_sets <- unlist(lapply(display_ids, function(cid) {
    idx <- which(cluster_id == cid)
    feature_sets[idx][correlation_order(masked[idx, , drop = FALSE])]
  }), use.names = FALSE)
  region_stats <- regions[, .(
    n_supported_regions = .N,
    supported_window_left = min(region_left),
    supported_window_right = max(region_right),
    strongest_region_effect = max(peak_ucell_effect),
    strongest_region_fwer_p = min(p_region)), by = feature_set]
  membership <- data.table(analysis_group = g, feature_set = feature_sets,
    set_name = feature_sets, cluster_id = cluster_id,
    region_count_stratum = ifelse(region_count_stratum == 3L, "3+",
      as.character(region_count_stratum)),
    row_position_top_to_bottom = match(feature_sets, ordered_sets),
    masked_signal_sum = rowSums(masked),
    masked_signal_max = apply(masked, 1L, max),
    supported_windows = rowSums(mask))
  membership <- merge(membership, region_stats, by = "feature_set",
    all.x = TRUE, sort = FALSE)
  setorder(membership, row_position_top_to_bottom)

  centroid_rows <- list(); summary_rows <- list()
  for (cid in display_ids) {
    idx <- which(cluster_id == cid)
    centroid <- colMeans(masked[idx, , drop = FALSE])
    centroid_rows[[cid]] <- data.table(analysis_group = g, cluster_id = cid,
      window_index = window_map$window_index,
      window_center = window_map$window_center,
      mean_significant_region_score = centroid,
      fraction_members_supported = colMeans(mask[idx, , drop = FALSE]))
    cors <- vapply(idx, function(i) safe_cor(masked[i, ], centroid), numeric(1))
    rank_index <- idx[order(-cors, -membership$masked_signal_sum[
      match(feature_sets[idx], membership$feature_set)], na.last = TRUE)]
    representative <- feature_sets[head(rank_index, 8L)]
    member_terms <- feature_sets[idx]
    themes <- annotate_themes(member_terms)
    membership[feature_set %chin% member_terms,
      lexical_theme := themes$assignments[match(feature_set, member_terms)]]
    peak_i <- which.max(centroid)
    centroid_weighted_pt <- weighted.mean(window_map$window_center,
      pmax(centroid, 0))
    phase <- if (window_map$window_center[[peak_i]] < bb$pt_median[[1L]])
      "before BB1" else if (window_map$window_center[[peak_i]] < bb$pt_median[[2L]])
      "between BB1 and BB2" else "after BB2"
    member_rows <- membership[feature_set %chin% member_terms]
    cluster_strata <- unique(member_rows$region_count_stratum)
    if (topology_aware) stopifnot(length(cluster_strata) == 1L)
    summary_rows[[cid]] <- data.table(analysis_group = g, cluster_id = cid,
      biological_summary = themes$primary,
      region_count_stratum = paste(cluster_strata, collapse = ","),
      n_gene_sets = length(idx),
      centroid_peak_pt = window_map$window_center[[peak_i]],
      centroid_weighted_pt = centroid_weighted_pt,
      centroid_peak_phase = phase,
      centroid_peak_score = centroid[[peak_i]],
      median_supported_windows = median(member_rows$supported_windows),
      median_supported_regions = median(member_rows$n_supported_regions),
      representative_gene_sets = paste(pretty_term(representative), collapse = "; "),
      dominant_lexical_themes = themes$display,
      recurring_name_keywords = common_keywords(member_terms))
  }
  centroids <- rbindlist(centroid_rows)
  summaries <- rbindlist(summary_rows)
  all_membership[[g]] <- membership
  all_centroids[[g]] <- centroids
  all_summaries[[g]] <- summaries

  plot_long <- as.data.table(as.table(masked))
  setnames(plot_long, c("feature_set", "window_column", "value"))
  plot_long[, feature_set := as.character(feature_set)]
  plot_long[, window_index := window_map$window_index[as.integer(window_column)]]
  plot_long[, window_center := window_map$window_center[as.integer(window_column)]]
  plot_long[, row_position_top_to_bottom := match(feature_set, ordered_sets)]
  plot_long[, y := length(ordered_sets) - row_position_top_to_bottom + 1L]
  ordered_cluster <- cluster_id[match(ordered_sets, feature_sets)]
  cluster_runs <- rle(ordered_cluster)
  run_end <- cumsum(cluster_runs$lengths)
  run_start <- c(1L, head(run_end, -1L) + 1L)
  y_breaks <- length(ordered_sets) - (run_start + run_end) / 2 + 1L
  boundary_y <- length(ordered_sets) - head(run_end, -1L) + 0.5
  tile_width <- min(diff(window_map$window_center))
  heatmap <- ggplot(plot_long, aes(window_center, y, fill = value)) +
    geom_tile(width = tile_width, height = 1) +
    geom_hline(yintercept = boundary_y, linewidth = .25, color = "#202020") +
    geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
      linewidth = .45, linetype = "dashed", color = "#374151") +
    scale_fill_gradientn(colors = c("#FFFFFF", "#FEE8C8", "#FDBB84",
      "#E34A33", "#7F0000"), limits = c(0, 1), name = "Residual\nUCell 0–1") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = y_breaks, labels = cluster_runs$values,
      expand = c(0, 0)) +
    labs(title = paste0(group_labels[[g]],
      ": supported-region residual-UCell k-means"),
      subtitle = paste0(length(selected), " GO sets; ",
        if (is.null(fixed_k)) "silhouette-selected" else "review-selected fine-resolution",
        " k=", chosen_k,
        "; white outside permutation-supported regions"),
      x = "Existing Bdev pseudotime (pt)", y = "K-means cluster") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold"))
  ggsave(file.path(stage_dir, paste0(g, "_significant_region_kmeans_heatmap.png")),
    heatmap, width = 12, height = max(8, min(20, length(selected) / 130)),
    dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, paste0(g, "_significant_region_kmeans_heatmap.pdf")),
    heatmap, width = 12, height = max(8, min(20, length(selected) / 130)),
    device = cairo_pdf, bg = "white")

  centroid_plot <- ggplot(centroids,
      aes(window_center, mean_significant_region_score, color = cluster_id)) +
    geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
      linewidth = .35, linetype = "dashed", color = "#6B7280") +
    geom_line(linewidth = .65, show.legend = FALSE) +
    facet_wrap(~cluster_id, ncol = 4) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
    labs(title = paste0(group_labels[[g]], ": cluster mean supported-region profiles"),
      x = "Existing Bdev pseudotime (pt)", y = "Mean residual-UCell 0–1") +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"))
  ggsave(file.path(stage_dir, paste0(g, "_significant_region_cluster_centroids.png")),
    centroid_plot, width = 12, height = max(5, ceiling(chosen_k / 4) * 2.2),
    dpi = 250, bg = "white")
  ggsave(file.path(stage_dir, paste0(g, "_significant_region_cluster_centroids.pdf")),
    centroid_plot, width = 12, height = max(5, ceiling(chosen_k / 4) * 2.2),
    device = cairo_pdf, bg = "white")

  report <- c(report, paste0("## ", group_labels[[g]]), "",
    paste0(if (is.null(fixed_k)) "Selected" else "Used fine-resolution",
      " **k = ", chosen_k, "**", if (is.null(fixed_k)) " by mean silhouette" else
      " after review of the coarse clustering", "; ", length(selected),
      " GO sets had at least one supported region (", nrow(regions),
      " total regions)."), "",
    vapply(seq_len(nrow(summaries)), function(i) paste0("- **",
      summaries$cluster_id[[i]], " — ", summaries$biological_summary[[i]],
      ".** ", summaries$n_gene_sets[[i]], " sets; supported-region count: ",
      summaries$region_count_stratum[[i]], "; centroid peaks at pt ",
      signif(summaries$centroid_peak_pt[[i]], 5), " (",
      summaries$centroid_peak_phase[[i]], "). Representative terms: ",
      summaries$representative_gene_sets[[i]], "."), character(1)), "",
    "| Cluster | Region-count stratum | Biological summary | Sets | Peak pt | Peak phase | Median supported windows | Dominant lexical themes |",
    "|---|---:|---|---:|---:|---|---:|---|",
    vapply(seq_len(nrow(summaries)), function(i) paste0("| ",
      summaries$cluster_id[[i]], " | ", summaries$region_count_stratum[[i]],
      " | ", summaries$biological_summary[[i]], " | ",
      summaries$n_gene_sets[[i]], " | ", signif(summaries$centroid_peak_pt[[i]], 5),
      " | ", summaries$centroid_peak_phase[[i]], " | ",
      summaries$median_supported_windows[[i]], " | ",
      summaries$dominant_lexical_themes[[i]], " |"), character(1)), "")
  completion_groups[[g]] <- list(eligible_gene_sets = length(selected),
    supported_regions = nrow(regions), selected_k = chosen_k,
    region_count_stratified = topology_aware,
    region_count_stratum_sizes = if (topology_aware)
      as.list(setNames(as.integer(stratum_sizes), c("one", "two", "three_plus"))) else NULL,
    region_count_stratum_cluster_allocation = if (topology_aware)
      as.list(setNames(as.integer(stratum_allocation), paste0(names(stratum_allocation),
        ifelse(names(stratum_allocation) == "3", "+", "")))) else NULL,
    silhouette_sample_size = length(sample_index),
    selected_mean_silhouette = silhouette_table[k == chosen_k,
      mean_silhouette][[1L]])
}

membership <- rbindlist(all_membership, use.names = TRUE, fill = TRUE)
summaries <- rbindlist(all_summaries, use.names = TRUE, fill = TRUE)
centroids <- rbindlist(all_centroids, use.names = TRUE, fill = TRUE)
silhouette_scan <- rbindlist(all_silhouette, use.names = TRUE, fill = TRUE)
fwrite(membership, file.path(stage_dir, "gene_set_cluster_membership.tsv"),
  sep = "\t", quote = FALSE, na = "NA")
fwrite(summaries, file.path(stage_dir, "cluster_summary.tsv"),
  sep = "\t", quote = FALSE, na = "NA")
fwrite(centroids, file.path(stage_dir, "cluster_centroid_profiles.tsv.gz"),
  sep = "\t", quote = FALSE, na = "NA")
fwrite(silhouette_scan, file.path(stage_dir, "silhouette_k_scan.tsv"),
  sep = "\t", quote = FALSE, na = "NA")
writeLines(report, file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", promoted = FALSE,
  source_report = report_root,
  score = paste("per-set 0-1 scaled RNA-feature-residualized UCell median",
    "inside permutation-supported positive regions; zero outside"),
  eligibility = "at least one permutation-supported positive region",
  clustering = paste0(if (topology_aware)
    "region-count-stratified " else "",
    "Euclidean k-means on significant-region-only 97-window profiles"),
  k_selection = if (topology_aware)
    paste0("review-selected fine total resolution: ", paste(names(fixed_k),
      fixed_k, sep = "=", collapse = ","),
      "; allocated across 1, 2, and 3+ region strata by sqrt(stratum size)") else if (is.null(fixed_k))
    "maximum mean silhouette on deterministic sample up to 1000 sets" else
    paste0("review-selected fine resolution: ", paste(names(fixed_k),
      fixed_k, sep = "=", collapse = ","),
      "; silhouette retained as diagnostic"),
  random_seed = seed, biological_summary_method = paste(
    "predeclared lexical theme annotation plus centroid-correlated",
    "representative GO terms; descriptive, not enrichment inference"),
  groups = completion_groups), file.path(stage_dir, "completion.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)
stopifnot(file.rename(stage_dir, output_dir))
message("Completed significant-region residual-UCell clustering at: ", output_dir)
