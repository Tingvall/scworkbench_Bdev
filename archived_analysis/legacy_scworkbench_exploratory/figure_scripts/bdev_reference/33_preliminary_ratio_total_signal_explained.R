#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 33_preliminary_ratio_total_signal_explained.R",
  "<smoothed_report> <turn_pattern_report> <depth_report> <output_dir>"))
smoothed_report <- normalizePath(args[[1L]], mustWork = TRUE)
turn_pattern_report <- normalizePath(args[[2L]], mustWork = TRUE)
depth_report <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
})
read_tsv <- function(path) read.delim(path, check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = "NA")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
r_squared <- function(observed, predicted) 1 -
  sum((observed - predicted)^2) / sum((observed - mean(observed))^2)

windows <- read_tsv(file.path(smoothed_report, "smoothed_slope_windows.tsv"))
totals <- read_tsv(file.path(depth_report,
  "smoothed_window_total_rna_atac.tsv"))
all_events <- read_tsv(file.path(turn_pattern_report,
  "turn_event_definition.tsv"))
events <- all_events[match(paste0("T", 1:3), all_events$region_id), ]
stopifnot(nrow(windows) == 2524L, nrow(totals) == nrow(windows),
  identical(windows$window_id, totals$window_id),
  identical(events$region_id, paste0("T", 1:3)),
  all(windows$smoothed_local_slope > 0))

model_data <- data.frame(window_id = windows$window_id,
  window_index = windows$window_index, pt = windows$pt_median,
  ratio = log2(windows$smoothed_local_slope),
  total_RNA = totals$log1p_total_RNA_mean,
  total_ATAC = totals$log1p_total_ATAC_mean)
turn_indices <- lapply(seq_len(3L), function(k)
  seq.int(events$left_window_index[[k]], events$right_window_index[[k]]))
scored_indices <- sort(unique(unlist(turn_indices)))
model_data$in_scored_region <- seq_len(nrow(model_data)) %in% scored_indices

full_fit <- lm(ratio ~ total_RNA + total_ATAC, data = model_data)
local_fit <- lm(ratio ~ total_RNA + total_ATAC,
  data = model_data[scored_indices, ])
model_data$full_fitted <- predict(full_fit, newdata = model_data)
model_data$full_residual <- model_data$ratio - model_data$full_fitted
model_data$local_fitted <- predict(local_fit, newdata = model_data)
model_data$local_residual <- model_data$ratio - model_data$local_fitted

# Contiguous folds prevent adjacent overlapping windows from being split
# randomly between train and test. This evaluates broad-path prediction.
fold <- cut(seq_len(nrow(model_data)), breaks = 10L, labels = FALSE)
blocked_prediction <- rep(NA_real_, nrow(model_data))
for (k in sort(unique(fold))) {
  test <- which(fold == k)
  train <- which(fold != k)
  fit <- lm(ratio ~ total_RNA + total_ATAC, data = model_data[train, ])
  blocked_prediction[test] <- predict(fit, newdata = model_data[test, ])
}
stopifnot(all(is.finite(blocked_prediction)))
model_data$blocked_cv_fitted <- blocked_prediction

# Central derivatives match the change-based comparisons used in the gene
# screens and ask whether total-signal changes explain local ratio changes.
derivative <- data.frame(window_index = 2:(nrow(model_data) - 1L),
  pt = model_data$pt[2:(nrow(model_data) - 1L)],
  ratio_derivative = (model_data$ratio[3:nrow(model_data)] -
    model_data$ratio[1:(nrow(model_data) - 2L)]) / 2,
  total_RNA_derivative = (model_data$total_RNA[3:nrow(model_data)] -
    model_data$total_RNA[1:(nrow(model_data) - 2L)]) / 2,
  total_ATAC_derivative = (model_data$total_ATAC[3:nrow(model_data)] -
    model_data$total_ATAC[1:(nrow(model_data) - 2L)]) / 2)
derivative$in_scored_region <- derivative$window_index %in% scored_indices
derivative_fit <- lm(ratio_derivative ~ total_RNA_derivative +
  total_ATAC_derivative, data = derivative)
derivative_local_fit <- lm(ratio_derivative ~ total_RNA_derivative +
  total_ATAC_derivative, data = derivative[derivative$in_scored_region, ])
derivative$fitted <- predict(derivative_fit, newdata = derivative)
derivative$local_fitted <- predict(derivative_local_fit,
  newdata = derivative)
derivative$residual <- derivative$ratio_derivative - derivative$local_fitted

per_turn <- do.call(rbind, lapply(seq_len(3L), function(k) {
  ix <- derivative$window_index %in% turn_indices[[k]]
  fit <- lm(ratio_derivative ~ total_RNA_derivative +
    total_ATAC_derivative, data = derivative[ix, ])
  data.frame(region_id = paste0("T", k), windows = sum(ix),
    derivative_R2 = r_squared(derivative$ratio_derivative[ix], fitted(fit)),
    beta_RNA = coef(fit)[["total_RNA_derivative"]],
    beta_ATAC = coef(fit)[["total_ATAC_derivative"]])
}))

# Leave one complete turn out. Center within the held-out turn before pooling,
# so between-turn offsets cannot masquerade as prediction of local shape.
loto_observed <- loto_predicted <- numeric()
loto_rows <- vector("list", 3L)
for (k in seq_len(3L)) {
  test <- derivative$window_index %in% turn_indices[[k]]
  train_turns <- unlist(turn_indices[-k])
  train <- derivative$window_index %in% train_turns
  fit <- lm(ratio_derivative ~ total_RNA_derivative +
    total_ATAC_derivative, data = derivative[train, ])
  predicted <- predict(fit, newdata = derivative[test, ])
  observed_centered <- derivative$ratio_derivative[test] -
    mean(derivative$ratio_derivative[test])
  predicted_centered <- predicted - mean(predicted)
  loto_observed <- c(loto_observed, observed_centered)
  loto_predicted <- c(loto_predicted, predicted_centered)
  loto_rows[[k]] <- data.frame(region_id = paste0("T", k),
    held_out_centered_derivative_R2 = r_squared(observed_centered,
      predicted_centered),
    held_out_centered_derivative_correlation = cor(observed_centered,
      predicted_centered))
}
loto <- do.call(rbind, loto_rows)
loto_pooled_r2 <- r_squared(loto_observed, loto_predicted)
loto_pooled_correlation <- cor(loto_observed, loto_predicted)

summary_table <- data.frame(
  analysis = c("level_full_in_sample", "level_scored_in_sample",
    "level_full_10_contiguous_fold_cv", "derivative_full_in_sample",
    "derivative_scored_in_sample", "derivative_leave_one_turn_out_centered"),
  R2 = c(r_squared(model_data$ratio, model_data$full_fitted),
    r_squared(model_data$ratio[scored_indices],
      model_data$local_fitted[scored_indices]),
    r_squared(model_data$ratio, model_data$blocked_cv_fitted),
    r_squared(derivative$ratio_derivative, derivative$fitted),
    r_squared(derivative$ratio_derivative[derivative$in_scored_region],
      derivative$local_fitted[derivative$in_scored_region]),
    loto_pooled_r2),
  correlation = c(cor(model_data$ratio, model_data$full_fitted),
    cor(model_data$ratio[scored_indices],
      model_data$local_fitted[scored_indices]),
    cor(model_data$ratio, model_data$blocked_cv_fitted),
    cor(derivative$ratio_derivative, derivative$fitted),
    cor(derivative$ratio_derivative[derivative$in_scored_region],
      derivative$local_fitted[derivative$in_scored_region]),
    loto_pooled_correlation))
write_tsv(summary_table, file.path(stage_dir,
  "total_signal_explained_summary.tsv"))
write_tsv(per_turn, file.path(stage_dir,
  "per_turn_derivative_total_signal_models.tsv"))
write_tsv(loto, file.path(stage_dir,
  "leave_one_turn_out_derivative_models.tsv"))
write_tsv(model_data, file.path(stage_dir,
  "ratio_total_signal_fitted_and_residual_windows.tsv"))
write_tsv(derivative, file.path(stage_dir,
  "ratio_derivative_total_signal_fitted_and_residual_windows.tsv"))

highlight_regions <- data.frame(region_id = events$region_id,
  xmin = windows$pt_median[events$left_window_index],
  xmax = windows$pt_median[events$right_window_index])
level_long <- rbind(
  data.frame(pt = model_data$pt, panel = "Observed and depth-fitted ratio",
    series = "Observed log2 RNA/ATAC", value = model_data$ratio),
  data.frame(pt = model_data$pt, panel = "Observed and depth-fitted ratio",
    series = "Fitted from total RNA + ATAC", value = model_data$full_fitted),
  data.frame(pt = model_data$pt, panel = "Residual after total-signal fit",
    series = "Residual ratio", value = model_data$full_residual))
level_long$panel <- factor(level_long$panel,
  levels = c("Observed and depth-fitted ratio",
    "Residual after total-signal fit"))
level_plot <- ggplot(level_long, aes(pt, value, color = series)) +
  geom_rect(data = highlight_regions,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = region_id),
    inherit.aes = FALSE, alpha = 0.10, color = NA) +
  geom_hline(yintercept = 0, color = "#888888", linetype = "dotted",
    linewidth = 0.3) +
  geom_line(linewidth = 0.8) +
  geom_vline(data = all_events, aes(xintercept = pt), inherit.aes = FALSE,
    color = "#555555", linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ panel, ncol = 1L, scales = "free_y") +
  scale_color_manual(values = c("Observed log2 RNA/ATAC" = "#111111",
    "Fitted from total RNA + ATAC" = "#7B3294",
    "Residual ratio" = "#008837"), name = NULL) +
  scale_fill_manual(values = c("T1" = "#E41A1C", "T2" = "#377EB8",
    "T3" = "#4DAF4A"), name = "Scored ±100 windows") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
    strip.text = element_text(face = "bold"), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15)) +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  labs(title = "RNA/ATAC ratio component associated with total signal",
    subtitle = paste0("Full-path linear R² = ",
      sprintf("%.3f", summary_table$R2[[1L]]),
      "; 10-block CV R² = ", sprintf("%.3f", summary_table$R2[[3L]]),
      "; scored-region derivative R² = ",
      sprintf("%.3f", summary_table$R2[[5L]])),
    x = "Bdev pseudotime", y = "log2(delta RNA / delta ATAC)")
ggsave(file.path(stage_dir, "ratio_total_signal_fit_and_residual.png"),
  level_plot, width = 11, height = 8, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "ratio_total_signal_fit_and_residual.pdf"),
  level_plot, width = 11, height = 8, device = cairo_pdf, bg = "white")

derivative_long <- rbind(
  data.frame(pt = derivative$pt, series = "Observed ratio derivative",
    value = derivative$ratio_derivative),
  data.frame(pt = derivative$pt, series = "Fitted from total-signal derivatives",
    value = derivative$local_fitted),
  data.frame(pt = derivative$pt, series = "Residual derivative",
    value = derivative$residual))
derivative_plot <- ggplot(derivative_long, aes(pt, value, color = series)) +
  geom_rect(data = highlight_regions,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = region_id),
    inherit.aes = FALSE, alpha = 0.10, color = NA) +
  geom_hline(yintercept = 0, color = "#888888", linetype = "dotted",
    linewidth = 0.3) +
  geom_line(linewidth = 0.65) +
  geom_vline(data = all_events, aes(xintercept = pt), inherit.aes = FALSE,
    color = "#555555", linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = c("Observed ratio derivative" = "#111111",
    "Fitted from total-signal derivatives" = "#7B3294",
    "Residual derivative" = "#008837"), name = NULL) +
  scale_fill_manual(values = c("T1" = "#E41A1C", "T2" = "#377EB8",
    "T3" = "#4DAF4A"), name = "Scored ±100 windows") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15)) +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  labs(title = "Local ratio changes after total-signal adjustment",
    subtitle = paste0("Model fitted in shaded T1-T3 regions; derivative R² = ",
      sprintf("%.3f", summary_table$R2[[5L]]),
      ". Full pseudotime is displayed."),
    x = "Bdev pseudotime", y = "Central change per window")
ggsave(file.path(stage_dir, "ratio_derivative_total_signal_residual.png"),
  derivative_plot, width = 11, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "ratio_derivative_total_signal_residual.pdf"),
  derivative_plot, width = 11, height = 5.5, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  target = "log2(smoothed delta-RNA/delta-ATAC)",
  predictors = list("smoothed log1p nCount_RNA", "smoothed log1p nCount_peaks"),
  model = "ordinary least squares with intercept",
  full_level_R2 = summary_table$R2[[1L]],
  full_level_10_contiguous_fold_cv_R2 = summary_table$R2[[3L]],
  scored_level_R2 = summary_table$R2[[2L]],
  full_derivative_R2 = summary_table$R2[[4L]],
  scored_derivative_R2 = summary_table$R2[[5L]],
  leave_one_turn_out_centered_derivative_R2 = summary_table$R2[[6L]],
  displayed_range = "entire pseudotime",
  highlighted_regions = "T1-T3 +/-100 windows",
  bundle_mutated = FALSE, capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
report <- c("# Total-signal contribution to the RNA/ATAC ratio", "",
  paste0("A linear model using smoothed log-total RNA and log-total ATAC",
    " explains ", sprintf("%.1f%%", 100 * summary_table$R2[[1L]]),
    " of the full-path ratio level in sample and ",
    sprintf("%.1f%%", 100 * summary_table$R2[[3L]]),
    " under 10 contiguous-fold prediction."),
  paste0("Within the scored T1-T3 regions, the in-sample level R2 is ",
    sprintf("%.1f%%", 100 * summary_table$R2[[2L]]), "."), "",
  paste0("For local changes, total-signal derivatives explain only ",
    sprintf("%.1f%%", 100 * summary_table$R2[[4L]]),
    " of the full-path ratio derivative and ",
    sprintf("%.1f%%", 100 * summary_table$R2[[5L]]),
    " inside the scored regions."),
  paste0("The turn-specific derivative R2 values are ",
    paste(sprintf("%s %.1f%%", per_turn$region_id,
      100 * per_turn$derivative_R2), collapse = ", "), "."),
  paste0("When each whole turn is held out and centered before pooling, R2 is ",
    sprintf("%.3f", summary_table$R2[[6L]]),
    " and correlation is ", sprintf("%.3f", summary_table$correlation[[6L]]),
    ", showing that a depth model learned at other turns does not generalize",
    " to the held-out local shape."), "",
  "## Interpretation and adjustment", "",
  paste("Total signal captures the broad developmental level of the ratio but",
    "does not explain most local T1-T3 changes. Directly subtracting the",
    "full-path fitted curve is available as a sensitivity analysis, but may",
    "remove biological cell-size/transcriptional changes together with depth."),
  paste("The preferred correction is upstream and conditional on pseudotime:",
    "for every RNA-PC and ATAC-LSI dimension, fit latent_dimension ~ spline(pt)",
    "+ log1p(modality depth), retain spline(pt) + residual, then recompute",
    "window centroids, cumulative distances, the ratio, and the gene ranking.",
    "This removes within-stage depth dependence while preserving the accepted",
    "Bdev pseudotime trend. Compare corrected and original ranks as a",
    "sensitivity analysis rather than treating either as ground truth."))
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed total-signal ratio diagnostic at", output_dir, "\n")
