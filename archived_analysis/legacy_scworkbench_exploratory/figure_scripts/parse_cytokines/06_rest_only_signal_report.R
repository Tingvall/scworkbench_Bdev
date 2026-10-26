#!/usr/bin/env Rscript
source("src/parse_rest_only_workflow.R")
main <- function() {
  stopifnot(file.exists(file.path(report, "scoring.json")), !file.exists(file.path(report, "report.json")))
  scoring <- read_json(file.path(report, "scoring.json"), simplifyVector = TRUE)
  status <- read_tsv(file.path(report, "all_signature_status.tsv"))
  samples <- read_tsv(file.path(report, "ucell_sample_medians.tsv"))
  profile <- read_tsv(file.path(report, "condition_profiles.tsv"))
  genes <- read_tsv(file.path(report, "selected_gene_evidence.tsv"))
  sets <- read_tsv(file.path(report, "rest_only_gene_sets.tsv"))
  old_path <- file.path(parent_cfg$signature_analysis$final$report, "all_signature_status.tsv")
  old_genes_path <- file.path(parent_cfg$signature_analysis$discovery$report, "selected_signature_members.tsv")
  old <- read_tsv(old_path); old_genes <- read_tsv(old_genes_path)
  passed <- status$final_status == "passes_rest_only_checks"
  ids <- status$signature_id[passed]
  stopifnot(nrow(status) == 52L, length(ids) == scoring$passed_sets,
    setequal(unique(sets$feature_set), ids), length(unique(samples$signature_id)) == scoring$scored_sets,
    nrow(samples) == scoring$scored_sets * 81L, nrow(profile) == scoring$scored_sets * 26L,
    !anyDuplicated(samples[c("signature_id", "sample_id")]),
    !anyDuplicated(profile[c("signature_id", "condition")]))
  for (id in unique(samples$signature_id)) {
    s <- samples[samples$signature_id == id & samples$treatment_category == "single", ]
    entry <- status[status$signature_id == id, ]
    stopifnot(nrow(s) == 78L, sum(s$condition == entry$cytokine) == 3L,
      isTRUE(all.equal(auc(s$oriented_sample_median, s$condition == entry$cytokine), entry$ucell_rest_auc, tolerance = 1e-12)))
  }
  stage <- wb_open_bundle(scoring$staged_bundle, read_only = TRUE, validate = FALSE)
  on.exit(wb_close_bundle(stage), add = TRUE)
  registered <- wb_feature_sets(stage, matrix_id = "RNA_counts", namespace = cfg$namespace)
  stopifnot(setequal(registered$feature_set, ids))
  for (id in ids) stopifnot(identical(wb_feature_set_get(stage, id, "RNA_counts"), sets$feature_id[sets$feature_set == id]))
  check <- function(x, cutoff) ifelse(is.na(x), "NOT_EVALUATED", ifelse(x >= cutoff, "PASS", "FAIL"))
  status$minimum_genes_check <- check(status$eligible_genes, cfg$analysis$minimum_genes)
  status$discovery_rest_check <- check(status$cv_rest_auc, cfg$analysis$discovery_rest_auc_minimum)
  status$fold_stability_check <- check(status$fold_gene_jaccard, cfg$analysis$fold_gene_jaccard_minimum)
  status$ucell_check[is.na(status$ucell_check)] <- "NOT_EVALUATED"
  for (stat in c("min", "max")) for (field in c("specificity_p_value", "specificity_FDR")) {
    status[[paste0("selected_gene_", field, "_", stat)]] <- vapply(status$signature_id, function(id) {
      v <- genes[genes$signature_id == id, field]
      if (!length(v)) NA_real_ else if (stat == "min") min(v) else max(v)
    }, numeric(1))
  }
  matching <- match(paste(status$cytokine, status$direction), paste(old$cytokine, old$direction))
  stopifnot(!anyNA(matching))
  status$old_signature_id <- old$signature_id[matching]
  status$old_status <- old$final_status[matching]
  status$old_selected_genes <- old$selected_genes[matching]
  status$gene_jaccard_with_old_selected_set <- vapply(seq_len(nrow(status)), function(i) {
    new <- genes$feature_id[genes$signature_id == status$signature_id[i]]
    prior <- old_genes$feature_id[old_genes$signature_id == status$old_signature_id[i]]
    if (!length(new) || !length(prior)) NA_real_ else length(intersect(new, prior)) / length(union(new, prior))
  }, numeric(1))
  gained <- passed & status$old_status != "accepted"
  lost <- !passed & status$old_status == "accepted"
  stopifnot(sum(passed) == 23L + sum(gained) - sum(lost))
  write_tsv(status, "signal_checks_and_comparison.tsv")
  write.csv(status, file.path(report, "signal_checks_and_comparison.csv"), row.names = FALSE, na = "NA")
  genes$final_status <- status$final_status[match(genes$signature_id, status$signature_id)]
  write_tsv(genes[genes$signature_id %in% ids, ], "rest_only_gene_evidence.tsv")

  colors <- c(passes_rest_only_checks = "#23835B", ucell_not_confirmed = "#D18B25",
    discovery_fail = "#9275B5", insufficient_genes = "#B7C1CA", insufficient_training_genes = "#B7C1CA")
  draw_inventory <- function() {
    treatments <- unique(status$cytokine); n <- length(treatments)
    par(mar = c(4, 7, 5, 2))
    plot(NA, xlim = c(.5, 2.5), ylim = c(.4, n + 2), xaxt = "n", yaxt = "n",
      xlab = "", ylab = "", main = "Parse signals without a treatment-control requirement", bty = "n")
    for (i in seq_len(nrow(status))) {
      x <- match(status$direction[i], c("up", "down")); y <- n + 1 - match(status$cytokine[i], treatments)
      rect(x - .45, y - .44, x + .45, y + .44, col = colors[[status$final_status[i]]], border = "white")
      text(x, y, if (status$selected_genes[i]) paste(status$selected_genes[i], "genes") else
        paste(status$eligible_genes[i], "eligible; no set"), cex = .8,
        col = if (status$final_status[i] == "insufficient_genes") "#193249" else "white")
    }
    axis(2, at = n:1, labels = treatments, las = 1, tick = FALSE, cex.axis = .8)
    axis(1, at = 1:2, labels = c("UP vs other stimuli", "DOWN vs other stimuli"), tick = FALSE)
    legend("top", c("Passes all new checks", "UCell fails", "Discovery fails", "Insufficient genes"),
      fill = colors[c("passes_rest_only_checks", "ucell_not_confirmed", "discovery_fail", "insufficient_genes")],
      bty = "n", cex = .7, horiz = TRUE)
    mtext("New results are staged for review; none have been promoted.", side = 1, line = 2, cex = .8)
  }
  draw_auc <- function() {
    s <- status[status$selected_genes > 0, ]; n <- nrow(s)
    values <- as.matrix(s[, c("cv_rest_auc", "fold_gene_jaccard", "ucell_rest_auc")])
    palette <- colorRampPalette(c("#FCF5EF", "#CCE3E8", "#28758A", "#0E354E"))(101)
    par(mar = c(5, 9, 4, 2))
    plot(NA, xlim = c(.5, 3.5), ylim = c(.5, n + .5), xaxt = "n", yaxt = "n",
      xlab = "", ylab = "", main = "Rest-only discovery, stability and UCell checks", bty = "n")
    for (i in seq_len(n)) for (j in 1:3) {
      v <- values[i, j]; y <- n + 1 - i
      rect(j - .47, y - .46, j + .47, y + .46, col = palette[1 + round(v * 100)], border = "white")
      text(j, y, sprintf("%.3f", v), cex = .68, col = if (v > .7) "white" else "#193249")
    }
    axis(2, at = n:1, labels = paste(s$cytokine, s$direction,
      ifelse(s$final_status == "passes_rest_only_checks", "", "*")), las = 1, tick = FALSE, cex.axis = .67)
    axis(1, at = 1:3, labels = c("CV rest AUC >=0.60", "Fold Jaccard >=0.10", "UCell rest AUC >=0.60"), tick = FALSE, cex.axis = .78)
    mtext("* UCell confirmation failed. AUC uses 3 target versus 75 other-single samples.", side = 1, line = 2.7, cex = .72)
  }
  draw_profile <- function() {
    s <- status[status$selected_genes > 0, ]; conditions <- sort(unique(profile$condition))
    mat <- matrix(NA_real_, nrow(s), length(conditions), dimnames = list(s$signature_id, conditions))
    mat[cbind(match(profile$signature_id, rownames(mat)), match(profile$condition, conditions))] <- profile$oriented_difference_from_other_single_median
    stopifnot(!anyNA(mat))
    limit <- max(abs(mat)); if (limit == 0) limit <- 1
    par(mar = c(8, 9, 4, 2))
    image(seq_along(conditions), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]),
      col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101), zlim = c(-limit, limit),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "Specificity profiles across individual stimuli", useRaster = TRUE)
    axis(1, at = seq_along(conditions), labels = conditions, las = 2, tick = FALSE, cex.axis = .7)
    axis(2, at = nrow(mat):1, labels = paste(s$cytokine, s$direction), las = 1, tick = FALSE, cex.axis = .67)
    points(match(s$cytokine, conditions), nrow(mat):1, pch = 0, cex = .85)
    mtext(sprintf("Blue: negative; white: zero; red: positive (range ±%.3f). Square = target.", limit), side = 3, line = .5, cex = .72)
    mtext("Oriented sample-median difference from that signature's other-single median; compare within rows.", side = 1, line = 5.8, cex = .7)
  }
  plots <- list(gene_set_inventory = draw_inventory, rest_only_checks = draw_auc, condition_profiles = draw_profile)
  for (name in names(plots)) {
    png(file.path(report, paste0(name, ".png")), width = 1800,
      height = if (name == "gene_set_inventory") 1850 else 2400, res = 160)
    plots[[name]](); dev.off()
  }
  pdf(file.path(report, "rest_only_figures.pdf"), width = 11.25, height = 15, onefile = TRUE)
  for (draw in plots) draw()
  dev.off()
  pdf(file.path(report, "sample_distributions.pdf"), width = 13, height = 8, onefile = TRUE)
  par(mar = c(9, 5, 4, 2))
  for (id in unique(samples$signature_id)) {
    s <- samples[samples$signature_id == id & samples$treatment_category == "single", ]
    entry <- status[status$signature_id == id, ]
    conditions <- sort(unique(s$condition))
    fill <- ifelse(conditions == entry$cytokine, "#D18B25", "#C9DDEA")
    boxplot(split(s$sample_median_ucell, factor(s$condition, levels = conditions)), col = fill,
      outline = FALSE, las = 2, ylab = "Raw sample-median UCell score", xlab = "", cex.axis = .8,
      main = paste(entry$cytokine, entry$direction, "—", entry$final_status),
      ylim = c(0, max(.01, max(s$sample_median_ucell) * 1.1)))
    points(match(s$condition, conditions) + (as.numeric(s$replicate) - 2) * .12,
      s$sample_median_ucell, pch = 21, bg = "white", cex = .65)
    mtext("Three biological samples per stimulus; boxes show quartile spread, points are sample medians.", side = 1, line = 6.6, cex = .8)
  }
  dev.off()

  fmt <- function(x) if (is.numeric(x)) ifelse(is.na(x), "NA", formatC(x, digits = 3, format = "fg")) else as.character(x)
  esc <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)
  columns <- c("cytokine", "direction", "eligible_genes", "selected_genes", "discovery_rest_check",
    "cv_rest_auc", "fold_stability_check", "fold_gene_jaccard", "ucell_check", "ucell_rest_auc", "final_status")
  display <- lapply(seq_len(nrow(status)), function(i) vapply(columns, function(f) fmt(status[[f]][i]), character(1)))
  labels <- c("Stimulus", "Direction", "Eligible genes", "Selected genes", "Discovery", "CV rest AUC",
    "Stability", "Fold Jaccard", "UCell", "UCell rest AUC", "Result")
  headline <- sprintf("%d sets pass all rest-only checks: %d up and %d down, spanning %d stimuli. All remain staged for review.",
    sum(passed), sum(passed & status$direction == "up"), sum(passed & status$direction == "down"), length(unique(status$cytokine[passed])))
  notes <- c(
    "Control comparisons are omitted from gene filtering, gene ranking, fold-direction requirements, candidate-size choice, discovery acceptance and UCell acceptance. No mean of control/rest AUCs or control-based dominance term is used.",
    "Genes require absolute, direction-consistent target-versus-rest log2FC ≥0.25, BH FDR ≤0.05 and the same sign in all three training folds. Technical feature exclusions and the original modelled gene universe are retained. At least five genes are required.",
    "Ranking uses signed specificity t plus specificity-only state concordance. Fold rankings omit state concordance, as in the original workflow. Candidate sizes are 5, 10, 15, 20, 25, 30, 40 and 50. Sizes within 0.02 of the best rest AUC are ranked by fold membership Jaccard, then smaller size. Training-fold raw-p cutoffs can relax from 0.05 to 0.20 when fewer than 50 genes qualify.",
    "Discovery acceptance requires rest AUC ≥0.60 and mean pairwise fold-gene Jaccard ≥0.10. UCell confirmation requires direction-oriented sample-median rest AUC ≥0.60. These retain the original individual rest thresholds; the original combined AUC ≥0.65 criterion is removed with the control comparison.",
    "Rest is the other 25 individual stimuli: three target versus 75 other-single biological samples. Combined treatments are excluded. All 154,209 source analysis cells are retained and scored; the 4,816 control cells do not enter the new AUCs or condition-comparison figures.",
    "The original full-data and three training-fold specificity fits are reused. Those original fits estimated variance using 81 samples including control; they are not new models fitted after deleting control samples. Existing Parse states also retain their original control-derived HVG construction. Only specificity effects enter state concordance here.",
    "UP and DOWN describe expression relative to the other individual stimuli, not induction/repression relative to unstimulated control. High raw down-set UCell means higher expression of those genes; the sign is inverted only for AUC and oriented differences. These signatures alone do not establish activation by a cytokine.",
    "P-values and FDR in the gene-evidence download are gene-level specificity tests, not signature-level p-values. Discovery CV is also used to choose size and is not nested validation. UCell confirmation uses the same source dataset, with only three target biological samples; it is not independent external validation. IQR is descriptive spread, not a confidence interval.",
    "The accepted Parse bundle and original 23 registered signatures are unchanged. A separate review bundle contains exact source counts and sample metadata; new scores and passing feature sets are held in a durable scWorkbench session. Zero result promotions were performed. Bdev integration was not changed."
  )
  changes <- sprintf("Compared with the original 23 accepted directions, %d additional directions now pass and %d formerly accepted directions do not pass. Gene memberships were reselected; these are not simply the old sets with a threshold removed.", sum(gained), sum(lost))
  md_lines <- c("# Parse cytokine signals without treatment-versus-control selection", "", headline, "",
    "[Passed gene sets (GMT)](rest_only_gene_sets.gmt) | [Gene membership (TSV)](rest_only_gene_sets.tsv) | [All checks and old/new comparison (CSV)](signal_checks_and_comparison.csv) | [Gene p-values/FDR](rest_only_gene_evidence.tsv)", "", changes, "",
    paste0("| ", paste(labels, collapse = " | "), " |"), paste0("|", paste(rep("---", length(labels)), collapse = "|"), "|"),
    vapply(display, function(x) paste0("| ", paste(x, collapse = " | "), " |"), character(1)), "",
    "![Gene set results](gene_set_inventory.png)", "", "![Discovery and UCell checks](rest_only_checks.png)", "",
    "![Other-stimulus profiles](condition_profiles.png)", "", "## Methods and interpretation", "", paste0(notes, "\n"),
    "## Review locations", "", paste("Input review bundle:", cfg$bundle_dir), "",
    paste("Durable session:", scoring$session_workspace), "", paste("Private staged bundle:", scoring$staged_bundle), "",
    "[Figure PDF](rest_only_figures.pdf) | [Per-sample distributions](sample_distributions.pdf) | [Scoring provenance](scoring.json)")
  writeLines(md_lines, file.path(report, "REPORT.md"))
  rows <- vapply(display, function(x) paste0("<tr>", paste0("<td>", esc(x), "</td>", collapse = ""), "</tr>"), character(1))
  writeLines(c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Parse rest-only cytokine signatures</title>',
    '<style>body{font:15px/1.5 system-ui,sans-serif;max-width:1450px;margin:25px auto;padding:0 20px;color:#193249;background:#f7f9fc}p{max-width:1150px}.scroll{overflow:auto;max-height:700px}table{border-collapse:collapse;width:100%;background:white;font-size:13px}th,td{padding:9px;border-bottom:1px solid #d8e0e8;text-align:left}thead th{position:sticky;top:0;background:#eaf0f6}a{color:#076690}img{max-width:100%;width:920px;background:white}input{font:inherit;padding:8px}code{overflow-wrap:anywhere}</style></head><body>',
    '<h1>Parse cytokine signals: treatment versus other stimuli</h1>', paste0('<p><b>', esc(headline), '</b></p>'),
    '<p><a href="rest_only_gene_sets.gmt">Passed gene sets: GMT</a> · <a href="rest_only_gene_sets.tsv">Gene membership: TSV</a> · <a href="signal_checks_and_comparison.csv">Checks and comparison: CSV</a> · <a href="rest_only_gene_evidence.tsv">Gene p-values and FDR</a> · <a href="discovery_gene_sets.gmt">All 47 discovery-passing sets</a></p>',
    paste0('<p>', esc(changes), '</p>'),
    '<p><b>UP/DOWN are relative to the other 25 individual stimuli. They do not demonstrate a change versus unstimulated control.</b></p>',
    '<label>Filter signals: <input type="search" id="search"></label><div class="scroll"><table id="signals"><thead><tr>',
    paste0('<th>', labels, '</th>', collapse = ""), '</tr></thead><tbody>', rows, '</tbody></table></div>',
    '<h2>Figures</h2><p><a href="rest_only_figures.pdf">Overview PDF</a> · <a href="sample_distributions.pdf">All 47 sample-distribution plots</a></p>',
    '<img src="gene_set_inventory.png" alt="Status and number of genes for all 52 directions">',
    '<img src="rest_only_checks.png" alt="Discovery AUC, fold stability and UCell AUC for all 47 scored sets">',
    '<img src="condition_profiles.png" alt="Signature profiles across all individual stimuli">',
    '<h2>Methods and interpretation</h2>', paste0('<p>', esc(notes), '</p>'),
    '<h2>Evidence and review state</h2><p><a href="ucell_sample_medians.tsv">Sample scores</a> · <a href="condition_profiles.tsv">Condition medians and IQR</a> · <a href="cell_distribution_quantiles.tsv">Cell-level spread</a> · <a href="selected_gene_evidence.tsv">All candidate gene evidence</a> · <a href="cross_validation_performance.tsv">All candidate sizes</a> · <a href="native_session_dag.json">Native operation history</a> · <a href="report.json">Report validation</a></p>',
    paste0('<p>Separate input bundle: <code>', esc(cfg$bundle_dir), '</code><br>Private staged results: <code>', esc(scoring$staged_bundle), '</code></p>'),
    '<script>document.getElementById("search").addEventListener("input",function(){const q=this.value.toLowerCase();document.querySelectorAll("#signals tbody tr").forEach(r=>{r.hidden=!r.textContent.toLowerCase().includes(q);});});</script></body></html>'), file.path(report, "REPORT.html"))
  stopifnot(identical(wb_source_sha256(scoring$score_rds), scoring$score_sha256))
  save_provenance("report", c(file.path(report, c("scoring.json", "all_signature_status.tsv", "ucell_sample_medians.tsv",
    "condition_profiles.tsv", "selected_gene_evidence.tsv", "rest_only_gene_sets.tsv")), old_path, old_genes_path,
    "figure_scripts/parse_cytokines/06_rest_only_signal_report.R"),
    list(passed_sets = sum(passed), gained_directions = status$signature_id[gained],
      lost_directions = status$signature_id[lost], registered_staged_memberships_exact = TRUE,
      all_sample_auc_values_recomputed = TRUE, score_payload_sha256_verified = TRUE,
      figures = c("gene_set_inventory.png", "rest_only_checks.png", "condition_profiles.png", "sample_distributions.pdf"),
      review_status = "awaiting_user_review", promotions = 0L))
  files <- list.files(report, full.names = TRUE)
  write_tsv(data.frame(file = basename(files), sha256 = vapply(files, wb_source_sha256, character(1))), "output_manifest.tsv")
  log_event(headline)
  log_event("Report:", file.path(report, "REPORT.html"))
}
main()
