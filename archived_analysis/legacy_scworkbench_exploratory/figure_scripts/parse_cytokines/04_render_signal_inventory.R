# Report-only rendering helpers; sourced by 04_cytokine_signal_inventory.R.
render_inventory <- function(inv, pairs, conditions, out, checked) {
  fmt <- function(x, digits = 3L) {
    if (is.numeric(x)) ifelse(is.na(x), "NA", formatC(x, format = "fg", digits = digits))
    else ifelse(is.na(x), "NA", as.character(x))
  }
  esc <- function(x) {
    x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub('"', "&quot;", x, fixed = TRUE)
  }
  short_id <- function(x) sub("^parse_cytokine\\.", "", x)
  labels <- paste(inv$cytokine, inv$direction)
  status_names <- c(accepted = "Accepted", ucell_transfer_not_confirmed = "UCell not confirmed",
    not_advanced_from_discovery = "Insufficient genes")
  status_colors <- c(accepted = "#23835B", ucell_transfer_not_confirmed = "#C77716",
    not_advanced_from_discovery = "#BBC3CD")
  scored <- inv[inv$selected_genes > 0L, ]
  score_labels <- paste(scored$cytokine, scored$direction,
    ifelse(scored$final_status == "accepted", "", "*"))
  genes_plot <- function() {
    par(mar = c(5, 7, 4, 2))
    cyto <- unique(inv$cytokine); n <- length(cyto)
    plot(NA, xlim = c(.5, 2.5), ylim = c(.3, n + 2.1), xaxt = "n", yaxt = "n",
      xlab = "", ylab = "", main = "All 52 Parse signature directions", bty = "n")
    for (i in seq_len(nrow(inv))) {
      x <- match(inv$direction[i], c("up", "down")); y <- n + 1 - match(inv$cytokine[i], cyto)
      rect(x - .45, y - .43, x + .45, y + .43,
        col = status_colors[[inv$final_status[i]]], border = "white")
      text(x, y, if (inv$selected_genes[i]) paste0(inv$selected_genes[i], " genes") else
        paste0("No set (", inv$full_eligible_genes[i], " eligible)"), cex = .82,
        col = if (inv$final_status[i] == "not_advanced_from_discovery") "#172B42" else "white")
    }
    axis(2, at = n:1, labels = cyto, las = 1, tick = FALSE, cex.axis = .82)
    axis(1, at = 1:2, labels = c("UP", "DOWN"), tick = FALSE)
    legend("top", legend = unname(status_names), fill = unname(status_colors), bty = "n", cex = .74, horiz = TRUE)
    mtext("Green: registered; orange: scored but not accepted; grey: no gene set created.", side = 1, line = 2.5, cex = .77)
  }
  auc_plot <- function() {
    cols <- c("discovery_cv_response_auc", "discovery_cv_specificity_auc", "ucell_response_auc", "ucell_specificity_auc")
    mat <- as.matrix(scored[, cols]); n <- nrow(mat)
    palette <- colorRampPalette(c("#FCF5EF", "#D5E6E8", "#28758A", "#0E354E"))(101)
    par(mar = c(6, 9, 4, 2))
    plot(NA, xlim = c(.5, 4.5), ylim = c(.5, n + .5), xaxt = "n", yaxt = "n",
      xlab = "", ylab = "", main = "Signature discrimination: sample-level AUC", bty = "n")
    for (i in seq_len(n)) for (j in 1:4) {
      value <- mat[i, j]; y <- n + 1 - i
      rect(j - .48, y - .47, j + .48, y + .47, col = palette[1L + round(value * 100)], border = "white")
      text(j, y, sprintf("%.3f", value), cex = .75, col = if (value > .7) "white" else "#172B42")
    }
    axis(2, at = n:1, labels = score_labels, las = 1, tick = FALSE, cex.axis = .73)
    axis(1, at = 1:4, labels = c("CV: control", "CV: rest", "UCell: control", "UCell: rest"), tick = FALSE, cex.axis = .83)
    mtext("AUC 0.5 = chance; 1 = perfect separation. Down scores are oriented toward repression.", side = 1, line = 2.5, cex = .72)
    mtext("* UCell transfer not confirmed. Three target and three control samples; rest = 75 samples.", side = 1, line = 3.6, cex = .72)
  }
  correlation_plot <- function() {
    ids <- scored$signature_id; n <- length(ids)
    mat <- diag(1, n); dimnames(mat) <- list(ids, ids)
    for (i in seq_len(nrow(pairs))) {
      mat[pairs$key_a[i], pairs$key_b[i]] <- pairs$spearman_cor[i]
      mat[pairs$key_b[i], pairs$key_a[i]] <- pairs$spearman_cor[i]
    }
    palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)
    layout(matrix(c(1, 2), nrow = 1), widths = c(12, 1))
    par(mar = c(10, 9, 4, 1))
    image(seq_len(n), seq_len(n), mat[, n:1], zlim = c(-1, 1), col = palette,
      xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "Native Parse UCell score correlations", useRaster = TRUE)
    axis(1, at = seq_len(n), labels = score_labels, las = 2, tick = FALSE, cex.axis = .64)
    axis(2, at = n:1, labels = score_labels, las = 1, tick = FALSE, cex.axis = .64)
    mtext("Spearman correlation across 154,209 cells; raw scores (down sets are not inverted).", side = 3, line = .5, cex = .7)
    par(mar = c(10, .2, 4, 3))
    image(1, seq(-1, 1, length.out = 201), matrix(seq(-1, 1, length.out = 201), nrow = 1),
      col = palette, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    axis(4, at = c(-1, -.5, 0, .5, 1), las = 1, cex.axis = .8)
    layout(1)
  }
  figures <- list(signature_inventory = genes_plot, discrimination_auc = auc_plot,
    signal_correlations = correlation_plot)
  for (name in names(figures)) {
    png(file.path(out, paste0(name, ".png")), width = 1600, height = 1750, res = 150)
    figures[[name]](); dev.off()
  }
  pdf(file.path(out, "inventory_figures.pdf"), width = 10.7, height = 11.7, onefile = TRUE)
  for (fun in figures) fun()
  dev.off()

  dictionary <- c(
    signature_id = "Exact Parse signature identifier; 26 stimuli x two directions.",
    cytokine = "Original source treatment label, including TGB-b1 and c-kit; labels are not silently corrected.",
    direction = "up: genes induced by target; down: genes repressed by target, in both discovery contrasts.",
    final_status = "accepted: registered in Parse; ucell_transfer_not_confirmed: scored candidate failed transfer thresholds; not_advanced_from_discovery: fewer than five eligible genes.",
    discovery_status = "Historical discovery status before UCell transfer confirmation.",
    registered_in_parse = "Whether this gene set is registered in the accepted Parse namespace.",
    modelled_genes = "Genes tested after CPM >=1 in at least three samples; includes technically excluded genes.",
    membership_eligible_modelled_genes = "Modelled genes excluding feature IDs starting mt-, Rpl, Rps.",
    de_genes_vs_control = "Eligible modelled genes with direction-signed log2FC >=0.5 and gene-level BH FDR <=0.05 versus control.",
    de_genes_vs_rest = "Eligible modelled genes with direction-signed log2FC >=0.25 and gene-level BH FDR <=0.05 versus the equal-weight mean of 25 other individual stimuli.",
    de_genes_both_contrasts = "Intersection of the two DE counts, before three-fold direction stability.",
    full_eligible_genes = "Genes passing both DE contrasts AND sign consistency in all three training folds.",
    selected_genes = "Genes in the actual scored candidate set. Zero means no set was created, not zero cytokine activity.",
    selected_gene_ids = "All selected feature IDs in saved rank order, separated by semicolons.",
    top_genes = "Leading gene IDs from the saved discovery/confirmation status table.",
    discovery_cv_response_auc = "Held-out sample pseudobulk rank-score AUC: target versus control; down scores are sign-inverted.",
    discovery_cv_specificity_auc = "Held-out sample pseudobulk rank-score AUC: target versus all 75 samples from the other 25 individual stimuli.",
    discovery_cv_balanced_auc = "Arithmetic mean of the two discovery AUCs; not balanced accuracy.",
    discovery_fold_gene_jaccard = "Saved mean pairwise Jaccard of selected-size training gene sets across three folds.",
    ucell_response_auc = "AUC of direction-oriented per-sample median UCell: three target versus three control samples.",
    ucell_specificity_auc = "AUC of direction-oriented per-sample median UCell: three target versus 75 other-single samples.",
    ucell_balanced_auc = "Arithmetic mean of the two UCell AUCs; not balanced accuracy.",
    program_scope = "broad_across_states when selected genes' median state concordance >=0.5; state_enriched otherwise. Supplementary saved Parse-state evidence.",
    selected_gene_median_response_logFC = "Median signed log2 fold-change of selected genes: target minus control; negative for down sets.",
    selected_gene_median_specificity_logFC = "Median signed log2 fold-change of selected genes: target minus equal-weight mean of other 25 single stimuli.",
    selected_gene_max_response_FDR = "Largest gene-level BH FDR among selected genes in target-versus-control contrast. NOT a signature-level FDR.",
    selected_gene_max_specificity_FDR = "Largest gene-level BH FDR among selected genes in target-versus-rest contrast. NOT a signature-level FDR.",
    selected_gene_median_state_concordance = "Median fraction of evaluable Parse states in which both gene effects have the expected sign. Descriptive composition sensitivity, not an adjusted causal estimate.",
    selected_gene_median_target_dominance_ratio = "Median abs(target response log2FC)/(maximum abs(other cytokine response log2FC)+0.25), across selected genes. Used in ranking, not a strict acceptance cutoff.",
    target_samples = "Number of target biological samples (3).", control_samples = "Number of control biological samples (3).",
    rest_samples = "Number of other-single biological samples (75); excludes control and combinations.",
    target_cells = "QC-accepted target cells.", control_cells = "QC-accepted control cells.", rest_cells = "QC-accepted cells under the other 25 individual stimuli.",
    scored_samples = "Samples with a saved median UCell value; NA if the direction was not scored.",
    zero_sample_medians = "Number of sample medians exactly zero (does not mean every cell has zero score).",
    nonzero_sample_medians = "Number of sample medians greater than zero.",
    minimum_sample_median = "Minimum raw sample median UCell across 81 samples.", maximum_sample_median = "Maximum raw sample median UCell across 81 samples.",
    control_group_median = "Median of the three control sample medians; raw UCell.", target_group_median = "Median of the three target sample medians; raw UCell.",
    other_single_group_median = "Median of all 75 other-single sample medians; raw UCell.",
    target_sample_median_iqr = "Q75-Q25 across three target sample medians; R quantile type 7. Descriptive spread, not a confidence interval.",
    control_sample_median_iqr = "Q75-Q25 across three control sample medians; type 7.",
    rest_sample_median_iqr = "Q75-Q25 across 75 other-single sample medians; type 7.",
    pooled_cell_score_iqr = "Native pooled cell-level UCell Q75-Q25 across all 154209 analysis cells. Mixes within-sample and between-treatment variation.",
    pooled_cell_score_sd = "Native pooled cell-level UCell standard deviation across all 154209 cells.",
    target_oriented_shift = "Direction sign x (target group median - control group median); positive means the expected cytokine response.",
    top_offtarget_condition = "Other stimulus with the largest oriented score shift from control; alphabetical choice if tied; can be zero or negative.",
    top_offtarget_oriented_shift = "Largest off-target oriented score shift; same signature, not a cross-signature score comparison.",
    top_offtarget_ties = "All conditions tied for largest off-target shift, absolute tolerance 1e-12.",
    target_minus_top_offtarget_shift = "Target shift minus strongest off-target shift; positive means the target exceeds every other condition median.",
    other_conditions_stronger = "Other single stimuli with oriented shift greater than target by >1e-12.",
    other_conditions_tied = "Other stimuli within absolute 1e-12 of the target shift.",
    target_response_rank = "1 + number of stronger other conditions, among 26 stimuli; shared ranks retain ties.",
    most_overlapping_signature = "Other scored candidate with greatest gene Jaccard; alphabetical tie-break. May be an unaccepted candidate.",
    maximum_gene_jaccard = "Maximum intersection/union of selected gene memberships over all other 30 scored candidates.",
    most_correlated_signature = "Other scored candidate with greatest absolute native raw-score Spearman correlation; alphabetical tie-break. May be unaccepted.",
    strongest_absolute_spearman_signed = "Signed Spearman correlation for the strongest absolute correlation; native cell-level descriptive metric, not replicate-level inference.",
    exclusion_stage = "Saved discovery/UCell stage and acceptance classification.", exclusion_reason = "Saved specific reason for exclusion, or acceptance explanation."
  )
  stopifnot(setequal(names(dictionary), names(inv)))
  write.table(data.frame(field = names(inv), definition = unname(dictionary[names(inv)])),
    file.path(out, "metric_dictionary.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
  methods <- c(
    "Scope: all 26 individual stimulus labels, 52 up/down directions. The pooled Parse analysis universe contains 154,209 QC-accepted cells from 81 biological samples (three per condition). Discovery and confirmation use biological samples as the unit; cell counts are not the replicate count. The bundle also contains five combined treatments, listed separately in all_conditions_scope.tsv; these have no individual signatures inferred here.",
    "Test versus control: target stimulus minus unstimulated control. Test versus rest: target minus the equal-weight mean of the other 25 individual stimuli. Rest excludes control and combinations. The specificity AUC separately compares the three target sample scores against all 75 other-single sample scores. These are not 25 separate pairwise differential-expression tests.",
    "Gene evidence: saved sample pseudobulk limma-voom contrasts. Eligible genes require direction-signed log2FC >=0.5 versus control and >=0.25 versus rest, BH FDR <=0.05 for both, and concordant directions in all three training folds; mitochondrial/ribosomal feature prefixes are excluded. Counts before and after each requirement are retained. Reported FDRs are gene-level, not signature-level.",
    "Discovery: three leave-one-replicate-label-out folds; candidate size is chosen using these same folds. Training-fold raw-p cutoffs could relax from 0.05 to 0.20 when fewer than 50 genes qualified; full-data gene FDR cutoffs remain 0.05. Discovery acceptance requires both AUCs >=0.60, their mean >=0.65, and fold-gene Jaccard >=0.10. The reported CV performance was also used for selection and is not a nested, independent test estimate.",
    "UCell transfer: saved scWorkbench scores on RNA_counts, max_rank=1500, without graph smoothing or cross-signature calibration. Confirmation uses per-biological-sample medians and requires both AUCs >=0.60 and their mean >=0.65. It checks transfer to cell-level scoring within the same source dataset; it is not external validation. AUC 0.5 indicates chance ordering and 1 perfect ordering; three target samples limit precision.",
    "Direction: raw UCell always measures expression of the genes in the set. A high raw DOWN score means higher expression of genes that the discovery stimulus normally represses, not greater repression. Only response-oriented shifts and AUCs invert down scores. No signed up-minus-down activity score is introduced.",
    "Comparisons: condition profiles compare each signature under all 26 stimuli; target rank and off-target gaps summarize those saved values. Native Pearson/Spearman correlations compare raw scores across all 154,209 pooled cells and can reflect treatment, cell state, composition and shared genes. Gene Jaccard describes membership overlap. These are descriptive metrics, not tests of unique or causal cytokine activity; raw score magnitudes are not calibrated across gene sets.",
    "Spread: IQR is Q75-Q25. Sample-level IQRs describe three target/control medians or 75 rest medians; pooled cell-level IQR and SD describe all analysis cells. Neither is a confidence interval. Missing metrics (NA) for unscored directions are unavailable, not evidence of no biological response.",
    "Dataset boundary: acceptance here refers to the Parse bundle. IFN-gamma up remains accepted in Parse even though its later Bdev gene-coverage check excluded it from that integration. Original source labels such as TGB-b1 are preserved. No new inference, bundle promotion, or Bdev analysis is performed by this report."
  )
  shared <- inv[inv$registered_in_parse & !is.na(inv$other_conditions_stronger) & inv$other_conditions_stronger > 0L, ]
  ap <- pairs[pairs$signature_a_status == "accepted" & pairs$signature_b_status == "accepted", ]
  ap <- ap[order(-ap$gene_jaccard, ap$key_a, ap$key_b), ]
  findings <- c(sprintf(
    "%d accepted signatures have a stronger median response under another stimulus: %s. Acceptance therefore does not imply a unique response to the named stimulus; specificity against the pooled rest can remain high.",
    nrow(shared), paste(paste(shared$cytokine, shared$direction), collapse = ", ")),
    sprintf("The largest gene-membership overlap among accepted sets is %s versus %s: %d shared genes, Jaccard %.3f. The paired raw-score Spearman correlation is %.3f.",
      short_id(ap$key_a[1]), short_id(ap$key_b[1]), ap$shared_genes[1], ap$gene_jaccard[1], ap$spearman_cor[1]))
  writeLines(c("# Parse cytokine signal inventory", "", sprintf(
    "26 stimuli; 52 directions; %d accepted; %d UCell failures; %d insufficient-gene directions.",
    sum(inv$registered_in_parse), sum(inv$final_status == "ucell_transfer_not_confirmed"),
    sum(inv$final_status == "not_advanced_from_discovery")), "",
    "[Interactive report](REPORT.html) | [Complete CSV](cytokine_signal_inventory.csv) | [Figures PDF](inventory_figures.pdf)", "",
    "| Stimulus | Direction | Selected genes | Stable eligible genes | UCell AUC vs control | UCell AUC vs rest | Status |",
    "|---|---|---:|---:|---:|---:|---|",
    vapply(seq_len(nrow(inv)), function(i) paste0("| ", paste(c(inv$cytokine[i], inv$direction[i],
      inv$selected_genes[i], inv$full_eligible_genes[i], fmt(inv$ucell_response_auc[i]),
      fmt(inv$ucell_specificity_auc[i]), status_names[[inv$final_status[i]]]), collapse = " | "), " |"), character(1)),
    "", "## Comparisons", "", paste0(findings, "\n"),
    "", "## Interpretation", "", paste0(methods, "\n"), "## Provenance", "",
    paste("Bundle:", checked$bundle_dir), "", paste("Accepted fingerprint:", checked$bundle_state_fingerprint), "",
    "All source files and report scripts are hashed in input_manifest.tsv. validation.json records alignment and membership checks; output_manifest.tsv hashes the final outputs."), file.path(out, "REPORT.md"))

  table_cols <- c("cytokine", "direction", "final_status", "selected_genes", "full_eligible_genes",
    "de_genes_vs_control", "de_genes_vs_rest", "ucell_response_auc", "ucell_specificity_auc",
    "discovery_cv_balanced_auc", "discovery_fold_gene_jaccard", "selected_gene_median_response_logFC",
    "selected_gene_median_specificity_logFC", "target_response_rank", "other_conditions_tied",
    "top_offtarget_condition", "maximum_gene_jaccard")
  table_labels <- c("Stimulus", "Direction", "Status", "Selected genes", "Stable eligible genes",
    "DE genes: control", "DE genes: rest", "UCell AUC: control", "UCell AUC: rest", "Discovery mean AUC",
    "Fold Jaccard", "Median log2FC: control", "Median log2FC: rest", "Target rank /26", "Other ties",
    "Top off-target", "Max gene overlap")
  header <- paste0('<th tabindex="0" title="', esc(dictionary[table_cols]), '">', table_labels, '</th>', collapse = "")
  rows <- vapply(seq_len(nrow(inv)), function(i) {
    cells <- vapply(table_cols, function(field) {
      value <- inv[[field]][i]
      if (field == "final_status") value <- status_names[[value]]
      if (field == "cytokine") return(paste0('<td><a href="#', inv$signature_id[i], '">', esc(value), '</a></td>'))
      paste0('<td>', esc(fmt(value)), '</td>')
    }, character(1))
    paste0('<tr data-status="', inv$final_status[i], '">', paste(cells, collapse = ""), '</tr>')
  }, character(1))
  details <- vapply(seq_len(nrow(inv)), function(i) {
    fields <- setdiff(names(inv), c("signature_id", "cytokine", "direction", "selected_gene_ids"))
    values <- vapply(fields, function(f) paste0('<tr><th title="', esc(dictionary[[f]]), '">', esc(f),
      '</th><td>', esc(fmt(inv[[f]][i], 5L)), '</td></tr>'), character(1))
    paste0('<details id="', inv$signature_id[i], '"><summary>', esc(labels[i]), ' — ',
      esc(status_names[[inv$final_status[i]]]), ' — ', inv$selected_genes[i], ' selected genes</summary><p>',
      esc(inv$exclusion_reason[i]), '</p><p><b>Selected genes:</b> ', esc(inv$selected_gene_ids[i]),
      '</p><table class="detail">', paste(values, collapse = ""), '</table></details>')
  }, character(1))
  downloads <- c(cytokine_signal_inventory.csv = "All 52 directions (CSV)", cytokine_signal_inventory.tsv = "All 52 directions (TSV)",
    accepted_signals.tsv = "23 accepted signatures", selected_gene_evidence.tsv = "Selected genes: effects, FDR and ranks",
    eligible_gene_evidence.tsv.gz = "All stable eligible genes, including unselected candidates", signal_pairwise_comparisons.tsv = "All 465 signal pairs: correlations and gene overlap",
    condition_score_profiles.tsv = "806 signature-by-condition score profiles", condition_expression_profiles.tsv = "806 signature-by-condition gene-effect profiles",
    sample_median_scores.tsv = "2,511 per-sample median scores", cell_score_quantiles.tsv = "Cell-level quantiles by target/control/rest",
    all_conditions_scope.tsv = "All 32 source conditions and inclusion scope", analysis_samples.tsv = "81 analysis samples and cell counts",
    metric_dictionary.tsv = "Metric definitions", inventory_figures.pdf = "Three-page figure PDF", validation.json = "Validation and provenance", input_manifest.tsv = "Input checksums")
  dl_html <- paste0('<li><a href="', names(downloads), '">', esc(downloads), '</a></li>', collapse = "\n")
  methods_html <- paste0('<p>', esc(methods), '</p>', collapse = "\n")
  dictionary_html <- paste0('<tr><th>', names(dictionary), '</th><td>', esc(dictionary), '</td></tr>', collapse = "\n")
  html <- c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
    '<title>Parse cytokine signal inventory</title><style>',
    'body{font:15px/1.5 system-ui,sans-serif;color:#183047;background:#f6f8fb;margin:0}main{max-width:1450px;margin:auto;padding:28px}h1{margin-bottom:6px}p{max-width:1100px}a{color:#096d97}.card,details{background:white;padding:16px;border:1px solid #dce3eb;border-radius:8px;margin:16px 0}.stats{font-size:20px;font-weight:650}.scroll{overflow:auto;max-height:720px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;padding:8px;border-bottom:1px solid #dce3eb}thead th{position:sticky;top:0;background:#eaf0f6;cursor:pointer;min-width:85px}td{font-variant-numeric:tabular-nums}tr[data-status="accepted"] td:first-child{border-left:5px solid #23835b}tr[data-status="ucell_transfer_not_confirmed"] td:first-child{border-left:5px solid #c77716}tr[data-status="not_advanced_from_discovery"] td:first-child{border-left:5px solid #bbc3cd}.controls{display:flex;gap:12px;flex-wrap:wrap;padding:12px 0}input,select{font:inherit;padding:9px;border:1px solid #aebdcc;border-radius:4px}figure{margin:20px 0}figure img{max-width:100%;width:850px;background:white}figcaption{max-width:950px}summary{cursor:pointer;font-weight:600}.detail th{width:330px;overflow-wrap:anywhere}details p{overflow-wrap:anywhere}.detail td{overflow-wrap:anywhere}.small{font-size:13px;color:#536579}code{overflow-wrap:anywhere}',
    '</style></head><body><main><h1>Parse cytokine signal inventory</h1>',
    '<p class="stats">26 individual stimuli · 52 directions · 23 accepted · 8 UCell failures · 21 insufficient-gene directions</p>',
    '<p>All available Parse discovery directions, including excluded candidates. Click column headings to sort, search to filter, and a stimulus name for its complete record. Hover over headings for metric definitions.</p>',
    '<p><b>Rest</b> = the other 25 individual stimuli. <b>High raw down-set score</b> = more expression of the repressed genes; only response-oriented metrics invert that sign. NA means unavailable.</p>',
    '<div class="card"><div class="controls"><input id="search" aria-label="Search signatures" placeholder="Search stimulus or direction"><select id="status" aria-label="Filter by status"><option value="">All statuses</option><option value="accepted">Accepted</option><option value="ucell_transfer_not_confirmed">UCell not confirmed</option><option value="not_advanced_from_discovery">Insufficient genes</option></select><span id="count">52 directions</span></div>',
    paste0('<div class="scroll"><table id="inventory"><thead><tr>', header, '</tr></thead><tbody>', paste(rows, collapse = "\n"), '</tbody></table></div></div>'),
    '<h2>Downloads</h2>', paste0('<ul>', dl_html, '</ul>'),
    '<h2>Comparison highlights</h2>', paste0('<p>', esc(findings), '</p>'),
    '<h2>Figures</h2>',
    '<figure><img src="signature_inventory.png" alt="Status and gene count for every up and down direction"><figcaption>Every individual stimulus is retained in the inventory. Orange gene counts refer to scored candidates that failed UCell confirmation; grey directions have no scored set.</figcaption></figure>',
    '<figure><img src="discrimination_auc.png" alt="Discovery and UCell AUC values for all 31 scored signatures"><figcaption>Sample-level discrimination against control and rest. Discovery selection and UCell confirmation use the same source dataset; neither is independent external validation.</figcaption></figure>',
    '<figure><img src="signal_correlations.png" alt="Pairwise native raw UCell Spearman correlations"><figcaption>Correlations across all pooled analysis cells; down scores are not inverted. Shared expression, cell composition and shared gene membership can contribute. * indicates an unaccepted candidate.</figcaption></figure>',
    '<h2>Interpretation and methods</h2>', methods_html,
    '<h2>Complete signature records</h2>', details,
    '<details><summary>Metric dictionary</summary><table class="detail">', dictionary_html, '</table></details>',
    '<h2>Source and validation</h2>', paste0('<p>Accepted Parse bundle: <code>', esc(checked$bundle_dir), '</code><br>Fingerprint: <code>', esc(checked$bundle_state_fingerprint), '</code></p>'),
    '<p class="small">Bundle opened read-only. Gene memberships, sample identities, cell counts, saved gene effects and score coverage are checked during rendering. Source and output checksums, dependency versions and validation results accompany this report.</p>',
    '<script>',
    'const body=document.querySelector("#inventory tbody");const search=document.querySelector("#search");const status=document.querySelector("#status");function filter(){let n=0;for(const row of body.rows){const keep=row.textContent.toLowerCase().includes(search.value.toLowerCase())&&(!status.value||row.dataset.status===status.value);row.hidden=!keep;if(keep)n++;}document.querySelector("#count").textContent=n+" directions";}search.addEventListener("input",filter);status.addEventListener("change",filter);document.querySelectorAll("#inventory th").forEach((th,index)=>{let ascending=true;const sort=()=>{const rows=Array.from(body.rows);rows.sort((a,b)=>{const x=a.cells[index].textContent.trim(),y=b.cells[index].textContent.trim();if(x==="NA"||y==="NA")return x===y?0:x==="NA"?1:-1;const result=x!==""&&y!==""&&Number.isFinite(Number(x))&&Number.isFinite(Number(y))?Number(x)-Number(y):x.localeCompare(y);return ascending?result:-result;});rows.forEach(row=>body.append(row));ascending=!ascending;};th.addEventListener("click",sort);th.addEventListener("keydown",e=>{if(e.key==="Enter")sort();});});document.querySelectorAll("a").forEach(a=>{if((a.getAttribute("href")||"").startsWith("#parse_cytokine"))a.addEventListener("click",()=>{document.getElementById(a.getAttribute("href").slice(1)).open=true;});});',
    '</script></main></body></html>')
  writeLines(html, file.path(out, "REPORT.html"))
}
