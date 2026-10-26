# Present saved acceptance decisions and gene-level p-values without new tests.
render_signal_checks <- function(inv, members, eligible, out, parameters, ucell_acceptance) {
  result <- function(value, cutoff) ifelse(is.na(value), "NOT_EVALUATED",
    ifelse(value >= cutoff, "PASS", "FAIL"))
  all_checks <- function(x) {
    if (any(x == "FAIL")) "FAIL" else if (any(x == "NOT_EVALUATED")) "NOT_EVALUATED" else "PASS"
  }
  # These are the existing rules in discovery stage 08 and confirmation stage 09.
  rules <- data.frame(
    check = c("minimum_eligible_genes", "discovery_vs_control", "discovery_vs_rest",
      "discovery_mean_auc", "discovery_fold_gene_stability", "ucell_vs_control",
      "ucell_vs_rest", "ucell_mean_auc"),
    source_field = c("full_eligible_genes", "discovery_cv_response_auc", "discovery_cv_specificity_auc",
      "discovery_cv_balanced_auc", "discovery_fold_gene_jaccard", "ucell_response_auc",
      "ucell_specificity_auc", "ucell_balanced_auc"),
    threshold = c(parameters$validation$minimum_signature_genes, .60, .60, .65, .10,
      ucell_acceptance$response_auc_minimum, ucell_acceptance$specificity_auc_minimum,
      ucell_acceptance$balanced_auc_minimum), stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(inv)), function(i) {
    values <- as.numeric(inv[i, rules$source_field])
    data.frame(signature_id = inv$signature_id[i], cytokine = inv$cytokine[i], direction = inv$direction[i],
      check = rules$check, result = result(values, rules$threshold), observed = values,
      minimum_required = rules$threshold, selected_genes = inv$selected_genes[i],
      signature_p_value = NA_real_, p_value_note = "Not computed: this is a threshold check, not a hypothesis test")
  })
  checks <- do.call(rbind, rows)
  summary <- inv[, c("signature_id", "cytokine", "direction", "de_genes_vs_control", "de_genes_vs_rest",
    "de_genes_both_contrasts", "full_eligible_genes", "selected_genes", "final_status")]
  check_matrix <- matrix(checks$result, nrow = nrow(inv), byrow = TRUE,
    dimnames = list(inv$signature_id, rules$check))
  for (name in colnames(check_matrix)) summary[[paste0(name, "_check")]] <- check_matrix[, name]
  summary$discovery_overall_check <- apply(check_matrix[, 2:5, drop = FALSE], 1, all_checks)
  summary$ucell_overall_check <- apply(check_matrix[, 6:8, drop = FALSE], 1, all_checks)
  stopifnot(identical(summary$discovery_overall_check == "PASS", inv$discovery_status == "validated"),
    identical(summary$minimum_eligible_genes_check == "FAIL", inv$discovery_status == "insufficient_genes"),
    identical(summary$ucell_overall_check == "PASS", inv$final_status == "accepted"))

  # The same selected membership is used for both p-value ranges. For directions
  # with no selected set, these values stay NA even if a few eligible genes exist.
  p_fields <- c("response_p_value", "specificity_p_value", "response_FDR", "specificity_FDR")
  for (field in p_fields) for (stat in c("min", "max")) {
    summary[[paste0("selected_gene_", field, "_", stat)]] <- vapply(summary$signature_id, function(id) {
      values <- members[members$signature_id == id, field]
      if (!length(values)) return(NA_real_)
      stopifnot(all(is.finite(values)), all(values >= 0 & values <= 1))
      if (stat == "min") min(values) else max(values)
    }, numeric(1))
  }
  exact <- members[, c("signature_id", "cytokine", "direction", "feature_id", "member_rank",
    "response_logFC", "response_p_value", "response_FDR", "specificity_logFC",
    "specificity_p_value", "specificity_FDR")]
  orientation <- ifelse(exact$direction == "up", 1, -1)
  exact$control_gene_check <- ifelse(orientation * exact$response_logFC >= parameters$modelling$response_absolute_log2fc_minimum &
    exact$response_FDR <= parameters$modelling$full_data_fdr_maximum_for_both_contrasts, "PASS", "FAIL")
  exact$rest_gene_check <- ifelse(orientation * exact$specificity_logFC >= parameters$modelling$specificity_absolute_log2fc_minimum &
    exact$specificity_FDR <= parameters$modelling$full_data_fdr_maximum_for_both_contrasts, "PASS", "FAIL")
  exact$final_signature_status <- inv$final_status[match(exact$signature_id, inv$signature_id)]
  stopifnot(all(exact$control_gene_check == "PASS"), all(exact$rest_gene_check == "PASS"),
    nrow(exact) == sum(summary$selected_genes))
  eligible_only <- eligible[eligible$signature_id %in% inv$signature_id[inv$selected_genes == 0L],
    c("signature_id", "cytokine", "direction", "feature_id", "response_logFC", "response_p_value",
      "response_FDR", "specificity_logFC", "specificity_p_value", "specificity_FDR")]
  write_out <- function(x, file) write.table(x, file.path(out, file), sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
  write_out(checks, "individual_check_results.tsv")
  write_out(summary, "check_summary.tsv")
  write.csv(summary, file.path(out, "check_summary.csv"), row.names = FALSE, na = "NA")
  write_out(exact, "selected_gene_pvalues.tsv")
  write_out(eligible_only, "unselected_direction_eligible_gene_pvalues.tsv")

  esc <- function(x) {
    x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }
  scientific <- function(x) ifelse(is.na(x), "NA", formatC(x, digits = 2, format = "e"))
  value_range <- function(prefix, i) paste(scientific(summary[[paste0(prefix, "_min")]][i]),
    scientific(summary[[paste0(prefix, "_max")]][i]), sep = " to ")
  headings <- c("Stimulus", "Direction", "Selected genes", "Eligible genes / ≥5 check", "Discovery checks",
    "UCell: control", "UCell: rest", "UCell: mean", "Gene p: control", "Gene p: rest",
    "Gene FDR: control", "Gene FDR: rest", "Final status")
  display <- lapply(seq_len(nrow(summary)), function(i) c(summary$cytokine[i], summary$direction[i],
    summary$selected_genes[i], paste(summary$full_eligible_genes[i], summary$minimum_eligible_genes_check[i], sep = " / "),
    summary$discovery_overall_check[i], summary$ucell_vs_control_check[i], summary$ucell_vs_rest_check[i],
    summary$ucell_mean_auc_check[i],
    if (summary$selected_genes[i]) value_range("selected_gene_response_p_value", i) else "NA",
    if (summary$selected_genes[i]) value_range("selected_gene_specificity_p_value", i) else "NA",
    if (summary$selected_genes[i]) value_range("selected_gene_response_FDR", i) else "NA",
    if (summary$selected_genes[i]) value_range("selected_gene_specificity_FDR", i) else "NA",
    summary$final_status[i]))
  explanation <- c(
    "PASS/FAIL applies to the saved acceptance criteria. NOT_EVALUATED means no signature reached that check; it is not a failed statistical test.",
    "The p-value and BH-adjusted FDR columns are ranges across the selected genes, not a single p-value for the signature. Exact gene values are downloadable. UCell AUC, cross-validation AUC and stability checks have no saved signature-level p-values, so none are invented here.",
    "Control = unstimulated control. Rest = equal-weight mean of the other 25 individual stimuli for gene tests; control and combined treatments are excluded. Score AUCs compare target sample scores with control or other-single sample scores.",
    "Each gene must have FDR ≤0.05 in both contrasts, direction-signed log2FC ≥0.5 versus control and ≥0.25 versus rest. Eligible genes also retain both effect signs in all three training folds; at least five are required. Selected genes are the smaller ranked subset actually scored.",
    "Discovery passes only if control AUC ≥0.60, rest AUC ≥0.60, their mean ≥0.65, and fold-gene Jaccard ≥0.10. UCell passes only if control AUC ≥0.60, rest AUC ≥0.60, and their mean ≥0.65. These are the original thresholds; no new tests or criteria were introduced.",
    "Zero selected genes means no set was created. Where one to four eligible genes exist, their actual p-values are retained in the separate unselected-direction download; the selected-set p-values remain NA. Gene-count stages below retain the separate contrast counts.",
    "Discovery and UCell confirmation use the same Parse source dataset, with three target and three control biological samples. A small selected-gene p-value does not establish that the resulting UCell signature passed its transfer checks. Overlap and correlation were descriptive comparisons, without acceptance tests."
  )
  md <- c("# Parse signals: checks, gene counts and p-values", "",
    "23 accepted; 8 failed UCell confirmation; 21 had fewer than five eligible genes.", "",
    "[Spreadsheet](check_summary.csv) | [Individual check outcomes](individual_check_results.tsv) | [Exact selected-gene p-values](selected_gene_pvalues.tsv)", "",
    paste0("| ", paste(headings, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(headings)), collapse = "|"), "|"),
    vapply(display, function(x) paste0("| ", paste(x, collapse = " | "), " |"), character(1)),
    "", paste0(explanation, "\n"))
  writeLines(md, file.path(out, "CHECKS.md"))
  cell <- function(x) {
    cls <- if (x == "PASS" || grepl(" / PASS$", x)) "pass" else if (x == "FAIL" || grepl(" / FAIL$", x)) "fail" else ""
    paste0('<td class="', cls, '">', esc(x), '</td>')
  }
  main_rows <- vapply(display, function(x) paste0("<tr>", paste(vapply(x, cell, character(1)), collapse = ""), "</tr>"), character(1))
  gene_counts <- inv[, c("cytokine", "direction", "de_genes_vs_control", "de_genes_vs_rest", "de_genes_both_contrasts", "full_eligible_genes", "selected_genes")]
  count_rows <- apply(gene_counts, 1, function(x) paste0("<tr>", paste0("<td>", esc(x), "</td>", collapse = ""), "</tr>"))
  detail_html <- vapply(inv$signature_id, function(id) {
    d <- checks[checks$signature_id == id, ]; info <- inv[inv$signature_id == id, ]
    details <- vapply(seq_len(nrow(d)), function(i) paste0("<tr><td>", d$check[i], "</td>", cell(d$result[i]),
      "<td>", if (is.na(d$observed[i])) "NA" else format(d$observed[i], digits = 5),
      "</td><td>≥ ", d$minimum_required[i], "</td><td>NA — not a p-value test</td></tr>"), character(1))
    paste0("<details><summary>", esc(paste(info$cytokine, info$direction)), " — individual checks</summary>",
      "<table><thead><tr><th>Check</th><th>Result</th><th>Observed</th><th>Required</th><th>Signature p-value</th></tr></thead><tbody>",
      paste(details, collapse = ""), "</tbody></table></details>")
  }, character(1))
  writeLines(c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
    '<title>Parse signal checks and p-values</title><style>',
    'body{font:15px/1.5 system-ui,sans-serif;color:#193249;background:#f7f9fc;margin:24px}h1{margin-bottom:5px}p{max-width:1150px}.scroll{overflow:auto;max-height:760px;border:1px solid #d8e0e8;background:white}table{border-collapse:collapse;font-size:13px;width:100%}th,td{padding:9px;text-align:left;border-bottom:1px solid #d8e0e8}thead th{position:sticky;top:0;background:#eaf0f6}td{font-variant-numeric:tabular-nums}.pass{color:#126238;font-weight:650}.fail{color:#a02d1c;font-weight:650}details{margin:12px 0;background:white;border:1px solid #d8e0e8;padding:12px}summary{cursor:pointer;font-weight:600}a{color:#076690}input{font:inherit;padding:8px;margin:10px 0}',
    '</style></head><body><h1>Parse signals: checks, gene counts and p-values</h1>',
    '<p><b>52 directions: 23 accepted · 8 UCell failures · 21 insufficient-gene directions.</b></p>',
    '<p><b>p and FDR are ranges across the selected genes.</b> Signature checks use thresholds and have no saved signature-level p-values. NOT_EVALUATED means the direction did not reach that check.</p>',
    '<p><a href="check_summary.csv">Summary spreadsheet</a> · <a href="individual_check_results.tsv">All eight checks per direction</a> · <a href="selected_gene_pvalues.tsv">Exact p-values for every selected gene</a> · <a href="unselected_direction_eligible_gene_pvalues.tsv">Eligible genes in directions with no selected set</a> · <a href="REPORT.html">Full metrics and figures</a></p>',
    '<label>Filter by stimulus, direction or result: <input id="filter" type="search"></label>',
    '<div class="scroll"><table id="summary"><thead><tr>', paste0('<th>', headings, '</th>', collapse = ""),
    '</tr></thead><tbody>', main_rows, '</tbody></table></div>',
    '<h2>How to read the checks</h2>', paste0('<p>', esc(explanation), '</p>'),
    '<details><summary>Gene counts at each filtering step</summary><div class="scroll"><table><thead><tr>',
    paste0('<th>', c("Stimulus", "Direction", "Pass control filter", "Pass rest filter", "Pass both", "Pass both + fold direction stability", "Selected"), '</th>', collapse = ""),
    '</tr></thead><tbody>', count_rows, '</tbody></table></div></details>',
    '<h2>Each check and its observed value</h2>', detail_html,
    '<p><a href="validation.json">Validation</a> · <a href="input_manifest.tsv">Input provenance</a>. Parse bundle opened read-only; no new fitting or scoring.</p>',
    '<script>document.getElementById("filter").addEventListener("input",function(){const q=this.value.toLowerCase();document.querySelectorAll("#summary tbody tr").forEach(row=>{row.hidden=!row.textContent.toLowerCase().includes(q);});});</script>',
    '</body></html>'), file.path(out, "CHECKS.html"))
}
