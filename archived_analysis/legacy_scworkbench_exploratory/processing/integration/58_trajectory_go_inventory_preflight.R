#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 58_trajectory_go_inventory_preflight.R",
  "<bundle_dir> <cytokine_inventory.tsv> <cytokine_membership.tsv> <output_dir>"))

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
cytokine_inventory_path <- normalizePath(args[[2L]], mustWork = TRUE)
cytokine_membership_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(jsonlite)
  library(msigdbr)
})

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE)
  write.table(as.data.frame(x), con, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ..., "\n")
  flush.console()
}

minimum_mapped_genes <- 15L
maximum_mapped_genes <- 300L
group_specs <- list(
  bp = list(label = "GO Biological Process", subcollection = "GO:BP",
    namespace = "msigdb_m5_go_bp_2026_1_mm"),
  mf = list(label = "GO Molecular Function", subcollection = "GO:MF",
    namespace = "msigdb_m5_go_mf_2026_1_mm"),
  cc = list(label = "GO Cellular Component", subcollection = "GO:CC",
    namespace = "msigdb_m5_go_cc_2026_1_mm"))

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells <- as.character(wb_cells(bundle, "all"))
metadata <- as.data.table(wb_meta(bundle, cells = cells))
matrix_features <- as.character(wb_feature_ids(bundle, "RNA_log1p_cp10k"))
stopifnot(length(cells) == 26723L,
  identical(as.character(metadata$cell_id), cells),
  sum(is.finite(metadata$pt)) == 12642L, sum(!is.finite(metadata$pt)) == 14081L,
  length(matrix_features) > 10000L, !anyDuplicated(matrix_features))

all_inventory <- list()
eligible_memberships <- list()
source_versions <- list()

for (group in names(group_specs)) {
  spec <- group_specs[[group]]
  group_dir <- file.path(stage_dir, group)
  dir.create(group_dir)
  log_event("Retrieving", spec$label)
  go <- as.data.table(msigdbr::msigdbr(db_species = "MM",
    species = "Mus musculus", collection = "M5",
    subcollection = spec$subcollection))
  required <- c("gene_symbol", "gs_name", "gs_description", "db_version",
    "gs_collection", "gs_subcollection")
  stopifnot(all(required %in% names(go)), nrow(go) > 0L,
    uniqueN(go$db_version) == 1L,
    identical(unique(go$gs_collection), "M5"),
    identical(unique(go$gs_subcollection), spec$subcollection))
  source_versions[[group]] <- unique(go$db_version)
  go <- unique(go[, ..required])
  setorder(go, gs_name, gene_symbol)
  source_inventory <- go[, .(
    n_genes_raw = uniqueN(gene_symbol),
    n_genes_mapped_to_bundle = uniqueN(
      gene_symbol[gene_symbol %chin% matrix_features]),
    description = first(gs_description)), by = .(feature_set = gs_name)]
  source_inventory[, mapped_fraction := n_genes_mapped_to_bundle / n_genes_raw]
  source_inventory[, preflight_eligible :=
    n_genes_mapped_to_bundle >= minimum_mapped_genes &
      n_genes_mapped_to_bundle <= maximum_mapped_genes]
  eligible_source_sets <- source_inventory[preflight_eligible == TRUE,
    feature_set]
  stopifnot(length(eligible_source_sets) > 0L)
  write_tsv_gz(source_inventory,
    file.path(group_dir, "source_mapping_inventory.tsv.gz"))

  # A three-term scWorkbench smoke test established compatibility of the
  # documented feature_set_import/feature_set_overlap_qc route. A full governed
  # BP import was benchmarked and stopped after registering 1,145/4,087 eligible
  # sets in ~13 minutes. The user-approved vectorized outside workflow therefore
  # performs the identical exact-ID intersection for this full inventory.
  inventory <- copy(source_inventory)
  inventory[, `:=`(analysis_group = group, collection_label = spec$label,
    set_name = feature_set,
    n_genes_detected_in_matrix = n_genes_mapped_to_bundle,
    eligible = preflight_eligible)]
  inventory[, exclusion_reason := fifelse(eligible, "eligible",
    fifelse(n_genes_mapped_to_bundle < minimum_mapped_genes,
      "fewer_than_15_mapped_genes",
      "more_than_300_mapped_genes"))]
  setcolorder(inventory, c("analysis_group", "collection_label", "feature_set",
    "set_name", "description", "n_genes_raw", "n_genes_mapped_to_bundle",
    "n_genes_detected_in_matrix", "mapped_fraction", "eligible",
    "exclusion_reason"))
  write_tsv_gz(inventory, file.path(group_dir, "set_inventory.tsv.gz"))
  all_inventory[[group]] <- inventory

  eligible_sets <- inventory[eligible == TRUE, feature_set]
  memberships <- go[gs_name %in% eligible_sets &
    gene_symbol %chin% matrix_features,
    .(analysis_group = group, feature_set = gs_name, feature_id = gene_symbol)]
  memberships <- unique(memberships)
  stopifnot(uniqueN(memberships$feature_set) == length(eligible_sets),
    all(memberships$feature_id %in% matrix_features))
  write_tsv_gz(memberships,
    file.path(group_dir, "eligible_mapped_memberships.tsv.gz"))
  eligible_memberships[[group]] <- memberships
  write_json(list(status = "complete", analysis_group = group,
    collection_label = spec$label, source_version = source_versions[[group]],
    source_sets = nrow(inventory),
    mapping_method = "exact gene-symbol intersection with scWorkbench RNA_log1p_cp10k feature IDs",
    eligible_sets = length(eligible_sets),
    excluded_sets = nrow(inventory) - length(eligible_sets),
    minimum_mapped_genes = minimum_mapped_genes,
    maximum_mapped_genes = maximum_mapped_genes,
    separate_group_multiple_testing = TRUE, promoted = FALSE),
    file.path(group_dir, "completion.json"))
  log_event("Completed", group, "preflight:", length(eligible_sets), "eligible")
}

log_event("Preparing reviewed cytokine inventory")
cytokine_inventory <- fread(cytokine_inventory_path)
cytokine_members <- fread(cytokine_membership_path)
stopifnot(nrow(cytokine_inventory) == 38L,
  !anyDuplicated(cytokine_inventory$signature_id),
  setequal(unique(cytokine_members$signature_id),
    cytokine_inventory$signature_id))
cytokine_counts <- cytokine_members[, .(
  n_genes_raw = uniqueN(feature_id),
  n_genes_mapped_to_bundle = uniqueN(feature_id[feature_id %in% matrix_features])),
  by = .(signature_id, cytokine, direction)]
cytokine_standard <- merge(cytokine_inventory[, .(signature_id, cytokine,
  direction)], cytokine_counts, by = c("signature_id", "cytokine", "direction"))
cytokine_standard[, `:=`(analysis_group = "cytokines",
  collection_label = "Reviewed cytokine signatures",
  feature_set = signature_id, set_name = paste(cytokine, direction),
  description = paste("Reviewed", cytokine, direction, "response signature"),
  n_genes_detected_in_matrix = n_genes_mapped_to_bundle,
  mapped_fraction = n_genes_mapped_to_bundle / n_genes_raw,
  eligible = n_genes_mapped_to_bundle >= minimum_mapped_genes &
    n_genes_mapped_to_bundle <= maximum_mapped_genes)]
cytokine_standard[, exclusion_reason := fifelse(eligible, "eligible",
  fifelse(n_genes_mapped_to_bundle < minimum_mapped_genes,
    "fewer_than_15_mapped_genes", "more_than_300_mapped_genes"))]
setcolorder(cytokine_standard, c("analysis_group", "collection_label",
  "feature_set", "set_name", "description", "n_genes_raw",
  "n_genes_mapped_to_bundle", "n_genes_detected_in_matrix", "mapped_fraction",
  "eligible", "exclusion_reason"))
all_inventory$cytokines <- cytokine_standard
cytokine_eligible <- cytokine_standard[eligible == TRUE, feature_set]
cytokine_mapped <- cytokine_members[signature_id %in% cytokine_eligible &
  feature_id %in% matrix_features,
  .(analysis_group = "cytokines", feature_set = signature_id, feature_id)]
eligible_memberships$cytokines <- cytokine_mapped

combined_inventory <- rbindlist(all_inventory, fill = TRUE, use.names = TRUE)
combined_memberships <- rbindlist(eligible_memberships, use.names = TRUE)
setorder(combined_inventory, analysis_group, feature_set)
setorder(combined_memberships, analysis_group, feature_set, feature_id)
write_tsv_gz(combined_inventory,
  file.path(stage_dir, "all_group_set_inventory.tsv.gz"))
write_tsv_gz(combined_memberships,
  file.path(stage_dir, "all_group_eligible_mapped_memberships.tsv.gz"))

group_summary <- combined_inventory[, .(
  source_sets = as.integer(.N),
  eligible_sets = as.integer(sum(eligible)),
  excluded_below_15 = as.integer(sum(
    n_genes_mapped_to_bundle < minimum_mapped_genes)),
  excluded_above_300 = as.integer(sum(
    n_genes_mapped_to_bundle > maximum_mapped_genes)),
  median_source_genes = as.numeric(median(n_genes_raw)),
  median_mapped_genes = as.numeric(median(n_genes_mapped_to_bundle)),
  median_mapped_fraction = as.numeric(median(mapped_fraction))),
  by = .(analysis_group, collection_label)]
group_summary[, group_order := match(analysis_group,
  c("bp", "mf", "cc", "cytokines"))]
setorder(group_summary, group_order)
group_summary[, group_order := NULL]
write_tsv(group_summary, file.path(stage_dir, "group_summary.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
completion <- list(status = "complete", purpose = "GO trajectory preflight",
  source = list(database = "Mouse MSigDB", version = unique(unlist(source_versions)),
    collection = "M5", msigdbr_version = as.character(packageVersion("msigdbr"))),
  eligibility = list(minimum_mapped_genes = minimum_mapped_genes,
    maximum_mapped_genes = maximum_mapped_genes,
    duplicate_mapped_memberships = "retained as distinct GO hypotheses"),
  planned_analysis = list(positive_only = TRUE,
    depth_adjustment = "score ~ ns(pt,df=8) + scaled log1p(nCount_RNA); subtract depth component only",
    window_width_fraction = 0.05, window_step_fraction = 0.01,
    relative_prominence_minimum = 0.10,
    merge_gap_fraction_strictly_less_than = 0.05,
    multiple_testing = "BH separately within BP, MF, CC, and cytokines",
    permutations = "adaptive; expensive stage not started"),
  implementation_boundary = list(
    scworkbench_documentation_consulted = TRUE,
    native_three_term_import_smoke_test = "complete",
    full_native_import_benchmark = "stopped after 1,145 of 4,087 eligible BP terms in approximately 13 minutes",
    full_inventory = "user-approved outside vectorized workflow using exact scWorkbench matrix feature IDs"),
  groups = split(as.data.frame(group_summary), group_summary$analysis_group),
  defined_pt_cells = sum(is.finite(metadata$pt)),
  missing_pt_cells = sum(!is.finite(metadata$pt)), promoted = FALSE,
  bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = fingerprint_after, bundle_unchanged = TRUE)
write_json(completion, file.path(stage_dir, "completion.json"))

table_lines <- c(
  "| Group | Source sets | Eligible | <15 mapped | >300 mapped | Median mapped | Median coverage |",
  "|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(group_summary)), function(i) paste0("| ",
    group_summary$collection_label[[i]], " | ", group_summary$source_sets[[i]],
    " | ", group_summary$eligible_sets[[i]], " | ",
    group_summary$excluded_below_15[[i]], " | ",
    group_summary$excluded_above_300[[i]], " | ",
    signif(group_summary$median_mapped_genes[[i]], 5), " | ",
    signif(group_summary$median_mapped_fraction[[i]], 4), " |"), character(1)))
report <- c("# GO trajectory enrichment preflight", "",
  "The expensive permutation stage has not been started.", "",
  table_lines, "", "## Fixed analysis decisions", "",
  "- Four separately reported groups: GO BP, GO MF, GO CC, and 38 reviewed cytokine signatures.",
  "- Eligibility: 15–300 genes mapped to `RNA_log1p_cp10k`.",
  "- Multiple testing: BH separately within each of the four groups.",
  "- Positive enrichment only; conditional RNA-depth adjustment using preserved Bdev `pt`.",
  "- Individual peak discovery precedes q-value filtering, 10% relative-prominence filtering, and the strict <5% gap merge.",
  "- Exact gene-symbol overlap was computed against feature IDs read from the scWorkbench `RNA_log1p_cp10k` matrix.",
  "- The documented native import/QC route passed a three-term smoke test, but its full BP benchmark was stopped after 1,145/4,087 eligible terms in about 13 minutes; the approved vectorized outside workflow is used at full scale.",
  "", "The bundle and exact `Bdev@meta.data$pt` values were not modified.")
writeLines(report, file.path(stage_dir, "REPORT.md"))

wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed GO trajectory preflight at", output_dir)
