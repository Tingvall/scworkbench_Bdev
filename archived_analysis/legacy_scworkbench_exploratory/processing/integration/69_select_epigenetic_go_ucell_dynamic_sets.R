#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(paste(
  "Usage: 69_select_epigenetic_go_ucell_dynamic_sets.R",
  "<go_dynamic_report> <output_dir>"))
report_root <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

completion <- read_json(file.path(report_root, "completion.json"),
  simplifyVector = TRUE)
stopifnot(identical(completion$status, "complete"),
  identical(completion$inference_mode, "dynamic_range"),
  completion$permutations == 1000L)

# Explicitly epigenetic/chromatin-related terms selected after excluding broad
# chromosome segregation, condensation, centromere-localization, and telomere
# maintenance terms that need not represent epigenetic regulation.
selection <- data.table(
  analysis_group = rep(c("bp", "mf", "cc"), each = 3L),
  feature_set = c(
    "GOBP_RANDOM_INACTIVATION_OF_X_CHROMOSOME",
    "GOBP_FACULTATIVE_HETEROCHROMATIN_FORMATION",
    "GOBP_NUCLEOSOME_ASSEMBLY",
    "GOMF_HISTONE_CHAPERONE_ACTIVITY",
    "GOMF_HISTONE_H3K9ME2_3_READER_ACTIVITY",
    "GOMF_NUCLEOSOME_BINDING",
    "GOCC_CHROMATIN_SILENCING_COMPLEX",
    "GOCC_CHROMOCENTER",
    "GOCC_PERICENTRIC_HETEROCHROMATIN"),
  selection_basis = c(
    "chromosome-wide epigenetic silencing",
    "repressive facultative heterochromatin formation",
    "nucleosome assembly and chromatin organization",
    "histone handling and nucleosome assembly",
    "reader of repressive H3K9me2/3 histone marks",
    "direct nucleosome interaction",
    "chromatin-silencing machinery",
    "heterochromatin-associated nuclear domain",
    "repressive pericentric heterochromatin"))

selected <- rbindlist(lapply(c("bp", "mf", "cc"), function(g) {
  x <- fread(file.path(report_root, g, "significant_gene_sets.tsv"))
  out <- merge(selection[analysis_group == g], x, by = "feature_set",
    all.x = TRUE, sort = FALSE)
  stopifnot(nrow(out) == 3L, !anyNA(out$q_dynamic),
    all(out$q_dynamic < 0.01))
  setorder(out, q_dynamic, p_dynamic, -observed_dynamic_range, feature_set)
  out[, rank_within_group := seq_len(.N)]
  out
}), use.names = TRUE, fill = TRUE)
setcolorder(selected, c("analysis_group", "rank_within_group",
  "feature_set", "set_name", "selection_basis", "q_dynamic", "p_dynamic",
  "observed_dynamic_range", "n_regions", "best_region_fwer_p",
  "strongest_ucell_effect", "peak_pseudotimes"))
fwrite(selected, file.path(stage_dir,
  "top3_epigenetic_dynamic_gene_sets_by_group.tsv"), sep = "\t",
  quote = FALSE, na = "NA")

fmt <- function(x) vapply(seq_len(nrow(x)), function(i) paste0("| ",
  x$rank_within_group[[i]], " | `", x$set_name[[i]], "` | ",
  x$selection_basis[[i]], " | ", signif(x$q_dynamic[[i]], 4), " | ",
  signif(x$p_dynamic[[i]], 4), " | ",
  signif(x$observed_dynamic_range[[i]], 5), " | ", x$n_regions[[i]],
  " | ", ifelse(is.na(x$peak_pseudotimes[[i]]), "—",
    x$peak_pseudotimes[[i]]), " |"), character(1))
report <- c("# Top epigenetic dynamic GO gene sets", "",
  "Up to three explicitly epigenetic/chromatin-related terms were selected independently within BP, MF, and CC.",
  "Broad chromosome segregation, chromosome condensation, centromere localization, and telomere-maintenance terms were excluded unless the GO term itself explicitly represented chromatin or epigenetic regulation.",
  "Within the curated relevant pool, ranking follows dynamic q, permutation p, and descending native-UCell peak-to-trough range.", "")
for (g in c("bp", "mf", "cc")) {
  x <- selected[analysis_group == g]
  report <- c(report, paste0("## ", toupper(g)), "",
    "| Rank | Gene set | Basis | q | p | UCell range | Regions | Region peak PT(s) |",
    "|---:|---|---|---:|---:|---:|---:|---|", fmt(x), "")
}
writeLines(report, file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", selected_sets = nrow(selected),
  selected_per_group = 3L,
  ranking = "dynamic q, permutation p, descending native-UCell range",
  broad_nonepigenetic_chromosome_terms_excluded = TRUE,
  source_report = report_root), file.path(stage_dir, "completion.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)

stopifnot(file.rename(stage_dir, output_dir))
message("Completed epigenetic GO selection at: ", output_dir)
