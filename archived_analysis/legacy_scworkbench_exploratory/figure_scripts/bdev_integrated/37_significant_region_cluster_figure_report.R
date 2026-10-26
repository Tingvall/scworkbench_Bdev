#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(paste(
  "Usage: 37_significant_region_cluster_figure_report.R",
  "<kmeans_report_dir> <output_dir>"))
source_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

completion <- read_json(file.path(source_dir, "completion.json"),
  simplifyVector = TRUE)
stopifnot(identical(completion$status, "complete"),
  identical(completion$promoted, FALSE),
  completion$groups$bp$selected_k == 25L,
  completion$groups$mf$selected_k == 20L,
  completion$groups$cc$selected_k == 20L)
summary_table <- fread(file.path(source_dir, "cluster_summary.tsv"))
membership <- fread(file.path(source_dir, "gene_set_cluster_membership.tsv"))
stopifnot(nrow(membership) == 1689L, !anyDuplicated(
  membership[, paste(analysis_group, feature_set)]),
  identical(sort(unique(summary_table$analysis_group)), c("bp", "cc", "mf")))

groups <- c("bp", "mf", "cc")
labels <- c(bp = "Biological Process", mf = "Molecular Function",
  cc = "Cellular Component")
k_values <- c(bp = 25L, mf = 20L, cc = 20L)
for (g in groups) {
  required <- file.path(source_dir, c(
    paste0(g, "_significant_region_kmeans_heatmap.png"),
    paste0(g, "_significant_region_cluster_centroids.png"),
    paste0(g, "_significant_region_kmeans_heatmap.pdf"),
    paste0(g, "_significant_region_cluster_centroids.pdf")))
  stopifnot(all(file.exists(required)))
}

overview <- c(
  "Across all three ontologies, the dominant trajectory is an early metabolic and receptor-associated state, followed by immune and biosynthetic remodeling near BB1, and then a late DNA-replication and mitotic program around and after BB2.",
  "The strongest cross-ontology agreement is among nucleotide metabolism, spliceosome organization, chromatin and DNA replication, and spindle and cytokinesis programs.")

bp_points <- c(
  "Very early clusters C01-C05, before BB1, emphasize lipid and hormone metabolism, biosynthesis, receptor signaling, and developmental differentiation.",
  "Near BB1, C06-C10 shift toward antigen processing and presentation, leukocyte migration, lymphocyte differentiation, immune signaling, and stress or redox responses.",
  "Between BB1 and BB2, C11-C19 show nucleotide biosynthesis, mitochondrial activity, amino-acid transport, protein folding, and developmental remodeling. C13 is especially coherent for pyrimidine metabolism, telomere maintenance, and mitochondrial protein import.",
  "C11, C13, and C15 have recurrent supported regions, indicating repeated or separated activation rather than one continuous peak.",
  "Approaching and following BB2, C20-C24 form a proliferative sequence: B-cell proliferation and isotype switching, DNA replication and repair, chromosome segregation, spindle organization, mitosis, and cytokinesis.",
  "C21 and C23 retain two supported regions, suggesting that components of the replication and mitotic programs recur at two trajectory positions. C25 is a small terminal differentiation and signaling group.")

mf_points <- c(
  "Early C01-C03 emphasize transporters, growth-factor and receptor binding, hormone binding, and receptor-kinase activity.",
  "C04-C06, approaching BB1, emphasize ligand-gated channels, cytokine and immune receptors, MHC and antigen binding, chemokine binding, and SMAD-associated functions.",
  "Around BB1, C07-C11 separate receptor and channel modes; C08 and C10 have two-region profiles, while C11 contains rarer three-region oxidoreductase, channel, and receptor functions.",
  "Between BB1 and BB2, C12-C16 progress through oxidoreductase and transporter functions, protein-folding chaperones and translation factors, RNA and polymerase binding, cyclic-nucleotide signaling, and transcriptional coregulation.",
  "After BB2, C17-C19 shift to chromatin and DNA binding, histone readers and helicases, kinase activation, ubiquitin machinery, cyclin-dependent kinase regulation, and microtubule motors.",
  "C20 is a small late extracellular-matrix and receptor-associated group.")

cc_points <- c(
  "Early C01-C03 emphasize ribosomal subunits, synapses, dense-core and secretory vesicles, and presynaptic membranes.",
  "Near BB1, C04-C05 highlight MHC-I loading, Golgi compartments, and ER or lumenal membranes, matching the Biological Process antigen-presentation wave.",
  "Between BB1 and BB2, C06-C07 contain extracellular-matrix, adhesion, connexin, and channel complexes, while C08-C11 show ER, mitochondrial, spliceosome, and snRNP compartments.",
  "The separate one-region and recurrent spliceosome clusters indicate distinct temporal modes of RNA-processing-complex activity rather than one uniform splicing signal.",
  "Later C15-C19 move into replication forks, DNA-damage sites, heterochromatin, DNA polymerase complexes, spindle poles, centrioles, cleavage furrows, kinetochore microtubules, APC/C, and cyclin-dependent kinase complexes.",
  "C20 is a small terminal synaptic, matrix, and contractile-compartment group.")

group_points <- list(bp = bp_points, mf = mf_points, cc = cc_points)
group_takeaway <- c(
  bp = "Early metabolism and development -> immune and antigen programs -> nucleotide and mitochondrial remodeling -> DNA replication -> mitosis and cytokinesis.",
  mf = "Receptor and transport functions -> immune and ion-channel functions -> chaperone and metabolic functions -> chromatin and DNA machinery -> kinase, ubiquitin, and cell-cycle regulation.",
  cc = "Ribosome and secretory compartments -> Golgi, MHC, and ER -> spliceosome and mitochondrial complexes -> replication and chromatin structures -> spindle and cytokinesis machinery.")

md <- c(
  "# Bdev significant-region residual-UCell trajectory clusters", "",
  "## Overview", "", paste0(overview, collapse = "\n\n"), "",
  "The dashed lines in every figure mark **BB1 (pt -19.16)** and **BB2 (pt -13.76)**.", "",
  "## Reading the figures", "",
  "- Heatmaps contain only GO sets with at least one permutation-supported positive region.",
  "- Each set was min-max scaled to 0-1 across its trajectory; values outside supported regions were set to zero.",
  "- K-means used these significant-region-only profiles with BP k=25, MF k=20, and CC k=20.",
  "- Centroid curves are cluster means of the same masked 0-1 profiles. They describe timing and shape, not absolute enrichment magnitude.",
  "- Cluster annotations summarize member GO names and representative terms; they are descriptive and are not an additional enrichment test.", "")

for (g in groups) {
  md <- c(md, paste0("## GO ", labels[[g]]), "",
    paste0("### Heatmap (k=", k_values[[g]], ")"), "",
    paste0("![GO ", labels[[g]], " significant-region heatmap](../",
      g, "_significant_region_kmeans_heatmap.png)"), "",
    paste0("[Vector PDF](../", g,
      "_significant_region_kmeans_heatmap.pdf)"), "",
    "### Cluster centroid curves", "",
    paste0("![GO ", labels[[g]], " centroid curves](../",
      g, "_significant_region_cluster_centroids.png)"), "",
    paste0("[Vector PDF](../", g,
      "_significant_region_cluster_centroids.pdf)"), "",
    "### Interpretation", "",
    paste0("- ", group_points[[g]]), "",
    paste0("**Trajectory summary:** ", group_takeaway[[g]]), "")
}
md <- c(md, "## Supporting tables", "",
  "- [Cluster summaries](../cluster_summary.tsv)",
  "- [Gene-set memberships](../gene_set_cluster_membership.tsv)",
  "- [Cluster centroid profiles](../cluster_centroid_profiles.tsv.gz)",
  "- [K diagnostic scan](../silhouette_k_scan.tsv)",
  "- [Analysis completion record](../completion.json)", "",
  "## Interpretation boundary", "",
  paste0("These results are descriptive trajectory-pattern clusters. Per-set scaling removes absolute score magnitude, ",
    "and GO redundancy means related terms can appear repeatedly. Biological conclusions should therefore emphasize ",
    "program-level timing and cross-ontology agreement rather than treating every GO term or cluster label as independent evidence."))
writeLines(md, file.path(stage_dir, "FIGURE_REPORT.md"))

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
li <- function(x) paste0("<li>", html_escape(x), "</li>")
section_html <- function(g) paste0(
  "<section><h2>GO ", labels[[g]], "</h2>",
  "<h3>Heatmap (k=", k_values[[g]], ")</h3>",
  "<a href='../", g, "_significant_region_kmeans_heatmap.pdf'>",
  "<img src='../", g, "_significant_region_kmeans_heatmap.png' alt='GO ",
  labels[[g]], " heatmap'></a>",
  "<p class='caption'>Significant-region-only residual-UCell profiles. Click for vector PDF.</p>",
  "<h3>Cluster centroid curves</h3>",
  "<a href='../", g, "_significant_region_cluster_centroids.pdf'>",
  "<img src='../", g, "_significant_region_cluster_centroids.png' alt='GO ",
  labels[[g]], " centroid curves'></a>",
  "<p class='caption'>Cluster means of the same masked 0-1 profiles. Click for vector PDF.</p>",
  "<h3>Interpretation</h3><ul>", paste(vapply(group_points[[g]], li,
    character(1)), collapse = ""), "</ul>",
  "<p class='takeaway'><strong>Trajectory summary:</strong> ",
  html_escape(group_takeaway[[g]]), "</p></section>")

html <- paste0("<!doctype html><html><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>Bdev significant-region residual-UCell trajectory clusters</title>",
  "<style>body{font-family:Arial,Helvetica,sans-serif;line-height:1.5;color:#172033;",
  "max-width:1500px;margin:0 auto;padding:32px}h1{font-size:2rem}h2{margin-top:3rem;",
  "border-bottom:2px solid #d6dbe5;padding-bottom:.4rem}h3{margin-top:1.6rem}",
  "img{display:block;max-width:100%;height:auto;border:1px solid #d6dbe5;background:white}",
  ".lede{font-size:1.08rem;background:#f3f6fa;padding:18px 22px;border-left:5px solid #334e68}",
  ".caption{color:#536171;font-size:.92rem}.takeaway{background:#f7f1e8;padding:12px 16px}",
  ".methods{background:#f8fafc;padding:16px 20px}li{margin:.42rem 0}",
  "a{color:#185fa5}footer{margin-top:3rem;border-top:1px solid #d6dbe5;padding-top:1rem;",
  "color:#536171}</style></head><body>",
  "<h1>Bdev significant-region residual-UCell trajectory clusters</h1>",
  "<div class='lede'><p>", html_escape(overview[[1L]]), "</p><p>",
  html_escape(overview[[2L]]), "</p><p>The dashed lines mark <strong>BB1 (pt -19.16)</strong> ",
  "and <strong>BB2 (pt -13.76)</strong>.</p></div>",
  "<div class='methods'><h2>Reading the figures</h2><ul>",
  "<li>Only GO sets with at least one permutation-supported positive region are shown.</li>",
  "<li>Each trajectory was min-max scaled to 0-1; unsupported windows were set to zero.</li>",
  "<li>K-means used the significant-region-only profiles: BP k=25, MF k=20, CC k=20.</li>",
  "<li>Centroid curves summarize timing and shape, not absolute enrichment magnitude.</li>",
  "<li>Biological annotations are descriptive summaries, not an additional enrichment test.</li>",
  "</ul></div>", paste(vapply(groups, section_html, character(1)), collapse = ""),
  "<section><h2>Supporting tables</h2><ul>",
  "<li><a href='../cluster_summary.tsv'>Cluster summaries</a></li>",
  "<li><a href='../gene_set_cluster_membership.tsv'>Gene-set memberships</a></li>",
  "<li><a href='../cluster_centroid_profiles.tsv.gz'>Cluster centroid profiles</a></li>",
  "<li><a href='../silhouette_k_scan.tsv'>K diagnostic scan</a></li>",
  "<li><a href='../completion.json'>Analysis completion record</a></li></ul></section>",
  "<footer><strong>Interpretation boundary.</strong> These are descriptive trajectory-pattern clusters. ",
  "Per-set scaling removes absolute magnitude, and GO redundancy means related terms are not independent. ",
  "Emphasis should be placed on program-level timing and cross-ontology agreement.</footer>",
  "</body></html>")
writeLines(html, file.path(stage_dir, "FIGURE_REPORT.html"))
write_json(list(status = "complete", promoted = FALSE,
  source_cluster_report = source_dir,
  source_completion = file.path(source_dir, "completion.json"),
  cluster_counts = as.list(k_values),
  figures = setNames(lapply(groups, function(g) list(
    heatmap_png = file.path(source_dir,
      paste0(g, "_significant_region_kmeans_heatmap.png")),
    centroids_png = file.path(source_dir,
      paste0(g, "_significant_region_cluster_centroids.png")))), groups),
  html = "FIGURE_REPORT.html", markdown = "FIGURE_REPORT.md"),
  file.path(stage_dir, "completion.json"), pretty = TRUE, auto_unbox = TRUE)
stopifnot(file.rename(stage_dir, output_dir))
message("Completed cluster figure report at: ", output_dir)
