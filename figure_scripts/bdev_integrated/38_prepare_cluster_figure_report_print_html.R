#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(paste(
  "Usage: 38_prepare_cluster_figure_report_print_html.R",
  "<source_html> <temporary_print_html>"))
source_html <- normalizePath(args[[1L]], mustWork = TRUE)
output_html <- normalizePath(args[[2L]], mustWork = FALSE)

html <- paste(readLines(source_html, warn = FALSE), collapse = "\n")
# Images were linked to vector PDFs in the interactive report. Retain the
# images but remove every anchor and remove the supporting-links section.
html <- gsub("<a href='[^']+'>(<img[^>]+>)</a>", "\\1", html, perl = TRUE)
html <- gsub("<section><h2>Supporting tables</h2>.*?</section>", "", html,
  perl = TRUE)
html <- gsub("<a href='[^']+'>(.*?)</a>", "\\1", html, perl = TRUE)
html <- gsub(" Click for vector PDF.", "", html, fixed = TRUE)
stopifnot(!grepl("<a ", html, fixed = TRUE),
  !grepl("Supporting tables", html, fixed = TRUE),
  !grepl("Click for vector PDF", html, fixed = TRUE))

base_uri <- paste0("file://", dirname(source_html), "/")
html <- sub("<head>", paste0("<head><base href='", base_uri, "'>"), html,
  fixed = TRUE)
print_css <- paste0(
  "<style>",
  "@page{size:A3 landscape;margin:12mm}",
  "body{max-width:none;padding:0;font-size:12pt;color:#111827}",
  "h1{font-size:24pt}h2{font-size:19pt;break-before:page}",
  "h3{font-size:15pt;break-after:avoid-page}",
  ".lede,.methods,.takeaway{break-inside:avoid-page}",
  "img{display:block;width:auto;max-width:100%;height:auto;max-height:205mm;",
  "object-fit:contain;margin:0 auto;break-inside:avoid-page;",
  "page-break-inside:avoid}",
  ".caption{break-after:page;margin-bottom:0}",
  "section ul{break-inside:auto}footer{break-before:page}",
  "</style>")
html <- sub("</head>", paste0(print_css, "</head>"), html, fixed = TRUE)
writeLines(html, output_html)
