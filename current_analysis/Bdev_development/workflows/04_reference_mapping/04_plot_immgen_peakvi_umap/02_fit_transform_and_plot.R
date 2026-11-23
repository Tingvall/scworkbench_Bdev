#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 02_fit_transform_and_plot.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbRepresentation)
  library(yaml)
  library(jsonlite)
  library(digest)
})

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
read_config <- function(path) {
  current <- yaml::read_yaml(path)
  base_name <- current$base_config
  current$base_config <- NULL
  if (is.null(base_name)) return(current)
  utils::modifyList(read_config(normalizePath(file.path(dirname(path), base_name))),
    current, keep.null = TRUE)
}
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
promoted_ref <- function(promoted, output, kind) {
  selected <- promoted$promoted$status == "ok" & promoted$promoted$output == output
  stopifnot(sum(selected) == 1L)
  wb_ref_artifact(promoted$promoted$artifact_id[selected][[1L]], kind)
}
with_devices <- function(stem, width, height, draw) {
  grDevices::pdf(paste0(stem, ".pdf"), width = width, height = height,
    useDingbats = FALSE)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(stem, ".png"), width = width, height = height,
    units = "in", res = 220, type = "cairo")
  draw()
  grDevices::dev.off()
}

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- read_config(cfg_path)
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "02_fit_transform_and_plot.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(sha256(file.path(cfg$inputs$mapping_bundle$path, "manifest.json")),
    cfg$inputs$mapping_bundle$manifest_sha256),
  identical(sha256(cfg$inputs$mapping_table$path), cfg$inputs$mapping_table$sha256),
  identical(sha256(cfg$inputs$authenticated_l4$path),
    cfg$inputs$authenticated_l4$sha256),
  identical(cfg$inputs$authenticated_l4$orientation, "native"),
  as.character(packageVersion("scworkbench")) == cfg$runtime$scworkbench$version,
  as.character(packageVersion("scwbRepresentation")) ==
    cfg$runtime$scwbRepresentation$version)
pack_manifest <- system.file("scworkbench", "operation-pack.json",
  package = "scwbRepresentation")
stopifnot(identical(sha256(pack_manifest),
  cfg$runtime$scwbRepresentation$manifest_sha256))
wb_operation_pack_activate(pack_manifest, lib.loc = .libPaths())
stopifnot(all(c("scwbrepresentation::umap_fit",
  "scwbrepresentation::umap_transform", "cell_set_derive") %in% wb_list_ops()))

source_manifest_before <- sha256(file.path(cfg$inputs$mapping_bundle$path,
  "manifest.json"))
derived_manifest_before <- sha256(file.path(cfg$storage$derived_bundle,
  "manifest.json"))

bundle <- wb_open_bundle(cfg$storage$derived_bundle, validate = FALSE,
  read_only = FALSE)
spaces <- wb_spaces(bundle)
umap_already_promoted <- cfg$umap$output_space %in% spaces$space_id
if (umap_already_promoted) {
  cat("Reusing already promoted UMAP space", cfg$umap$output_space, "\n")
} else stopifnot(identical(source_manifest_before, derived_manifest_before))
required_spaces <- c(cfg$inputs$mapping_bundle$reference_space,
  cfg$inputs$mapping_bundle$query_space,
  cfg$inputs$mapping_bundle$authenticated_union_space)
stopifnot(all(required_spaces %in% spaces$space_id))
reference_latent <- wb_space(bundle, cfg$inputs$mapping_bundle$reference_space)
query_latent <- wb_space(bundle, cfg$inputs$mapping_bundle$query_space)
stopifnot(nrow(reference_latent) == cfg$inputs$mapping_bundle$expected_reference_cells,
  nrow(query_latent) == cfg$inputs$mapping_bundle$expected_query_libraries,
  ncol(reference_latent) == cfg$inputs$mapping_bundle$expected_dimensions,
  ncol(query_latent) == cfg$inputs$mapping_bundle$expected_dimensions,
  !length(intersect(rownames(reference_latent), rownames(query_latent))))
design <- as.character(cfg$umap$design)[[1L]]
if (identical(design, "authenticated_union_fit")) {
  fit_space <- cfg$inputs$mapping_bundle$authenticated_union_space
  fit_latent <- wb_space(bundle, fit_space)
  stopifnot(nrow(fit_latent) == nrow(reference_latent) + nrow(query_latent),
    ncol(fit_latent) == cfg$inputs$mapping_bundle$expected_dimensions,
    setequal(rownames(fit_latent), c(rownames(reference_latent), rownames(query_latent))))
  fit_note <- "UMAP fit on authenticated model-bound Bdev/ImmGen 10D PeakVI union"
  fit_prefix <- "step4_immgen_peakvi_joint_umap"
  fit_description <- "Joint UMAP model for the authenticated Bdev and ImmGen PeakVI union"
  if (!umap_already_promoted)
    cat("Fitting UMAP on authenticated union of", nrow(fit_latent), "objects\n")
} else if (identical(design,
    "reference_only_umap_fit_then_frozen_model_query_transform")) {
  cat("Promoting explicit query cell set for the UMAP transform contract\n")
  cell_set_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_query_cell_set"))
  cell_set_run <- wb_session_run_op(cell_set_session, bundle, "cell_set_derive",
    inputs = list(cell_ids = rownames(query_latent),
      target_cell_set = cfg$umap$query_cell_set, role = "analysis_universe",
      status = "accepted", source = "authenticated ImmGen PeakVI query space",
      created_by = "scworkbench_Bdev Step 4 derived ImmGen UMAP",
      description = "The 28 ImmGen ATAC libraries eligible for frozen UMAP transform",
      min_cells = 1L),
    notes = "Provenance-bound query cell set for frozen UMAP transform")
  cell_set_run$session <- wb_session_resume(bundle,
    cell_set_run$session$durable$workspace)
  cell_set_promoted <- wb_session_promote_current(bundle, cell_set_run$session,
    outputs = c("manifest", "rule_table"),
    name_prefix = "step4_immgen_umap_query_cell_set",
    description = "Exact 28-library ImmGen UMAP query cell set",
    created_by = "scworkbench_Bdev Step 4 derived ImmGen UMAP", fail_fast = TRUE)
  cell_set_promoted$session <- wb_session_discard(cell_set_promoted$session)
  stopifnot(identical(wb_cells(bundle, cell_set = cfg$umap$query_cell_set),
    rownames(query_latent)))
  fit_space <- cfg$inputs$mapping_bundle$reference_space
  fit_latent <- reference_latent
  fit_note <- "Reference-only UMAP fit on qualified 10D Bdev PeakVI latent"
  fit_prefix <- "step4_immgen_peakvi_reference_umap"
  fit_description <- "Reference-only UMAP model for the qualified Bdev PeakVI latent"
  if (!umap_already_promoted)
    cat("Fitting reference-only UMAP on", nrow(fit_latent), "Bdev cells\n")
} else stop("Unsupported UMAP design: ", design)

if (!umap_already_promoted) {
  fit_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_umap_fit"))
  fit <- wb_session_run_op(fit_session, bundle, "scwbrepresentation::umap_fit",
    inputs = list(space = fit_space,
      dims = as.integer(unlist(cfg$umap$dimensions)),
      space_id_out = cfg$umap$output_space,
      n_neighbors = as.integer(cfg$umap$n_neighbors),
      min_dist = as.numeric(cfg$umap$min_dist), metric = cfg$umap$metric,
      nn_method = cfg$umap$nn_method, spread = as.numeric(cfg$umap$spread),
      seed = as.integer(cfg$umap$seed)), notes = fit_note)
  stopifnot(nrow(fit$run$outputs$space_2d$values) == nrow(fit_latent),
    ncol(fit$run$outputs$space_2d$values) == 2L)
  fit$session <- wb_session_resume(bundle, fit$session$durable$workspace)
  fit_promoted <- wb_session_promote_current(bundle, fit$session, outputs = "model",
    name_prefix = fit_prefix, description = fit_description,
    created_by = "scworkbench_Bdev Step 4 derived ImmGen UMAP", fail_fast = TRUE)
  model_ref <- promoted_ref(fit_promoted, "model", "model")
  fit_promoted$session <- wb_session_discard(fit_promoted$session)
}

if (!umap_already_promoted && identical(design,
    "reference_only_umap_fit_then_frozen_model_query_transform")) {
  cat("Transforming", nrow(query_latent), "ImmGen libraries through frozen UMAP model\n")
  transform_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_query_umap_transform"))
  transform <- wb_session_run_op(transform_session, bundle,
    "scwbrepresentation::umap_transform", inputs = list(model = model_ref,
      in_space = cfg$inputs$mapping_bundle$query_space,
      dims = as.integer(unlist(cfg$umap$dimensions)),
      space_id_out = cfg$umap$output_space, overwrite = TRUE,
      cell_set = cfg$umap$query_cell_set, min_cells = 1L),
    notes = "Frozen-model UMAP transform of 28 ImmGen ATAC libraries")
  stopifnot(nrow(transform$run$outputs$space_2d$values) == nrow(query_latent),
    ncol(transform$run$outputs$space_2d$values) == 2L)
  transform$session <- wb_session_resume(bundle,
    transform$session$durable$workspace)
  transform_promoted <- wb_session_promote_current(bundle, transform$session,
    outputs = "projection_membership",
    name_prefix = "step4_immgen_peakvi_umap_projection",
    description = "Frozen UMAP projection of 28 ImmGen ATAC libraries",
    created_by = "scworkbench_Bdev Step 4 derived ImmGen UMAP", fail_fast = TRUE)
  transform_promoted$session <- wb_session_discard(transform_promoted$session)
}
wb_close_bundle(bundle)

bundle <- wb_open_bundle(cfg$storage$derived_bundle, validate = FALSE,
  read_only = TRUE)
coordinates_raw <- wb_space(bundle, cfg$umap$output_space)
wb_close_bundle(bundle)
if (is.data.frame(coordinates_raw) &&
    all(c("cell_id", "x", "y") %in% colnames(coordinates_raw))) {
  coordinates <- as.matrix(coordinates_raw[, c("x", "y"), drop = FALSE])
  rownames(coordinates) <- as.character(coordinates_raw$cell_id)
} else coordinates <- as.matrix(coordinates_raw)
expected_cells <- c(rownames(reference_latent), rownames(query_latent))
stopifnot(nrow(coordinates) == length(expected_cells),
  setequal(rownames(coordinates), expected_cells), ncol(coordinates) == 2L,
  all(is.finite(coordinates)))
coordinates <- coordinates[expected_cells, , drop = FALSE]

mapping <- read.delim(cfg$inputs$mapping_table$path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
stopifnot(nrow(mapping) == nrow(query_latent),
  setequal(mapping$gsm, rownames(query_latent)))
mapping <- mapping[match(rownames(query_latent), mapping$gsm), , drop = FALSE]
l4_family <- readRDS(cfg$inputs$authenticated_l4$path)
l4_selected <- as.character(l4_family$key) == cfg$inputs$authenticated_l4$key &
  is.finite(as.numeric(l4_family$value))
reference_l4 <- stats::setNames(as.numeric(l4_family$value[l4_selected]),
  as.character(l4_family$cell_id[l4_selected]))
stopifnot(!anyDuplicated(names(reference_l4)),
  all(rownames(reference_latent) %in% names(reference_l4)),
  all(is.finite(mapping$accepted_native_l4)))

coordinate_table <- data.frame(
  cell_id = rownames(coordinates),
  source = c(rep("Bdev_reference", nrow(reference_latent)),
    rep("ImmGen_query_library", nrow(query_latent))),
  UMAP_1 = as.numeric(coordinates[, 1L]),
  UMAP_2 = as.numeric(coordinates[, 2L]),
  population = c(rep(NA_character_, nrow(reference_latent)), mapping$population),
  support_call = c(rep(NA_character_, nrow(reference_latent)), mapping$support_call),
  native_l4 = c(unname(reference_l4[rownames(reference_latent)]),
    mapping$accepted_native_l4),
  stringsAsFactors = FALSE)
coordinate_path <- file.path(dirs[["tables"]], "01_bdev_immgen_peakvi_umap_coordinates.tsv")
write_tsv(coordinate_table, coordinate_path)

reference_plot <- coordinate_table[coordinate_table$source == "Bdev_reference", ]
query_plot <- coordinate_table[coordinate_table$source == "ImmGen_query_library", ]
populations <- sort(unique(query_plot$population), method = "radix")
population_colors <- stats::setNames(grDevices::palette.colors(length(populations),
  "Alphabet"), populations)
query_plot$fill <- unname(population_colors[query_plot$population])
query_plot$pch <- ifelse(query_plot$support_call == "borderline", 24L, 21L)
query_plot$border <- ifelse(query_plot$support_call == "borderline",
  "#D55E00", "#202020")
borderline_index <- which(query_plot$support_call == "borderline")
borderline_label_position <- c(3L, 2L, 4L)

draw_population <- function() {
  old <- par(mar = c(4.8, 5.0, 4.2, 1.2), xaxs = "i", yaxs = "i")
  on.exit(par(old))
  xpad <- diff(range(coordinate_table$UMAP_1)) * 0.05
  ypad <- diff(range(coordinate_table$UMAP_2)) * 0.07
  plot(reference_plot$UMAP_1, reference_plot$UMAP_2, pch = 16,
    col = grDevices::adjustcolor("#A8A8A8", alpha.f = 0.48), cex = 0.62,
    xlab = "PeakVI UMAP 1", ylab = "PeakVI UMAP 2",
    xlim = range(coordinate_table$UMAP_1) + c(-xpad, xpad),
    ylim = range(coordinate_table$UMAP_2) + c(-ypad, ypad),
    main = "ImmGen ATAC libraries in the joint PeakVI UMAP")
  mtext("Authenticated 10D union: 1,382 Bdev cells + 28 ImmGen libraries; borderline libraries labeled",
    side = 3, line = 0.5, cex = 0.82)
  points(query_plot$UMAP_1, query_plot$UMAP_2, pch = query_plot$pch,
    bg = query_plot$fill, col = query_plot$border, cex = 1.8, lwd = 1.7)
  text(query_plot$UMAP_1[borderline_index], query_plot$UMAP_2[borderline_index],
    labels = query_plot$cell_id[borderline_index],
    pos = borderline_label_position, offset = 0.55, cex = 0.64, font = 2)
  legend("topleft", legend = populations, pt.bg = population_colors,
    col = "#202020", pch = 21, pt.cex = 1.25, cex = 0.64, ncol = 2,
    bty = "n", title = "ImmGen population")
  legend("topright", legend = c("supported", "borderline"),
    pch = c(21, 24), pt.bg = "white", col = c("#202020", "#D55E00"),
    pt.lwd = 1.7, bty = "n", cex = 0.75, title = "Reference support")
}
with_devices(file.path(dirs[["figures"]],
  "01_immgen_population_peakvi_umap"), 11.5, 8.5, draw_population)

l4_range <- range(coordinate_table$native_l4, finite = TRUE)
l4_palette <- grDevices::hcl.colors(256L, "Inferno")
l4_index <- function(x) {
  scaled <- (x - l4_range[[1L]]) / diff(l4_range)
  pmax(1L, pmin(256L, 1L + floor(scaled * 255L)))
}
reference_plot$l4_color <- l4_palette[l4_index(reference_plot$native_l4)]
query_plot$l4_color <- l4_palette[l4_index(query_plot$native_l4)]
l4_ticks <- pretty(l4_range, n = 5L)
l4_ticks <- l4_ticks[l4_ticks >= l4_range[[1L]] & l4_ticks <= l4_range[[2L]]]

draw_l4 <- function() {
  old <- par(mar = c(4.8, 5.0, 4.2, 1.2), xaxs = "i", yaxs = "i")
  on.exit(par(old))
  xpad <- diff(range(coordinate_table$UMAP_1)) * 0.05
  ypad <- diff(range(coordinate_table$UMAP_2)) * 0.07
  plot(reference_plot$UMAP_1, reference_plot$UMAP_2, pch = 16,
    col = grDevices::adjustcolor(reference_plot$l4_color, alpha.f = 0.68),
    cex = 0.68, xlab = "PeakVI UMAP 1", ylab = "PeakVI UMAP 2",
    xlim = range(coordinate_table$UMAP_1) + c(-xpad, xpad),
    ylim = range(coordinate_table$UMAP_2) + c(-ypad, ypad),
    main = "Native fitted Slingshot L4 on the joint PeakVI UMAP")
  mtext("Large symbols are ImmGen libraries; only borderline libraries are labeled",
    side = 3, line = 0.5, cex = 0.82)
  points(query_plot$UMAP_1, query_plot$UMAP_2, pch = 21,
    bg = query_plot$l4_color, col = query_plot$border, cex = 2.0,
    lwd = ifelse(query_plot$support_call == "borderline", 2.8, 1.25))
  text(query_plot$UMAP_1[borderline_index], query_plot$UMAP_2[borderline_index],
    labels = query_plot$cell_id[borderline_index],
    pos = borderline_label_position, offset = 0.58, cex = 0.64, font = 2)
  legend("topleft", legend = format(l4_ticks, trim = TRUE), pch = 16,
    col = l4_palette[l4_index(l4_ticks)], pt.cex = 1.25, bty = "n",
    cex = 0.78, title = "Native L4")
  legend("topright", legend = c("ImmGen library", "borderline support"),
    pch = 21, pt.bg = "white", col = c("#202020", "#D55E00"),
    pt.lwd = c(1.25, 2.8), pt.cex = 1.35, bty = "n", cex = 0.75)
}
with_devices(file.path(dirs[["figures"]],
  "02_native_l4_peakvi_umap"), 11.5, 8.5, draw_l4)

source_manifest_after <- sha256(file.path(cfg$inputs$mapping_bundle$path,
  "manifest.json"))
derived_manifest_after <- sha256(file.path(cfg$storage$derived_bundle,
  "manifest.json"))
stopifnot(identical(source_manifest_before, source_manifest_after),
  umap_already_promoted || !identical(derived_manifest_before, derived_manifest_after),
  sum(query_plot$support_call == "supported") == 25L,
  sum(query_plot$support_call == "borderline") == 3L)

summary <- list(schema_version = "bdev.immgen_peakvi_umap_summary.v1",
  status = "complete", run_identity = cfg$run_identity,
  design = design,
  resumed_promoted_umap = umap_already_promoted,
  reference_cells = nrow(reference_plot), query_libraries = nrow(query_plot),
  populations = length(populations), supported = sum(query_plot$support_call == "supported"),
  borderline = sum(query_plot$support_call == "borderline"),
  out_of_reference = sum(query_plot$support_call == "out_of_reference"),
  umap = cfg$umap, native_l4_orientation = "native",
  source_manifest_sha256 = source_manifest_after,
  derived_manifest_sha256 = derived_manifest_after,
  coordinate_table = coordinate_path,
  figures = list(population = file.path(dirs[["figures"]],
    "01_immgen_population_peakvi_umap.png"), native_l4 =
    file.path(dirs[["figures"]], "02_native_l4_peakvi_umap.png")),
  paper_analysis = FALSE, implementation_unlocked = FALSE, sif_built = FALSE)
write_json(summary, file.path(dirs[["reports"]],
  "01_immgen_peakvi_umap_summary.json"))
receipt <- c(summary, list(config_path = cfg_path, config_sha256 = sha256(cfg_path),
  scwbRepresentation_manifest_sha256 = sha256(pack_manifest),
  coordinates_sha256 = sha256(coordinate_path), source_mutated = FALSE))
write_json(receipt, file.path(dirs[["manifests"]], "02_execution_receipt.json"))
cat("Completed", design, "UMAP with", nrow(coordinate_table),
  "coordinates\n")
