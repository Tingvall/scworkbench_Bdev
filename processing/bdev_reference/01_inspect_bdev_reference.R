# Read-only source inspection; does not create or modify a scWorkbench bundle.
# Usage: Rscript --vanilla processing/bdev_reference/01_inspect_bdev_reference.R SOURCE OUTPUT_DIR
# Run in the supplied scWorkbench image with SOURCE's directory mounted :ro.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 01_inspect_bdev_reference.R SOURCE OUTPUT_DIR")
source_path <- args[[1L]]
out_dir <- args[[2L]]
if (!file.exists(source_path)) stop("Source file does not exist: ", source_path)
if (file.exists(file.path(out_dir, "inspection.json"))) stop("Use a new output directory; inspection.json already exists")
options(width = 160, warn = 1)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Signac)
  library(slingshot)
  library(SingleCellExperiment)
  library(jsonlite)
})
log_event <- function(...) { cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n"); flush.console() }
safe <- function(expr) tryCatch(expr, error = function(e) list(error = conditionMessage(e)))
is_error <- function(x) is.list(x) && identical(names(x), "error")
info_before <- file.info(source_path)[, c("size", "mtime", "ctime")]
con <- gzfile(source_path, open = "rb")
workspace <- new.env(parent = baseenv())
log_event("Loading source workspace through a read-only connection")
loaded <- load(con, envir = workspace)
close(con)
log_event("Loaded objects:", paste(loaded, collapse = ", "))
inventory <- lapply(loaded, function(nm) {
  x <- get(nm, envir = workspace, inherits = FALSE)
  list(name = nm, class = class(x), dimensions = safe(dim(x)))
})

id_sets <- list()
id_records <- list()
register_ids <- function(label, ids, expected_count = NULL, axis = "cells") {
  present <- !is.null(ids)
  if (present) ids <- as.character(ids)
  id_sets[[label]] <<- if (present) ids else NA_character_
  id_records[[label]] <<- list(
    label = label, axis = axis, present = present, n = if (present) length(ids) else NULL,
    expected_count = expected_count,
    duplicated = if (present) sum(duplicated(ids)) else NULL,
    missing_or_empty = if (present) sum(is.na(ids) | !nzchar(ids)) else NULL,
    first = if (present) head(ids, 3L) else NULL,
    last = if (present) tail(ids, 3L) else NULL
  )
  label
}
matrix_summary <- function(x, label, cell_axis = 2L) {
  d <- dim(x)
  cells <- if (cell_axis == 2L) colnames(x) else rownames(x)
  list(class = class(x), dim = d,
       feature_count = if (cell_axis == 2L) d[[1L]] else NULL,
       feature_ids_present = if (cell_axis == 2L) !is.null(rownames(x)) else NULL,
       duplicated_feature_ids = if (cell_axis == 2L && !is.null(rownames(x))) sum(duplicated(rownames(x))) else NULL,
       cell_id_record = register_ids(label, cells, d[[cell_axis]]))
}
seurat_objects <- list()
candidates <- list()
seurat_refs <- list()
candidate_refs <- list()
name_hits <- list()
scan_omissions <- list()
scanned_s4_slots <- list()

inspect_trajectory <- function(x, path) {
  log_event("Inspecting trajectory candidate:", path)
  candidate_refs[[path]] <<- x
  result <- list(path = path, class = class(x))
  accessors <- c(pseudotime = "slingPseudotime", weights = "slingCurveWeights", reduced_dim = "slingReducedDim", cluster_labels = "slingClusterLabels")
  for (nm in names(accessors)) {
    v <- safe(getExportedValue("slingshot", accessors[[nm]])(x))
    if (is_error(v)) result[[nm]] <- v else {
      result[[nm]] <- matrix_summary(v, paste0(path, "::", nm), 1L)
      result[[nm]]$column_names <- colnames(v)
    }
  }
  result$lineages <- safe(slingshot::slingLineages(x))
  result$params <- safe(slingshot::slingParams(x))
  curves <- safe(slingshot::slingCurves(x))
  if (is_error(curves)) result$curves <- curves else {
    result$curves <- lapply(seq_along(curves), function(i) {
      cv <- curves[[i]]
      prefix <- paste0(path, "::curve[", i, "]")
      list(name = names(curves)[i], class = class(cv), components = names(cv),
           s = if (!is.null(cv$s)) matrix_summary(cv$s, paste0(prefix, "$s"), 1L) else NULL,
           lambda_length = length(cv$lambda), lambda_named = !is.null(names(cv$lambda)),
           lambda_id_record = if (!is.null(names(cv$lambda))) register_ids(paste0(prefix, "$lambda"), names(cv$lambda), length(cv$lambda)) else NULL,
           weights_length = length(cv$w),
           weights_id_record = if (!is.null(names(cv$w))) register_ids(paste0(prefix, "$w"), names(cv$w), length(cv$w)) else NULL,
           order_length = length(cv$ord), order_is_permutation = if (!is.null(cv$s)) identical(sort(as.integer(cv$ord)), seq_len(nrow(cv$s))) else NULL)
    })
  }
  candidates[[path]] <<- result
  invisible(NULL)
}

scan <- function(x, path, depth = 0L) {
  if (depth > 20L) { scan_omissions[[path]] <<- "Depth limit (20)"; return(invisible(NULL)) }
  cls <- class(x)
  if (inherits(x, "SlingshotDataSet") || inherits(x, "PseudotimeOrdering")) {
    inspect_trajectory(x, path)
    return(invisible(NULL))
  }
  if (inherits(x, "Seurat")) {
    inspect_seurat(x, path)
    return(invisible(NULL))
  }
  if (inherits(x, "SingleCellExperiment")) {
    if ("slingshot" %in% names(SummarizedExperiment::colData(x))) inspect_trajectory(x, path)
    scan(S4Vectors::metadata(x), paste0(path, "@metadata"), depth + 1L)
    cd <- SummarizedExperiment::colData(x)
    for (nm in names(cd)) if (isS4(cd[[nm]]) || is.list(cd[[nm]])) scan(cd[[nm]], paste0(path, "::colData[[", encodeString(nm, quote = '"'), "]]"), depth + 1L)
    return(invisible(NULL))
  }
  if (is.environment(x) || is.function(x) || is.language(x)) {
    scan_omissions[[path]] <<- paste("Not traversed:", paste(cls, collapse = "/"))
    return(invisible(NULL))
  }
  if (is.data.frame(x)) {
    for (nm in names(x)) if (is.list(x[[nm]]) || isS4(x[[nm]])) scan(x[[nm]], paste0(path, "[[", encodeString(nm, quote = '"'), "]]"), depth + 1L)
    return(invisible(NULL))
  }
  if (inherits(x, "Matrix") || is.matrix(x) || is.atomic(x) || is.null(x)) return(invisible(NULL))
  if (is.list(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      nm <- if (!is.null(nms) && !is.na(nms[i]) && nzchar(nms[i])) nms[i] else as.character(i)
      child_path <- paste0(path, "[[", if (!is.null(nms) && !is.na(nms[i]) && nzchar(nms[i])) encodeString(nms[i], quote = '"') else i, "]]")
      if (grepl("sling|lineage|pseudo|trajectory|^sds", nm, ignore.case = TRUE)) name_hits[[child_path]] <<- list(class = class(x[[i]]), dim = safe(dim(x[[i]])))
      scan(x[[i]], child_path, depth + 1L)
    }
  } else if (isS4(x)) {
    scanned_s4_slots[[path]] <<- slotNames(x)
    for (nm in slotNames(x)) scan(slot(x, nm), paste0(path, "@", nm), depth + 1L)
  }
  invisible(NULL)
}

inspect_seurat <- function(x, path) {
  log_event("Inspecting Seurat object:", path)
  seurat_refs[[path]] <<- x
  scanned_s4_slots[[path]] <<- slotNames(x)
  global_cells <- colnames(x)
  obj <- list(path = path, class = class(x), object_version = as.character(x@version),
              default_assay = SeuratObject::DefaultAssay(x), n_cells = length(global_cells),
              global_cell_ids = register_ids(paste0(path, "::cells"), global_cells, length(global_cells)))
  md <- x[[]]
  obj$metadata <- list(n_rows = nrow(md), n_columns = ncol(md),
                       cell_id_record = register_ids(paste0(path, "::metadata"), rownames(md), nrow(md)),
                       columns = lapply(names(md), function(nm) list(name = nm, class = class(md[[nm]]), missing = sum(is.na(md[[nm]])))))
  obj$active_ident <- list(n = length(SeuratObject::Idents(x)), levels = levels(SeuratObject::Idents(x)), cell_id_record = register_ids(paste0(path, "::active_ident"), names(SeuratObject::Idents(x)), length(SeuratObject::Idents(x))))
  obj$assays <- list()
  log_event("Cells:", length(global_cells), "Assays:", paste(SeuratObject::Assays(x), collapse = ", "))
  for (an in SeuratObject::Assays(x)) {
    a <- x[[an]]
    ar <- list(class = class(a), n_features = nrow(a), n_cells = ncol(a), feature_id_example = head(rownames(a), 5L),
               cell_id_record = register_ids(paste0(path, "::assay:", an), colnames(a), ncol(a)),
               layers = list())
    layer_names <- SeuratObject::Layers(a, search = NA)
    for (ln in layer_names) ar$layers[[ln]] <- safe(matrix_summary(SeuratObject::LayerData(a, layer = ln), paste0(path, "::assay:", an, ":", ln), 2L))
    ar$feature_metadata_columns <- safe(names(a[[]]))
    if (inherits(a, "ChromatinAssay")) {
      ar$genome <- safe(as.character(unique(GenomeInfoDb::genome(Signac::granges(a)))))
      ar$fragments <- safe(lapply(Signac::Fragments(a), function(f) list(path = Signac::GetFragmentData(f, slot = "path"), n_cells = length(Signac::GetFragmentData(f, slot = "cells")))))
    }
    obj$assays[[an]] <- ar
    assay_path <- paste0(path, "@assays[[", encodeString(an, quote = '"'), "]]")
    scanned_s4_slots[[assay_path]] <<- slotNames(a)
    for (sn in setdiff(slotNames(a), c("counts", "data", "scale.data", "layers", "cells", "features"))) {
      scan(slot(a, sn), paste0(assay_path, "@", sn), 1L)
    }
  }
  obj$reductions <- list()
  for (rn in SeuratObject::Reductions(x)) {
    r <- x[[rn]]
    em <- SeuratObject::Embeddings(r)
    obj$reductions[[rn]] <- c(matrix_summary(em, paste0(path, "::reduction:", rn), 1L),
                              list(assay = SeuratObject::DefaultAssay(r), key = SeuratObject::Key(r), columns = colnames(em),
                                   loadings_dim = safe(dim(SeuratObject::Loadings(r))), stdev_length = length(SeuratObject::Stdev(r)),
                                   nonfinite_values = sum(!is.finite(em))))
    reduction_path <- paste0(path, "@reductions[[", encodeString(rn, quote = '"'), "]]")
    scanned_s4_slots[[reduction_path]] <<- slotNames(r)
    for (sn in setdiff(slotNames(r), c("cell.embeddings", "feature.loadings", "feature.loadings.projected"))) {
      scan(slot(r, sn), paste0(reduction_path, "@", sn), 1L)
    }
  }
  obj$misc_entries <- as.character(names(x@misc))
  obj$tools_entries <- as.character(names(x@tools))
  obj$graphs <- names(x@graphs)
  obj$commands <- names(x@commands)
  seurat_objects[[path]] <<- obj
  # Save partial inventory before traversing potentially complex source extras.
  jsonlite::write_json(list(workspace_objects = inventory, seurat_objects = seurat_objects),
                       file.path(out_dir, "inventory.json"), pretty = TRUE, auto_unbox = TRUE,
                       null = "null", digits = 16, na = "null")
  for (sn in setdiff(slotNames(x), c("assays", "reductions", "active.ident"))) {
    scan(slot(x, sn), paste0(path, "@", sn), 1L)
  }
  invisible(NULL)
}

for (nm in loaded) scan(get(nm, envir = workspace, inherits = FALSE), paste0("workspace[[", encodeString(nm, quote = '"'), "]]"))

comparisons <- list()
for (sp in names(seurat_objects)) {
  base_label <- seurat_objects[[sp]]$global_cell_ids
  baseline <- id_sets[[base_label]]
  comparisons[[sp]] <- lapply(names(id_records), function(label) {
    rec <- id_records[[label]]
    ids <- id_sets[[label]]
    valid <- isTRUE(rec$present)
    list(target = label, n = rec$n, ids_present = valid,
         same_set = if (valid) setequal(baseline, ids) else NULL,
         same_order = if (valid) identical(baseline, ids) else NULL,
         missing_from_target = if (valid) length(setdiff(baseline, ids)) else NULL,
         extra_in_target = if (valid) length(setdiff(ids, baseline)) else NULL,
         missing_examples = if (valid) head(setdiff(baseline, ids), 5L) else NULL,
         extra_examples = if (valid) head(setdiff(ids, baseline), 5L) else NULL)
  })
}

coordinate_matches <- list()
for (cp in names(candidates)) {
  candidate <- candidate_refs[[cp]]
  reduced <- safe(slingshot::slingReducedDim(candidate))
  if (is_error(reduced) || is.null(dim(reduced))) next
  for (sp in names(seurat_objects)) {
    seurat <- seurat_refs[[sp]]
    for (rn in SeuratObject::Reductions(seurat)) {
      em <- SeuratObject::Embeddings(seurat[[rn]])
      can_align <- !is.null(rownames(reduced)) && !is.null(rownames(em)) && !anyDuplicated(rownames(reduced)) && all(rownames(reduced) %in% rownames(em))
      rec <- list(candidate = cp, seurat_object = sp, reduction = rn,
                  trajectory_dim = dim(reduced), reduction_dim = dim(em), can_align_cells = can_align)
      if (can_align && ncol(reduced) == ncol(em)) {
        aligned <- em[rownames(reduced), , drop = FALSE]
        rec$max_absolute_difference <- max(abs(reduced - aligned))
        rec$coordinates_equal_1e8 <- is.finite(rec$max_absolute_difference) && rec$max_absolute_difference <= 1e-8
      }
      coordinate_matches[[length(coordinate_matches) + 1L]] <- rec
    }
  }
}
info_after <- file.info(source_path)[, c("size", "mtime", "ctime")]
report <- list(source = list(path = source_path, format = "R workspace (RDX3), loaded into isolated environment", before = info_before, after = info_after, stat_unchanged = identical(info_before, info_after)),
               runtime = list(R = R.version.string, packages = lapply(c("SeuratObject", "Signac", "slingshot", "SingleCellExperiment"), function(p) list(name = p, version = as.character(packageVersion(p))))),
               workspace_objects = inventory, seurat_objects = seurat_objects,
               trajectory_candidates = candidates, trajectory_name_hits = name_hits,
               cell_identity_records = id_records, cell_comparisons = comparisons,
               trajectory_coordinate_comparisons = coordinate_matches,
               scanned_s4_slots = scanned_s4_slots, scan_omissions = scan_omissions)
jsonlite::write_json(report, file.path(out_dir, "inspection.json"), pretty = TRUE, auto_unbox = TRUE, null = "null", digits = 16, na = "null")
capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))
log_event("Inspection complete; report:", file.path(out_dir, "inspection.json"))
