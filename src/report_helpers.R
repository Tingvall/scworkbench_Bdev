format_cell_count <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)

select_table_columns <- function(x, columns) {
  if (!is.data.frame(x)) stop("x must be a data.frame")
  x[, intersect(columns, names(x)), drop = FALSE]
}

report_run_directory <- function(reports_root, bundle_id, run_id, create = FALSE) {
  values <- c(bundle_id = bundle_id, run_id = run_id)
  if (any(!grepl("^[A-Za-z0-9_.-]+$", values))) {
    stop("bundle_id and run_id may contain only letters, numbers, dot, underscore, and hyphen")
  }
  path <- file.path(reports_root, bundle_id, run_id)
  if (isTRUE(create)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}
