# Project bundle-registry helpers. These describe repository configuration;
# scWorkbench remains authoritative for bundle contents and operations.

read_project_yaml <- function(path) {
  if (!file.exists(path)) stop("Configuration file does not exist: ", path)
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The yaml package is required to read project configuration")
  }
  value <- yaml::read_yaml(path)
  if (!is.list(value)) stop("Expected a YAML mapping in: ", path)
  value
}

load_bundle_registry <- function(repository_root = ".") {
  path <- file.path(repository_root, "config", "bundles.yml")
  registry <- read_project_yaml(path)
  if (!identical(as.integer(registry$schema_version), 1L)) {
    stop("Unsupported bundle registry schema_version")
  }
  if (!is.list(registry$bundles) || !length(registry$bundles)) {
    stop("Bundle registry contains no bundle entries")
  }
  registry
}

bundle_registry_entry <- function(registry, bundle_id) {
  bundle_id <- as.character(bundle_id)[1L]
  entry <- registry$bundles[[bundle_id]]
  if (is.null(entry)) stop("Unknown bundle_id: ", bundle_id)
  entry$bundle_id <- bundle_id
  entry
}

validate_bundle_registry <- function(registry, repository_root = ".") {
  ids <- names(registry$bundles)
  if (is.null(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Bundle registry IDs must be non-empty and unique")
  }
  allowed_status <- c("planned", "working", "accepted")
  for (bundle_id in ids) {
    entry <- bundle_registry_entry(registry, bundle_id)
    if (!entry$status %in% allowed_status) {
      stop("Invalid status for ", bundle_id, ": ", entry$status)
    }
    if (!is.character(entry$bundle_dir) || length(entry$bundle_dir) != 1L ||
        !nzchar(entry$bundle_dir)) {
      stop("bundle_dir must be a non-empty string for ", bundle_id)
    }
    parents <- as.character(entry$parents %||% character(0))
    if (bundle_id %in% parents) stop(bundle_id, " cannot be its own parent")
    unknown <- setdiff(parents, ids)
    if (length(unknown)) {
      stop(bundle_id, " has unknown parents: ", paste(unknown, collapse = ", "))
    }
    config_path <- file.path(repository_root, entry$config)
    if (!file.exists(config_path)) stop("Missing bundle config: ", config_path)
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
