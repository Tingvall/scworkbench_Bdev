# Shared configuration and public scWorkbench operations for Parse-to-Bdev.
suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
  library(yaml)
})
source(file.path(getwd(), "src", "bundle_config.R"))
source(file.path(getwd(), "src", "bundle_validation.R"))
source(file.path(getwd(), "src", "reference_validation.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, grepl("^[A-Za-z0-9_]+$", args[[1]]))
integration_run <- args[[1]]
registry <- load_bundle_registry()
cfg <- read_project_yaml("config/bdev_integrated.yml")
parse_cfg <- read_project_yaml("config/parse_cytokines.yml")
ref_cfg <- read_project_yaml("config/bdev_reference.yml")
root <- cfg$storage_root
report_dir <- file.path(root, "reports", integration_run)
bundle_dir <- cfg$bundle$path
session_root <- file.path(root, "sessions", integration_run)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_root, recursive = TRUE, showWarnings = FALSE)
integration_runtime <- registry$bundles$bdev_integrated$runtime %||% registry$bundles$parse_cytokines$runtime
stopifnot(identical(as.character(packageVersion("scworkbench")), integration_runtime$scworkbench_version))
options(warn = 1, width = 180)
log_event <- function(...) { cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz="UTC"), ..., "\n"); flush.console() }
write_json <- function(x, name) jsonlite::write_json(x, file.path(report_dir, name),
  pretty=TRUE, auto_unbox=TRUE, null="null", na="null", digits=NA)
write_tsv <- function(x, name) write.table(x, file.path(report_dir, name),
  sep="\t", quote=TRUE, row.names=FALSE, na="NA")
read_json <- function(name) jsonlite::read_json(file.path(report_dir, name), simplifyVector=TRUE)
open_parent <- function(id) wb_open_bundle(registry$bundles[[id]]$bundle_dir, read_only=TRUE, validate=FALSE)
check_parents <- function() {
  out <- lapply(c("bdev_reference", "parse_cytokines"), function(id) {
    checked <- validate_registered_bundle(registry, id, run_recipient_validation=FALSE)
    path <- if (id == "bdev_reference") ref_cfg$validation$report else registry$bundles[[id]]$validation_report
    validation <- jsonlite::read_json(path, simplifyVector=FALSE)
    stopifnot(isTRUE(validation$ok))
    checked$accepted_validation_report <- path
    checked$accepted_validation_sha256 <- wb_source_sha256(path)
    checked
  })
  names(out) <- c("bdev_reference", "parse_cytokines")
  out
}
artifact_id <- function(promotion, output) {
  rows <- promotion$promoted[promotion$promoted$output == output & promotion$promoted$status == "ok", , drop=FALSE]
  stopifnot(nrow(rows) == 1L)
  as.character(rows$artifact_id[[1]])
}
# Each operation has its own durable session. Completed promotions are reused
# on restart; an incomplete operation resumes the last verified checkpoint.
run_operation <- function(bundle, label, op, inputs, outputs, notes="") {
  receipt_path <- file.path(report_dir, paste0(label, "_promotion.rds"))
  if (file.exists(receipt_path)) {
    receipt <- readRDS(receipt_path)
    log_event("Reusing promoted operation", label)
    return(receipt)
  }
  interface <- wb_op_interface(op, bundle=bundle)
  write_json(interface, paste0(label, "_interface.json"))
  workspace <- file.path(session_root, label)
  session <- if (dir.exists(workspace)) wb_session_resume(bundle, workspace) else
    wb_session_create(bundle, workspace_root=session_root, session_id=label)
  log_event("Running", label, op)
  if (!is.null(session$durable$promotion_receipt)) {
    stop("Operation is already promoted but its project receipt is missing: ", label,
         ". Inspect the durable promotion artifact IDs before restoring the receipt; do not rerun the operation.")
  } else {
    if (length(session$runs)) {
      stopifnot(length(session$runs)==1L)
      completed <- session$runs[[session$head_run_id]]
      stopifnot(all(outputs %in% names(completed$outputs)))
      log_event("Promoting the existing completed checkpoint", label)
    } else {
      rr <- wb_session_run_op(session, bundle, op, inputs=inputs, notes=notes)
      session <- rr$session
    }
    promotion <- wb_session_promote_current(bundle, session, outputs=outputs,
      name_prefix=label, description=notes, created_by="jtg", fail_fast=TRUE)
  }
  staged <- wb_session_view(promotion$session)
  stopifnot(!identical(normalizePath(staged$dir),normalizePath(bundle$dir)))
  wb_close_bundle(staged)
  promotion <- list(promoted=promotion$promoted,session_workspace=workspace,op=op)
  saveRDS(promotion, receipt_path)
  write_json(promotion$promoted, paste0(label, "_promotion.json"))
  promotion
}
load_output <- function(bundle, promotion, output) wb_art_load(bundle, artifact_id(promotion, output))
pt_family <- function(bundle) {
  cells <- wb_cells(bundle, "all")
  md <- wb_meta(bundle, cells=cells)
  use <- is.finite(md$pt)
  out <- wb_prim_cell_vector_family(rep("bdev_pt", sum(use)), cells[use], md$pt[use],
    meta=list(source_bundle="bdev_reference", source_field="pt",
      source_artifact_id=ref_cfg$reference_state$pseudotime$artifact_id,
      authoritative_values="exact source metadata doubles; legacy scalar artifact rounds values",
      source_fingerprint=registry$bundles$bdev_reference$accepted_state$bundle_state_fingerprint,
      direction="as_is", fitted_model=FALSE))
  stopifnot(identical(out$value, as.numeric(md$pt[use])), sum(use)==12642L,
            sum(is.na(md$pt))==14081L)
  out
}

# Payload read evidence contains bundle-specific storage paths. Compare every
# coordinate, dimension and identifier without that non-coordinate attribute.
integration_space_values_equal <- function(a,b) {
  attr(a,"scworkbench_payload_read") <- NULL
  attr(b,"scworkbench_payload_read") <- NULL
  identical(a,b)
}
