#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(scworkbench))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
source_root <- normalizePath(args[1], mustWork = TRUE)
report_dir <- normalizePath(args[2], mustWork = TRUE)
contract_test_skip <- function(reason) stop("Required recovery contract skipped: ", reason)
run_contract <- function(test) {
  if (!"bundle_dir" %in% names(formals(test))) return(test())
  root <- tempfile("scworkbench-recovery-contract-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  bundle <- wb_init_bundle(root)
  cells <- paste0("cell", seq_len(12L))
  features <- paste0("gene", seq_len(8L))
  wb_define_cells(bundle, cells, dataset_id = "recovery_test")
  wb_attach_matrix(bundle, "RNA_counts", Matrix::Matrix(
    matrix(seq_len(96L), nrow = 8L, dimnames = list(features, cells)), sparse = TRUE
  ), feature_type = "gene")
  wb_close_bundle(bundle)
  test(bundle_dir = root)
}
for (name in c("atomic_authoritative_manifest", "cell_vector_durable_identity_roundtrip", "deep_bundle_integrity_validation",
               "op_feature_set_import_qc")) {
  environment <- new.env(parent = globalenv())
  sys.source(file.path(source_root, "tests", paste0("test_", name, "_contract.R")), environment)
  tests <- ls(environment, pattern = "^test_.*contract$")
  stopifnot(length(tests) > 0L)
  for (test in tests) {
    cat("RUN", test, "\n")
    run_contract(get(test, envir = environment))
    cat("PASS", test, "\n")
  }
}
audit <- wb_refresh_semantic_catalogs(docs_dir = file.path(report_dir, "semantic_catalogs"))
jsonlite::write_json(audit, file.path(report_dir, "semantic_audit.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null")
cat("PASS: semantic catalog audit\n")
