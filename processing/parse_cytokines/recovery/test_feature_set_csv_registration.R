#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(scworkbench))
main <- function() {
  root <- tempfile("scworkbench-csv-registration-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cells <- paste0("c", seq_len(12L))
  genes <- paste0("g", seq_len(8L))
  bundle_dir <- file.path(root, "bundle")
  b <- wb_init_bundle(bundle_dir)
  wb_define_cells(b, cells, dataset_id = "csv_fixture")
  wb_attach_matrix(b, "RNA_counts", Matrix::Matrix(
    matrix(seq_len(96L), nrow = 8L, dimnames = list(genes, cells)), sparse = TRUE
  ), feature_type = "gene")
  wb_close_bundle(b)
  b <- wb_open_bundle(bundle_dir)
  on.exit(try(wb_close_bundle(b), silent = TRUE), add = TRUE)
  input <- data.frame(feature_set = rep(c("signal.up", "signal.down"), each = 3L),
                       set_name = rep(c("signal.up", "signal.down"), each = 3L),
                       feature_id = genes[c(3L, 1L, 2L, 8L, 6L, 7L)],
                       description = 'Description with a comma, and "quotes".')
  csv <- file.path(root, "sets.csv")
  write.csv(input, csv, row.names = FALSE)
  stopifnot(identical(read.csv(csv, stringsAsFactors = FALSE), input))
  run <- wb_session_run_op(wb_session_new(), b, "feature_set_import", inputs = list(
    matrix = "RNA_counts", source_path = csv, source_format = "table", namespace = "csv_test",
    set_id_col = "feature_set", set_name_col = "set_name", gene_col = "feature_id",
    description_col = "description", overwrite = FALSE, strict = TRUE, is_readonly = TRUE
  ))
  manifest <- run$run$outputs$import_manifest$table
  stopifnot(nrow(manifest) == 2L, all(manifest$status == "imported"),
            all(manifest$n_genes_raw == manifest$n_genes_mapped_to_bundle))
  wb_session_promote_current(b, run$session, outputs = list(
    wb_ref_run(run$run$run_id, "import_summary"), wb_ref_run(run$run$run_id, "import_manifest"),
    wb_ref_run(run$run$run_id, "rejected_sets")
  ), fail_fast = TRUE)
  registered <- wb_feature_sets(b, matrix_id = "RNA_counts", namespace = "csv_test")
  stopifnot(nrow(registered) == 2L)
  for (id in unique(input$feature_set)) {
    stopifnot(identical(wb_feature_set_get(b, id, "RNA_counts"), input$feature_id[input$feature_set == id]))
  }
  stopifnot(all(registered$is_readonly == 1L))
  stopifnot(isTRUE(wb_finalize_recipient_bundle(b)$ready))
  wb_close_bundle(b)
  validation <- wb_bundle_validation_report(bundle_dir, level = "recipient")
  stopifnot(isTRUE(validation$ok))
  cat("PASS: CSV registration preserves exact ordered gene membership and read-only flags, with recipient validation.\n")
}
main()
