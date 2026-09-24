#!/usr/bin/env Rscript
source("src/integration_workflow.R")
main <- function() {
  log_event("Checking accepted parent identities and retained validation evidence")
  parents <- check_parents()
  bdev <- open_parent("bdev_reference"); on.exit(wb_close_bundle(bdev), add=TRUE)
  parse <- open_parent("parse_cytokines"); on.exit(wb_close_bundle(parse), add=TRUE)
  cells <- wb_cells(bdev, "all")
  md <- wb_meta(bdev, cells=cells)
  stopifnot(length(cells)==26723L, identical(md$cell_id, cells), !anyDuplicated(cells))
  pt <- pt_family(bdev)
  features <- wb_feature_ids(bdev, "RNA_counts")
  stopifnot(length(features)==32285L, !anyDuplicated(features))
  sets <- read.delim(parse_cfg$signature_analysis$final$feature_set_tsv, check.names=FALSE)
  stopifnot(length(unique(sets$feature_set))==23L, nrow(sets)==720L)
  registered <- wb_feature_sets(parse, matrix_id="RNA_counts", namespace="parse_cytokine_signatures")
  stopifnot(setequal(registered$feature_set, unique(sets$feature_set)))
  coverage <- do.call(rbind, lapply(unique(sets$feature_set), function(id) {
    genes <- sets$feature_id[sets$feature_set==id]
    stopifnot(identical(as.character(wb_feature_set_get(parse,id,"RNA_counts")), as.character(genes)))
    mapped <- intersect(genes, features)
    data.frame(feature_set=id, cytokine=sub("\\.(up|down)$", "", sub("^parse_cytokine\\.","",id)),
      direction=sub(".*\\.","",id), original_genes=length(genes), mapped_genes=length(mapped),
      coverage=length(mapped)/length(genes), eligible=length(mapped)>=5 && length(mapped)/length(genes)>=0.8,
      missing_genes=paste(setdiff(genes,features),collapse=";"))
  }))
  sample_coverage <- do.call(rbind, lapply(sort(unique(md$orig.ident)), function(id) {
    rows <- md$orig.ident==id; values <- md$pt[rows & is.finite(md$pt)]
    data.frame(sample_id=id,cells=sum(rows),pt_cells=length(values),
      pt_min=if(length(values)) min(values) else NA_real_,pt_max=if(length(values)) max(values) else NA_real_)
  }))
  write_tsv(coverage,"feature_coverage.tsv")
  write_tsv(sample_coverage,"sample_coverage.tsv")
  write_tsv(sets,"accepted_source_gene_sets.tsv")
  ops <- c("feature_set_import","feature_set_overlap_qc","ucell","trajectory_signal_dynamics",
    "cell_vector_family_qc","cell_vector_family_calibrate","metadata_cell_vector",
    "matrix_log_normalize","feature_set_namespace_transfer","trajectory_geneset_dynamics")
  for(op in ops) write_json(wb_op_interface(op,bundle=bdev),paste0("preflight_",op,"_interface.json"))
  write_json(list(status="complete",parents=parents,cells=length(cells),defined_pt=length(pt$value),
    missing_pt=sum(is.na(md$pt)),pseudotime_exact=TRUE,source_sets=23L,eligible_sets=sum(coverage$eligible),
    source_feature_set_sha256=wb_source_sha256(parse_cfg$signature_analysis$final$feature_set_tsv),
    runtime=registry$bundles$parse_cytokines$runtime,
    sample_replication_status="four source sample IDs; independence not established",
    genome_policy="gene-symbol transfer only; no peak/genome-coordinate integration"),"preflight.json")
  capture.output(sessionInfo(),file=file.path(report_dir,"session_info.txt"))
  log_event("Preflight complete",sum(coverage$eligible),"eligible signatures")
  print(coverage); print(sample_coverage)
}
main()
