#!/usr/bin/env Rscript
source("src/integration_workflow.R")
main <- function() {
  preflight <- read_json("preflight.json")
  stopifnot(preflight$status=="complete", !file.exists(file.path(report_dir,"creation.json")))
  parents <- check_parents()
  bdev <- open_parent("bdev_reference"); on.exit(wb_close_bundle(bdev),add=TRUE)
  cells <- wb_cells(bdev,"all")
  md <- wb_meta(bdev,cells=cells)
  metadata <- md[,setdiff(names(md),"dataset_id"),drop=FALSE]
  metadata$reference_source_bundle <- "bdev_reference"
  metadata$reference_source_cell_id <- cells
  space_ids <- ref_cfg$bundle$spaces
  spaces <- setNames(lapply(space_ids, function(id) {
    x <- wb_space(bdev,id,cells=cells)
    if(is.data.frame(x)) { m <- as.matrix(x[,c("x","y")]); rownames(m) <- x$cell_id; m } else x
  }),space_ids)
  source_file <- file.path(report_dir,"preflight.json")
  source_identity <- list(dataset_name="bdev_integrated",phase="parse_cytokine_pseudotime",
    reference_parent="bdev_reference",parents=parents,
    inputs=list(wb_source_input(path=source_file, source_id="accepted_parent_preflight", accession=integration_run,
      role="parent_bundle_identity_and_signature_source",format="json",scope="external_file",
      sha256=wb_source_sha256(source_file))),
    extraction="scWorkbench public getters from accepted read-only Bdev; selected RNA_counts, all metadata and five spaces",
    pseudotime=list(field="pt",policy="preserve exact doubles and missing entries",refit=FALSE),
    parse_input=list(namespace="parse_cytokine_signatures",source_sha256=preflight$source_feature_set_sha256))
  wb_source_preflight(source_identity)
  log_event("Creating Bdev-centered integration bundle")
  b <- if (dir.exists(bundle_dir)) {
    existing <- wb_open_bundle(bundle_dir,validate=FALSE)
    stopifnot(wb_manifest(existing)$bundle_id=="bdev_integrated")
    if(length(wb_cells(existing,"all"))) stopifnot(file.exists(file.path(report_dir,"integration_ingest_report.json")))
    existing
  } else wb_init_bundle(bundle_dir,bundle_id="bdev_integrated",bundle_name="Bdev cytokine response along existing pseudotime",overwrite=FALSE)
  on.exit(try(wb_close_bundle(b),silent=TRUE),add=TRUE)
  if (!length(wb_cells(b,"all"))) {
    ingest <- wb_ingest_apply(b,set_name="bdev_reference",cell_ids=cells,metadata=metadata,
    spaces=spaces, matrices=list(RNA_counts=list(mat=methods::as(wb_matrix(bdev,"RNA_counts"),"dgCMatrix"),layer="counts",
      feature_type="gene",description="Exact accepted Bdev reference RNA counts; no cross-dataset cell merge")),
    source="accepted bdev_reference and parse_cytokines; Parse gene-set transfer",
    source_identity=source_identity,created_by="jtg",apply_mode="atomic",emit_report_artifact=TRUE)
    write_json(attr(ingest,"ingest_report"),"integration_ingest_report.json")
  }
  out_md <- wb_meta(b,cells=cells)
  stopifnot(identical(wb_cells(b,"all"),cells))
  metadata_checks <- vapply(setdiff(names(md),"dataset_id"),function(nm) reference_values_equal(md[[nm]],out_md[[nm]]),logical(1))
  stopifnot(all(metadata_checks))
  for(id in space_ids) {
    out <- wb_space(b,id,cells=cells)
    if(is.data.frame(out)) out <- as.matrix(out[,c("x","y")])
    stopifnot(identical(as.numeric(out),as.numeric(spaces[[id]])))
  }
  features <- wb_feature_ids(bdev,"RNA_counts")
  stopifnot(identical(wb_feature_ids(b,"RNA_counts"),features))
  log_event("Verifying every imported RNA count in bounded feature slices")
  for(start in seq(1L,length(features),by=512L)) {
    genes <- features[start:min(start+511L,length(features))]
    old <- as.matrix(wb_matrix_read_features(bdev,"RNA_counts",genes,cells,missing="NA"))
    new <- as.matrix(wb_matrix_read_features(b,"RNA_counts",genes,cells,missing="NA"))
    stopifnot(identical(dim(old),dim(new)),identical(as.numeric(old),as.numeric(new)))
    if(start %% 4096L == 1L) log_event("RNA features checked through",min(start+511L,length(features)))
  }
  make_set <- function(id, selected, description) run_operation(b,paste0("set_",id),"cell_set_derive",list(
    target_cell_set=id,cell_ids=selected,min_cells=1L,role="analysis_universe",status="accepted",
    source="exact Bdev reference metadata",created_by="jtg",description=description,overwrite=FALSE),
    c("manifest","rule_table"),description)
  make_set("bdev_pt_defined",cells[is.finite(md$pt)],"Bdev cells with finite authoritative pt; no imputation")
  for(id in sort(unique(md$orig.ident))) make_set(paste0("bdev_sample_",id),cells[md$orig.ident==id],paste("Reference orig.ident",id))
  for(phase in sort(unique(md$Phase))) make_set(paste0("bdev_phase_",phase),cells[md$Phase==phase],paste("Reference cell-cycle Phase",phase))
  saveRDS(pt_family(b),file.path(report_dir,"exact_pt_family.rds"))
  write_json(list(status="complete",bundle_dir=bundle_dir,cells=length(cells),matrices="RNA_counts",spaces=space_ids,
    exact_metadata=as.list(metadata_checks),all_rna_values_exact=TRUE,all_coordinates_exact=TRUE,
    pt_exact=TRUE,pt_missingness_exact=TRUE,parents=parents),"creation.json")
  log_event("Derived bundle created with exact RNA counts, metadata and coordinates")
}
main()
