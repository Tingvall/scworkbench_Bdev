# Native scWorkbench distribution summaries for pooled pseudotime bins.
# The grouping repeats the engine's equal-width bin assignment; it is verified
# against the retained native trajectory counts and medians before use.
pooled_signature_summary <- function(scores, pt, retained_curve, n_bins=30L, min_cells=30L) {
  stopifnot(!anyDuplicated(pt$cell_id),all(is.finite(pt$value)),length(pt$value)==12642L)
  bins <- as.integer(cut(pt$value,breaks=n_bins,include.lowest=TRUE,labels=FALSE))
  stopifnot(all(is.finite(bins)),all(bins>=1L & bins<=n_bins))
  grouping <- wb_prim_cell_vector(setNames(sprintf("bin_%02d",bins),pt$cell_id),value_type="categorical",
    meta=list(source_field="Bdev@meta.data$pt",population="pooled",n_bins=n_bins,
      binning="cut(pt, breaks=n_bins, include.lowest=TRUE); native equal-width trajectory convention"))
  take <- scores$cell_id %in% pt$cell_id
  signal <- wb_prim_cell_vector_family(scores$key[take],scores$cell_id[take],scores$value[take],
    meta=list(source="retained Bdev UCell scores",population="all cells with defined source pt"))
  stopifnot(length(signal$value)==22L*length(pt$value))
  for(key in unique(signal$key)) stopifnot(identical(signal$cell_id[signal$key==key],pt$cell_id))
  # Exported native operations are pure computations on primitives. No session
  # promotion, artifact registration or bundle modification occurs here.
  native <- wb_op_cell_vector_family_qc(signal=signal,grouping=grouping,correlation_methods="spearman")
  native_pt <- wb_op_cell_vector_family_qc(signal=pt,grouping=grouping,correlation_methods="spearman")
  spread <- native$outputs$group_summary$table
  pt_summary <- native_pt$outputs$group_summary$table
  spread$bin <- as.integer(sub("bin_","",spread$group_id))
  spread$feature_set <- spread$key
  ix <- match(spread$group_id,pt_summary$group_id)
  stopifnot(!anyNA(ix),all(spread$n_finite==pt_summary$n_finite[ix]))
  spread$pseudotime <- pt_summary$score_median[ix]
  spread$pseudotime_min <- pt_summary$score_min[ix]
  spread$pseudotime_max <- pt_summary$score_max[ix]
  spread$iqr <- spread$score_q75-spread$score_q25
  spread$eligible_bin <- spread$n_finite>=min_cells
  spread <- spread[order(spread$feature_set,spread$bin),]
  primary <- spread[spread$eligible_bin,]
  old_key <- paste(retained_curve$feature_set,retained_curve$bin)
  new_key <- paste(primary$feature_set,primary$bin)
  stopifnot(!anyDuplicated(new_key),setequal(old_key,new_key))
  old <- retained_curve[match(new_key,old_key),]
  # Earlier promoted tables round numeric values to four decimal places.
  stopifnot(identical(as.integer(primary$n_finite),as.integer(old$n_cells)),
    all(abs(primary$score_median-old$value)<=0.00005001),
    all(abs(primary$pseudotime-old$pseudotime)<=0.00005001),
    all(spread$score_q25<=spread$score_median),all(spread$score_median<=spread$score_q75),
    all(spread$iqr>=0),all(is.finite(primary$score_sd)),all(primary$score_sd>=0),
    all(primary$score_q25>=0 & primary$score_q75<=1))
  counts <- tabulate(bins,nbins=n_bins)
  limits <- range(pt$value)
  edges <- seq(limits[[1]],limits[[2]],length.out=n_bins+1L)
  support <- data.frame(bin=seq_len(n_bins),pt_left=head(edges,-1L),pt_right=tail(edges,-1L),
    n_cells=counts,eligible_bin=counts>=min_cells)
  stopifnot(sum(counts)==length(pt$value),sum(counts[counts>=min_cells])==sum(primary$n_finite[primary$key==primary$key[[1]]]))
  list(native=native,native_pt=native_pt,grouping=grouping,spread=spread,primary=primary,support=support,
    checks=list(pooled_cell_ids_exact=TRUE,all_defined_pt_cells_accounted_for=TRUE,
      native_bin_counts_match=TRUE,native_medians_match_within_recorded_rounding=TRUE,
      quartiles_ordered_and_bounded=TRUE,standard_deviations_finite=TRUE),
    parameters=list(population="pooled",n_bins=n_bins,min_cells_per_bin=min_cells,
      center="median",band="25th to 75th percentiles (cell-level IQR)",quantile_type=7L,
      sd="sample SD of cell scores within each bin",line="connect all retained bin centers",
      smoothing="none",bootstrap_n=0L,inference="descriptive cell-to-cell spread, not a confidence interval"))
}
