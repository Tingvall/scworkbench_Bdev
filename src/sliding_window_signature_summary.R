# Overlapping subsets are selected on the unchanged source pt coordinate.
# scWorkbench computes every median, quartile and SD from the selected cells.
sliding_window_signature_summary <- function(scores,pt,width_fractions,step_fraction,min_cells=30L) {
  stopifnot(!anyDuplicated(pt$cell_id),all(is.finite(pt$value)),
    all(width_fractions>0 & width_fractions<=1),step_fraction>0,
    step_fraction<min(width_fractions),min_cells>=2L)
  keys <- unique(scores$key)
  values <- vapply(keys,function(key) {
    use <- scores$key==key
    stopifnot(!anyDuplicated(scores$cell_id[use]))
    ix <- match(pt$cell_id,scores$cell_id[use]); stopifnot(!anyNA(ix))
    scores$value[use][ix]
  },numeric(length(pt$value)))
  rownames(values) <- pt$cell_id
  stopifnot(all(is.finite(values)),all(values>=0 & values<=1))
  limits <- range(pt$value); span <- diff(limits)
  stopifnot(span>0)
  grid_fractions <- seq(0,1,by=step_fraction)
  if(tail(grid_fractions,1)<1) grid_fractions <- c(grid_fractions,1)
  centers <- limits[[1]]+grid_fractions*span
  centers[[1]] <- limits[[1]]; centers[[length(centers)]] <- limits[[2]]
  tables <- list(); support <- list(); membership <- list(); native_summaries <- list(); cursor <- 0L
  for(fraction in width_fractions) {
    width <- fraction*span
    covered <- rep(FALSE,length(pt$value))
    for(i in seq_along(centers)) {
      left <- max(limits[[1]],centers[[i]]-width/2)
      right <- min(limits[[2]],centers[[i]]+width/2)
      ix <- which(pt$value>=left & pt$value<=right)
      id <- sprintf("width_%03d_window_%03d",round(fraction*1000),i)
      covered[ix] <- TRUE
      membership[[id]] <- pt$cell_id[ix]
      stopifnot(!anyDuplicated(membership[[id]]))
      cursor <- cursor+1L
      support[[cursor]] <- data.frame(window_id=id,width_fraction=fraction,window_index=i,
        pseudotime=centers[[i]],window_left=left,window_right=right,
        nominal_width=width,actual_width=right-left,edge_truncated=(right-left)<width-1e-12,
        n_cells=length(ix),eligible_window=length(ix)>=min_cells)
      if(length(ix)) {
        signal <- wb_prim_cell_vector_family(rep(keys,each=length(ix)),rep(pt$cell_id[ix],times=length(keys)),
          as.numeric(values[ix,,drop=FALSE]),meta=list(window_id=id,source_field="Bdev@meta.data$pt",population="pooled"))
        native <- wb_op_cell_vector_family_qc(signal,correlation_methods="pearson")
        native_summaries[[id]] <- native$outputs$key_summary
        tab <- native$outputs$key_summary$table
        stopifnot(setequal(tab$key,keys),all(tab$n_finite==length(ix)))
        tab <- tab[match(keys,tab$key),c("key","score_min","score_q25","score_median","score_mean","score_q75","score_max","score_sd")]
        names(tab)[[1]] <- "feature_set"
      } else {
        tab <- data.frame(feature_set=keys,score_min=NA_real_,score_q25=NA_real_,score_median=NA_real_,
          score_mean=NA_real_,score_q75=NA_real_,score_max=NA_real_,score_sd=NA_real_)
      }
      tables[[cursor]] <- cbind(support[[cursor]][rep(1L,length(keys)),],tab)
      if(i %% 25L==0L) log_event("Sliding windows:",paste0(100*fraction,"%"),i,"of",length(centers))
    }
    stopifnot(all(covered))
  }
  summaries <- do.call(rbind,tables); windows <- do.call(rbind,support)
  rownames(summaries) <- rownames(windows) <- NULL
  summaries$iqr <- summaries$score_q75-summaries$score_q25
  retained <- summaries[summaries$eligible_window,]
  stopifnot(!anyDuplicated(summaries[c("window_id","feature_set")]),
    all(is.finite(retained$score_median)),all(is.finite(retained$score_sd)),
    all(retained$score_q25<=retained$score_median & retained$score_median<=retained$score_q75),
    all(retained$score_q25>=0 & retained$score_q75<=1),all(retained$iqr>=0),
    all(vapply(seq_len(nrow(windows)),function(i) length(membership[[windows$window_id[[i]]]])==windows$n_cells[[i]],logical(1))))
  list(summaries=summaries,windows=windows,membership=membership,native_summaries=native_summaries,
    parameters=list(population="pooled",source_pt_range=limits,width_fractions=width_fractions,
      widths_pt_units=width_fractions*span,step_fraction=step_fraction,step_pt_units=step_fraction*span,
      centers_per_width=length(centers),min_cells_per_window=min_cells,selection="left <= original pt <= right",
      edge_policy="clip windows to the observed pt range",center="median",band="25th to 75th percentiles (IQR)",
      quantile_type=7L,extra_smoothing="none",bootstrap_n=0L,random_sampling=FALSE,
      native_operation="wb_op_cell_vector_family_qc; only key distribution summaries are used",
      interpretation="cell-to-cell spread; overlapping windows are not independent; band is not a confidence interval"),
    checks=list(cell_score_alignment_exact=TRUE,window_members_unique=TRUE,window_counts_verified=TRUE,
      every_defined_pt_cell_covered_for_each_width=TRUE,quartiles_ordered_and_bounded=TRUE,
      retained_values_finite=TRUE,unique_signature_window_rows=TRUE))
}
