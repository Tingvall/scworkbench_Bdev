# Display existing native window summaries; no statistics are recomputed.
sliding_window_iqr_overview <- function(summaries,support,coverage,fraction,parameters) {
  d <- summaries[summaries$width_fraction==fraction & summaries$eligible_window,]
  windows <- support[support$width_fraction==fraction,]
  stopifnot(nrow(windows)>0,setequal(unique(d$feature_set),coverage$feature_set[coverage$eligible]),
    !anyDuplicated(d[c("window_id","feature_set")]),all(is.finite(d$score_median)),
    all(d$score_q25<=d$score_median & d$score_median<=d$score_q75))
  d$direction <- factor(coverage$direction[match(d$feature_set,coverage$feature_set)],levels=c("up","down"))
  low <- windows[!windows$eligible_window,]
  limits <- as.numeric(parameters$source_pt_range)
  half_step <- parameters$step_pt_units/2
  low$left <- pmax(limits[[1]],low$pseudotime-half_step)
  low$right <- pmin(limits[[2]],low$pseudotime+half_step)
  palette <- c(up="#B2182B",down="#2166AC")
  p <- ggplot(d,aes(pseudotime,score_median,group=feature_set,color=direction,fill=direction))+
    geom_rect(data=low,aes(xmin=left,xmax=right,ymin=-Inf,ymax=Inf),inherit.aes=FALSE,fill="#64748B",alpha=.13)+
    geom_ribbon(aes(ymin=score_q25,ymax=score_q75),alpha=.20,colour=NA)+
    geom_line(linewidth=.85)+scale_color_manual(values=palette)+scale_fill_manual(values=palette)+
    facet_wrap(~feature_set,scales="free_y",ncol=4,labeller=labeller(feature_set=function(x) sub("^parse_cytokine[.]","",x)))+
    labs(title=sprintf("Pooled cytokine signatures: %.0f%% sliding-window median and IQR",100*fraction),
      subtitle=sprintf("All %d scored gene sets; %d of %d window centers retained",length(unique(d$feature_set)),sum(windows$eligible_window),nrow(windows)),
      x="Existing Bdev pseudotime (pt)",y="UCell score",color="Parse direction",fill="Parse direction",
      caption=sprintf("Window: %.0f%% of pt range (%.3f units); step: %.0f%% (%.3f units). Line: median; band: 25th-75th percentiles.\nGray marks windows below %d cells. Connecting segments across these regions are for display; endpoint windows are clipped.",
        100*fraction,fraction*diff(limits),100*parameters$step_fraction,parameters$step_pt_units,parameters$min_cells_per_window))+
    theme(legend.position="bottom")
  line <- ggplot_build(p)$data[[3]]
  stopifnot(length(unique(line$group))==length(unique(d$feature_set)),nrow(line)==nrow(d),
    all(is.finite(line$x)),all(is.finite(line$y)))
  p
}
