#!/usr/bin/env Rscript

# Launch the reviewed implementation as a top-level Rscript process. This is
# required because its bundle cleanup handlers intentionally live for the full
# process lifetime and must not be evaluated through source().
args <- commandArgs(trailingOnly = TRUE)
status <- system2(file.path(R.home("bin"), "Rscript"),
  c("--vanilla",
    "processing/parse_cytokines/17_flt3_il7_competitor_refinement_pilot.R",
    args), env = "SCWB_REFINEMENT_SCOPE=all_cytokines")
quit(save = "no", status = status, runLast = FALSE)
