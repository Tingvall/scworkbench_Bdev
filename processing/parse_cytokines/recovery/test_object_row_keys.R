#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(scworkbench))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  env <- new.env(parent = asNamespace("scworkbench"))
  sys.source(args[1], env)
  fast <- env$.wb_object_character_row_key_tokens
} else fast <- getFromNamespace(".wb_object_character_row_key_tokens", "scworkbench")
legacy <- function(keys) vapply(seq_len(nrow(keys)), function(i) {
  jsonlite::serializeJSON(as.list(keys[i, , drop = FALSE]), digits = NA)
}, character(1))
values <- c("", "plain", '"', "\\", '\",\"', '</script>/', "αβ 🧬 café", "NA", "null",
            vapply(1:31, intToUtf8, character(1)), 'ends\\', 'ends"', 'starts"', '\\",\"\\')
latin <- iconv("café", from = "UTF-8", to = "latin1")
values <- c(values, latin)
for (cols in list(list(cell_id = values), list(key = rev(values), cell_id = values),
                  setNames(list(values, rev(values)), c('name"\\', "α\n")))) {
  keys <- as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE)
  expected <- legacy(keys)
  for (chunk in c(1L, 2L, 7L, 10000L)) stopifnot(identical(fast(keys, chunk), expected))
}
set.seed(702)
keys <- data.frame(key = sample(values, 10000L, replace = TRUE),
                   cell_id = sample(values, 10000L, replace = TRUE))
old_time <- system.time(expected <- legacy(keys))
new_time <- system.time(observed <- fast(keys))
stopifnot(identical(expected, observed), identical(fast(keys[FALSE, ]), character()))
stopifnot(identical(fast(data.frame(cell_id = rep("", 5L))),
                    legacy(data.frame(cell_id = rep("", 5L)))))
if (!length(args)) {
  dispatch <- getFromNamespace(".wb_object_row_key_tokens", "scworkbench")
  for (keys in list(data.frame(cell_id = c("a", NA_character_)),
                    data.frame(key = c("a", "b"), point = c(1L, 2L)),
                    data.frame(key = factor(c("a", "b"))),
                    data.frame(cell_id = I(c("a", "b"))))) {
    stopifnot(identical(dispatch(keys), legacy(keys)))
  }
}
cat("PASS: character row tokens match the original serializer, including escapes, encodings, and empty strings.\n")
cat("Legacy elapsed:", old_time[["elapsed"]], "seconds; vectorized elapsed:", new_time[["elapsed"]], "seconds.\n")
