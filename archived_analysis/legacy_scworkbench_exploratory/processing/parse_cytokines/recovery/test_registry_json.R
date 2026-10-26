#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(scworkbench))
hash_rows <- getFromNamespace(".wb_manifest_rows_sha256", "scworkbench")
identity <- getFromNamespace(".wb_manifest_table_identity", "scworkbench")
sha_object <- getFromNamespace(".wb_prov_sha256_object", "scworkbench")
reference_hash <- function(rows) digest::digest(as.character(jsonlite::toJSON(
  rows, dataframe = "rows", auto_unbox = TRUE, null = "null", na = "null"
)), algo = "sha256", serialize = FALSE)
check <- function(rows) {
  expected <- reference_hash(rows)
  stopifnot(identical(hash_rows(rows), expected))
  for (size in c(1L, 2L, 7L, 1024L)) {
    stopifnot(identical(hash_rows(rows, stream_threshold_bytes = -1,
                                 string_chunk_chars = size), expected))
  }
}
fixture <- data.frame(
  id = 1:8,
  text = c('quotes" and backslash\\', "tabs\tand\nnewlines\r", "αβγ 🧬 café",
           "</script>/", "\b\f", NA_character_, "", "normal"),
  numeric = c(0, -1.23456789012, 1e20, 1e-20, NA_real_, Inf, -Inf, NaN),
  flag = c(TRUE, FALSE, NA, TRUE, FALSE, TRUE, NA, FALSE),
  stringsAsFactors = FALSE
)
names(fixture)[2L] <- 'column"\\\n'
check(fixture)
check(fixture[FALSE, , drop = FALSE])
check(data.frame(a = character(), b = numeric()))
check(data.frame(a = c(NA_character_, NA_character_), b = c(NA_real_, NA_real_)))
check(data.frame(a = factor(c("a", "b", NA)), b = as.Date(c("2026-01-01", NA, "2026-01-03"))))
set.seed(101)
random <- data.frame(
  id = seq_len(100),
  text = vapply(seq_len(100), function(i) paste(sample(fixture[[2]][!is.na(fixture[[2]])],
                                                   4, replace = TRUE), collapse = ""), character(1)),
  value = rnorm(100)
)
check(random)
con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
DBI::dbWriteTable(con, "fixture", random[rep(seq_len(100), 51), ])
legacy_identity <- function(con, table) {
  fields <- DBI::dbGetQuery(con, paste0("PRAGMA table_info(", DBI::dbQuoteIdentifier(con, table), ");"))
  header <- list(schema = "scworkbench.manifest.registry-table.v1", table = table,
                 columns = fields[, c("name", "type", "notnull", "pk"), drop = FALSE])
  state <- sha_object(header)
  rows <- DBI::dbGetQuery(con, paste0("SELECT * FROM ", table, " ORDER BY rowid;"))
  for (start in seq.int(1L, nrow(rows), by = 5000L)) {
    chunk <- rows[start:min(nrow(rows), start + 4999L), , drop = FALSE]
    state <- sha_object(list(previous = state, rows_sha256 = reference_hash(chunk)))
  }
  list(row_count = as.integer(nrow(rows)), content_sha256 = state)
}
stopifnot(identical(identity(con, "fixture"), legacy_identity(con, "fixture")))
before <- identity(con, "fixture")
DBI::dbExecute(con, "UPDATE fixture SET text = 'changed' WHERE rowid = 5001;")
stopifnot(!identical(before, identity(con, "fixture")))
DBI::dbDisconnect(con)
# Exercise the default streaming threshold on one large provenance-like field.
large <- data.frame(id = "large", provenance_json = strrep('{"key":"cell\\name"}', 4000000L))
stopifnot(identical(hash_rows(large), reference_hash(large)))
cat("PASS: streaming JSON matches original jsonlite hashes; chunk boundaries and tamper detection preserved.\n")
