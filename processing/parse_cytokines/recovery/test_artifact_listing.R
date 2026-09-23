#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(scworkbench))
# Disposable SQL fixtures compare the public getter against its original query.
con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
bundle <- list(con = con)
rows <- data.frame(
  artifact_id = paste0("a", 1:9),
  type = rep(c("obs_table", "cell_vector"), length.out = 9),
  created_at = c("2026-09-23", "2026-09-22", "2026-09-23", NA, "", "z", "é", "a", "a"),
  provenance_json = c(rep("{}", 8), paste(rep("x", 1000000), collapse = "")),
  stringsAsFactors = FALSE
)
DBI::dbWriteTable(con, "artifacts", rows)
for (lifecycle in c(FALSE, TRUE)) {
  if (lifecycle) DBI::dbWriteTable(con, "entity_lifecycle", data.frame(
    entity_kind = "artifact", entity_id = rows$artifact_id,
    lifecycle_state = rep(c("current", "retracted", "superseded"), 3),
    provenance_bound = rep(0:1, length.out = 9)
  ))
  for (history in c(FALSE, TRUE)) for (type in c(NA_character_, "obs_table", "absent")) {
    select <- if (lifecycle) paste0(
      "SELECT a.*, COALESCE(el.lifecycle_state, 'current') AS lifecycle_state, ",
      "COALESCE(el.provenance_bound, 0) AS provenance_bound FROM artifacts a ",
      "LEFT JOIN entity_lifecycle el ON el.entity_kind = 'artifact' AND el.entity_id = a.artifact_id"
    ) else "SELECT a.*, 'current' AS lifecycle_state, 0 AS provenance_bound FROM artifacts a"
    conditions <- character()
    if (!is.na(type)) conditions <- c(conditions, paste0("a.type = '", type, "'"))
    if (lifecycle && !history) conditions <- c(conditions, "COALESCE(el.lifecycle_state, 'current') = 'current'")
    query <- paste0(select, if (length(conditions)) paste0(" WHERE ", paste(conditions, collapse = " AND ")) else "",
                    " ORDER BY a.created_at DESC;")
    expected <- DBI::dbGetQuery(con, query)
    observed <- wb_art_list(bundle, type = if (is.na(type)) NULL else type, include_history = history)
    stopifnot(identical(observed, expected))
  }
}
DBI::dbDisconnect(con)
cat("PASS: artifact-list ordering, ties, missing values, filtering, lifecycle, and payload retention\n")
