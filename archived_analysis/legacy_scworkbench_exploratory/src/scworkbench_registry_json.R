# Installed into a separate scWorkbench recovery build by the stage-10 runner.
# This is engine code, not a client-side validation override. Keep the existing
# registry-table identity scheme and exact jsonlite bytes, including 5,000-row
# hash boundaries, while avoiding a single oversized character string.
.wb_manifest_rows_sha256 <- function(rows, stream_threshold_bytes = 64 * 1024^2,
                                     string_chunk_chars = 32 * 1024^2) {
  rows <- .wb_json_sanitize(rows)
  text_bytes <- sum(vapply(rows, function(column) {
    if (is.character(column)) sum(nchar(column, type = "bytes"), na.rm = TRUE) else 0
  }, numeric(1)))
  encode <- function(x) as.character(jsonlite::toJSON(
    x, dataframe = "rows", auto_unbox = TRUE, null = "null", na = "null"
  ))
  if (text_bytes <= stream_threshold_bytes) {
    return(digest::digest(encode(rows), algo = "sha256", serialize = FALSE))
  }
  path <- tempfile("scworkbench-registry-json-", fileext = ".json")
  on.exit(unlink(path), add = TRUE)
  output <- file(path, open = "wb")
  closed <- FALSE
  on.exit(if (!closed) close(output), add = TRUE)
  emit <- function(text) writeChar(enc2utf8(text), output, eos = NULL, useBytes = TRUE)
  emit("[")
  for (i in seq_len(nrow(rows))) {
    if (i > 1L) emit(",")
    emit("{")
    for (j in seq_along(rows)) {
      if (j > 1L) emit(",")
      value <- rows[[j]][i]
      if (is.character(value) && length(value) == 1L && !is.na(value) &&
          nchar(value, type = "bytes") > string_chunk_chars) {
        emit(encode(names(rows)[j]))
        emit(":\"")
        nchars <- nchar(value, type = "chars")
        for (start in seq.int(1, nchars, by = string_chunk_chars)) {
          fragment <- encode(substr(value, start, min(nchars, start + string_chunk_chars - 1)))
          # jsonlite escapes / only when preceded by <. Preserve that context
          # when the two characters lie on opposite sides of a chunk boundary.
          if (start > 1 && identical(substr(value, start - 1, start), "</")) emit("\\")
          emit(substr(fragment, 2L, nchar(fragment, type = "chars") - 1L))
        }
        emit("\"")
      } else {
        # Delegate scalar types, missing values, and column-name escaping to
        # jsonlite, removing only the outer [{ and }] of this one-column row.
        fragment <- encode(rows[i, j, drop = FALSE])
        emit(substr(fragment, 3L, nchar(fragment, type = "chars") - 2L))
      }
    }
    emit("}")
  }
  emit("]")
  close(output)
  closed <- TRUE
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}
