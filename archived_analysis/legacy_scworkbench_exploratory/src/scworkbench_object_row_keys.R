# Internal engine fast path for plain character row-key columns. The original
# serializer remains the fallback for missing values, attributed columns, and
# non-character keys. Return exactly the legacy serializeJSON row tokens.
.wb_object_character_row_key_tokens <- function(keys, chunk_size = 10000L) {
  n <- nrow(keys)
  if (!n) return(character())
  header <- paste0('{"type":"list","attributes":{"names":',
                   jsonlite::serializeJSON(names(keys), digits = NA), '},"value":[')
  out <- character(n)
  for (start in seq.int(1L, n, by = chunk_size)) {
    index <- start:min(n, start + chunk_size - 1L)
    rows <- rep(header, length(index))
    for (j in seq_along(keys)) {
      # A trailing empty sentinel compensates for strsplit dropping its final
      # empty field. Escaped string contents cannot contain the unescaped
      # array-element delimiter ","; literal quotes always have a backslash.
      encoded <- as.character(jsonlite::toJSON(c(keys[[j]][index], ""),
                                               auto_unbox = FALSE, na = "null"))
      middle <- substr(encoded, 3L, nchar(encoded, type = "chars") - 2L)
      values <- strsplit(middle, '\",\"', fixed = TRUE)[[1L]]
      stopifnot(length(values) == length(index))
      rows <- paste0(rows, if (j > 1L) "," else "",
                     '{"type":"character","attributes":{},"value":["',
                     values, '"]}')
    }
    out[index] <- paste0(rows, "]}")
  }
  out
}
