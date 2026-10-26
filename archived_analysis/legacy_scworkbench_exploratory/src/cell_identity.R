# Reusable cell-identity checks for source inspection and integration preflight.

cell_id_diagnostics <- function(cell_ids, label = "cells") {
  ids <- as.character(cell_ids)
  list(
    label = label,
    n = length(ids),
    duplicated = sum(duplicated(ids)),
    missing_or_empty = sum(is.na(ids) | !nzchar(ids)),
    valid = length(ids) > 0L && !anyDuplicated(ids) && !any(is.na(ids) | !nzchar(ids))
  )
}

compare_ordered_cell_ids <- function(reference_ids, candidate_ids,
                                     reference_label = "reference",
                                     candidate_label = "candidate") {
  reference_ids <- as.character(reference_ids)
  candidate_ids <- as.character(candidate_ids)
  list(
    reference = cell_id_diagnostics(reference_ids, reference_label),
    candidate = cell_id_diagnostics(candidate_ids, candidate_label),
    same_set = setequal(reference_ids, candidate_ids),
    same_order = identical(reference_ids, candidate_ids),
    missing_from_candidate = setdiff(reference_ids, candidate_ids),
    extra_in_candidate = setdiff(candidate_ids, reference_ids)
  )
}

assert_disjoint_cell_id_collections <- function(cell_ids_by_bundle) {
  if (!is.list(cell_ids_by_bundle) || length(cell_ids_by_bundle) < 2L) {
    stop("cell_ids_by_bundle must be a named list with at least two bundles")
  }
  if (is.null(names(cell_ids_by_bundle)) || any(!nzchar(names(cell_ids_by_bundle)))) {
    stop("cell_ids_by_bundle must have non-empty bundle names")
  }
  for (name in names(cell_ids_by_bundle)) {
    check <- cell_id_diagnostics(cell_ids_by_bundle[[name]], name)
    if (!isTRUE(check$valid)) stop("Invalid cell IDs for bundle: ", name)
  }
  pairs <- utils::combn(names(cell_ids_by_bundle), 2L, simplify = FALSE)
  overlaps <- lapply(pairs, function(pair) {
    shared <- intersect(
      as.character(cell_ids_by_bundle[[pair[[1L]]]]),
      as.character(cell_ids_by_bundle[[pair[[2L]]]])
    )
    list(left = pair[[1L]], right = pair[[2L]], n_shared = length(shared),
         examples = head(shared, 10L))
  })
  if (any(vapply(overlaps, function(x) x$n_shared > 0L, logical(1)))) {
    stop("Cross-bundle cell-ID collisions detected")
  }
  invisible(overlaps)
}
