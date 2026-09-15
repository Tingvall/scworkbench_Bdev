# Validate a registered accepted bundle through scWorkbench public interfaces.

validate_registered_bundle <- function(registry, bundle_id,
                                       run_recipient_validation = TRUE,
                                       progress = TRUE) {
  if (!requireNamespace("scworkbench", quietly = TRUE)) {
    stop("Run validate_registered_bundle() inside the supplied scWorkbench image")
  }
  entry <- bundle_registry_entry(registry, bundle_id)
  if (!identical(entry$status, "accepted")) {
    stop("Bundle is not registered as accepted: ", bundle_id)
  }
  if (!dir.exists(entry$bundle_dir)) stop("Bundle directory does not exist: ", entry$bundle_dir)

  bundle <- scworkbench::wb_open_bundle(entry$bundle_dir, read_only = TRUE, validate = FALSE)
  manifest <- tryCatch(
    scworkbench::wb_manifest(bundle),
    finally = scworkbench::wb_close_bundle(bundle)
  )
  expected <- entry$accepted_state$bundle_state_fingerprint
  observed <- manifest$current_state$fingerprint
  identity_ok <- identical(manifest$bundle_id, bundle_id)
  fingerprint_ok <- identical(observed, expected)
  if (!identity_ok) stop("Bundle ID does not match registry entry: ", bundle_id)
  if (!fingerprint_ok) stop("Bundle state fingerprint does not match accepted registry state: ", bundle_id)

  recipient <- NULL
  if (isTRUE(run_recipient_validation)) {
    recipient <- scworkbench::wb_bundle_validation_report(
      entry$bundle_dir, level = "recipient", progress = progress
    )
    if (!isTRUE(recipient$ok)) stop("Recipient validation failed for: ", bundle_id)
  }

  list(
    bundle_id = bundle_id,
    bundle_dir = entry$bundle_dir,
    identity_ok = identity_ok,
    bundle_state_fingerprint = observed,
    fingerprint_ok = fingerprint_ok,
    recipient_validation = recipient
  )
}
