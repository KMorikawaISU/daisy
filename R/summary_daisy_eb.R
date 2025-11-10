# summary_daisy_eb.R
# S3 summary() and print() for objects returned by EB_est() (class "daisy_eb").

#' Summarize EB_est results (daisy)
#'
#' Prints the point estimate and diagnostics. If `auto = TRUE`, it shows
#' the selected first-step model and a leaderboard with Entropy2, D1, and D2
#' for all candidates. If `auto = FALSE`, it shows the fixed model with its
#' Entropy2, D1, and D2.
#'
#' @param object An object returned by EB_est() with class "daisy_eb".
#' @param digits Number of digits for printing.
#' @param ... Unused.
#' @return An object of class "summary.daisy_eb".
#' @export
#' @method summary daisy_eb
summary.daisy_eb <- function(object, digits = 4, ...) {
  # 1) point estimate
  est <- tryCatch(unname(object$result[["est"]]), error = function(e) NA_real_)

  # 2) diagnostics table
  lb <- object$leaderboard
  if (is.null(lb) && !is.null(object$all_results)) {
    # rebuild from candidates (preferable: preserves per-candidate D2)
    lb <- do.call(
      rbind,
      lapply(object$all_results, function(z) {
        data.frame(
          divergence = z$model,
          r          = if (identical(z$model, "KL")) NA_real_ else z$r,
          Entropy2   = z$Entropy2,
          D1         = z$D1,
          D2         = z$D2,
          stringsAsFactors = FALSE
        )
      })
    )
  }
  if (is.null(lb)) {
    # fixed-model fallback (single row)
    lb <- data.frame(
      divergence = object$model,
      r          = if (identical(object$model, "KL")) NA_real_ else object$r,
      Entropy2   = object$Entropy2,
      D1         = object$D1,
      D2         = object$D2,
      stringsAsFactors = FALSE
    )
  } else {
    # normalize column name if needed
    if ("model" %in% names(lb) && !("divergence" %in% names(lb))) {
      names(lb)[names(lb) == "model"] <- "divergence"
    }
    # complete columns
    if (!"D1" %in% names(lb)) lb$D1 <- object$D1
    if (!"D2" %in% names(lb)) lb$D2 <- NA_real_
    keep <- intersect(c("divergence", "r", "Entropy2", "D1", "D2"), names(lb))
    lb   <- lb[, keep, drop = FALSE]
  }
  if ("Entropy2" %in% names(lb)) {
    lb <- lb[order(lb$Entropy2, decreasing = TRUE), , drop = FALSE]
  }

  # 3) selected model
  auto_flag <- isTRUE(object$auto)
  if (auto_flag && !is.null(object$best_model)) {
    sel <- object$best_model
    if (!is.null(sel$model) && is.null(sel$divergence)) sel$divergence <- sel$model
  } else {
    sel <- list(
      divergence = object$model,
      r          = if (identical(object$model, "KL")) NA_real_ else object$r
    )
  }

  out <- list(
    estimate       = est,
    auto           = auto_flag,
    selected_model = sel,
    leaderboard    = lb
  )
  class(out) <- "summary.daisy_eb"
  out
}

#' @export
print.summary.daisy_eb <- function(x, digits = 4, ...) {
  cat("Point estimate: ", format(x$estimate, digits = digits), "\n", sep = "")

  if (isTRUE(x$auto)) {
    cat("Selected model (auto): ", x$selected_model$divergence, sep = "")
  } else {
    cat("Model: ", x$selected_model$divergence, sep = "")
  }
  if (!is.null(x$selected_model$r) && !is.na(x$selected_model$r)) {
    cat(sprintf(" (r = %s)", format(signif(x$selected_model$r, digits = digits))))
  }
  cat("\n\n")

  cat("Diagnostics by model (higher Entropy2 is better):\n")
  lb <- x$leaderboard
  numcols <- vapply(lb, is.numeric, logical(1))
  lb[numcols] <- lapply(lb[numcols], function(v) signif(v, digits))
  print(lb, row.names = FALSE)
  invisible(x)
}
