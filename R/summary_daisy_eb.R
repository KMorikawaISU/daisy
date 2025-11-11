# R/summary_daisy_eb.R
# S3 summary() and print() for objects returned by EB_est() (class "daisy_eb").
# This version uses the name "estimate" consistently for:
#   - the top-level point estimate in the summary object, and
#   - the per-model point estimate column in the leaderboard table.

#' Summarize EB_est results (daisy)
#'
#' Creates a summary object containing:
#'   - the overall point estimate (for the returned/selected run),
#'   - the selected model information (if auto = TRUE),
#'   - and a leaderboard-like table that includes the per-model point
#'     estimate `estimate` as the leftmost column.
#'
#' If `auto = TRUE`, the table includes *all* candidate models; otherwise, it
#' contains only the fixed model. Higher Entropy2 is considered better.
#'
#' @param object An object returned by EB_est() (class "daisy_eb").
#' @param digits Number of digits to carry for numeric fields in the returned object.
#'   (The print method controls on-screen rounding.)
#' @param ... Unused.
#' @return An object of class "summary.daisy_eb" with elements:
#'   \itemize{
#'     \item \code{estimate} numeric point estimate for the *selected* (or only) run.
#'     \item \code{selected} list(divergence, r) describing the selected model
#'           when `auto = TRUE`, or the fixed model when `auto = FALSE`.
#'     \item \code{table} data.frame with columns
#'           \code{estimate}, \code{divergence}, \code{r}, \code{Entropy2}, \code{D1}, \code{D2}.
#'     \item \code{auto} logical indicating whether auto model search was used.
#'   }
#' @export
#' @method summary daisy_eb
summary.daisy_eb <- function(object, digits = getOption("digits"), ...) {
  # Extract the point estimate for the returned object (selected run or fixed run)
  estimate <- NA_real_
  if (!is.null(object$result) && !is.null(object$result["est"])) {
    # Note: EB_est_one() stores the scalar under name "est"; convert to "estimate" here.
    estimate <- as.numeric(object$result["est"])
  }

  auto_flag <- isTRUE(object$auto)

  # Build a table from a list of candidate results (all_results), including per-model estimate
  from_all_with_estimate <- function(lst) {
    do.call(
      rbind,
      lapply(lst, function(z) {
        # Per-candidate point estimate if available (stored as "est" in EB_est_one())
        z_est <- NA_real_
        if (!is.null(z$result) && !is.null(z$result["est"])) {
          z_est <- as.numeric(z$result["est"])
        }
        data.frame(
          estimate   = z_est,
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

  if (auto_flag) {
    if (!is.null(object$all_results)) {
      # Preferred path: reconstruct full table from all_results
      tbl <- from_all_with_estimate(object$all_results)
    } else if (!is.null(object$leaderboard)) {
      # Fallback: leaderboard may be missing per-model estimates
      tbl <- object$leaderboard
      if ("model" %in% names(tbl) && !("divergence" %in% names(tbl))) {
        names(tbl)[names(tbl) == "model"] <- "divergence"
      }
      if (!("estimate" %in% names(tbl))) {
        tbl$estimate <- NA_real_
      }
      if (!("D1" %in% names(tbl))) tbl$D1 <- object$D1
      if (!("D2" %in% names(tbl))) tbl$D2 <- object$D2
      have <- names(tbl)
      tbl  <- tbl[, intersect(c("estimate","divergence","r","Entropy2","D1","D2"), have), drop = FALSE]
    } else {
      # Minimal fallback (rare)
      tbl <- data.frame(
        estimate   = estimate,
        divergence = object$model,
        r          = if (identical(object$model, "KL")) NA_real_ else object$r,
        Entropy2   = object$Entropy2,
        D1         = object$D1,
        D2         = object$D2,
        stringsAsFactors = FALSE
      )
    }

    # Sort by Entropy2 (descending) when available
    if ("Entropy2" %in% names(tbl)) {
      tbl <- tbl[order(tbl$Entropy2, decreasing = TRUE), , drop = FALSE]
    }

    # Selected model info (harmonize field name)
    sel <- object$best_model
    if (!is.null(sel) && "model" %in% names(sel) && is.null(sel$divergence)) {
      sel$divergence <- sel$model
    }

  } else {
    # Fixed-model path: single row with estimate included
    tbl <- data.frame(
      estimate   = estimate,
      divergence = object$model,
      r          = if (identical(object$model, "KL")) NA_real_ else object$r,
      Entropy2   = object$Entropy2,
      D1         = object$D1,
      D2         = object$D2,
      stringsAsFactors = FALSE
    )
    sel <- list(
      divergence = object$model,
      r          = if (identical(object$model, "KL")) NA_real_ else object$r
    )
  }

  # Keep only the intended columns (in the required order)
  wanted <- intersect(c("estimate","divergence","r","Entropy2","D1","D2"), names(tbl))
  tbl <- tbl[, wanted, drop = FALSE]

  out <- list(
    estimate = estimate,  # selected run's point estimate
    selected = sel,       # selected or fixed model
    table    = tbl,       # table with per-model estimate at the leftmost column
    auto     = auto_flag
  )
  attr(out, "digits") <- digits
  class(out) <- "summary.daisy_eb"
  out
}

#' Print method for summary.daisy_eb
#'
#' Nicely prints:
#'   - the point estimate for the selected (or fixed) run,
#'   - the selected model (auto) or the fixed model,
#'   - a table with per-model point estimates (estimate) and diagnostics.
#'
#' @param x A summary.daisy_eb object (returned by summary()).
#' @param digits Number of digits for printing numeric values.
#' @param ... Unused.
#' @export
#' @method print summary.daisy_eb
print.summary.daisy_eb <- function(x, digits = getOption("digits"), ...) {
  # Allow digits passed via summary()
  d_attr <- attr(x, "digits", exact = TRUE)
  if (is.numeric(d_attr) && length(d_attr) == 1L) digits <- d_attr

  cat("daisy::EB_est summary\n")

  # Selected (or fixed) run's point estimate
  if (!is.null(x$estimate) && is.finite(x$estimate)) {
    cat(sprintf("  Point estimate (selected): %.*f\n", max(0, digits - 2L), x$estimate))
  } else {
    cat("  Point estimate (selected): NA\n")
  }

  # Model header
  if (isTRUE(x$auto)) {
    cat("  Selected model (auto): ", x$selected$divergence, sep = "")
  } else {
    cat("  Model (fixed): ", x$selected$divergence, sep = "")
  }
  if (!is.null(x$selected$r) && !is.na(x$selected$r)) {
    cat(sprintf(" (r = %s)", format(signif(x$selected$r, digits = max(3, digits - 2L)))))
  }
  cat("\n\n")

  # Leaderboard / diagnostics table with per-model estimate
  cat("Diagnostics by model (higher Entropy2 is better):\n")
  df <- x$table
  num_cols <- intersect(names(df), c("estimate","r","Entropy2","D1","D2"))
  for (cc in num_cols) df[[cc]] <- round(df[[cc]], max(0, digits - 2L))
  print(df, row.names = FALSE)

  invisible(x)
}
