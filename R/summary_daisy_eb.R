# R/summary_daisy_eb.R
# S3 summary() and print() for objects returned by EB_est() (class "daisy_eb").
# This ensures the printed summary always shows:
#  - the point estimate,
#  - the selected model (if auto = TRUE),
#  - and a table with Entropy2, D1, D2 for all candidates (auto=TRUE) or for the fixed model.

#' Summarize EB_est results (daisy)
#'
#' Builds a compact summary object for printing. If `auto = TRUE`, it
#' includes the selected first-step model and a leaderboard with
#' Entropy2, D1, and D2 for all candidates. If `auto = FALSE`, it shows
#' the fixed model with those diagnostics.
#'
#' @param object An object returned by EB_est() (class "daisy_eb").
#' @param digits Number of digits to carry for numeric fields in the returned object.
#'   (The print method controls printed rounding.)
#' @param ... Unused.
#' @return An object of class "summary.daisy_eb" with elements:
#'   \itemize{
#'     \item \code{est}: numeric point estimate (mean of w2 * y).
#'     \item \code{selected}: list(divergence, r) for the chosen model.
#'     \item \code{table}: data.frame with (divergence, r, Entropy2, D1, D2).
#'     \item \code{auto}: logical flag (TRUE if auto search was used).
#'   }
#' @export
#' @method summary daisy_eb
summary.daisy_eb <- function(object, digits = getOption("digits"), ...) {
  # Extract point estimate from the EB_est result:
  est <- NA_real_
  if (!is.null(object$result) && !is.null(object$result["est"])) {
    est <- as.numeric(object$result["est"])
  }

  auto_flag <- isTRUE(object$auto)

  # Helper: build a leaderboard-like table from all_results
  from_all <- function(lst) {
    do.call(
      rbind,
      lapply(lst, function(z) {
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

  if (auto_flag) {
    # Prefer the per-candidate list if present
    if (!is.null(object$all_results)) {
      tbl <- from_all(object$all_results)
    } else if (!is.null(object$leaderboard)) {
      # Fallback to provided leaderboard (harmonize column names)
      tbl <- object$leaderboard
      if ("model" %in% names(tbl) && !("divergence" %in% names(tbl))) {
        names(tbl)[names(tbl) == "model"] <- "divergence"
      }
      if (!"D1" %in% names(tbl)) tbl$D1 <- object$D1
      if (!"D2" %in% names(tbl)) tbl$D2 <- NA_real_
      keep <- intersect(c("divergence","r","Entropy2","D1","D2"), names(tbl))
      tbl  <- tbl[, keep, drop = FALSE]
    } else {
      # Minimal fallback (should rarely happen)
      tbl <- data.frame(
        divergence = object$model,
        r          = if (identical(object$model, "KL")) NA_real_ else object$r,
        Entropy2   = object$Entropy2,
        D1         = object$D1,
        D2         = object$D2,
        stringsAsFactors = FALSE
      )
    }
    # Sort by Entropy2 (descending)
    if ("Entropy2" %in% names(tbl)) {
      tbl <- tbl[order(tbl$Entropy2, decreasing = TRUE), , drop = FALSE]
    }
    sel <- object$best_model
    if (!is.null(sel) && "model" %in% names(sel) && is.null(sel$divergence)) {
      sel$divergence <- sel$model
    }
  } else {
    # Fixed-model path
    tbl <- data.frame(
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

  # Keep only the intended columns (order as documented)
  wanted <- intersect(c("divergence","r","Entropy2","D1","D2"), names(tbl))
  tbl <- tbl[, wanted, drop = FALSE]

  out <- list(
    est      = est,
    selected = sel,
    table    = tbl,
    auto     = auto_flag
  )
  class(out) <- "summary.daisy_eb"
  out
}

#' Print method for summary.daisy_eb
#'
#' Always prints the point estimate and the model choice (auto or fixed),
#' followed by a table of Entropy2/D1/D2. Higher Entropy2 is better.
#'
#' @param x A summary.daisy_eb object (returned by summary()).
#' @param digits Number of digits for printing numeric values.
#' @param ... Unused.
#' @export
#' @method print summary.daisy_eb
print.summary.daisy_eb <- function(x, digits = 4, ...) {
  cat("daisy::EB_est summary\n")

  # Point estimate always shown:
  if (!is.null(x$est) && is.finite(x$est)) {
    cat(sprintf("  Point estimate: %.*f\n", digits, x$est))
  } else {
    cat("  Point estimate: NA\n")
  }

  # Selected model (auto or fixed):
  if (isTRUE(x$auto)) {
    cat("  Selected model (auto): ", x$selected$divergence, sep = "")
  } else {
    cat("  Model (fixed): ", x$selected$divergence, sep = "")
  }
  if (!is.null(x$selected$r) && !is.na(x$selected$r)) {
    cat(sprintf(" (r = %s)", format(signif(x$selected$r, digits = digits))))
  }
  cat("\n\n")

  # Leaderboard / diagnostics table:
  cat("Diagnostics by model (higher Entropy2 is better):\n")
  df <- x$table
  num_cols <- intersect(names(df), c("r","Entropy2","D1","D2"))
  for (cc in num_cols) df[[cc]] <- round(df[[cc]], digits)
  print(df, row.names = FALSE)
  invisible(x)
}
