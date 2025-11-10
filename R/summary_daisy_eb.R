# R/summary_daisy_eb.R
# S3 summary() and print() for objects returned by EB_est() (class "daisy_eb").

#' Summarize EB_est results (daisy)
#'
#' Prints the point estimate and diagnostics. If `auto = TRUE`, it shows the
#' selected first-step model and a leaderboard with Entropy2, D1, and D2 for
#' all candidates. If `auto = FALSE`, it shows the fixed model with its
#' Entropy2, D1, and D2.
#'
#' @param object An object returned by EB_est() (class "daisy_eb").
#' @param digits Number of digits to print.
#' @param ... Unused.
#' @return An object of class "summary.daisy_eb".
#' @export
#' @method summary daisy_eb
summary.daisy_eb <- function(object, digits = getOption("digits"), ...) {
  est <- NA_real_
  if (!is.null(object$result) && !is.null(object$result["est"])) {
    est <- as.numeric(object$result["est"])
  }

  auto_flag <- isTRUE(object$auto)

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
    if (!is.null(object$all_results)) {
      tbl <- from_all(object$all_results)
    } else if (!is.null(object$leaderboard)) {
      tbl <- object$leaderboard
      if ("model" %in% names(tbl) && !("divergence" %in% names(tbl))) {
        names(tbl)[names(tbl) == "model"] <- "divergence"
      }
      if (!"D1" %in% names(tbl)) tbl$D1 <- object$D1
      if (!"D2" %in% names(tbl)) tbl$D2 <- NA_real_
      keep <- intersect(c("divergence","r","Entropy2","D1","D2"), names(tbl))
      tbl  <- tbl[, keep, drop = FALSE]
    } else {
      tbl <- data.frame(
        divergence = object$model,
        r          = if (identical(object$model, "KL")) NA_real_ else object$r,
        Entropy2   = object$Entropy2,
        D1         = object$D1,
        D2         = object$D2,
        stringsAsFactors = FALSE
      )
    }
    if ("Entropy2" %in% names(tbl)) {
      tbl <- tbl[order(tbl$Entropy2, decreasing = TRUE), , drop = FALSE]
    }
    sel <- object$best_model
    if (!is.null(sel) && "model" %in% names(sel) && is.null(sel$divergence)) {
      sel$divergence <- sel$model
    }
  } else {
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

#' @export
print.summary_daisy_eb <- function(x, digits = 4, ...) {
  cat("daisy::EB_est summary\n")
  if (!is.null(x$est) && is.finite(x$est)) {
    cat(sprintf("  Point estimate: %.*f\n", digits, x$est))
  } else {
    cat("  Point estimate: NA\n")
  }

  if (isTRUE(x$auto)) {
    cat("  Selected model (auto): ", x$selected$divergence, sep = "")
  } else {
    cat("  Model (fixed): ", x$selected$divergence, sep = "")
  }
  if (!is.null(x$selected$r) && !is.na(x$selected$r)) {
    cat(sprintf(" (r = %s)", format(signif(x$selected$r, digits = digits))))
  }
  cat("\n\n")

  cat("Diagnostics by model (higher Entropy2 is better):\n")
  df <- x$table
  num_cols <- intersect(names(df), c("r","Entropy2","D1","D2"))
  for (cc in num_cols) df[[cc]] <- round(df[[cc]], digits)
  print(df, row.names = FALSE)
  invisible(x)
}
