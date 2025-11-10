# R/summary_daisy_eb.R
# S3 summary() and print() for objects returned by EB_est() (class "daisy_eb").

#' Summarize EB_est results (daisy)
#'
#' Prints the point estimate and diagnostics. If `auto = TRUE`, it shows
#' the selected first-step model and a leaderboard with Entropy2, D1, and D2
#' for all candidates. If `auto = FALSE`, it shows the fixed model with its
#' Entropy2, D1, and D2.
#'
#' @param object An object returned by EB_est() (class "daisy_eb").
#' @param digits Integer; number of digits to print.
#' @param ... Unused.
#' @return An object of class "summary.daisy_eb" that prints nicely.
#' @export
#' @method summary daisy_eb
summary.daisy_eb <- function(object, digits = getOption("digits"), ...) {
  # point estimate
  est <- NA_real_
  if (!is.null(object$result) && !is.null(object$result["est"])) {
    est <- as.numeric(object$result["est"])
  }

  auto_flag <- isTRUE(object$auto)

  # helper: rebuild diagnostics table from all_results when available
  from_all_results <- function(lst) {
    do.call(
      rbind,
      lapply(lst, function(z) {
        data.frame(
          divergence = if (!is.null(z$model)) z$model else NA_character_,
          r          = if (!is.null(z$r))     z$r     else NA_real_,
          Entropy2   = if (!is.null(z$Entropy2)) z$Entropy2 else NA_real_,
          D1         = if (!is.null(z$D1))    z$D1    else if (!is.null(object$D1)) object$D1 else NA_real_,
          D2         = if (!is.null(z$D2))    z$D2    else NA_real_,
          check.names = FALSE
        )
      })
    )
  }

  if (auto_flag) {
    if (!is.null(object$all_results)) {
      tbl <- from_all_results(object$all_results)
    } else if (!is.null(object$leaderboard)) {
      tbl <- object$leaderboard
      if ("model" %in% names(tbl) && !("divergence" %in% names(tbl))) {
        names(tbl)[names(tbl) == "model"] <- "divergence"
      }
      if (!"D1" %in% names(tbl)) tbl$D1 <- if (!is.null(object$D1)) object$D1 else NA_real_
      if (!"D2" %in% names(tbl)) tbl$D2 <- NA_real_
      keep <- intersect(c("divergence","r","Entropy2","D1","D2"), names(tbl))
      tbl  <- tbl[, keep, drop = FALSE]
    } else {
      tbl <- data.frame(
        divergence = if (!is.null(object$model)) object$model else NA_character_,
        r          = if (!is.null(object$r))     object$r     else NA_real_,
        Entropy2   = if (!is.null(object$Entropy2)) object$Entropy2 else NA_real_,
        D1         = if (!is.null(object$D1)) object$D1 else NA_real_,
        D2         = if (!is.null(object$D2)) object$D2 else NA_real_,
        check.names = FALSE
      )
    }

    if ("Entropy2" %in% names(tbl)) {
      tbl <- tbl[order(tbl$Entropy2, decreasing = TRUE), , drop = FALSE]
    }

    sel <- object$best_model
    if (!is.null(sel) && "model" %in% names(sel) && !("divergence" %in% names(sel))) {
      sel$divergence <- sel$model
    }

  } else {
    tbl <- data.frame(
      divergence = if (!is.null(object$model)) object$model else NA_character_,
      r          = if (!is.null(object$r))     object$r     else NA_real_,
      Entropy2   = if (!is.null(object$Entropy2)) object$Entropy2 else NA_real_,
      D1         = if (!is.null(object$D1)) object$D1 else NA_real_,
      D2         = if (!is.null(object$D2)) object$D2 else NA_real_,
      check.names = FALSE
    )
    sel <- list(
      divergence = if (!is.null(object$model)) object$model else NA_character_,
      r          = if (!is.null(object$r))     object$r     else NA_real_
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
print.summary.daisy_eb <- function(x, digits = 4, ...) {
  cat("daisy::EB_est summary\n")
  if (!is.null(x$est) && is.finite(x$est)) {
    cat(sprintf("  Point estimate: %.*f\n", digits, x$est))
  } else {
    cat("  Point estimate: NA\n")
  }

  if (isTRUE(x$auto)) {
    div <- if (!is.null(x$selected$divergence)) x$selected$divergence else
      if (!is.null(x$selected$model))       x$selected$model      else "NA"
    rr  <- if (!is.null(x$selected$r))           x$selected$r          else NA_real_
    cat("  Selected model (auto): ", div, sep = "")
    if (is.finite(rr)) cat(sprintf(" (r = %s)", format(rr)))
    cat("\n")
  } else {
    div <- if (!is.null(x$selected$divergence)) x$selected$divergence else
      if (!is.null(x$selected$model))       x$selected$model      else "NA"
    rr  <- if (!is.null(x$selected$r))           x$selected$r          else NA_real_
    cat("  Model (fixed): ", div, sep = "")
    if (is.finite(rr)) cat(sprintf(" (r = %s)", format(rr)))
    cat("\n")
  }

  if (!is.null(x$table) && nrow(x$table) > 0) {
    cat("\nDiagnostics by model:\n")
    df <- x$table
    num_cols <- intersect(names(df), c("r","Entropy2","D1","D2"))
    for (cc in num_cols) df[[cc]] <- round(df[[cc]], digits)
    print(df, row.names = FALSE)
  } else {
    cat("\n(no diagnostics table available)\n")
  }
  invisible(x)
}
