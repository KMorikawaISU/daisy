#' Entropy-Balancing Estimator (with optional auto selection)
#'
#' This version enforces **no constant columns (e.g., 1)** in `MU_int` and **no
#' explicit intercept (1)** in `MU_ext` from the user. Any constant/near-constant
#' columns in `MU_int` are dropped with a warning. If `MU_ext` appears to contain
#' an intercept (leading value ~ 1) it is dropped with a warning.
#' An intercept `(Intercept) = 1` is then **added internally** to both
#' `MU_int` and `MU_ext` before estimation so that the core logic can continue to
#' assume an intercept is present in the first column/element.
#'
#' Added features:
#' - Clear input validation and informative errors/warnings.
#' - S3 `summary()` method (lm-like) for the returned object.
#'   (In the summary table, the `r` column is hidden for KL rows.)
#'
#' @param dat_int Data frame with column `y_int` and predictors.
#' @param MU_int Matrix (n x p) of internal moments (**do not include a column of 1s**).
#' @param MU_ext Numeric vector of length p (**do not include a leading 1**).
#' @param eta Scalar external target for step-2.
#' @param divergence "KL","LW","QLS","TS".
#' @param r Tuning parameter for LW/QLS/TS (ignored for KL).
#' @param w_type Logical; TRUE → weighted second-step regression.
#' @param second_covariate Optional extra covariates for step-2 constraints
#'   (n x k numeric). If provided, the first k elements of `MU_ext` are assumed
#'   to correspond to the external moments of these covariates (same order).
#' @param link "identity","logit","probit","gamma".
#' @param M Threshold for step-2 Lagrange multiplier.
#' @param auto If TRUE, searches across KL and the set of (LW/QLS/TS) with r in `r_set`,
#'   by Entropy2, and picks the best.
#' @param r_set Numeric vector of candidate r's for `auto=TRUE`.
#' @return A list (class `"eb_est"`) with the same fields as before plus metadata
#'   for `summary()`. For `auto=TRUE`, it also contains `leaderboard`, `all_results`,
#'   and `best_model`. A precomputed `summary(fit)` is stored in `$summary` for convenience.
#' @export
#' @importFrom stats lm glm binomial cov optim predict quantile var
EB_est <- function(
    dat_int, MU_int, MU_ext, eta,
    divergence = "KL", r = 1,
    w_type = FALSE, second_covariate = NULL,
    link = "identity", M = 10,
    auto = FALSE, r_set = c(0.01, 0.1, 0.5, 1)
) {
  # ---------- local helpers (not exported) ----------
  .near <- function(x, target = 1, tol = 1e-8) {
    is.finite(x) && is.finite(target) && abs(as.numeric(x) - target) < tol
  }

  .validate_mu_data <- function(MU, name = "MU", const_tol = 1e-12) {
    # Accept numeric matrix or data.frame (numeric-only); no NA/NaN/Inf.
    if (is.null(MU)) stop(sprintf("`%s` cannot be NULL.", name), call. = FALSE)
    if (is.data.frame(MU)) {
      is_num <- vapply(MU, is.numeric, logical(1))
      if (!all(is_num)) {
        bad <- names(MU)[!is_num]
        stop(sprintf("`%s` must be numeric. Non-numeric columns: %s.",
                     name, paste(bad, collapse = ", ")), call. = FALSE)
      }
      MU <- as.matrix(MU)
    } else if (!is.matrix(MU)) {
      stop(sprintf("`%s` must be a matrix or data.frame.", name), call. = FALSE)
    }
    if (!is.numeric(MU)) stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
    if (anyNA(MU) || any(!is.finite(MU))) {
      stop(sprintf("`%s` contains NA/NaN/Inf.", name), call. = FALSE)
    }
    if (nrow(MU) < 1L || ncol(MU) < 1L) {
      stop(sprintf("`%s` must have positive dimensions.", name), call. = FALSE)
    }

    # Detect constant / near-constant columns; drop them with warnings.
    cn <- colnames(MU); if (is.null(cn)) cn <- paste0("V", seq_len(ncol(MU)))
    rng <- apply(MU, 2, function(col) diff(range(col)))
    sds <- suppressWarnings(apply(MU, 2, stats::sd))
    const <- (is.na(sds) & (rng < const_tol)) | (!is.na(sds) & (sds < const_tol))

    if (any(const)) {
      mns <- colMeans(MU)
      near1 <- const & vapply(mns, .near, logical(1), target = 1, tol = 1e-8)
      if (any(near1)) {
        warning(sprintf("`%s`: constant column(s) equal to 1 detected and dropped: %s.",
                        name, paste(cn[near1], collapse = ", ")), call. = FALSE)
        warning("Do not include constants (e.g., 1) in MU_int / MU_ext; an intercept is added internally.",
                call. = FALSE)
      }
      other_const <- const & !near1
      if (any(other_const)) {
        warning(sprintf("`%s`: constant/near-constant column(s) detected and dropped: %s.",
                        name, paste(cn[other_const], collapse = ", ")), call. = FALSE)
      }
      MU <- MU[, !const, drop = FALSE]
      cn <- cn[!const]
    }

    if (ncol(MU) == 0L) {
      stop(sprintf("All columns of `%s` were constant and dropped. Provide at least one non-constant covariate.",
                   name), call. = FALSE)
    }

    colnames(MU) <- cn

    # Rank deficiency warning (without intercept)
    if (ncol(MU) > 1L) {
      rnk <- qr(MU)$rank
      if (!is.na(rnk) && rnk < ncol(MU)) {
        warning(sprintf("`%s` appears rank-deficient (rank %d < %d).",
                        name, rnk, ncol(MU)), call. = FALSE)
      }
    }
    MU
  }

  .sanitize_mu_ext <- function(MU_ext, p_expected, tol = 1e-8) {
    # Accept numeric vector or 1-row/1-col matrix; coerce to numeric vector.
    if (is.null(MU_ext)) stop("`MU_ext` cannot be NULL.", call. = FALSE)

    if (is.matrix(MU_ext)) {
      if (min(nrow(MU_ext), ncol(MU_ext)) != 1L) {
        stop("`MU_ext` must be a numeric vector or a 1D matrix (1 x p or p x 1).", call. = FALSE)
      }
      MU_ext <- as.numeric(MU_ext)
    } else if (is.data.frame(MU_ext)) {
      MU_ext <- as.numeric(as.matrix(MU_ext))
    } else if (!is.numeric(MU_ext)) {
      stop("`MU_ext` must be numeric.", call. = FALSE)
    }
    if (anyNA(MU_ext) || any(!is.finite(MU_ext))) {
      stop("`MU_ext` contains NA/NaN/Inf.", call. = FALSE)
    }

    len <- length(MU_ext)
    if (len == p_expected + 1L && .near(MU_ext[1], 1, tol)) {
      # Drop user-supplied intercept and warn.
      warning("`MU_ext` appears to include an intercept (=1). It will be ignored; an intercept is added internally.",
              call. = FALSE)
      MU_ext <- MU_ext[-1L]
      len <- length(MU_ext)
    }
    if (len != p_expected) {
      stop(sprintf("Length mismatch in `MU_ext`: expected %d (no intercept), got %d.",
                   p_expected, len), call. = FALSE)
    }
    MU_ext
  }

  # ---------- input validation (dat_int / basic args) ----------
  if (!is.data.frame(dat_int)) stop("`dat_int` must be a data.frame.", call. = FALSE)
  if (!("y_int" %in% names(dat_int))) stop("`dat_int` must contain a column named `y_int`.", call. = FALSE)
  if (!is.numeric(dat_int$y_int)) stop("`dat_int$y_int` must be numeric.", call. = FALSE)

  if (!is.logical(w_type) || length(w_type) != 1L || is.na(w_type)) {
    stop("`w_type` must be TRUE/FALSE.", call. = FALSE)
  }
  link <- match.arg(link, c("identity","logit","probit","gamma"))
  if (!is.numeric(M) || length(M) != 1L || !is.finite(M)) {
    stop("`M` must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.logical(auto) || length(auto) != 1L || is.na(auto)) {
    stop("`auto` must be TRUE/FALSE.", call. = FALSE)
  }
  if (isTRUE(auto)) {
    if (!is.numeric(r_set) || length(r_set) < 1L) stop("`r_set` must be a numeric vector.", call. = FALSE)
  } else {
    divergence <- match.arg(divergence, c("KL","LW","QLS","TS"))
    if (!is.numeric(r) || length(r) != 1L || !is.finite(r)) {
      stop("`r` must be a finite numeric scalar.", call. = FALSE)
    }
  }

  # ---------- validate MU_int (drop constants; no user intercept) ----------
  MU_int <- .validate_mu_data(MU_int, name = "MU_int", const_tol = 1e-12)

  # Alignment to n
  n <- nrow(dat_int)
  if (nrow(MU_int) != n) {
    stop(sprintf("Row mismatch: `nrow(MU_int)` (%d) must equal `nrow(dat_int)` (%d).",
                 nrow(MU_int), n), call. = FALSE)
  }

  # second_covariate checks
  if (!is.null(second_covariate)) {
    if (is.vector(second_covariate)) second_covariate <- matrix(second_covariate, ncol = 1L)
    if (is.data.frame(second_covariate)) {
      is_num <- vapply(second_covariate, is.numeric, logical(1))
      if (!all(is_num)) {
        bad <- names(second_covariate)[!is_num]
        stop(sprintf("`second_covariate` must be numeric. Non-numeric columns: %s.",
                     paste(bad, collapse = ", ")), call. = FALSE)
      }
      second_covariate <- as.matrix(second_covariate)
    } else if (!is.matrix(second_covariate)) {
      stop("`second_covariate` must be a numeric matrix/data.frame or NULL.", call. = FALSE)
    }
    if (nrow(second_covariate) != n) {
      stop("`second_covariate` must have the same number of rows as `dat_int`.", call. = FALSE)
    }
  }

  # ---------- sanitize MU_ext (no user intercept; length must match) ----------
  p <- ncol(MU_int)  # user-supplied moments without intercept (after dropping constants)
  MU_ext <- .sanitize_mu_ext(MU_ext, p_expected = p, tol = 1e-8)

  # ---------- internally add intercept = 1 ----------
  MU_int2 <- cbind(`(Intercept)` = 1, MU_int)
  MU_ext2 <- c(1, as.numeric(MU_ext))

  # ---------- dispatch (auto vs single) ----------
  if (!isTRUE(auto)) {
    res <- EB_est_one(dat_int, MU_int2, MU_ext2, eta,
                      divergence = divergence, r = r,
                      w_type = w_type, second_covariate = second_covariate,
                      link = link, M = M)
    # Attach metadata and class for summary()
    res$.__eb_meta__ <- list(
      call            = match.call(),
      auto            = FALSE,
      intercept_added = TRUE,
      mu_cols         = colnames(MU_int), # user-supplied names (without intercept)
      mu_ncol         = ncol(MU_int),
      validation_ok   = TRUE
    )
    class(res) <- unique(c("eb_est", class(res)))
    # Precompute summary object for convenience
    res$summary <- tryCatch(summary(res), error = function(e) NULL)
    return(res)
  }

  # auto = TRUE: enumerate candidates
  cand <- list(list(div = "KL", r = 1))
  for (dv in c("LW","QLS","TS")) for (rr in r_set) cand[[length(cand) + 1L]] <- list(div = dv, r = rr)

  all_results <- list(); ent_vec <- numeric(length(cand)); keys <- character(length(cand))
  for (i in seq_along(cand)) {
    dv <- cand[[i]]$div; rr <- cand[[i]]$r
    key <- if (dv == "KL") "KL" else paste0(dv, "_r=", format(rr, trim = TRUE, scientific = FALSE))
    keys[i] <- key
    res_i <- EB_est_one(dat_int, MU_int2, MU_ext2, eta,
                        divergence = dv, r = rr,
                        w_type = w_type, second_covariate = second_covariate,
                        link = link, M = M)
    all_results[[key]] <- res_i
    ent_vec[i] <- ifelse(is.na(res_i$Entropy2), -Inf, res_i$Entropy2)
  }

  best_idx <- which.max(ent_vec)
  best <- all_results[[keys[best_idx]]]
  leaderboard <- data.frame(
    candidate  = keys,
    divergence = vapply(cand, function(z) z$div, ""),
    r          = vapply(cand, function(z) z$r, 0),
    Entropy2   = ent_vec,
    stringsAsFactors = FALSE
  )
  leaderboard <- leaderboard[order(-leaderboard$Entropy2), ]

  best$best_model  <- list(divergence = cand[[best_idx]]$div, r = cand[[best_idx]]$r)
  best$leaderboard <- leaderboard
  best$all_results <- all_results
  best$auto        <- TRUE

  # Attach metadata and class for summary()
  best$.__eb_meta__ <- list(
    call            = match.call(),
    auto            = TRUE,
    intercept_added = TRUE,
    mu_cols         = colnames(MU_int), # user-supplied names (without intercept)
    mu_ncol         = ncol(MU_int),
    validation_ok   = TRUE
  )
  class(best) <- unique(c("eb_est", class(best)))
  best$summary <- tryCatch(summary(best), error = function(e) NULL)
  best
}

#' @keywords internal
EB_est_one <- function(
    dat_int, MU_int, MU_ext, eta,
    divergence = "KL", r = 1,
    w_type = FALSE, second_covariate = NULL,
    link = "identity", M = 10
) {
  # `MU_int` is assumed to include an intercept in the first column (added upstream).
  # `MU_ext` is assumed to include an intercept as its first element (added upstream).
  y_int <- dat_int$y_int

  # D1 (Mahalanobis-like) based on non-intercept columns
  mu_diff <- drop(colMeans(MU_int) - MU_ext)
  if (ncol(MU_int) <= 1L) {
    D1 <- NA_real_
    Sx <- matrix(NA_real_, 0, 0)
  } else {
    Sx <- stats::cov(MU_int[, -1, drop = FALSE])
    eig <- eigen(Sx, symmetric = TRUE, only.values = TRUE)$values
    if (any(eig < .Machine$double.eps)) Sx <- Sx + diag(1e-8, nrow(Sx))
    D1 <- c(sqrt(t(mu_diff[-1]) %*% solve(Sx) %*% mu_diff[-1]))
  }

  # First-step weights via lambda optimization
  k1 <- 0; ok1 <- 0
  while ((k1 < 1000) & (ok1 == 0)) {
    lmd.est00 <- lapply(1:10, function(.) try(
      stats::optim(stats::runif(ncol(MU_int), -0.1, 0.1),
                   lmd.fun(MU_int, MU_ext, divergence, r),
                   method = "BFGS"), silent = TRUE))
    vals <- sapply(lmd.est00, function(o) if (inherits(o, "try-error")) Inf else o$value)
    lmd.est0 <- lmd.est00[[which.min(vals)]]
    if (inherits(lmd.est0, "try-error")) { k1 <- k1 + 1; next }
    ok1 <- as.numeric(lmd.est0$convergence == 0)
    k1  <- k1 + 1
  }
  LMD.est <- c(MU_int %*% lmd.est0$par)
  w.hat   <- w.hat.fun(LMD.est, divergence, r)
  D2      <- sqrt(mean(w.hat^2) - 1)

  # Second-step regression (weighted or unweighted)
  ww <- if (isFALSE(w_type)) rep(1, length(w.hat)) else w.hat
  reg_int <- switch(link,
                    identity = try(stats::lm(y_int ~ ., data = dat_int, weights = ww), silent = TRUE),
                    logit    = try(stats::glm(y_int ~ ., data = dat_int,
                                              family = stats::binomial("logit"),  weights = ww), silent = TRUE),
                    probit   = try(stats::glm(y_int ~ ., data = dat_int,
                                              family = stats::binomial("probit"), weights = ww), silent = TRUE),
                    gamma    = try(stats::glm(y_int ~ ., data = dat_int,
                                              family = stats::Gamma(link = "log"), weights = ww), silent = TRUE),
                    stop("Unknown link.")
  )

  entropy2 <- NA_real_
  if (!inherits(reg_int, "try-error")) {
    h <- stats::predict(reg_int); calH <- w.hat * h

    if (is.null(second_covariate)) {
      # Only intercept and h
      eta_int     <- as.matrix(cbind(1, calH))
      eta_ext_vec <- c(1, eta)
    } else {
      # Intercept + (w.hat * second_covariate) + h
      eta_int <- as.matrix(cbind(1, second_covariate * w.hat, calH))
      # Assumption retained: the first k entries of MU_ext (without intercept upstream)
      # correspond to external moments of second_covariate (same order).
      k_sc <- ncol(second_covariate)
      eta_ext_vec <- c(MU_ext[1:(k_sc + 1L)], eta)  # MU_ext includes intercept at [1]
    }

    # Second-step lambda optimization (KL with box constraint by M)
    k2 <- 0; ok2 <- 0
    while ((k2 < 1000) & (ok2 == 0)) {
      lmd2.est00 <- lapply(1:10, function(.) try(
        stats::optim(stats::runif(ncol(eta_int), -0.1, 0.1),
                     lmd.fun(eta_int, eta_ext_vec, divergence = "KL", r = 0),
                     method = "L-BFGS-B", lower = -5, upper = 5), silent = TRUE))
      vals2 <- sapply(lmd2.est00, function(o) if (inherits(o, "try-error")) Inf else o$value)
      lmd2.est0 <- lmd2.est00[[which.min(vals2)]]
      if (inherits(lmd2.est0, "try-error")) { k2 <- k2 + 1; next }
      ok2 <- as.numeric(lmd2.est0$convergence == 0 & sum(lmd2.est0$par >= M) < 1)
      k2  <- k2 + 1
    }

    # Second-step weights built from the SECOND-step linear predictor
    # (use KL for mapping to weights; this is consistent with Entropy2 below)
    LMD2.est <- c(eta_int %*% lmd2.est0$par)
    w2.hat   <- w.hat.fun(LMD2.est, divergence = "KL", r = 1)

    entropy2 <- ent.fun(w2.hat, divergence = "KL", r = 1)
    est      <- if ((k1 == 1000) | (k2 == 1000)) NA_real_ else mean(w2.hat * y_int)

    return(list(model = divergence, r = r, w_type = w_type,
                Entropy2 = entropy2, res.w1 = lmd.est0, res.w2 = lmd2.est0,
                result = c(est = est, lmd2.est = lmd2.est0$par), w2.hat = w2.hat,
                D1 = D1, D2 = D2))
  } else {
    return(list(model = divergence, r = r, w_type = w_type,
                Entropy2 = entropy2, res.w1 = reg_int, res.w2 = reg_int,
                result = c(est = NA_real_, lmd2.est = NA_real_), w2.hat = NA,
                D1 = D1, D2 = D2))
  }
}

# ------------------------------
# S3 methods: summary() and print()
# ------------------------------

#' @export
summary.eb_est <- function(object, ...) {
  x <- object
  meta <- x$.__eb_meta__
  auto <- isTRUE(meta$auto)

  # Gather candidate models for auto; else single-model fallback.
  models <- NULL
  if (!is.null(x$all_results) && is.list(x$all_results) && length(x$all_results) > 0L) {
    models <- unname(x$all_results)
  } else {
    models <- list(x)
  }

  # Helper to fetch numeric scalars safely
  get_num <- function(m, keys) {
    for (k in keys) {
      if (!is.null(m[[k]]) && is.numeric(m[[k]]) && length(m[[k]]) == 1L) return(as.numeric(m[[k]]))
      if (grepl("\\.", k, fixed = TRUE)) {
        parts <- strsplit(k, ".", fixed = TRUE)[[1]]
        z <- m; ok <- TRUE
        for (p in parts) { if (is.null(z[[p]])) { ok <- FALSE; break } else z <- z[[p]] }
        if (ok && is.numeric(z) && length(z) == 1L) return(as.numeric(z))
      }
    }
    NA_real_
  }

  # Build model comparison table (hide r for KL by setting it to NA for display)
  mt <- lapply(seq_along(models), function(i) {
    m <- models[[i]]

    div <- if (!is.null(m$model)) as.character(m$model) else NA_character_
    r_raw <- if (!is.null(m$r)) as.numeric(m$r) else NA_real_
    r_disp <- if (!is.na(div) && div == "KL") NA_real_ else r_raw

    D1  <- get_num(m, c("D1", "metrics.D1", "diagnostics.D1"))
    D2  <- get_num(m, c("D2", "metrics.D2", "diagnostics.D2"))
    ent <- get_num(m, c("Entropy2", "metrics.Entropy2"))
    est <- NA_real_; if (!is.null(m$result) && is.numeric(m$result["est"])) est <- as.numeric(m$result["est"])

    data.frame(
      model      = i,
      divergence = div,
      r          = r_disp,
      D1         = D1,
      D1_lt_1    = ifelse(is.finite(D1), D1 < 1, NA),
      D2         = D2,
      D2_lt_1    = ifelse(is.finite(D2), D2 < 1, NA),
      Entropy2   = ent,
      est        = est,
      stringsAsFactors = FALSE
    )
  })
  model_table <- do.call(rbind, mt)

  # Select best model (maximize Entropy2 if available, otherwise minimize D1 + D2)
  best_idx <- NA_integer_
  if (any(is.finite(model_table$Entropy2))) {
    score <- ifelse(is.finite(model_table$Entropy2), model_table$Entropy2, -Inf)
    best_idx <- which.max(score)
  } else {
    s <- rowSums(model_table[, c("D1", "D2")], na.rm = TRUE)
    s[!is.finite(model_table$D1) & !is.finite(model_table$D2)] <- Inf
    best_idx <- which.min(s)
  }
  best_idx <- as.integer(best_idx)

  # Coefficients (if available). Many EB pipelines do not store them; NA is fine.
  best_model <- models[[best_idx]]
  coef_table <- NULL
  if (!is.null(best_model$coef) || !is.null(best_model$coefficients)) {
    bet <- best_model$coef; if (is.null(bet)) bet <- best_model$coefficients
    nm <- names(bet); if (is.null(nm)) nm <- paste0("b", seq_along(bet))
    V <- NULL; for (nmV in c("vcov", "V", "Sigma")) if (!is.null(best_model[[nmV]])) { V <- best_model[[nmV]]; break }
    se <- rep(NA_real_, length(bet)); z <- p <- rep(NA_real_, length(bet))
    if (is.matrix(V) && nrow(V) == length(bet) && ncol(V) == length(bet)) {
      se <- sqrt(pmax(diag(V), 0))
      if (all(is.finite(se)) && all(se > 0)) {
        z <- bet / se
        p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
      }
    }
    coef_table <- data.frame(
      term = nm, estimate = as.numeric(bet),
      std.error = as.numeric(se), statistic = as.numeric(z), p.value = as.numeric(p),
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    call         = if (!is.null(meta$call)) meta$call else x$call,
    auto         = auto,
    model_table  = model_table,
    best_index   = best_idx,
    coef_table   = coef_table,
    intercept    = "(Intercept) added internally",
    mu_cols      = if (!is.null(meta$mu_cols)) meta$mu_cols else NULL
  )
  class(out) <- "summary.eb_est"
  out
}

#' @export
print.summary.eb_est <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Summary of eb_est object\n")
  if (!is.null(x$call)) { cat("Call: "); print(x$call) }
  cat("\nIntercept: ", x$intercept, "\n", sep = "")
  if (!is.null(x$mu_cols)) {
    cat("Moment columns (user-supplied, without intercept): ",
        paste(x$mu_cols, collapse = ", "), "\n", sep = "")
  }
  cat("\nModel comparison:\n")
  tbl <- utils::head(x$model_table, n = nrow(x$model_table))

  # Pretty-print: show blank for r where it's NA (i.e., KL rows)
  if ("r" %in% names(tbl)) {
    tbl$r <- ifelse(is.na(tbl$r), "", format(tbl$r, trim = TRUE))
  }

  print(tbl, digits = digits, row.names = FALSE)
  cat("\nSelected (best) model index: ", x$best_index, "\n", sep = "")
  if (!is.null(x$coef_table)) {
    cat("\nCoefficients of the selected model:\n")
    print(x$coef_table, digits = digits, row.names = FALSE)
  } else {
    cat("\n(No coefficient table available in the object.)\n")
  }
  invisible(x)
}
