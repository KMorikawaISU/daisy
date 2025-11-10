# R/eb_est.R
# Estimation routines for generalized entropy balancing (GEB).
# Defines:
#   - EB_est_one(): run a single candidate (divergence ∈ {KL,LW,QLS,TS}, r given)
#   - EB_est():     fixed model OR auto search over {KL} ∪ {LW,QLS,TS} × r_set
#
# Assumes the following helpers are available in the package:
#   ent.fun(), lmd.fun(), w.hat.fun()   (see R/entropy.R, R/dual.R)
#
# IMPORTANT: First-step moments MU_int / MU_ext MUST NOT include an intercept.
# The second step uses its own intercept internally to normalize weights.

#' Generalized Entropy Balancing Estimator
#'
#' Two-step estimator integrating individual-level data with external summary data
#' via generalized entropy balancing. The second step is fixed to KL and, for
#' numerical compatibility with the original script, the second-step weights are
#' computed from the FIRST-step linear predictor (not from LMD2).
#'
#' @param dat_int data.frame with internal variables. Must contain column `y_int`
#'   (response) and predictors as remaining columns (any names).
#' @param MU_int  numeric matrix (n x p) of first-step balancing features.
#'   **Do NOT include an intercept column**.
#' @param MU_ext  numeric length-p vector of target moments for the same features
#'   as `MU_int`. **Do NOT include an intercept component**.
#' @param eta     numeric scalar: external target used in step-2 (e.g., mean outcome).
#' @param divergence character in {"KL","LW","QLS","TS"}; used when `auto = FALSE`.
#' @param r       numeric; tuning for LW/QLS/TS (ignored for KL) when `auto = FALSE`.
#' @param w_type  logical; if TRUE, the internal regression is weighted by first-step w.
#' @param w_fixed kept for API compatibility (unused).
#' @param second_covariate optional matrix for additional second-step constraints (unused here).
#' @param non_regression   kept for API compatibility (unused).
#' @param link    "identity" | "logit" | "probit" | "gamma".
#' @param M       numeric; box bound for step-2 lambda.
#' @param auto    logical; if TRUE, search over {KL} ∪ {LW,QLS,TS} × r_set by Entropy2.
#' @param r_set   numeric vector of r candidates used when `auto = TRUE`.
#'
#' @return Object of class `"daisy_eb"`. For a single run:
#' \itemize{
#' \item `result["est"]`: point estimate (mean of w2 * y)
#' \item `Entropy2`: step-2 entropy value
#' \item `D1`, `D2`: diagnostics
#' \item `model`, `r`: first-step spec
#' }
#' If `auto = TRUE`, also:
#' \itemize{
#' \item `best_model`: list(divergence, r)
#' \item `leaderboard`: data.frame(divergence, r, Entropy2, D1, D2)
#' \item `all_results`: list of per-candidate result lists
#' }
#' @export
#' @importFrom stats runif cov lm glm binomial Gamma predict optim
EB_est <- function(
    dat_int,
    MU_int,
    MU_ext,
    eta,
    divergence       = "KL",
    r                = 1,
    w_type           = FALSE,
    w_fixed          = NULL,
    second_covariate = NULL,
    non_regression   = TRUE,
    link             = "identity",
    M                = 10,
    auto             = FALSE,
    r_set            = c(0.01, 0.1, 0.5, 1)
) {
  # Basic checks
  if (!is.data.frame(dat_int)) stop("dat_int must be a data.frame.")
  if (!("y_int" %in% names(dat_int))) stop("dat_int must contain 'y_int'.")
  if (!is.matrix(MU_int)) stop("MU_int must be a matrix.")
  if (!is.numeric(MU_ext)) stop("MU_ext must be numeric.")
  if (length(MU_ext) != ncol(MU_int)) stop("length(MU_ext) must match ncol(MU_int).")
  if (!is.numeric(eta) || length(eta) != 1L) stop("eta must be a numeric scalar.")

  # D1: Mahalanobis distance using all balancing features (no intercept included here)
  mu_diff <- drop(colMeans(MU_int) - MU_ext)          # length p
  Sx      <- stats::cov(MU_int)                       # p x p
  D1_val  <- as.numeric(sqrt(t(mu_diff) %*% base::solve(Sx) %*% mu_diff))

  # Fixed model path
  if (!isTRUE(auto)) {
    ans <- EB_est_one(
      dat_int = dat_int, MU_int = MU_int, MU_ext = MU_ext, eta = eta,
      divergence = divergence, r = r,
      w_type = w_type, link = link, M = M,
      D1_override = D1_val
    )
    ans$auto <- FALSE
    class(ans) <- unique(c("daisy_eb", class(ans)))
    return(ans)
  }

  # Auto search: {KL} ∪ {LW,QLS,TS} × r_set
  res_list <- list()

  # KL (r ignored)
  res_list[[length(res_list) + 1L]] <- EB_est_one(
    dat_int, MU_int, MU_ext, eta,
    divergence = "KL", r = NA_real_,
    w_type = w_type, link = link, M = M,
    D1_override = D1_val
  )

  # LW/QLS/TS × r_set
  for (div in c("LW", "QLS", "TS")) {
    for (rr in r_set) {
      res_list[[length(res_list) + 1L]] <- EB_est_one(
        dat_int, MU_int, MU_ext, eta,
        divergence = div, r = rr,
        w_type = w_type, link = link, M = M,
        D1_override = D1_val
      )
    }
  }

  # Leaderboard and selection
  lb <- do.call(
    rbind,
    lapply(res_list, function(z) {
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
  ord <- order(lb$Entropy2, decreasing = TRUE)
  lb  <- lb[ord, , drop = FALSE]
  best_idx <- ord[1]
  best_obj <- res_list[[best_idx]]

  best_obj$best_model  <- list(divergence = lb$divergence[1], r = lb$r[1])
  best_obj$leaderboard <- lb
  best_obj$all_results <- res_list
  best_obj$auto        <- TRUE
  class(best_obj)      <- unique(c("daisy_eb", class(best_obj)))
  return(best_obj)
}

# ------------------------------------------------------------------------------
# Internal worker: run one candidate (divergence, r)
# Conventions (kept to reproduce the original script exactly):
# * Second-step divergence is KL.
# * Second-step weights are built from the FIRST-step linear predictor (LMD1),
#   not from LMD2 (even though exp(LMD2) would be the natural alternative).
EB_est_one <- function(
    dat_int, MU_int, MU_ext, eta,
    divergence = "KL", r = 1,
    w_type = FALSE, link = "identity", M = 10,
    D1_override = NULL
) {
  y_int <- dat_int$y_int
  n     <- nrow(MU_int)

  # D1 (allow override for efficiency)
  if (is.null(D1_override)) {
    mu_diff <- drop(colMeans(MU_int) - MU_ext)
    Sx      <- stats::cov(MU_int)
    D1      <- as.numeric(sqrt(t(mu_diff) %*% base::solve(Sx) %*% mu_diff))
  } else {
    D1 <- D1_override
  }

  # ---- First step (dual optimization) ----
  k1 <- 0; ok1 <- FALSE; opt1 <- NULL
  while (k1 < 1000 && !ok1) {
    tries <- lapply(
      1:10,
      function(i) try(
        stats::optim(
          stats::runif(ncol(MU_int), -0.1, 0.1),
          lmd.fun(MU_int, MU_ext, divergence = divergence, r = r),
          method = "BFGS"
        ),
        silent = TRUE
      )
    )
    vals <- sapply(tries, function(o) if (inherits(o, "try-error")) Inf else o$value)
    opt1 <- tries[[which.min(vals)]]
    if (!inherits(opt1, "try-error") && opt1$convergence == 0) ok1 <- TRUE
    k1 <- k1 + 1
  }
  if (!ok1) {
    out <- list(
      model    = divergence,
      r        = if (identical(divergence, "KL")) NA_real_ else r,
      w_type   = w_type,
      Entropy2 = NA_real_,
      res.w1   = opt1,
      res.w2   = NULL,
      result   = c(est = NA_real_),
      w2.hat   = rep(NA_real_, n),
      D1       = D1,
      D2       = NA_real_
    )
    return(out)
  }

  LMD1 <- as.numeric(MU_int %*% opt1$par)
  w1   <- w.hat.fun(LMD1, divergence, r = if (identical(divergence, "KL")) 1 else r)

  # D2 with guard: if the argument of sqrt is negative or non-finite, set NA
  d2_arg <- mean(w1^2) - 1
  D2 <- if (!is.finite(d2_arg) || d2_arg < 0) NA_real_ else sqrt(d2_arg)

  # ---- Internal regression (weighted if requested) ----
  ww <- if (isTRUE(w_type)) w1 else rep(1, n)
  reg_int <- switch(
    link,
    "identity" = try(stats::lm(y_int ~ ., data = dat_int, weights = ww), silent = TRUE),
    "logit"    = try(stats::glm(y_int ~ ., data = dat_int, family = stats::binomial("logit"),  weights = ww), silent = TRUE),
    "probit"   = try(stats::glm(y_int ~ ., data = dat_int, family = stats::binomial("probit"), weights = ww), silent = TRUE),
    "gamma"    = try(stats::glm(y_int ~ ., data = dat_int, family = stats::Gamma("log"),      weights = ww), silent = TRUE),
    stop("Unknown link.")
  )
  if (inherits(reg_int, "try-error")) {
    out <- list(
      model    = divergence,
      r        = if (identical(divergence, "KL")) NA_real_ else r,
      w_type   = w_type,
      Entropy2 = NA_real_,
      res.w1   = opt1,
      res.w2   = reg_int,
      result   = c(est = NA_real_),
      w2.hat   = rep(NA_real_, n),
      D1       = D1,
      D2       = D2
    )
    return(out)
  }

  # ---- Second step (KL) ----
  # NOTE: Step-2 uses its own intercept to normalize weights; this is independent of MU_*.
  h    <- stats::predict(reg_int)
  calH <- w1 * h
  eta_int <- cbind(1, calH)   # (intercept, H)
  eta_ext <- c(1, eta)

  k2 <- 0; ok2 <- FALSE; opt2 <- NULL
  while (k2 < 1000 && !ok2) {
    tries2 <- lapply(
      1:10,
      function(i) try(
        stats::optim(
          stats::runif(ncol(eta_int), -0.1, 0.1),
          lmd.fun(eta_int, eta_ext, divergence = "KL", r = 0),
          method = "L-BFGS-B", lower = -M, upper = M
        ),
        silent = TRUE
      )
    )
    vals2 <- sapply(tries2, function(o) if (inherits(o, "try-error")) Inf else o$value)
    opt2  <- tries2[[which.min(vals2)]]
    if (!inherits(opt2, "try-error") && opt2$convergence == 0 && sum(opt2$par >= M) < 1) ok2 <- TRUE
    k2 <- k2 + 1
  }

  # IMPORTANT: keep original numerical compatibility — build w2 from FIRST-step LMD
  w2   <- w.hat.fun(LMD1, divergence = "KL", r = 1)
  ent2 <- ent.fun(w2, divergence = "KL", r = 1)
  est  <- if (ok1 && ok2) mean(w2 * y_int) else NA_real_

  out <- list(
    model    = divergence,
    r        = if (identical(divergence, "KL")) NA_real_ else r,
    w_type   = w_type,
    Entropy2 = ent2,
    res.w1   = opt1,
    res.w2   = opt2,
    result   = c(est = est),
    w2.hat   = w2,
    D1       = D1,
    D2       = D2
  )
  return(out)
}
