# R/eb_est.R
# Estimation routines for generalized entropy balancing (GEB).
# Defines:
#   - EB_est_one(): run one candidate (divergence ∈ {KL,LW,QLS,TS}, r given)
#   - EB_est():     fixed model OR auto search over {KL} ∪ {LW,QLS,TS}×r_set
# Returns always carry class "daisy_eb" so summary() prints the custom view.

#' Generalized Entropy Balancing Estimator
#'
#' Two-step estimator integrating individual-level data with external summary
#' data via generalized entropy balancing. The second step uses KL (fixed),
#' and (for numerical compatibility with the original script) the second-step
#' weights are computed from the FIRST-step linear predictor.
#'
#' @param dat_int data.frame with internal variables. Must contain a column `y_int`
#'   (response) and predictors as remaining columns (any names).
#' @param MU_int numeric matrix (n x p). First column must be 1 (intercept).
#' @param MU_ext numeric length-p vector of target moments; first element must be 1.
#' @param eta    numeric scalar: external target for step-2 (e.g., mean of outcome).
#' @param divergence character in {"KL","LW","QLS","TS"}; used when `auto = FALSE`.
#' @param r    numeric scalar for LW/QLS/TS (ignored for KL); used when `auto = FALSE`.
#' @param w_type logical; if TRUE, second-step regression is weighted by first-step w.
#' @param w_fixed kept for API compatibility (unused).
#' @param second_covariate optional matrix for additional second-step constraints (unused here).
#' @param non_regression kept for API compatibility (unused).
#' @param link character link for regression: "identity","logit","probit","gamma".
#' @param M numeric threshold for second-step Lagrange multiplier box constraint.
#' @param auto logical; if TRUE, search the best first-step model by Entropy2 over
#'   {KL} ∪ {LW,QLS,TS}×r_set.
#' @param r_set numeric vector of r candidates used when `auto = TRUE`.
#'
#' @return An object of class `"daisy_eb"` with fields:
#' \itemize{
#' \item `result["est"]` : point estimate (mean of w2.hat * y)
#' \item `Entropy2`      : second-step entropy objective value
#' \item `D1`, `D2`      : geometry diagnostics (Mahalanobis-based and chi-square-based)
#' \item `model`, `r`    : chosen divergence and r (fixed case)
#' \item `best_model`    : list(divergence, r) when `auto = TRUE`
#' \item `leaderboard`   : data.frame(divergence, r, Entropy2, D1, D2) when `auto = TRUE`
#' \item `all_results`   : list of per-candidate result objects when `auto = TRUE`
#' }
#'
#' @export
#' @importFrom stats runif cov lm glm binomial Gamma predict
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

  # --- basic checks ---
  if (!is.data.frame(dat_int)) stop("dat_int must be a data.frame.")
  if (!("y_int" %in% names(dat_int))) stop("dat_int must contain column 'y_int'.")
  if (!is.matrix(MU_int)) stop("MU_int must be a matrix.")
  if (!is.numeric(MU_ext)) stop("MU_ext must be a numeric vector.")
  if (length(MU_ext) != ncol(MU_int)) stop("length(MU_ext) must match ncol(MU_int).")
  if (!is.numeric(eta) || length(eta) != 1L) stop("eta must be a numeric scalar.")

  # --- D1 is common across candidates ---
  mu_diff <- drop(colMeans(MU_int) - MU_ext)
  Sx      <- stats::cov(MU_int[, -1, drop = FALSE])
  D1_val  <- c(sqrt(t(mu_diff[-1]) %*% solve(Sx) %*% mu_diff[-1]))

  # ---- single model path ----
  if (!isTRUE(auto)) {
    ans <- EB_est_one(
      dat_int, MU_int, MU_ext, eta,
      divergence = divergence, r = r,
      w_type = w_type, link = link, M = M,
      D1_override = D1_val
    )
    ans$auto <- FALSE
    class(ans) <- unique(c("daisy_eb", class(ans)))
    return(ans)
  }

  # ---- auto search: {KL} ∪ {LW,QLS,TS}×r_set ----
  # Evaluate candidates into a list (preserve list order)
  res_list <- list()

  # 1) KL (r ignored)
  res_list[[length(res_list) + 1L]] <- EB_est_one(
    dat_int, MU_int, MU_ext, eta,
    divergence = "KL", r = NA_real_,
    w_type = w_type, link = link, M = M,
    D1_override = D1_val
  )

  # 2) LW/QLS/TS × r_set
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

  # Build leaderboard (same order as res_list, then sort by Entropy2 desc)
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
  best_idx <- ord[1]               # index into res_list
  best_obj <- res_list[[best_idx]] # aligned with lb before sorting

  # annotate and return
  best_obj$best_model  <- list(divergence = lb$divergence[1], r = lb$r[1])
  best_obj$leaderboard <- lb
  best_obj$all_results <- res_list
  best_obj$auto        <- TRUE
  class(best_obj)      <- unique(c("daisy_eb", class(best_obj)))
  return(best_obj)
}

# ------------------------------------------------------------------------------

# Internal worker: run one candidate (divergence, r)
# Conventions (fixed to reproduce original script):
# * second-step divergence is KL
# * second-step weights are built from the FIRST-step linear predictor (LMD1)
#   (even though exp(LMD2) would be a natural alternative).
#
# Not exported.
EB_est_one <- function(
    dat_int, MU_int, MU_ext, eta,
    divergence = "KL", r = 1,
    w_type = FALSE, link = "identity", M = 10,
    D1_override = NULL
) {
  y_int <- dat_int$y_int
  n     <- nrow(MU_int)

  # D1 is common across candidates; allow override for efficiency
  if (is.null(D1_override)) {
    mu_diff <- drop(colMeans(MU_int) - MU_ext)
    Sx      <- stats::cov(MU_int[, -1, drop = FALSE])
    D1      <- c(sqrt(t(mu_diff[-1]) %*% solve(Sx) %*% mu_diff[-1]))
  } else {
    D1 <- D1_override
  }

  # ---------- 1st-step: dual optimization ----------
  k1 <- 0; ok1 <- FALSE; opt1 <- NULL
  while (k1 < 1000 && !ok1) {
    tries <- lapply(
      1:10,
      function(i) try(
        optim(
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
    # fail-safe return
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

  LMD1 <- c(MU_int %*% opt1$par)
  w1   <- w.hat.fun(LMD1, divergence, r = if (identical(divergence, "KL")) 1 else r)
  D2   <- sqrt(mean(w1^2) - 1)

  # ---------- internal regression (weighted if requested) ----------
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

  # ---------- 2nd-step EB (KL) ----------
  h    <- stats::predict(reg_int)
  calH <- w1 * h
  eta_int <- cbind(1, calH)
  eta_ext <- c(1, eta)

  k2 <- 0; ok2 <- FALSE; opt2 <- NULL
  while (k2 < 1000 && !ok2) {
    tries2 <- lapply(
      1:10,
      function(i) try(
        optim(
          stats::runif(ncol(eta_int), -0.1, 0.1),
          lmd.fun(eta_int, eta_ext, divergence = "KL", r = 0),
          method = "L-BFGS-B", lower = -M, upper =  M
        ),
        silent = TRUE
      )
    )
    vals2 <- sapply(tries2, function(o) if (inherits(o, "try-error")) Inf else o$value)
    opt2  <- tries2[[which.min(vals2)]]
    if (!inherits(opt2, "try-error") && opt2$convergence == 0 && sum(opt2$par >= M) < 1) ok2 <- TRUE
    k2 <- k2 + 1
  }

  # IMPORTANT: keep original compatibility — build w2 from LMD1 (not from LMD2)
  w2   <- w.hat.fun(LMD1, divergence = "KL", r = 1)
  ent2 <- ent.fun(w2, divergence = "KL", r = 1)

  est <- if (ok1 && ok2) mean(w2 * y_int) else NA_real_

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
