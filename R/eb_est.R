# R/eb_est.R
# Estimation routines for generalized entropy balancing (GEB).
# Defines:
#   - EB_est_one(): run a single candidate (divergence ∈ {KL,LW,QLS,TS}, r given)
#   - EB_est():     fixed model OR auto search over {KL} ∪ {LW,QLS,TS} × r_set
#
# IMPORTANT: The user should pass MU_int and MU_ext WITHOUT an intercept.
# Internally, the first-step EB augments an intercept column (1) to the design
# and a leading 1 to MU_ext. The second step adds its own intercept as (1, H).
#
# Assumes the following helpers are available in the package:
#   ent.fun(), lmd.fun(), w.hat.fun()   (see R/entropy.R, R/dual.R)

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
#'   **Do NOT include an intercept column** (the function will add it internally).
#' @param MU_ext  numeric length-p vector of target moments for the same features
#'   as `MU_int`. **Do NOT include an intercept component** (the function will add it internally).
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
#' \item `D1`, `D2`: diagnostics (D1 from features WITHOUT intercept)
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

  # D1: Mahalanobis distance using features WITHOUT intercept (by design)
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
# Internal worker: run one candidate (divergence, r).
# USER INPUTS MUST NOT include intercept in MU_int/MU_ext. Here we ADD the
# intercept internally ONLY for the first-step EB dual.
#
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
  # Internal outcome
  y_int <- dat_int$y_int

  # Internal sample size: read automatically from dat_int
  n <- nrow(dat_int)
  if (nrow(MU_int) != n) {
    stop("dat_int and MU_int must have the same number of rows.")
  }

  # ---- D1 (Mahalanobis distance, using features WITHOUT intercept) ----
  if (is.null(D1_override)) {
    mu_diff <- drop(colMeans(MU_int) - MU_ext)
    Sx      <- stats::cov(MU_int)
    D1      <- as.numeric(sqrt(t(mu_diff) %*% base::solve(Sx) %*% mu_diff))
  } else {
    D1 <- D1_override
  }

  # ---- First step (dual optimization) ----
  # Internally add intercept to both the design and target moments:
  X1  <- cbind(1, MU_int)         # n x (p+1)
  MU1 <- c(1, MU_ext)             # length p+1

  k1 <- 0; ok1 <- FALSE; opt1 <- NULL
  while (k1 < 1000 && !ok1) {
    tries <- lapply(
      1:10,
      function(i) try(
        stats::optim(
          stats::runif(ncol(X1), -0.1, 0.1),                 # dimension matches (p+1)
          lmd.fun(X1, MU1, divergence = divergence, r = r),  # pass augmented X/MU
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

  # FIRST-step linear predictor uses the augmented design:
  LMD1 <- as.numeric(X1 %*% opt1$par)
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
  # Step-2 uses its own intercept to normalize weights; independent of MU_*.
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
    if (!inherits(opt2, "try-error") &&
        opt2$convergence == 0 &&
        sum(opt2$par >= M) < 1) {
      ok2 <- TRUE
    }
    k2 <- k2 + 1
  }

  if (inherits(opt2, "try-error")) {
    w2   <- rep(NA_real_, n)
    ent2 <- NA_real_
    est  <- NA_real_
  } else {
    LMD2 <- as.numeric(eta_int %*% opt2$par)
    w2   <- w.hat.fun(LMD2, "KL", r = 1)
    ent2 <- ent.fun(w2, divergence = "KL", r = 1)
    est  <- if (ok1 && ok2) mean(w2 * y_int) else NA_real_
  }

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



#' Bootstrap Variance (Appendix B; no auto inside bootstrap)
#'
#' @inheritParams EB_est
#' @param n_ext External sample size (n1).
#' @param BB Bootstrap repetitions (B1 = B2 = BB).
#' @param external.boot TRUE resamples external summaries; FALSE fixes them.
#' @param seed Optional RNG seed.
#' @param max_redraws Maximum number of resampling attempts per bootstrap
#'   iteration when the estimator fails to produce a finite estimate.
#'   If this limit is exceeded, the corresponding bootstrap replicate is stored
#'   as NA.
#' @return List with elements:
#'   \itemize{
#'     \item \code{point_estimate} Point estimate from the original data.
#'     \item \code{bootstrap_se}, \code{bootstrap_var} Bootstrap standard error/variance.
#'     \item \code{ci_normal_95}, \code{ci_percentile_95} 95\% confidence intervals.
#'     \item \code{SigmaW_hat} Estimated covariance of (mu_x, eta).
#'     \item \code{theta_boot} Vector of bootstrap estimates.
#'     \item \code{bootstrap_redraws_total} Total number of redraws due to non-convergence.
#'     \item \code{bootstrap_redraws_per_iteration} Integer vector of redraw counts per bootstrap iteration.
#'   }
#' @export
EB_bootstrap_var <- function(
    dat_int, MU_int, MU_ext, eta,
    n_ext, BB = 200,
    divergence = "KL", r = 1,
    w_type = FALSE, second_covariate = NULL,
    link = "identity", M = 10,
    external.boot = TRUE, seed = NULL,
    max_redraws = 50L
) {
  if (!is.null(seed)) set.seed(seed)

  stopifnot(
    is.matrix(MU_int),
    is.numeric(MU_ext),
    ncol(MU_int) == length(MU_ext),
    is.data.frame(dat_int),
    "y_int" %in% names(dat_int)
  )

  max_redraws <- as.integer(max_redraws)
  if (is.na(max_redraws) || max_redraws < 0L) {
    stop("'max_redraws' must be a non-negative integer.")
  }

  # Internal sample size taken automatically from dat_int
  n <- nrow(dat_int)
  if (nrow(MU_int) != n) {
    stop("dat_int and MU_int must have the same number of rows.")
  }

  p  <- ncol(MU_int)
  dW <- p + 1L  # (mu_x[1:p], eta)

  # Point estimate on the original data
  point <- EB_est(dat_int, MU_int, MU_ext, eta,
                  divergence = divergence, r = r,
                  w_type = w_type, second_covariate = second_covariate,
                  link = link, M = M, auto = FALSE)
  theta_hat <- as.numeric(point$result["est"])

  # ------------------------------------------------------------------
  # First-level bootstrap: estimate Sigma_W for (mu_x, eta)
  # ------------------------------------------------------------------
  W_mat <- matrix(NA_real_, nrow = BB, ncol = dW)
  colnames(W_mat) <- c(paste0("mu_x[", 1:p, "]"), "eta")

  for (b1 in seq_len(BB)) {
    idx <- sample.int(n, n, TRUE)
    MU_b  <- MU_int[idx, , drop = FALSE]
    di_b  <- dat_int[idx, , drop = FALSE]
    W_mat[b1, ] <- c(colMeans(MU_b), mean(di_b$y_int))
  }

  SigmaW_hat <- stats::cov(W_mat, use = "complete.obs")
  eig <- eigen(SigmaW_hat, symmetric = TRUE, only.values = TRUE)$values
  if (any(eig < .Machine$double.eps)) {
    SigmaW_hat <- SigmaW_hat + diag(1e-8, nrow(SigmaW_hat))
  }

  mu_ext_obs  <- as.numeric(MU_ext)
  eta_ext_obs <- as.numeric(eta)
  mean_vec    <- c(mu_ext_obs, eta_ext_obs)
  # length(mean_vec) == p + 1 == ncol(SigmaW_hat)

  Sigma_scaled <- (n / n_ext) * SigmaW_hat

  theta_boot    <- numeric(BB)
  redraw_counts <- integer(BB)

  .rmvnorm_ <- function(n, mean, sigma) {
    d <- length(mean)
    if (requireNamespace("MASS", quietly = TRUE)) {
      return(MASS::mvrnorm(n = n, mu = mean, Sigma = sigma))
    } else if (requireNamespace("mvtnorm", quietly = TRUE)) {
      return(mvtnorm::rmvnorm(n = n, mean = mean, sigma = sigma))
    } else {
      ev <- eigen(sigma, symmetric = TRUE)
      ev$values[ev$values < 0] <- 0
      L <- ev$vectors %*% diag(sqrt(ev$values), d, d)
      Z <- matrix(stats::rnorm(n * d), n, d)
      return(sweep(Z %*% t(L), 2, mean, FUN = "+"))
    }
  }

  # ------------------------------------------------------------------
  # Second-level bootstrap: resample until finite estimate or max_redraws
  # ------------------------------------------------------------------
  for (b2 in seq_len(BB)) {
    redraws_b2 <- 0L
    theta_b2   <- NA_real_
    converged  <- FALSE

    while (!converged) {
      idx2  <- sample.int(n, n, TRUE)
      MU_b2 <- MU_int[idx2, , drop = FALSE]
      di_b2 <- dat_int[idx2, , drop = FALSE]

      if (isTRUE(external.boot)) {
        draw <- as.numeric(.rmvnorm_(1, mean = mean_vec, sigma = Sigma_scaled))
        mu_ext_b2  <- draw[1:p]
        eta_ext_b2 <- draw[p + 1L]
      } else {
        mu_ext_b2  <- mu_ext_obs
        eta_ext_b2 <- eta_ext_obs
      }

      # MU_ext_b2 is passed WITHOUT intercept; EB_est() will add it internally
      MU_ext_b2 <- mu_ext_b2

      res <- try(
        EB_est(di_b2, MU_b2, MU_ext_b2, eta_ext_b2,
               divergence = divergence, r = r,
               w_type = w_type, second_covariate = second_covariate,
               link = link, M = M, auto = FALSE),
        silent = TRUE
      )

      theta_b2 <- NA_real_
      if (!inherits(res, "try-error") &&
          !is.null(res$result) &&
          !is.null(res$result["est"])) {
        theta_b2 <- as.numeric(res$result["est"])
      }

      converged <- is.finite(theta_b2)

      if (!converged) {
        redraws_b2 <- redraws_b2 + 1L
        if (redraws_b2 > max_redraws) {
          warning(sprintf(
            "Bootstrap iteration %d: failed to obtain a finite estimate after %d redraws; storing NA.",
            b2, max_redraws
          ))
          break
        }
      }
    }

    theta_boot[b2]    <- theta_b2
    redraw_counts[b2] <- redraws_b2
  }

  # ------------------------------------------------------------------
  # Bootstrap variance, SE, and CIs
  # ------------------------------------------------------------------
  boot_var <- stats::var(theta_boot, na.rm = TRUE)
  boot_se  <- sqrt(boot_var)

  ci_norm <- c(
    lower = theta_hat - 1.96 * boot_se,
    upper = theta_hat + 1.96 * boot_se
  )

  qs <- stats::quantile(theta_boot, c(0.025, 0.975),
                        na.rm = TRUE, names = FALSE)
  ci_perc <- c(lower = qs[1], upper = qs[2])

  total_redraws <- sum(redraw_counts)
  if (total_redraws > 0L) {
    n_iters_redrawn <- sum(redraw_counts > 0L)
    message(sprintf(
      "Bootstrap: %d out of %d iterations required redraws (total redraws = %d).",
      n_iters_redrawn, BB, total_redraws
    ))
  }

  list(
    point_estimate                  = theta_hat,
    bootstrap_se                    = boot_se,
    bootstrap_var                   = boot_var,
    ci_normal_95                    = ci_norm,
    ci_percentile_95                = ci_perc,
    B1                              = BB,
    B2                              = BB,
    n_int                           = n,
    n_ext                           = n_ext,
    SigmaW_hat                      = SigmaW_hat,
    theta_boot                      = theta_boot,
    call_divergence                 = divergence,
    call_r                          = r,
    external.boot                   = external.boot,
    max_redraws                     = max_redraws,
    bootstrap_redraws_total         = total_redraws,
    bootstrap_redraws_per_iteration = redraw_counts
  )
}
