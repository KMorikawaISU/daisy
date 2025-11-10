# eb_est.R
# Generalized Entropy Balancing estimator (two-step EB)
# This file defines:
#   - EB_est_one(): run EB once for a specified divergence and r
#   - EB_est():     run EB once (auto=FALSE) or do auto model search (auto=TRUE)
# The function relies on helper functions defined elsewhere:
#   ent.fun(), lmd.fun(), w.hat.fun()  (see your R/entropy.R, R/dual.R, etc.)

#' One-run EB estimator for a specified first-step divergence
#'
#' @param dat_int data.frame; must contain a numeric column `y_int` and predictors
#' @param MU_int  n x p matrix; first column must be 1 (intercept)
#' @param MU_ext  length-p target moments; MU_ext[1] must be 1 (intercept)
#' @param eta     scalar; external target for step-2 (e.g., external mean outcome)
#' @param divergence character; one of "KL","LW","QLS","TS"
#' @param r       numeric; tuning parameter for LW/QLS/TS (ignored for KL)
#' @param w_type  logical; if TRUE, weight internal regression by first-step w
#' @param second_covariate optional matrix used in step-2 constraints (default NULL)
#' @param link    "identity" | "logit" | "probit" | "gamma" (for internal regression)
#' @param M       numeric; upper bound threshold for step-2 lambda
#' @return list with fields:
#'   model, r, w_type, Entropy2, res.w1, res.w2, result (est, lmd2.est...),
#'   w2.hat, D1, D2
#' @keywords internal
#' @importFrom stats lm glm binomial Gamma predict cov optim runif
EB_est_one <- function(
    dat_int,
    MU_int,
    MU_ext,
    eta,
    divergence = "KL",
    r = 1,
    w_type = FALSE,
    second_covariate = NULL,
    link = "identity",
    M = 10
) {
  # --- Input objects ---
  y_int <- dat_int$y_int
  n     <- nrow(MU_int)

  # --- D1: Mahalanobis distance between internal and external X means ---
  mu_diff <- drop(colMeans(MU_int) - MU_ext)                 # length p
  Sx      <- stats::cov(MU_int[, -1, drop = FALSE])          # exclude intercept
  D1      <- as.numeric(sqrt(t(mu_diff[-1]) %*% solve(Sx) %*% mu_diff[-1]))

  # ---------- First-step EB (dual optimization) ----------
  k1 <- 0
  ok1 <- 0
  lmd1_best <- NULL

  while ((k1 < 1000) && (ok1 == 0)) {
    tries <- lapply(
      1:10,
      function(i) {
        set.seed(i + 123)
        init <- stats::runif(ncol(MU_int), -0.1, 0.1)
        try(
          stats::optim(
            par    = init,
            fn     = lmd.fun(MU_int, MU_ext, divergence = divergence, r = r),
            method = "BFGS"
          ),
          silent = TRUE
        )
      }
    )
    vals <- sapply(tries, function(o) if (inherits(o, "try-error")) Inf else o$value)
    pos  <- which.min(vals)
    lmd1_best <- tries[[pos]]
    if (!inherits(lmd1_best, "try-error")) {
      ok1 <- as.numeric(lmd1_best$convergence == 0)
    }
    k1 <- k1 + 1
  }

  # Linear predictor and first-step weights
  LMD1   <- as.numeric(MU_int %*% lmd1_best$par)
  w_hat1 <- w.hat.fun(LMD1, divergence = divergence, r = r)
  D2     <- sqrt(mean(w_hat1^2) - 1)

  # ---------- Internal regression (weighted or not) ----------
  ww <- if (isTRUE(w_type)) w_hat1 else rep(1, length(w_hat1))
  reg_int <- switch(
    link,
    "identity" = try(stats::lm(y_int ~ ., data = dat_int, weights = ww), silent = TRUE),
    "logit"    = try(stats::glm(y_int ~ ., data = dat_int, family = stats::binomial(link = "logit"), weights = ww), silent = TRUE),
    "probit"   = try(stats::glm(y_int ~ ., data = dat_int, family = stats::binomial(link = "probit"), weights = ww), silent = TRUE),
    "gamma"    = try(stats::glm(y_int ~ ., data = dat_int, family = stats::Gamma(link = "log"), weights = ww), silent = TRUE),
    stop("Unknown link: ", link)
  )

  # Default outputs in case regression fails
  entropy2 <- NA_real_
  if (inherits(reg_int, "try-error")) {
    res.est <- c(est = NA_real_, lmd2.est = NA_real_)
    out <- list(
      model    = divergence,
      r        = r,
      w_type   = w_type,
      Entropy2 = entropy2,
      res.w1   = reg_int,     # carry the error object
      res.w2   = reg_int,
      result   = res.est,
      w2.hat   = NA_real_,
      D1       = D1,
      D2       = D2
    )
    return(out)
  }

  # ---------- Second-step EB (KL fixed) ----------
  h    <- stats::predict(reg_int)      # internal linear predictor
  calH <- w_hat1 * h                   # H_i in the paper notation

  if (is.null(second_covariate)) {
    eta_int <- cbind(1, calH)
    eta_ext <- c(1, eta)
  } else {
    eta_int <- cbind(1, second_covariate * w_hat1, calH)
    eta_ext <- c(MU_ext[1:(ncol(second_covariate) + 1)], eta)
  }

  k2 <- 0
  ok2 <- 0
  lmd2_best <- NULL
  while ((k2 < 1000) && (ok2 == 0)) {
    tries2 <- lapply(
      1:10,
      function(i) {
        set.seed(i + 456)
        init <- stats::runif(ncol(eta_int), -0.1, 0.1)
        try(
          stats::optim(
            par    = init,
            fn     = lmd.fun(eta_int, eta_ext, divergence = "KL", r = 0),
            method = "L-BFGS-B",
            lower  = -5,
            upper  =  5
          ),
          silent = TRUE
        )
      }
    )
    vals2 <- sapply(tries2, function(o) if (inherits(o, "try-error")) Inf else o$value)
    pos2  <- which.min(vals2)
    lmd2_best <- tries2[[pos2]]
    if (!inherits(lmd2_best, "try-error")) {
      ok2 <- as.numeric(lmd2_best$convergence == 0 & sum(lmd2_best$par >= M) < 1)
    }
    k2 <- k2 + 1
  }

  # IMPORTANT: to reproduce the original implementation exactly,
  # we compute the second-step weight from FIRST-step LMD (not from LMD2):
  w2_hat <- w.hat.fun(LMD1, divergence = "KL", r = 1)
  entropy2 <- ent.fun(w2_hat, divergence = "KL", r = 1)

  est <- if ((k1 == 1000) || (k2 == 1000)) NA_real_ else mean(w2_hat * y_int)

  # Collect results
  res.est <- c(
    est       = est,
    lmd2.est  = if (!inherits(lmd2_best, "try-error")) lmd2_best$par else NA_real_
  )

  out <- list(
    model    = divergence,
    r        = r,
    w_type   = w_type,
    Entropy2 = entropy2,
    res.w1   = lmd1_best,
    res.w2   = lmd2_best,
    result   = res.est,
    w2.hat   = w2_hat,
    D1       = D1,
    D2       = D2
  )
  out
}

#' Generalized Entropy Balancing estimator
#'
#' If `auto = TRUE`, the function computes \{KL\} ∪ \{LW,QLS,TS\}×`r_set`
#' and selects the model that maximizes Entropy2. If `auto = FALSE`,
#' it runs the single specified model (`divergence`, `r`) and returns diagnostics.
#'
#' The second step (integration of outcome) is fixed to KL to reproduce
#' the original implementation. For numerical compatibility, the second-step
#' weight is computed from the FIRST-step linear predictor (as in the script).
#'
#' @param dat_int data.frame; must contain a numeric column `y_int` and predictors
#' @param MU_int  n x p matrix; first column must be 1 (intercept)
#' @param MU_ext  length-p target moments; MU_ext[1] must be 1 (intercept)
#' @param eta     scalar; external target for step-2 (e.g., external mean outcome)
#' @param divergence character; one of "KL","LW","QLS","TS" (used when auto=FALSE)
#' @param r       numeric; tuning parameter for LW/QLS/TS (ignored for KL; used when auto=FALSE)
#' @param w_type  logical; if TRUE, weight internal regression by first-step w
#' @param second_covariate optional matrix used in step-2 constraints (default NULL)
#' @param non_regression kept for API compatibility; not used here (eta is supplied)
#' @param link    "identity" | "logit" | "probit" | "gamma"
#' @param M       numeric; upper bound threshold for step-2 lambda
#' @param auto    logical; if TRUE, run model selection over KL and {LW,QLS,TS}×r_set
#' @param r_set   numeric vector of candidate r values for LW/QLS/TS (when auto=TRUE)
#' @return a list with class "daisy_eb" containing (for single run):
#'   \item{model, r, w_type}{first-step model spec}
#'   \item{Entropy2}{step-2 entropy value}
#'   \item{result}{named numeric; includes `est`}
#'   \item{D1, D2}{diagnostics}
#'   If auto=TRUE, it also contains:
#'   \item{best_model}{list(divergence, r)}
#'   \item{leaderboard}{data.frame of all candidates (model, r, Entropy2, D1, D2)}
#'   \item{all_results}{list of per-candidate result objects}
#' @export
#' @importFrom stats cov
EB_est <- function(
    dat_int,
    MU_int,
    MU_ext,
    eta,
    divergence       = "KL",
    r                = 1,
    w_type           = FALSE,
    second_covariate = NULL,
    non_regression   = TRUE,
    link             = "identity",
    M                = 10,
    auto             = FALSE,
    r_set            = c(0.01, 0.1, 0.5, 1)
) {
  if (!isTRUE(auto)) {
    # ---- single model run ----
    ans <- EB_est_one(
      dat_int = dat_int, MU_int = MU_int, MU_ext = MU_ext, eta = eta,
      divergence = divergence, r = r,
      w_type = w_type, second_covariate = second_covariate,
      link = link, M = M
    )
    ans$auto <- FALSE
    class(ans) <- unique(c("daisy_eb", class(ans)))
    return(ans)
  }

  # ---- auto model search: KL ∪ {LW,QLS,TS} × r_set ----
  cand <- list()
  # KL (r ignored)
  cand[["KL"]] <- list(model = "KL", r = NA_real_)
  for (d in c("LW", "QLS", "TS")) {
    for (rr in r_set) {
      key <- paste0(d, "_", rr)
      cand[[key]] <- list(model = d, r = rr)
    }
  }

  # Evaluate candidates
  all_results <- lapply(
    cand,
    function(sp) {
      EB_est_one(
        dat_int = dat_int, MU_int = MU_int, MU_ext = MU_ext, eta = eta,
        divergence = sp$model, r = if (is.na(sp$r)) 1 else sp$r,
        w_type = w_type, second_covariate = second_covariate,
        link = link, M = M
      )
    }
  )

  # Build leaderboard (include D1 & D2 per candidate)
  lb <- do.call(
    rbind,
    Map(
      function(k, res) {
        data.frame(
          model    = res$model,
          r        = if (identical(res$model, "KL")) NA_real_ else res$r,
          Entropy2 = res$Entropy2,
          D1       = res$D1,
          D2       = res$D2,
          stringsAsFactors = FALSE,
          row.names = k
        )
      },
      names(all_results),
      all_results
    )
  )
  lb <- lb[order(lb$Entropy2, decreasing = TRUE), , drop = FALSE]

  # Select best by Entropy2 (largest)
  best_i <- which.max(lb$Entropy2)
  best_model <- list(divergence = lb$model[best_i], r = lb$r[best_i])

  # Return best result with metadata
  best_res <- all_results[[best_i]]
  best_res$best_model  <- best_model
  best_res$leaderboard <- lb
  best_res$all_results <- all_results
  best_res$auto        <- TRUE

  class(best_res) <- unique(c("daisy_eb", class(best_res)))
  best_res
}
