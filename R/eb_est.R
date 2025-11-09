#' Entropy-Balancing Estimator (with optional auto selection)
#'
#' @param dat_int Data frame with column y_int and predictors.
#' @param MU_int Matrix (n x p) of internal moments (1 in first col).
#' @param MU_ext Length-p target external moments (first is 1).
#' @param eta Scalar external target for step-2.
#' @param divergence "KL","LW","QLS","TS".
#' @param r Tuning parameter for LW/QLS/TS (ignored for KL).
#' @param w_type Logical; TRUE → weighted second-step regression.
#' @param second_covariate Optional extra covariates for step-2 constraints.
#' @param link "identity","logit","probit","gamma".
#' @param M Threshold for step-2 Lagrange multiplier.
#' @param auto If TRUE, searches across KL and the set of (LW/QLS/TS) with r in \code{r_set},
#' by Entropy2, and picks the best.
#' @param r_set Numeric vector of candidate r's for auto=TRUE.
#' @return A list with \code{result[["est"]]} (point estimate), \code{Entropy2},
#' selected model info (when \code{auto=TRUE}), and diagnostics.
#' @export
#' @importFrom stats lm glm binomial cov optim predict quantile var
EB_est <- function(
    dat_int, MU_int, MU_ext, eta,
    divergence = "KL", r = 1,
    w_type = FALSE, second_covariate = NULL,
    link = "identity", M = 10,
    auto = FALSE, r_set = c(0.01, 0.1, 0.5, 1)
) {
  stopifnot(is.matrix(MU_int), is.numeric(MU_ext), ncol(MU_int) == length(MU_ext),
            is.data.frame(dat_int), "y_int" %in% names(dat_int))

  if (!isTRUE(auto)) {
    return(
      EB_est_one(dat_int, MU_int, MU_ext, eta,
                 divergence, r, w_type, second_covariate, link, M)
    )
  }

  cand <- list(list(div="KL", r=1))
  for (dv in c("LW","QLS","TS")) for (rr in r_set) cand[[length(cand)+1]] <- list(div=dv, r=rr)

  all_results <- list(); ent_vec <- numeric(length(cand)); keys <- character(length(cand))
  for (i in seq_along(cand)) {
    dv <- cand[[i]]$div; rr <- cand[[i]]$r
    key <- if (dv=="KL") "KL" else paste0(dv, "_r=", format(rr, trim=TRUE, scientific=FALSE))
    keys[i] <- key
    res <- EB_est_one(dat_int, MU_int, MU_ext, eta, dv, rr, w_type, second_covariate, link, M)
    all_results[[key]] <- res
    ent_vec[i] <- ifelse(is.na(res$Entropy2), -Inf, res$Entropy2)
  }

  best_idx <- which.max(ent_vec)
  best <- all_results[[keys[best_idx]]]
  leaderboard <- data.frame(candidate=keys,
                            divergence=vapply(cand, function(z) z$div, ""),
                            r=vapply(cand, function(z) z$r, 0),
                            Entropy2=ent_vec, stringsAsFactors = FALSE)
  leaderboard <- leaderboard[order(-leaderboard$Entropy2), ]

  best$best_model  <- list(divergence = cand[[best_idx]]$div, r = cand[[best_idx]]$r)
  best$leaderboard <- leaderboard
  best$all_results <- all_results
  best$auto        <- TRUE
  best
}

#' @keywords internal
EB_est_one <- function(
    dat_int, MU_int, MU_ext, eta,
    divergence = "KL", r = 1,
    w_type = FALSE, second_covariate = NULL,
    link = "identity", M = 10
) {
  y_int <- dat_int$y_int

  mu_diff <- drop(colMeans(MU_int) - MU_ext)
  Sx <- stats::cov(MU_int[, -1, drop = FALSE])
  eig <- eigen(Sx, symmetric = TRUE, only.values = TRUE)$values
  if (any(eig < .Machine$double.eps)) Sx <- Sx + diag(1e-8, nrow(Sx))
  D1 <- c(sqrt(t(mu_diff[-1]) %*% solve(Sx) %*% mu_diff[-1]))

  k1 <- 0; ok1 <- 0
  while ((k1 < 1000) & (ok1 == 0)) {
    lmd.est00 <- lapply(1:10, function(.) try(
      stats::optim(stats::runif(ncol(MU_int), -0.1, 0.1),
                   lmd.fun(MU_int, MU_ext, divergence, r),
                   method = "BFGS"), silent = TRUE))
    vals <- sapply(lmd.est00, function(o) if (inherits(o,"try-error")) Inf else o$value)
    lmd.est0 <- lmd.est00[[which.min(vals)]]
    if (inherits(lmd.est0,"try-error")) { k1 <- k1 + 1; next }
    ok1 <- as.numeric(lmd.est0$convergence == 0)
    k1  <- k1 + 1
  }
  LMD.est <- c(MU_int %*% lmd.est0$par)
  w.hat   <- w.hat.fun(LMD.est, divergence, r)
  D2      <- sqrt(mean(w.hat^2) - 1)

  ww <- if (isFALSE(w_type)) rep(1, length(w.hat)) else w.hat
  reg_int <- switch(link,
                    identity = try(stats::lm(y_int ~ ., data = dat_int, weights = ww), silent = TRUE),
                    logit    = try(stats::glm(y_int ~ ., data = dat_int, family = stats::binomial("logit"),  weights = ww), silent = TRUE),
                    probit   = try(stats::glm(y_int ~ ., data = dat_int, family = stats::binomial("probit"), weights = ww), silent = TRUE),
                    gamma    = try(stats::glm(y_int ~ ., data = dat_int, family = stats::Gamma(link = "log"), weights = ww), silent = TRUE),
                    stop("Unknown link.")
  )
  entropy2 <- NA_real_
  if (!inherits(reg_int, "try-error")) {
    h <- stats::predict(reg_int); calH <- w.hat * h
    if (is.null(second_covariate)) {
      eta_int <- as.matrix(cbind(1, calH));  eta_ext_vec <- c(1, eta)
    } else {
      eta_int <- as.matrix(cbind(1, second_covariate * w.hat, calH))
      eta_ext_vec <- c(MU_ext[1:(ncol(second_covariate) + 1)], eta)
    }

    k2 <- 0; ok2 <- 0
    while ((k2 < 1000) & (ok2 == 0)) {
      lmd2.est00 <- lapply(1:10, function(.) try(
        stats::optim(stats::runif(ncol(eta_int), -0.1, 0.1),
                     lmd.fun(eta_int, eta_ext_vec, divergence = "KL", r = 0),
                     method = "L-BFGS-B", lower = -5, upper = 5), silent = TRUE))
      vals2 <- sapply(lmd2.est00, function(o) if (inherits(o,"try-error")) Inf else o$value)
      lmd2.est0 <- lmd2.est00[[which.min(vals2)]]
      if (inherits(lmd2.est0,"try-error")) { k2 <- k2 + 1; next }
      ok2 <- as.numeric(lmd2.est0$convergence == 0 & sum(lmd2.est0$par >= M) < 1)
      k2  <- k2 + 1
    }

    # 互換仕様：第二段の重みは第一段の LMD から作成
    w2.hat   <- w.hat.fun(LMD.est, divergence = "KL", r = 1)
    entropy2 <- ent.fun(w2.hat, divergence = "KL", r = 1)
    est      <- if ((k1==1000) | (k2==1000)) NA_real_ else mean(w2.hat * y_int)

    return(list(model=divergence, r=r, w_type=w_type,
                Entropy2=entropy2, res.w1=lmd.est0, res.w2=lmd2.est0,
                result=c(est=est, lmd2.est=lmd2.est0$par), w2.hat=w2.hat,
                D1=D1, D2=D2))
  } else {
    return(list(model=divergence, r=r, w_type=w_type,
                Entropy2=entropy2, res.w1=reg_int, res.w2=reg_int,
                result=c(est=NA_real_, lmd2.est=NA_real_), w2.hat=NA,
                D1=D1, D2=D2))
  }
}
