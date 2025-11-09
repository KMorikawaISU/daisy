#' Bootstrap Variance (Appendix B; no auto inside bootstrap)
#'
#' @inheritParams EB_est
#' @param n_ext External sample size (n1).
#' @param BB Bootstrap repetitions (B1=B2=BB).
#' @param external.boot TRUE resamples external summaries; FALSE fixes them.
#' @param seed Optional RNG seed.
#' @return List: point_estimate, bootstrap_se/var, 95% CIs, SigmaW_hat, theta_boot.
#' @export
EB_bootstrap_var <- function(
    dat_int, MU_int, MU_ext, eta,
    n_ext, BB = 200,
    divergence = "KL", r = 1,
    w_type = FALSE, second_covariate = NULL,
    link = "identity", M = 10,
    external.boot = TRUE, seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(is.matrix(MU_int), is.numeric(MU_ext), ncol(MU_int) == length(MU_ext),
            is.data.frame(dat_int), "y_int" %in% names(dat_int))

  n <- nrow(MU_int); p <- ncol(MU_int); dW <- (p - 1) + 1

  point <- EB_est(dat_int, MU_int, MU_ext, eta,
                  divergence=divergence, r=r,
                  w_type=w_type, second_covariate=second_covariate,
                  link=link, M=M, auto=FALSE)
  theta_hat <- as.numeric(point$result["est"])

  W_mat <- matrix(NA_real_, nrow = BB, ncol = dW)
  colnames(W_mat) <- c(paste0("mu_x[", 1:(p - 1), "]"), "eta")
  for (b1 in 1:BB) {
    idx <- sample.int(n, n, TRUE)
    MU_b  <- MU_int[idx,,drop=FALSE]
    di_b  <- dat_int[idx,,drop=FALSE]
    W_mat[b1,] <- c(colMeans(MU_b[,-1,drop=FALSE]), mean(di_b$y_int))
  }
  SigmaW_hat <- stats::cov(W_mat, use = "complete.obs")
  eig <- eigen(SigmaW_hat, symmetric = TRUE, only.values = TRUE)$values
  if (any(eig < .Machine$double.eps)) SigmaW_hat <- SigmaW_hat + diag(1e-8, nrow(SigmaW_hat))

  mu_ext_obs <- as.numeric(MU_ext[-1]); eta_ext_obs <- as.numeric(eta)
  mean_vec <- c(mu_ext_obs, eta_ext_obs)
  Sigma_scaled <- (n / n_ext) * SigmaW_hat

  theta_boot <- numeric(BB)
  .rmvnorm_ <- function(n, mean, sigma) {
    d <- length(mean)
    if (requireNamespace("MASS", quietly = TRUE)) {
      return(MASS::mvrnorm(n = n, mu = mean, Sigma = sigma))
    } else if (requireNamespace("mvtnorm", quietly = TRUE)) {
      return(mvtnorm::rmvnorm(n = n, mean = mean, sigma = sigma))
    } else {
      ev <- eigen(sigma, symmetric = TRUE); ev$values[ev$values < 0] <- 0
      L <- ev$vectors %*% diag(sqrt(ev$values), d, d)
      Z <- matrix(stats::rnorm(n * d), n, d)
      return(sweep(Z %*% t(L), 2, mean, FUN = "+"))
    }
  }

  for (b2 in 1:BB) {
    idx2 <- sample.int(n, n, TRUE)
    MU_b2 <- MU_int[idx2,,drop=FALSE]
    di_b2 <- dat_int[idx2,,drop=FALSE]

    if (isTRUE(external.boot)) {
      draw <- as.numeric(.rmvnorm_(1, mean=mean_vec, sigma=Sigma_scaled))
      mu_ext_b2 <- draw[1:(p-1)]; eta_ext_b2 <- draw[p]
    } else {
      mu_ext_b2 <- mu_ext_obs;    eta_ext_b2 <- eta_ext_obs
    }
    MU_ext_b2 <- c(1, mu_ext_b2)

    res <- EB_est(di_b2, MU_b2, MU_ext_b2, eta_ext_b2,
                  divergence=divergence, r=r,
                  w_type=w_type, second_covariate=second_covariate,
                  link=link, M=M, auto=FALSE)
    theta_boot[b2] <- as.numeric(res$result["est"])
  }

  boot_var <- stats::var(theta_boot, na.rm=TRUE)
  boot_se  <- sqrt(boot_var)
  ci_norm  <- c(lower = theta_hat - 1.96*boot_se,
                upper = theta_hat + 1.96*boot_se)
  qs <- stats::quantile(theta_boot, c(0.025,0.975), na.rm=TRUE, names=FALSE)
  ci_perc <- c(lower=qs[1], upper=qs[2])

  list(point_estimate=theta_hat, bootstrap_se=boot_se, bootstrap_var=boot_var,
       ci_normal_95=ci_norm, ci_percentile_95=ci_perc,
       B1=BB, B2=BB, n_int=n, n_ext=n_ext,
       SigmaW_hat=SigmaW_hat, theta_boot=theta_boot,
       call_divergence=divergence, call_r=r, external.boot=external.boot)
}
