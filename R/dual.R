# R/dual.R
# Dual objective and gradient for generalized entropy balancing (GEB).
# NOTE: X and MU MUST NOT include an intercept column/component. If an intercept
# is needed (e.g., in step-2), it is added explicitly within that step.

#' Dual objective factory for generalized EB
#'
#' @param X  numeric matrix n x p of balancing features (no intercept column).
#' @param MU numeric length-p vector of target feature means (no intercept component).
#' @param divergence character in {"KL","LW","QLS","TS"}.
#' @param r numeric tuning for LW/QLS/TS (ignored for KL).
#' @return A function of `lambda` returning the dual objective:
#'   `- t(lambda) %*% MU + mean( F(X %*% lambda) )`, where F depends on the divergence.
#' @keywords internal
lmd.fun <- function(X, MU, divergence, r = 1) {
  function(lmd) {
    LMD <- as.numeric(as.matrix(X) %*% lmd)
    if (divergence == "KL") {
      return(- drop(crossprod(lmd, MU)) + mean(exp(LMD)))
    } else if (divergence == "LW") {
      # principal branch of Lambert W from {lamW}; must be called with namespace
      lw <- lamW::lambertW0(r * exp(LMD))
      return(- drop(crossprod(lmd, MU)) + mean((lw^2 + 2 * lw) / (2 * r)))
    } else if (divergence == "QLS") {
      qls <- sqrt(1 + 4 * r * exp(LMD))
      return(- drop(crossprod(lmd, MU)) + mean((qls - log((1 + qls) / 2)) / r))
    } else if (divergence == "TS") {
      # dilogarithm from {gsl}; also call with namespace
      ts <- gsl::dilog(-exp(r * LMD))
      return(- drop(crossprod(lmd, MU)) + mean(- ts / r^2 - pi^2 / (6 * r^2)))
    } else {
      stop("Unknown 'divergence'. Use 'KL', 'LW', 'QLS', or 'TS'.")
    }
  }
}

#' Dual gradient factory for generalized EB
#'
#' @param X  numeric matrix n x p of balancing features (no intercept column).
#' @param MU numeric length-p vector of target feature means (no intercept component).
#' @param divergence character in {"KL","LW","QLS","TS"}.
#' @param r numeric tuning for LW/QLS/TS (ignored for KL).
#' @return A function of `lambda` returning the gradient
#'   `E[rho(X lambda) X] - MU`, where rho depends on the divergence.
#' @keywords internal
grad.fun <- function(X, MU, divergence, r = 1) {
  n <- nrow(X)
  function(lmd) {
    LMD <- as.numeric(as.matrix(X) %*% lmd)
    if (divergence == "KL") {
      return(drop(crossprod(X, exp(LMD)) / n) - MU)
    } else if (divergence == "LW") {
      lw <- lamW::lambertW0(r * exp(LMD))
      return(drop(crossprod(X, lw / r) / n) - MU)
    } else if (divergence == "QLS") {
      qls <- sqrt(1 + 4 * r * exp(LMD))
      return(drop(crossprod(X, (qls - 1) / (2 * r)) / n) - MU)
    } else if (divergence == "TS") {
      return(drop(crossprod(X, log(1 + exp(r * LMD)) / r) / n) - MU)
    } else {
      stop("Unknown 'divergence'. Use 'KL', 'LW', 'QLS', or 'TS'.")
    }
  }
}
