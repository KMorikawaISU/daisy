# R/entropy.R
# Special functions and links used by the GEB estimator.
# NOTE: All comments and docs are in English. No intercept should be included
# in MU_int / MU_ext. Intercepts are handled internally where needed.

#' Primal entropy objective: returns -mean G(w)
#'
#' @param w numeric vector of weights.
#' @param divergence character in {"KL","LW","QLS","TS"}.
#' @param r numeric tuning for LW/QLS/TS (ignored for KL).
#' @keywords internal
ent.fun <- function(w, divergence, r = 1) {
  if (divergence == "KL") {
    # G(w) = w log w - w
    return(-mean(w * log(w) - w))
  } else if (divergence == "LW") {
    # G_L(w) = (w log w - w) + (r/2) w^2
    return(-mean(w * log(w) - w + (r / 2) * w^2))
  } else if (divergence == "QLS") {
    # G_Q(w) = (w log w - w) + ( (1+rw) log(1+rw) - (1+rw) ) / r
    return(-mean(w * log(w) - w + ((1 + r * w) * log(1 + r * w) - (1 + r * w)) / r))
  } else if (divergence == "TS") {
    # G_S(w) = (1/2) w^2 + (1/r^2) Li2( e^{- r w} )
    return(-mean(w^2 / 2 + gsl::dilog(exp(- r * w)) / r^2))
  } else {
    stop("Unknown 'divergence'. Use 'KL', 'LW', 'QLS', or 'TS'.")
  }
}

#' Link rho: map linear predictor -> primal weight
#'
#' @param eta numeric vector (linear predictor).
#' @param divergence character in {"KL","LW","QLS","TS"}.
#' @param r numeric tuning for LW/QLS/TS (ignored for KL).
#' @keywords internal
w.hat.fun <- function(eta, divergence, r = 1) {
  if (divergence == "KL") {
    return(exp(eta))
  } else if (divergence == "LW") {
    # principal branch of Lambert W from {lamW}; must be called with namespace
    lw <- lamW::lambertW0(r * exp(eta))
    return(lw / r)
  } else if (divergence == "QLS") {
    qls <- sqrt(1 + 4 * r * exp(eta))
    return((qls - 1) / (2 * r))
  } else if (divergence == "TS") {
    return(log(1 + exp(r * eta)) / r)
  } else {
    stop("Unknown 'divergence'. Use 'KL', 'LW', 'QLS', or 'TS'.")
  }
}
