#' @keywords internal
ent.fun <- function(w, divergence, r = 1) {
  if (divergence == "KL") return(-mean(w * log(w) - w))
  if (divergence == "LW") return(-mean(w * log(w) - w + r/2 * w^2))
  if (divergence == "QLS") return(-mean(w * log(w) - w + ((1 + r*w)*log(1 + r*w) - (1 + r*w))/r))
  if (divergence == "TS") return(-mean(w^2/2 + gsl::dilog(exp(-r*w))/r^2))
  stop("Unknown 'divergence'. Use 'KL','LW','QLS','TS'.")
}

#' @keywords internal
w.hat.fun <- function(eta, divergence, r = 1) {
  if (divergence == "KL") return(exp(eta))
  if (divergence == "LW") return(lamW::lambertW0(r * exp(eta)) / r)
  if (divergence == "QLS") return((sqrt(1 + 4*r*exp(eta)) - 1) / (2*r))
  if (divergence == "TS") return(log(1 + exp(r*eta)) / r)
  stop("Unknown 'divergence'.")
}
