#' @keywords internal
lmd.fun <- function(X, MU, divergence, r = 1) {
  function(lmd) {
    LMD <- c(as.matrix(X) %*% lmd)
    if (divergence == "KL")  return(-c(crossprod(lmd, MU)) + mean(exp(LMD)))
    if (divergence == "LW")  { lw <- lamW::lambertW0(r*exp(LMD)); return(-c(crossprod(lmd, MU)) + mean((lw^2 + 2*lw)/(2*r))) }
    if (divergence == "QLS") { q  <- sqrt(1 + 4*r*exp(LMD));     return(-c(crossprod(lmd, MU)) + mean((q - log((1 + q)/2))/r)) }
    if (divergence == "TS")  { ts <- gsl::dilog(-exp(r*LMD));     return(-c(crossprod(lmd, MU)) + mean(-ts/r^2 - pi^2/(6*r^2))) }
    stop("Unknown 'divergence'.")
  }
}
