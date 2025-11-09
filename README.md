
<!-- README.md is generated from README.Rmd. Please edit that file -->

# daisy

**daisy** provides estimation and bootstrap tools for **Data-Adaptive
Integration with Summary Data** using **generalized entropy balancing
(GEB)**.

Main functions

- `EB_est()` — estimator (supports `auto = TRUE` model selection)  
- `EB_bootstrap_var()` — bootstrap variance (no auto inside bootstrap;
  specify the model explicitly)

## Installation

``` r
# using pak (recommended)
# install.packages("pak")
pak::pak("KMorikawaISU/daisy")

# or using remotes
# install.packages("remotes")
remotes::install_github("KMorikawaISU/daisy")
```

## Quick start (paper-style toy data)

> Simplest setting in the paper: **x1 ~ N(0,1)**, **x2 ~ Bern(0.5)**,  
> \*\*y \| x ~ N(0.5\*x1 - x2, 1)**, internal **n = 200**, external
> summaries **n1 = 2000**.  
> To keep the README light, this chunk is **not executed\*\*
> (`eval = FALSE`).

``` r
library(daisy)

## -------------------------------
## 1) Generate toy data (simple)
## -------------------------------
set.seed(1)
n  <- 200
n1 <- 2000

# Internal (no labels)
x1_int <- rnorm(n, 0, 1)
x2_int <- rbinom(n, 1, 0.5)
y_int  <- rnorm(n, 0.5 * x1_int - x2_int, 1)

# External (summary only; same DGP, mean of x1 = 0)
x1_ext <- rnorm(n1, 0, 1)
x2_ext <- rbinom(n1, 1, 0.5)
y_ext  <- rnorm(n1, 0.5 * x1_ext - x2_ext, 1)

# dat_int: data.frame with y_int and predictors
dat_int <- data.frame(
  x_int.1 = x1_int,
  x_int.2 = x2_int,
  y_int   = y_int
)

# MU_int: balance over (1, x1, x1^2, x2)
MU_int <- cbind(1, x1_int, x1_int^2, x2_int)

# MU_ext: target external moments (1, E[x1], E[x1^2], E[x2])
MU_ext <- c(1, mean(x1_ext), mean(x1_ext^2), mean(x2_ext))

# eta: external target for step-2 (mean of outcome)
eta <- mean(y_ext)

## -------------------------------
## 2) Estimation (auto model search)
## -------------------------------
fit <- EB_est(
  dat_int, MU_int, MU_ext, eta,
  auto  = TRUE,
  r_set = c(0.01, 0.1, 0.5, 1),
  link  = "identity"    # or "logit", "probit", "gamma"
)

fit$best_model      # list(divergence=..., r=...)
fit$result["est"]   # point estimate
fit$Entropy2

## -------------------------------
## 3) Bootstrap variance (no auto)
##    fix the selected model from step 2
## -------------------------------
bt <- EB_bootstrap_var(
  dat_int, MU_int, MU_ext, eta,
  n_ext = n1,                    
  BB    = 200,
  divergence    = fit$best_model$divergence,
  r             = fit$best_model$r,
  external.boot = TRUE           # TRUE: resample external summaries; FALSE: fix them
)

bt$bootstrap_se
bt$ci_percentile_95
```

## Notes

- Step-2 uses KL (fixed). For numerical compatibility with the original
  implementation, the second-step weights are built from the
  **first-step** linear predictor.
