

# spdgp

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/spdgp.png)](https://CRAN.R-project.org/package=spdgp)
[![CRAN
Release](https://www.r-pkg.org/badges/last-release/spdgp.png)](https://CRAN.R-project.org/package=spdgp)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-20b2aa.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/JosiahParry/spdgp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JosiahParry/spdgp/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

*Spatial Data Generation Processes*

`{spdgp}` is an R port of the [`pysal`](https://pysal.org/) module
[`dgp`](https://pysal.org/spreg/api.html#dgp) within
[`spreg`](https://pysal.org/spreg/) library.

`spdgp` is designed around the
[`spdep`](https://r-spatial.github.io/spdep/) package’s `listw` object
for representing spatial weights matrices.

## Overview

Use `spdgp` to generate data for the following models:

- OLS: `sim_ols()`
- Spatial Error Model (SEM): `sim_sem()`
- Spatial Lag Model (SAR): `sim_sar()`
- Spatially Lagged X Model (SLX): `sim_slx()`
- Spatial Lagged X Error Model (SLX Error): `sim_slx_error()`
- Spatial Autoregressive Model with Autoregressive Errors (SARAR / SAC /
  “Combo” model): `sim_sarar()`
- Spatial Durbin Model: `sim_durbin()`
- General Nested Model (GNM): `sim_gns()`
- Matrix Exponential Spatial Lag Model (MESS): `sim_mess()`

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=spdgp) with:

``` r
install.packages("spdgp")
```

- Install from source code on
  [GitHub](https://github.com/JosiahParry/spdgp) with:

``` r
if (!requireNamespace("pak")) {
    install.packages("pak")
}

pak::pak("josiahparry/spdgp")
```

## Basic usage:

We first need to create a spatial weights matrix to simulate based off
of:

``` r
library(spdgp)

set.seed(42)

n <- 50
listw <- sim_grid_listw(10, 5)
```

Next we can simulate our error term, x from our betas.

``` r
# simulate error 
u <- make_error(n, method = "normal")

# simulate x values based on uniform distribution
x <- make_x(n, method = "uniform")

# create x's according to an intercept and beta value
x_beta <- make_xb(x, c(1, 5))
```

Next, we’ll simulate the y and specify the autoregrssive parameter
$\rho = 0.5$.

``` r
# simulate y from error and the x_beta
y <- sim_sar(u, x_beta, listw, rho = 0.5)
```

Fit an SAR model using simulated data.

``` r
library(spatialreg)

sar_mod <- lagsarlm(y ~ x$x_1, listw = listw)

summary(sar_mod)
```

    #> 
    #> Call:lagsarlm(formula = y ~ x$x_1, listw = listw)
    #> 
    #> Residuals:
    #>       Min        1Q    Median        3Q       Max 
    #> -2.531747 -0.611036 -0.043396  0.739112  2.200584 
    #> 
    #> Type: lag 
    #> Coefficients: (asymptotic standard errors) 
    #>             Estimate Std. Error z value Pr(>|z|)
    #> (Intercept)  0.86583    1.30543  0.6632   0.5072
    #> x$x_1        5.13379    0.16814 30.5326   <2e-16
    #> 
    #> Rho: 0.49234, LR test value: 38.315, p-value: 6.0202e-10
    #> Asymptotic standard error: 0.059849
    #>     z-value: 8.2265, p-value: 2.2204e-16
    #> Wald statistic: 67.675, p-value: 2.2204e-16
    #> 
    #> Log likelihood: -78.41619 for lag model
    #> ML residual variance (sigma squared): 1.2853, (sigma: 1.1337)
    #> Number of observations: 50 
    #> Number of parameters estimated: 4 
    #> AIC: 164.83, (AIC for lm: 201.15)
    #> LM test for residual autocorrelation
    #> test value: 0.11376, p-value: 0.7359

In the model we can see that the estimate of `rho` is quite close to the
specified value of `0.5`.
