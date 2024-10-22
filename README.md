# Spatial Data Generation Processes

`{spdgp}` is an R port of the [`pysal`](https://pysal.org/) module
[`dgp`](https://pysal.org/spreg/api.html#dgp) within
[`spreg`](https://pysal.org/spreg) library.

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

`spdgp` can be installed from github using:

``` r
if (!requireNamespace("pak")) {
    install.packages("pak")
}

pak::pak("josiahparry/spgdp")
```

## Basic usage:

We first need to create a spatial weights matrix to simulate based off
of:

``` r
library(spdgp)

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
    #> -3.100635 -0.377549  0.042234  0.703459  2.810884 
    #> 
    #> Type: lag 
    #> Coefficients: (asymptotic standard errors) 
    #>             Estimate Std. Error z value Pr(>|z|)
    #> (Intercept)  0.49641    0.97505  0.5091   0.6107
    #> x$x_1        4.98587    0.19082 26.1281   <2e-16
    #> 
    #> Rho: 0.52947, LR test value: 54.755, p-value: 1.3656e-13
    #> Asymptotic standard error: 0.050408
    #>     z-value: 10.504, p-value: < 2.22e-16
    #> Wald statistic: 110.33, p-value: < 2.22e-16
    #> 
    #> Log likelihood: -83.83562 for lag model
    #> ML residual variance (sigma squared): 1.5825, (sigma: 1.258)
    #> Number of observations: 50 
    #> Number of parameters estimated: 4 
    #> AIC: 175.67, (AIC for lm: 228.43)
    #> LM test for residual autocorrelation
    #> test value: 1.8893, p-value: 0.16928

In the model we can see that the estimate of `rho` is quite close to the
specified value of `0.5`.
