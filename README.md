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

# simulate y from error and the x_beta
y <- sim_sar(u, x_beta, listw)
```

Fit an SAR model using simulated data.

``` r
library(spatialreg)

sar_mod <- lagsarlm(y ~ x$x_1, listw = listw)

summary(sar_mod)
```


    Call:lagsarlm(formula = y ~ x$x_1, listw = listw)

    Residuals:
          Min        1Q    Median        3Q       Max 
    -2.756735 -0.617262  0.058604  0.642419  1.916998 

    Type: lag 
    Coefficients: (asymptotic standard errors) 
                Estimate Std. Error z value Pr(>|z|)
    (Intercept) -0.92650    1.29918 -0.7131   0.4758
    x$x_1        5.06428    0.15051 33.6475   <2e-16

    Rho: 0.57344, LR test value: 52.782, p-value: 3.7259e-13
    Asymptotic standard error: 0.057748
        z-value: 9.9301, p-value: < 2.22e-16
    Wald statistic: 98.607, p-value: < 2.22e-16

    Log likelihood: -70.13125 for lag model
    ML residual variance (sigma squared): 0.90401, (sigma: 0.9508)
    Number of observations: 50 
    Number of parameters estimated: 4 
    AIC: 148.26, (AIC for lm: 199.05)
    LM test for residual autocorrelation
    test value: 1.5878, p-value: 0.20764
