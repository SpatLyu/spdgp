#' Generate spatial weights matrix for a grid
#' 
#' Create a spatial weights matrix based on a square grid structure.
#' 
#' @param nrow the number of rows in the grid.
#' @param ncol defaults to `nrow`. The number columns in the grid.
#' @param style the spatial weights style. Defaults to row standardized. See [spdep::nb2listw()] for more.
#' @param type default `"queen"`. Can also be `"rook"`.
#' @export
#' @examples
#' sim_grid_listw(10, 5)
sim_grid_listw <- function(nrow, ncol = nrow, style = "W", type = c("queen", "rook")) {
  check_number_whole(nrow)
  check_number_whole(ncol)
  type <- rlang::arg_match(type)
  spdep::nb2listw(spdep::cell2nb(nrow, ncol, type), style = style)
}

#' Create a square grid
#' 
#' Creates a square grid with `ncol` and `nrow` dimensions.
#' 
#' @inheritParams sim_grid_listw
#' @export
#' @examples
#' make_square_grid(3, 2)
make_square_grid <- function(nrow, ncol = nrow) {
  sf::st_make_grid(
    cellsize = c(1, 1),
    offset = c(0, 0),
    n = c(ncol, nrow)
  )
}

#' Simulate an error term
#' 
#' @param n the number of values to simulate
#' @param mu the sample average.
#' @param var the sample variance. The `sqrt(var)` is passed to `rnorm()` and `rlnorm()` for normal and laplace distributions. `sqrt(var / 2)` is used for `laplace()`
#' @param method must be one of `"normal"`, `"laplace"`, `"cauchy"`, or `"lognormal"`.
#' 
#' @details
#' 
#' - `"normal"`: fit with `rnorm()`
#' - `"laplace"`: fit with `smoothmest::rdoublex()`
#' - `"cauchy"`: fit with `rcauchy()`
#' - `"lognormal"`: fit with `rlnorm()`
#' 
#' @export
make_error <- function(
    n = 10,
    mu = 0,
    var = 1,
    method = c("normal", "laplace", "cauchy", "lognormal")
) {

  method <- match.arg(method, method, FALSE)
  switch(
    method, 
    "normal" = rnorm(n, mu, sqrt(var)),
    "laplace" = smoothmest::rdoublex(n, mu, sqrt(var / 2)),
    "cauchy" = rcauchy(n),
    "lognormal" = rlnorm(n, mu, sqrt(var)),
    stop("method not supported")
  )
}

#' @param cor correlation between bivariate normal
#' @inheritParams make_error
#' @rdname make_x
#' @export
make_x_bivariate <- function(n = 5, mu = 1, cor = 0.25, var = c(1, 1)) {
  # check the length of mu
  if (!length(mu) %in% 1:2) {
    cli:::cli_abort("Bivariate normal data requires 2 values of {.arg mu}")
  } 
  if (rlang::is_scalar_atomic(mu)) {
    mu <- rep(mu, 2)
  } 

  if (!length(var) %in% 1:2) {
    cli:::cli_abort("Bivariate normal data requires 2 values of {.arg var}")
  }
  
  if (rlang::is_scalar_atomic(var)) {
    var <- rep(var, 2)
  } 

  ucov <- cor * sqrt(var[1] * var[2])

  mcov <- matrix(
    c(c(var[1], ucov), c(ucov, var[2])),
    byrow = TRUE,
    ncol = 2
  )
  MASS::mvrnorm(n, mu, mcov)
}


#' @rdname make_x
#' @export
make_x_uniform <- function(n = 5, var = 1) {
  check_number_decimal(var)
  check_number_whole(n)
  runif(n, 0, sqrt(12 * var))
}

#' @rdname make_x
#' @export
make_x_normal <- function(n = 5, mu = 0, var = 1) {
  check_number_decimal(var)
  check_number_whole(n)
  check_number_decimal(mu)
  rnorm(n, mu, var)
}

#' Simulate X variables
#' 
#' Simulates independent variables. 
#' 
#' @param method must be one of `"uniform"` (default), `"normal"`, or `"bivnormal"` (bivariate normal)
#' @rdname make_x
#' @export
#' @examples
#' make_x(10, mu = c(0.5, 1.2), var = c(1, 0.5)) 
make_x <- function(n = 5, mu = 0, var = 1, cor = 0, method = c("uniform", "normal", "bivnormal")) {
  check_number_whole(n)

  # check the method 
  method <- rlang::arg_match0(method, c("uniform", "normal", "bivnormal"))

  if (!rlang::is_bare_numeric(mu)) {
    rlang::abort("`mu` must be a numeric vector")
  }

  if (!rlang::is_bare_numeric(var)) {
    rlang::abort("`var` must be a numeric vector")
  }

  if (!rlang::is_bare_numeric(cor)) {
    rlang::abort("`var` must be a numeric vector")
  }

  fmls <- rlang::fn_fmls()
  param_sizes <- c(length(mu), length(var))
  # determine how many columns
  k <- max(param_sizes)

  if (any(k == 0) || !all(param_sizes == 1 | param_sizes == k) ) {
    rlang::abort("`mu` and `var` must either be scalars or the same length")
  }

  if (identical(method, "bivnormal")) {
    res <- as.data.frame(make_x_bivariate(n, mu, cor, var))
    colnames(res) <- c("x_1", "x_2")
    return(res)
  }

  res <- switch(
    method,
    uniform = Map(make_x_uniform, n, var),
    normal = Map(make_x_normal, n, mu, var)
  )

  res <- as.data.frame(res) 
  colnames(res) <- paste0("x_", 1:k)
  res
}

#' Create Spatial Lags of variables
#' 
#' Given a dataframe of numeric values and a spatial weights matrix, calculate the spatial lag of each variable.
#' 
#' @returns 
#' 
#' A `data.frame` of the spatially lagged variables.
#' 
#' @param order unused. 
#' @export
#' @examples
#' listw <- sim_grid_listw(10, 10)
#' x_vars <- make_x(100, mu = c(0.5, 1.2), var = c(1, 0.5)) 
#' res <- make_wx(x_vars, listw)
#' head(res)
make_wx <- function(x, listw, order = NULL) {
  lagged <- lapply(x, function(.x) spdep::lag.listw(listw, .x))
  names(lagged) <- paste0(names(lagged), "_lag")
  data.frame(lagged)
}


#' Calculate predicted X values based on coefficients
#' 
#' This function calculates predicted x values based on regression coefficients.
#' The results of this function can be passed to other `sim_*()` functions. 
#' 
#' @param x a data frame of x variables generated with `make_x()`
#' @param beta a vector of the beta coefficients for each of the variables. There must be `ncol(x) + 1` values. The first element of the vector is the intercept.
#' 
#' @examples
#' grid <- make_square_grid(5)
#' listw <- sim_grid_listw(5, 5)
#' n <- 25
#' x <- make_x(25, c(0,1), c(1,4))
#' 
#' betas <- c(1, 1.5, -2)
#' make_xb(x, betas)
#' @export
make_xb <- function(x, beta) {
  n <- nrow(x)
  k <- ncol(x)
  
  if ((k + 1) != length(beta)) {
    stop("Incompatible dimensions")
  } else {
    b <- matrix(beta, ncol = 1)
    x1 <- cbind(1, as.matrix(x)) # Add a column of ones for the constant term
    xb <- x1 %*% b
    return(as.numeric(xb))
  }
}




#' Calculate the effect of spatially lagged X variables
#' 
#' This function computes the contribution of spatially lagged X variables based on 
#' provided coefficients. The function takes the spatially lagged variables (`wx`, see [make_wx()])
#' and multiplies them by their corresponding regression coefficients (`gamma`), returning
#' the predicted influence of the spatial lags. Only spatial lags are considered; 
#' the original X variables are not included in this calculation.
#'
#' @param wx a matrix of spatially lagged x variables.
#' @param gamma a vector of coefficients for the spatially lagged x variables. Its length must match the number of columns in wx.
#' 
#' @examples
#' grid <- make_square_grid(5)
#' listw <- spdep::nb2listw(spdep::poly2nb(grid))
#' n <- 25
#' x <- make_x(25, c(0,1), c(1,4))
#' wx <- make_wx(x, listw)
#' gamma <- c(1.75, 0.4)
#' make_wxg(wx, gamma) 
#' @export
make_wxg <- function(wx, gamma) {
  wx <- as.matrix(wx)  # Convert wx to a matrix
  k <- ncol(wx)
  
  if (k > 1) {
    if (k != length(gamma)) {
      stop("Incompatible dimensions")
    } else {
      g <- matrix(gamma, ncol = 1)
      wxg <- wx %*% g  # Matrix multiplication
    }
  } else {
    wxg <- wx * gamma  # Scalar multiplication if wx has 1 column
  }
  
  return(as.numeric(wxg))
}

# calculate the inverse product 
inverse_prod <- function(listw, x, scalar) {
  w <- as(listw, "CsparseMatrix")
  n <- vctrs::vec_size(x)
  Matrix::solve(Matrix::Diagonal(n) - scalar * w, as.matrix(x))
}

#' Simulate Spatial Error Process
#' 
#' This function generates a pure spatial error process, which is useful when 
#' you only want to simulate the error structure without including any deterministic
#' part (i.e., no xb term). This can be used to analyze or simulate the behavior
#' of spatially dependent errors in isolation.
#' 
#' @inheritParams sim_slx
#' @inheritParams make_wx
#' @param lambda a value value between -1 and 1. The spatial autoregressive coefficient for the error term.
#' @param model default `"sar"`. Which model should be simulated. Provide `"ma"` for the moving average.
#' 
#' @references See [`spreg.dgp.dgp_errproc`](https://pysal.org/spreg/generated/spreg.dgp.dgp_errproc.html#spreg.dgp.dgp_errproc)
#' @export
#' @examples
#' listw <- sim_grid_listw(5) 
#' u <- make_error(25)
#' sim_error(u, listw)
sim_error <- function(u, listw, lambda = 0.5, model = c("sar", "ma")) {
  check_number_decimal(lambda, min = -1, max = 1)
  if (!rlang::is_bare_numeric(u)) {
    rlang::abort("`u` must be a numeric vector")
  }
  model <- match.arg(model)
  n0 <- length(u)  # Get the number of elements in the random error vector

  if (length(listw$weights) != n0) {
    stop("Error: incompatible weights dimensions")
  }

  if (model == 'sar') {
    # Use the inverse_prod function for SAR model
    y <- inverse_prod(listw, u, lambda)
  } else if (model == 'ma') {
    # Use the lagged values for MA model
    y <- u + lambda * spdep::lag.listw(listw, u)  # Using lag.listw from spdep
  } else {
    stop("Error: unsupported model type")
  }
  as.vector(y)
}


#' Simulate OLS
#' 
#' Simulate a y variable for an Ordinary Least Squares (OLS) regression.
#' 
#' @examples
#' u <- make_error(50, method = "normal")
#' x <- make_x(50)
#' xb <- make_xb(x, c(1,2))
#' y <- sim_ols(u, xb)
#' lm(y ~ x[[1]])
#'
#' @references [`spreg.dgp.dgp_ols`](https://pysal.org/spreg/generated/spreg.dgp.dgp_ols.html#spreg.dgp.dgp_ols)
#' @export
#' @inheritParams sim_slx
sim_ols <- function(u, xb) {
  # Check if the lengths of u and xb are the same
  if (vctrs::vec_size(u) != vctrs::vec_size(xb)) {
    stop("Error: dimension mismatch")
  }
  
  # Compute y as the sum of xb and u
  as.numeric((xb + u)[[1]])
}

#' Simulate Spatially Lagged X (SLX) model
#' 
#' This function simulates the y values of an SLX model, where the dependent 
#' variable is influenced by both the original and spatially lagged x variables. 
#' 
#' @param u an error vector 
#' @param xb predicted x values as calculated by `make_xb()`
#' @param wxg predicted spatial lag effect as calculated by `make_wxg()`
#' @examples
#' ncol <- 20
#' n <- ncol^2
#' listw <- sim_grid_listw(ncol, ncol)  # Create spatial weights for a grid
#' u <- make_error(n, method = "normal")  # Simulate random errors
#' x <- make_x(n, method = "uniform")  # Generate x variables
#' xb <- make_xb(x, c(1, 2))  # Calculate xb using the original x and coefficients
#' wx <- make_wx(x, listw)  # Generate spatially lagged x variables
#' wxg <- make_wxg(wx, 0.5)  # Calculate the effect of the spatial lags
#' y <- sim_slx(u, xb, wxg)  # Simulate the SLX model outcome
#' df <- data.frame(y, x)
#' spatialreg::lmSLX(y ~ ., data = df, listw = listw)  # Estimate the SLX model
#' @export
#' @references [`spreg.dgp.dgp_slx`](https://pysal.org/spreg/generated/spreg.dgp.dgp_slx.html#spreg.dgp.dgp_slx)
sim_slx <- function(u, xb, wxg) {
  # Check if the sizes of u, xb, and wxg are the same
  if (vctrs::vec_size(u) != vctrs::vec_size(xb)) {
    stop("Error: dimension mismatch between u and xb")
  } else if (vctrs::vec_size(xb) != vctrs::vec_size(wxg)) {
    stop("Error: dimension mismatch between xb and wxg")
  }
  
  # Compute y as the sum of xb, wxg, and u
  xb + wxg + u
}


#' Simulate Spatial Error Model (SEM)
#' 
#' Simulate the y values for an SEM model.
#' 
#' @references [`spreg.dgp.dgp_sperror`](https://pysal.org/spreg/generated/spreg.dgp.dgp_sperror.html#spreg.dgp.dgp_sperror)
#' @export
#' @inheritParams sim_error 
#' @inheritParams sim_slx
#' @examples
#' ncol <- 10
#' n <- ncol^2
#' listw <- sim_grid_listw(ncol, ncol)  # Create spatial weights for a grid
#' u <- make_error(n)  # Simulate random errors
#' x <- make_x(
#'   n,
#'   mu = c(0.25, 5),
#'   var = c(1, 0.75),
#'   method = "normal"
#' )  # Generate x variables
#' 
#' # create xb with intercept = 1, beta1 = 2, beta2 = -3
#' xb <- make_xb(x, c(1, 2, -3))
#' y <- sim_sem(u, xb, listw)
#' 
#' # combine data 
#' df <- cbind(y = y, x)
#' 
#' # fit SEM model
#' # Note lambda, x_1, and x_2 estimates.
#' spatialreg::errorsarlm(y ~ ., df, listw)
sim_sem <- function(u, xb, listw, lambda = 0.5, model = c("sar", "ma")) {
  n_lw <- length(listw$neighbours)
  n_u <- vctrs::vec_size(u)
  n_xb <- vctrs::vec_size(xb)
  all_n <- c(n_lw, n_u, n_xb)
  check_number_decimal(lambda)

  model <- rlang::arg_match(model)
  
  if (!max(all_n) == min(all_n)) {
    rlang::abort("`u`, `xb`, and `listw` must have the same number of features")
  }

  u1 <- switch(
    model,
    "sar" = inverse_prod(listw, u, lambda),
    "ma" = u + lambda * spdep::lag.listw(listw, u)
  )
  as.numeric(xb + u1)
}

#' Simulate Spatially Lagged X Error Model
#' @references [`spreg.dgp.dgp_slxerror`](https://pysal.org/spreg/generated/spreg.dgp.dgp_slxerror.html#spreg.dgp.dgp_slxerror)
#' @export
#' @inheritParams sim_sem
#' @inheritParams sim_slx
sim_slx_error <- function(u, xb, wxg, lambda = 0.5, model = c("sar", "ma")) {
  n_lw <- length(listw$neighbours)
  n_u <- vctrs::vec_size(u)
  n_xb <- vctrs::vec_size(xb)
  n_wxg <- vctrs::vec_size(wxg)
  all_n <- c(n_lw, n_u, n_xb, n_wxg)
  check_number_decimal(lambda)
  model <- rlang::arg_match(model)
  if (!max(all_n) == min(all_n)) {
    rlang::abort("`u`, `xb`, `wxg`,  and `listw` must have the same number of features")
  }

  u1 <- switch(
    model,
    "sar" = inverse_prod(listw, u, lambda),
    "ma" = u + lambda * spdep::lag.listw(listw, u)
  )
  xb + wxg + u1
}

#' Simulate Spatial Lag Model (SAR)
#' 
#' Simulate y for a SAR model. 
#' 
#' @inheritParams sim_slx
#' @param rho the spatial autoregressive coefficient for the spatially lagged dependent variable.
#' @references [`spreg.dgp.dgp_lag`](https://pysal.org/spreg/generated/spreg.dgp.dgp_lag.html#spreg.dgp.dgp_lag)
#' @export
#' ncol <- 20
#' n <- ncol^2
#' listw <- sim_grid_listw(ncol, ncol)  # Create spatial weights for a grid
#' u <- make_error(n)  # Simulate random errors
#' x <- make_x(
#'   n,
#'   mu = c(0.25, 5),
#'   var = c(1, 0.75),
#'   method = "normal"
#' )  # Generate x variables
#' 
#' # create xb with intercept = 1, beta1 = 2, beta2 = -3
#' xb <- make_xb(x, c(1, 2, -3))
#' y <- sim_sar(u, xb, listw)
#' 
#' # combine data 
#' df <- cbind(y = y, x)
#' 
#' # fit SAR model
#' # Note lambda, x_1, and x_2 estimates.
#' spatialreg::stsls(y ~ ., df, listw)
sim_sar <- function(u, xb, listw, rho = 0.5) {
  n_lw <- length(listw$neighbours)
  n_u <- vctrs::vec_size(u)
  n_xb <- vctrs::vec_size(xb)
  all_n <- c(n_lw, n_u, n_xb)
  check_number_decimal(rho)
  if (!max(all_n) == min(all_n)) {
    rlang::abort("`u`, `xb`, and `listw` must have the same number of features")
  }
  y1 <- xb + u
  as.numeric(inverse_prod(listw, y1, rho))
}

#' Simulate the Spatial Durbin Model
#' 
#' @references [`spreg.dgp.dgp_spdurbin`](https://pysal.org/spreg/generated/spreg.dgp.dgp_spdurbin.html#spreg.dgp.dgp_spdurbin)
#' @export
#' ncol <- 20
#' n <- ncol^2
#' listw <- sim_grid_listw(ncol, ncol)  # Create spatial weights for a grid
#' u <- make_error(n)  # Simulate random errors
#' x <- make_x(
#'   n,
#'   mu = c(0.25, 5),
#'   var = c(1, 0.75),
#'   method = "normal"
#' )  # Generate x variables
#' 
#' # create xb with intercept = 1, beta1 = 2, beta2 = -3
#' xb <- make_xb(x, c(1, 2, -3))
#' wx <- make_wx(x, listw)
#' wxg <- make_wxg(wx, c(-2, 1.5))
#' y <- sim_durbin(u, xb, wxg, listw, rho = 0.5)
#' 
#' # combine data 
#' df <- cbind(y = y, x)
#' 
#' # fit SDM
#' spatialreg::lagsarlm(y ~ ., df, listw, Durbin = TRUE)
sim_durbin <- function(u, xb, wxg, listw, rho = 0.5) {
  n_lw <- length(listw$neighbours)
  n_u <- vctrs::vec_size(u)
  n_wxg <- vctrs::vec_size(wxg)
  n_xb <- vctrs::vec_size(xb)
  all_n <- c(n_lw, n_u, n_xb, n_wxg)
  check_number_decimal(rho)
  if (!max(all_n) == min(all_n)) {
    rlang::abort("`u`, `xb`, `wxg`, and `listw` must have the same number of features")
  }
  y1 <- xb + wxg + u
  as.numeric(inverse_prod(listw, y1, rho))
}

#' Simulate the Spatial Autoregressive Model with Autoregressive Errors 
#' 
#' Generate `y` values for the "combo" / SARAR / SAC model. 
#' @references [`spreg.dgp.dgp_lagerr`](https://pysal.org/spreg/generated/spreg.dgp.dgp_lagerr.html#spreg.dgp.dgp_lagerr)
#' @export
#' @inheritParams sim_sem
#' @inheritParams sim_sar
sim_sarar <- function(u, xb, listw, rho = 0.5, lambda = 0.2, model = c("sar", "ma")) {
  n_lw <- length(listw$neighbours)
  n_u <- vctrs::vec_size(u)
  n_xb <- vctrs::vec_size(xb)
  all_n <- c(n_lw, n_u, n_xb)
  check_number_decimal(lambda)
  check_number_decimal(rho)

  model <- rlang::arg_match(model)
  
  if (!max(all_n) == min(all_n)) {
    rlang::abort("`u`, `xb`, and `listw` must have the same number of features")
  }

  u1 <- switch(
    model,
    "sar" = inverse_prod(listw, u, lambda),
    "ma" = u + lambda * spdep::lag.listw(listw, u)
  )

  y1 <- xb + u1
  inverse_prod(listw, y1, rho)

}


#' Simulate General Nested Model
#' 
#' @references [`spreg.dgp.dgp_gns`](https://pysal.org/spreg/generated/spreg.dgp.dgp_gns.html#spreg.dgp.dgp_gns)
#' @export
sim_gns <- function(u, xb, wxg, listw, rho = 0.5, lambda = 0.2, model = c("sar", "ma")) {
  n_lw <- length(listw$neighbours)
  n_u <- vctrs::vec_size(u)
  n_xb <- vctrs::vec_size(xb)
  n_wxg <- vctrs::vec_size(wxg)
  all_n <- c(n_lw, n_u, n_xb, n_wxg)
  check_number_decimal(lambda)
  check_number_decimal(rho)

  model <- rlang::arg_match(model)
  
  if (!max(all_n) == min(all_n)) {
    rlang::abort("`u`, `xb`, `wxg`, and `listw` must have the same number of features")
  }

  u1 <- switch(
    model,
    "sar" = inverse_prod(listw, u, lambda),
    "ma" = u + lambda * spdep::lag.listw(listw, u)
  )

  y1 <- xb + wxg + u1 
  inverse_prod(listw, y1, rho)
}

#' Simiulate Matrix Exponential Spatial Lag Model
#' @importFrom spatialreg as_dgRMatrix_listw 
#' @export
#' @references [`dgp_mess`](https://pysal.org/spreg/_modules/spreg/dgp.html#dgp_mess)
sim_mess <- function(u, xb, listw, rho = 0.5) {
  n_lw <- length(listw$neighbours)
  n_u <- vctrs::vec_size(u)
  n_xb <- vctrs::vec_size(xb)
  all_n <- c(n_lw, n_u, n_xb)
  if (!max(all_n) == min(all_n)) {
    rlang::abort("`u`, `xb`, `wxg`, and `listw` must have the same number of features")
  }
  check_number_decimal(rho)
  
  alpha <- log(1 - rho)
  w <- as(listw, "CsparseMatrix")
  aw <- w * alpha * -1
  xbu <- xb + u
  as.numeric(Matrix::expm(aw) %*% xbu)
}
