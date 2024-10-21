# Helper function to generate data
#' @export
make_grid <- function(ncol, nrow = ncol) {
  sf::st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(ncol, nrow))
}

# grid <- make_grid(5)
# listw <- spdep::nb2listw(spdep::poly2nb(grid))
# n = 25

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
    "lognormal" = rlnorm(n, mu, sqrt(var))
  )
}


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
  MASS::mvrnorm(10, mu, mcov)
}


#' @export
make_x_uniform <- function(n = 5, var = 1) {
  check_number_decimal(var)
  check_number_whole(n)
  runif(n, 0, sqrt(12 * var))
}


#' @export
make_x_normal <- function(n = 5, mu = 0, var = 1) {
  check_number_decimal(var)
  check_number_whole(n)
  check_number_decimal(mu)
  rnorm(n, mu, var)
}


#' @export
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

#' @param order unused. 
#' @export
make_wx <- function(x, listw, order = NULL) {
  lagged <- lapply(x, function(.x) spdep::lag.listw(listw, .x))
  names(lagged) <- paste0(names(lagged), "_lag")
  cbind(x, data.frame(lagged))
}

#' @examples
#' grid <- make_grid(5)
#' listw <- spdep::nb2listw(spdep::poly2nb(grid))
#' n <- 25
#' x <- make_x(25, c(0,1), c(1,4))
#' betas <- runif(3, max = 2)
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


#' @examples
# grid <- make_grid(5)
# listw <- spdep::nb2listw(spdep::poly2nb(grid))
# n <- 25
# x <- make_x(25, c(0,1), c(1,4))
# wx <- make_wx(x, listw)
# gamma <- runif(4, max = 2)
# make_wxg(wx, gamma) 
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



# n <- 100
# grid <- make_grid(n)
# listw <- spdep::nb2listw(spdep::poly2nb(grid))
# u <- make_x(n^2, mu = c(0.5, 0.1), var = rep(1, 2), method = "normal")
# dgp_errorproc(u, listw)

#' Simulate Spatial Error Model
#' 
#' @references See [`spreg.dgp.dgp_errproc`](https://pysal.org/spreg/generated/spreg.dgp.dgp_errproc.html#spreg.dgp.dgp_errproc)
sim_error <- function(u, listw, lambda = 0.5, model = c("sar", "ma")) {
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
#' @examples
#' u = make_x(50, method = "normal")
#' x = make_x(50)
#' xb = make_xb(x, c(1,2))
#' y <- dgp_ols(u, xb)
#' mod <- lm(y ~ x[[1]])
#' summary(mod)
#' @references [`spreg.dgp.dgp_ols`](https://pysal.org/spreg/generated/spreg.dgp.dgp_ols.html#spreg.dgp.dgp_ols)
#' @export
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
#' @examples
#' n <- 100
#' grid <- make_grid(sqrt(n), sqrt(n))
#' listw <- spdep::nb2listw(spdep::poly2nb(grid))
#' u <- make_error(n, method = "normal")
#' x <- make_x(n, method = "uniform")
#' xb <- make_xb(x, c(1, 2))
#' wx <- make_wx(x, listw)
#' wxg <- make_wxg(wx, c(2, 2))
#' y <- dgp_slx(u, xb, wxg)
#' df <- data.frame(y, x)
#' spatialreg::lmSLX(y ~ ., data = df, listw = listw) 
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
#' @references [`spreg.dgp.dgp_sperror`](https://pysal.org/spreg/generated/spreg.dgp.dgp_sperror.html#spreg.dgp.dgp_sperror)
#' @export
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
#' @references [`spreg.dgp.dgp_lag`](https://pysal.org/spreg/generated/spreg.dgp.dgp_lag.html#spreg.dgp.dgp_lag)
#' @export
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
  inverse_prod(listw, y1, rho)
}

#' Simulate the Spatial Durbin Model
#' 
#' @references [`spreg.dgp.dgp_spdurbin`](https://pysal.org/spreg/generated/spreg.dgp.dgp_spdurbin.html#spreg.dgp.dgp_spdurbin)
#' @export
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
  inverse_prod(listw, y1, rho)
}

#' Simulate the Spatial Autoregressive Model with Autoregressive Errors 
#' 
#' Generate `y` values for the "combo" / SARAR / SAC model. 
#' @references [`spreg.dgp.dgp_lagerr`](https://pysal.org/spreg/generated/spreg.dgp.dgp_lagerr.html#spreg.dgp.dgp_lagerr)
#' @export
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
#' @export
#' @referece [`dgp_mess`](https://pysal.org/spreg/_modules/spreg/dgp.html#dgp_mess)
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
