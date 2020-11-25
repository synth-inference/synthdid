#' Jackknife standard error of function `theta` at samples `x`.
#' @param x vector of samples
#' @param theta a function which returns a scalar estimate
#' @importFrom stats var
#' @keywords internal
jackknife = function(x, theta) {
  n = length(x)
  u = rep(0, n)
  for (i in 1:n) {
    u[i] = theta(x[-i])
  }
  jack.se = sqrt(((n - 1) / n) * (n - 1) * var(u))

  jack.se
}

# collapse Y to an N0+1 x T0+1 vector by averaging the last N1=nrow(Y)-N0 rows and T1=ncol(Y)-T0 columns
collapsed.form = function(Y, N0, T0) {
  N = nrow(Y); T = ncol(Y)
  rbind(cbind(Y[1:N0, 1:T0, drop = FALSE], rowMeans(Y[1:N0, (T0 + 1):T, drop = FALSE])),
    cbind(t(colMeans(Y[(N0 + 1):N, 1:T0, drop = FALSE])), mean(Y[(N0 + 1):N, (T0 + 1):T, drop = FALSE])))
}

# return the component-wise sum of decreasing vectors in which NA is taken to mean that the vector has stopped decreasing
# and we can use the last non-na element. Where both are NA, leave as NA.
pairwise.sum.decreasing = function(x, y) {
  na.x = is.na(x)
  na.y = is.na(y)
  x[is.na(x)] = min(x[!na.x])
  y[is.na(y)] = min(y[!na.y])
  pairwise.sum = x + y
  pairwise.sum[na.x & na.y] = NA
  pairwise.sum
}

#' Simulate correlated mean zero normal variables
#'
#' Simulate correlated mean zero normal variables with covariance matrix
#' `sigma` using the Cholesky decomposition. A simple convenience
#' function to simulate synthetic data with.
#'
#' @param n The number of samples.
#' @param sigma The covariance matrix.
#' @return A n X ncol(sigma) matrix of draws.
#' @keywords internal
sdid_rmvnorm = function(n, sigma) {
  K = ncol(sigma)
  A = chol(sigma)
  z = matrix(rnorm(n * K), n, K)
  z %*% A
}
