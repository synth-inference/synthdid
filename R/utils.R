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

# pass column names / indices. For now, assume indices. Defaults correspond to our built-in datasets
make.panel = function(panel, unit=2, time=3, outcome=4, treatment=5) {
    panel = panel[order(panel[, treatment], panel[, unit], panel[, time]), ]
    T = length(unique(panel[, time]))
    N = length(unique(panel[, unit]))
    Y = matrix(panel[,outcome], N, T, byrow = TRUE,
               dimnames = list(unique(panel[,unit]), unique(panel[,time])))
    W = matrix(panel[,treatment], N, T, byrow = TRUE,
               dimnames = list(unique(panel[,unit]), unique(panel[,time])))
    N0 = N-sum(W[,T])
    T0 = T-sum(W[N,])
    if(! (all(W[1:N0,] == 0) && all(W[,1:T0] == 0) && all(W[N0+1,T0+1]==1))) { 
	error('The package cannot use this data. Treatment adoption is not simultaneous')
    }
    list(Y=Y, N0=N0, T0=T0)
}

# A convenience function for generating data for unit tests.
random.low.rank = function() {
  n_0 <- 100
  n_1 <- 10
  T_0 <- 120
  T_1 <- 20
  n <- n_0 + n_1
  T <- T_0 + T_1
  tau <- 1
  sigma <- .5
  rank <- 2
  rho <- 0.7
  var <- outer(1:T, 1:T, FUN=function(x, y) rho^(abs(x-y)))

  W <- (1:n > n_0) %*% t(1:T > T_0)
  U <- matrix(rpois(rank * n, sqrt(sample(1:n)) / sqrt(n)), n, rank)
  V <- matrix(rpois(rank * T, sqrt(1:T) / sqrt(T)), T, rank)
  alpha <- outer(10*sample(1:n)/n, rep(1,T))
  beta <-  outer(rep(1,n), 10*(1:T)/T)
  mu <- U %*% t(V) + alpha + beta
  error <- mvtnorm::rmvnorm(n, sigma = var, method = "chol")
  Y <- mu + tau * W  + sigma * error
  rownames(Y) = 1:n
  colnames(Y) = 1:T
  list(Y=Y, L=mu, N0=n_0, T0=T_0)
}
