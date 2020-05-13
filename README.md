[![Build Status](https://travis-ci.com/synth-inference/synthdid.svg?branch=master)](https://travis-ci.com/synth-inference/synthdid)

# synthdid: Synthetic Difference in Differences Estimation

This package implements the synthetic difference in difference estimator (SDID) for the
average treatment effect in panel data, as proposed in Arkhangelsky et al (2019).
We consider a setting in which we observe a matrix `Y = L + tau W + noise` where W
is a matrix of indicators for treatment.  All treated units must begin treatment simultaneously,
so W indicates a treated block, i.e. `W[i,j] = 1` for `i > N_0, j > T_0` and is zero otherwise.
This applies, in particular, to the case of a single treated unit.

This package is currently in beta and the functionality and interface is subject to change.

To install this package in R, run the following commands:
```R
library(devtools)
install_github("synth-inference/synthdid")
```
Example usage:

```R
library(synthdid)
library(mvtnorm)

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
U <- matrix(rpois(rank * n, sqrt(1:n) / sqrt(n)), n, rank)
V <- matrix(rpois(rank * T, sqrt(1:T) / sqrt(T)), T, rank)
alpha <- outer(10*(1:n)/n, rep(1,T))
beta <-  outer(rep(1,n), 10*(1:T)/T)
mu <- U %*% t(V) + alpha + beta
error <- rmvnorm(n, sigma = var, method = "chol")
Y <- mu + tau * W  + sigma * error
rownames(Y) = 1:n
colnames(Y) = 1:T

tau.hat = synthdid_estimate(Y,n_0,T_0)
se = synthdid_se(tau.hat)

print(paste("true tau:", tau))
print(paste0("point estimate: ", round(tau.hat, 2)))
print(paste0("95% CI for tau: (", round(tau.hat - 1.96 * se, 2), ", ", round(tau.hat + 1.96 * se, 2), ")"))
plot(tau.hat)
```

#### References
Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
<b>Synthetic Difference in Differences</b>
2019.
[<a href="https://arxiv.org/abs/1812.09970">arxiv</a>]
