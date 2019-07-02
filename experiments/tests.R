library(mvtnorm)
library(synthdid)

N_0 <- 100
N_1 <- 20 
T_0 <- 120
T_1 <- 5
N <- N_0 + N_1
T <- T_0 + T_1
tau <- 1
sigma <- 2
rank <- 2
rho <- 0.7
var <- outer(1:T, 1:T, FUN=function(x, y) rho^(abs(x-y)))
W <- (1:N > N_0) %*% t(1:T > T_0)

U <- matrix(rpois(rank * N, sqrt(1:N) / sqrt(N)), N, rank)
V <- matrix(rpois(rank * T, sqrt(1:T) / sqrt(T)), T, rank)
mu <- U %*% t(V)
error <- rmvnorm(N, sigma = var, method = "chol")
Y <- mu + tau * W  + sigma * error 

## These should be small if we've correctly formulated our mathematical programs in sc_weight
max(abs(attr(synthdid_estimate(Y,N_0,T_0, solver='OSQP'), 'lambda') - 
        attr(synthdid_estimate(Y,N_0,T_0, solver='ECOS'), 'lambda')))
max(abs(attr(synthdid_estimate(Y,N_0,T_0, lambda.intercept=TRUE, solver='OSQP'), 'lambda') - 
        attr(synthdid_estimate(Y,N_0,T_0, lambda.intercept=TRUE, solver='ECOS'), 'lambda')))


