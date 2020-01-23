library(mvtnorm)
library(osqp)
source('../R/sdid_lib.R')

### Outdated: Check for equivalence between solutions using ECOS and OSQP Solvers. OSQP removed because it too often failed.
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

## and these should be near one
cosine = function(x,y){ sum(x*y)/sqrt(sum(x^2)*sum(y^2)) }
cosine(attr(synthdid_estimate(Y,N_0,T_0, solver='OSQP'), 'lambda'),
        attr(synthdid_estimate(Y,N_0,T_0, solver='ECOS'), 'lambda'))
cosine(attr(synthdid_estimate(Y,N_0,T_0, lambda.intercept=TRUE, solver='OSQP'), 'lambda'),
        attr(synthdid_estimate(Y,N_0,T_0, lambda.intercept=TRUE, solver='ECOS'), 'lambda'))

### Check for equivalence of ECOS and FW (Frank-Wolfe) solvers. Seems to work.

source('../R/sdid_lib.R')

M = rnorm(25*50)
dim(M) = c(25,50)
target = rnorm(25)

max.iter=500
weights = sc_weight.fw(M, target, zeta=0, max.iter=max.iter)
weights.ls = sc_weight.fw(M, target, zeta=0, max.iter=max.iter, step='linesearch')

## correctness
for(zeta in c(0,1)) {
 for(intercept in c(FALSE, TRUE)) { 
    weights.ref = sc_weight(M,target,    intercept=intercept, zeta=zeta)
    weights     = sc_weight.fw(M,target, intercept=intercept, zeta=zeta)
    weights.ls  = sc_weight.fw(M,target, intercept=intercept, zeta=zeta, step='linesearch')
    print(c(max(abs(weights - weights.ref)), max(abs(weights - weights.ref))))
    }
}

## convergence rate: line search not very useful
plot.data = data.frame(iter=c(1:length(attr(weights,'objective')), 1:length(attr(weights.ls, 'objective'))),
                       type=rep(c('decreasing', 'ls'), c(length(attr(weights,'objective')), length(attr(weights.ls, 'objective')))),
                       val=c(attr(weights,'objective'), attr(weights.ls, 'objective')))
ggplot(plot.data[plot.data$iter > 50,]) + geom_line(aes(x=iter, y=val, color=type))


## larger example with signal
# line search converges much faster
M = rnorm(2500*5000)
dim(M) = c(2500,5000)
x = 1/(1:ncol(M))^2
x = x / sum(x)
target = M %*% x
max.iter=500
weights = sc_weight.fw(M,target, intercept=TRUE, zeta=0, max.iter=max.iter, step='decreasing')
weights.ls = sc_weight.fw(M,target, intercept=TRUE, zeta=0, max.iter=max.iter, step='linesearch')
sum((M%*% weights.ls - target)^2) - sum((M%*% weights - target)^2) 
sum((M%*% weights.ls - target)^2) - sum((M%*% x - target)^2) 
c(max(abs(weights.ls - x)), max(abs(weights - x)))

plot.data = data.frame(iter=c(1:length(attr(weights,'objective')), 1:length(attr(weights.ls, 'objective'))),
                       type=rep(c('decreasing', 'ls'), c(length(attr(weights,'objective')), length(attr(weights.ls, 'objective')))),
                       val=c(attr(weights,'objective'), attr(weights.ls, 'objective')))
ggplot(plot.data) + geom_line(aes(x=iter, y=val, color=type)) + ylim(0,.01)




