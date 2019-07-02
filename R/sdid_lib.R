library(osqp)

sc_weight_no_intercept = function(M, target, zeta = 1) {
    if (nrow(M) != length(target)) {
        stop("invalid dimensions")
    }
    
    # solve.QP cannot have 0 penalty for quadratic term
    if (zeta == 0) { zeta = 1e-06 }
    
    # we solve a QP with parameters [weights, imbalance] 
    # where we use an equality constraint to impose that
    # imbalance = M * weights - target
    # our objective is encoded as zeta*||weights||^2 + || imbalance ||^2 / length(target)
    # = [ weights, imbalance]' * [zeta*I, 0; 0, (1/length(target)) I] * [weights, imbalance] 
    # in our call to solve.QP, the parameter Dmat is this block-diagonal matrix, 
    # and we pass dvec=0 because we have no linear term
    #Dmat = diag(c(rep(zeta, ncol(M)), rep(1 / length(target), nrow(M))))
    #dvec = rep(0, ncol(M) + nrow(M))
    # our first nrow(M)+1 constraints are equality constraints
    # the first nrow(M) impose that M*weights - imbalance = target
    # the next imposes that sum(weights)=1
    # and the remaining constraints impose the positivity of our weights
    #AT = rbind(cbind(M, diag(1, nrow(M))),
    #           c(rep(1, ncol(M)), rep(0, nrow(M))),
    #           cbind(diag(1, ncol(M)), matrix(0, ncol(M), nrow(M))))
    #lb = c(target, 1, rep(0, ncol(M)))
    #ub = c(target, 1, rep(Inf, ncol(M)))
    #soln = solve_osqp(Dmat, -dvec, AT, lb, ub, osqpSettings(verbose=FALSE))
    #gamma = soln$x[1:ncol(M)]
    
    Dmat = t(M) %*% M + zeta * diag(ncol(M))
    dvec = t(target) %*% M
    meq = 1
    AT = rbind(rep(1, ncol(M)), 
               diag(1, ncol(M)))
    lb = c(1, rep(0, ncol(M)))
    ub = c(1, rep(Inf, ncol(M)))
    soln = solve_osqp(Dmat, -dvec, AT, lb, ub, osqpSettings(verbose=FALSE))
    gamma = soln$x
    return(gamma)
}

sc_weight_intercept = function(M, target, zeta = 1) {
    if (nrow(M) != length(target)) {
        stop("invalid dimensions")
    }
    
    # solve.QP cannot have 0 penalty for quadratic term
    if (zeta == 0) { zeta = 1e-06 }
    
    # parametrization is [weights, imbalance, c] 
    Dmat = diag(c(rep(zeta, ncol(M)), rep(1 / length(target), nrow(M)), 1e-06))
    dvec = rep(0, ncol(M) + nrow(M) + 1)
    AT = rbind(cbind(M, diag(1, nrow(M)), -1),
               c(rep(1, ncol(M)), rep(0, nrow(M)), 0),
               cbind(diag(1, ncol(M)), matrix(0, ncol(M), nrow(M)), 0))
    lb = c(target, 1, rep(0, ncol(M)))
    ub = c(target, 1, rep(Inf, ncol(M)))
    
    soln = solve_osqp(Dmat, -dvec, AT, lb, ub, osqpSettings(verbose=FALSE))
    gamma = soln$x[1:ncol(M)]
    return(gamma)
}

#' Computes regularized synthetic control weights
#  if intercept = FALSE, these weights gamma solve
#' argmin_{gamma}{||M gamma - target||_2^2 / length(target) + zeta * ||gamma||_2^2 : sum gamma = 1, gamma >= 0}.
#' if intercept = TRUE, these weights gamma solve
#' argmin_{gamma,c}{||M gamma - target - c||_2^2 / length(target) + zeta * ||gamma||_2^2 : sum gamma = 1, gamma >= 0}
#' Allowing an intercept may be useful when there are global time trends; see (2.7) in the paper.
#' @param M
#' @param target
#' @param zeta. Defaults to 1.
#' @param intercept. Defaults to FALSE
#' @return gamma
#' @export sc_weight
sc_weight = function(M, target, zeta=1, intercept=FALSE) {
    if(intercept) {
        sc_weight_intercept(M,target,zeta)
    } else {
        sc_weight_no_intercept(M,target,zeta)
    }
}

#' synthetic diff-in-diff and synthetic control predictions 
#' for a single missing observation, element (n, T) in an n x T matrix.
#' @param Y, the observation matrix
#' @param zeta, the weight on an L2 penalty on the weights. See (6.1) in the paper. Defaults to var(Y).
#' @return a 2-vector of estimates, synthetic diff-in-diff followed by synthetic control
#' @export sdid_predict
sdid_predict = function(Y, zeta = var(as.numeric(Y)), lambda.intercept=FALSE, omega.intercept=FALSE) {
    NN = nrow(Y)
    TT = ncol(Y)
    lambda.weight = sc_weight(Y[-NN, -TT], Y[-NN, TT], zeta = zeta, intercept = lambda.intercept)
    omega.weight = sc_weight(t(Y[-NN, -TT]), Y[NN, -TT], zeta = zeta, intercept = omega.intercept)
    SC.transpose.est = sum(lambda.weight * Y[NN, -TT])
    SC.est = sum(omega.weight * Y[-NN, TT])
    interact.est = omega.weight %*% Y[-NN, -TT] %*% lambda.weight
    sdid.est = SC.est + SC.transpose.est - interact.est
    attr(sdid.est, 'sc.estimate') = SC.est
    attr(sdid.est, 'lambda') = lambda.weight
    attr(sdid.est, 'omega') = omega.weight
    return(sdid.est)
}

#' Computes synthetic diff-in-diff estimate for an average treatment effect with a treated block. 
#' See Section 4.1 of the paper. Also computes a jacknife estimate of its standard error.
#' This uses the method described in the footnote in Section 5 of the paper if fast.var=FALSE.
#' If fast.var=T, we use a variant which is less computationally demanding,
#' in which we do not re-estimate omega in our jackknife replications.
#' @param Y, the observation matrix. 
#' @param n_0, the number of control units. Rows 1-n_0 of Y correspond to the control units.
#' @param T_0, the number of pre-treatment time steps. Columns 1-T_0 of Y correspond to pre-treatment time steps.
#' @param zeta, the weight on an L2 penalty on the weights. See (6.1) in the paper. Defaults to var(Y).
#' @param fast.var. Defaults to TRUE.
#' @return An average treatment effect estimate, with a standard error estimate attached as the attribute 'se'
#' @export sdid_est
sdid_est <- function(Y, n_0, T_0, zeta=var(as.numeric(Y)), fast.var=T, lambda.intercept=FALSE, omega.intercept=FALSE){
    n = nrow(Y)
    T = ncol(Y)
    Y_00 <- Y[1:n_0,1:T_0]
    Y_10 <- Y[(n_0+1):n,1:T_0]
    Y_01 <- Y[1:n_0,(T_0+1):T]
    Y_11 <- Y[(n_0+1):n,(T_0+1):T]
    omega.weight <- sc_weight(t(Y_00), colMeans(Y_10), zeta = zeta, intercept = omega.intercept)
    lambda.weight <- sc_weight(Y_00, rowMeans(Y_01),   zeta = zeta, intercept = lambda.intercept)
        
    est <- sdid_simple(omega.weight, lambda.weight, Y_00, Y_10, Y_01, Y_11)
    
    est_jk <- rep(0,n-n_0)
    for(i in 1:(n-n_0)) {
        if(fast.var) {
            est_jk[i] <- sdid_simple(omega.weight, lambda.weight, Y_00, Y_10[-i,], Y_01, Y_11[-i,])
        } else {
            est_jk[i] <- sdid_est(Y[-(n_0+i),], n_0, T_0, zeta=zeta, fast.var=T, 
                                  lambda.intercept = lambda.intercept, omega.intercept = omega.intercept)
        }
    }
    V.hat <- (n - n_0 - 1) * mean((est_jk - est)^2)
    attr(est, 'se')=sqrt(V.hat)
    attr(est, 'lambda') = lambda.weight
    attr(est, 'omega') = omega.weight
    return(est)
}
sdid_simple <- function(omega.weight, lambda.weight, Y_00, Y_10, Y_01, Y_11){
    interact_est <- omega.weight %*% Y_00 %*% lambda.weight
    time_est <- sum(lambda.weight* colMeans(Y_10))
    unit_est <- sum(omega.weight*rowMeans(Y_01))
    est <- as.vector(mean(Y_11) - unit_est - time_est+ interact_est)
    return(est)
}
        

#' Returns a diff-in-diff estimate for treatment effect with a treated block. 
#' Also returns, as an attribute of the return value, a jacknife estimate of its standard error.
#' @param Y, the observation matrix. 
#' @param n_0, the number of control units. Rows 1-n_0 of Y correspond to the control units.
#' @param T_0, the number of pre-treatment time steps. Columns 1-T_0 of Y correspond to pre-treatment time steps.
#' @return An average treatment effect estimate, with a standard error estimate attached as the attribute 'se'
#' @export did_est
did_est <- function(Y, n_0, T_0){
    n = nrow(Y)
    T = ncol(Y)
    Y_00 <- Y[1:n_0,1:T_0]
    Y_10 <- Y[(n_0+1):n,1:T_0]
    Y_01 <- Y[1:n_0,(T_0+1):T]
    Y_11 <- Y[(n_0+1):n,(T_0+1):T]
        
    est <- did_simple(Y_00, Y_10, Y_01, Y_11)
    est_jk <- rep(0, nrow=n)
    for(i in 1:n) {
        if (i <= n_0){
            Y_00_jk <- Y_00[-i,] 
            Y_10_jk <- Y_10
            Y_01_jk <- Y_01[-i,]
            Y_11_jk <- Y_11
            est_jk[i] <- did_simple(Y_00_jk, Y_10_jk, Y_01_jk, Y_11_jk)
        } else {
            Y_00_jk <- Y_00
            Y_10_jk <- Y_10[-(i - n_0 +1),]
            Y_01_jk <- Y_01
            Y_11_jk <- Y_11[-(i - n_0 +1),]
            est_jk[i] <- did_simple(Y_00_jk, Y_10_jk, Y_01_jk, Y_11_jk)
        }
    }
    V.hat <- (n - 1) * mean((est_jk - est)^2)
    attr(est,'se')=sqrt(V.hat)
    return(est)
}
did_simple <- function(Y_00, Y_10, Y_01, Y_11){
    return(as.vector(mean(Y_11) - mean(Y_01) - mean(Y_10) + mean(Y_00)))
}


#' Imputes control potential outcomes for each element of a treated block, Y[i,j] for i > n_0, j > T_0. 
#' Each such entry is imputed separately using synthetic diff-in-diff based on control and pretreatment observations in Y.
#' Use sdid_est instead of this when estimating average treatment effects with multiple missing units,
#' as it is better to impute the mean of the missing entries than to impute each missing entry individually and average them.
#' @param Y, the observation matrix. 
#' @param n_0, the number of control units. Rows 1-n_0 of Y correspond to the control units.
#' @param T_0, the number of pre-treatment time steps. Columns 1-T_0 of Y correspond to pre-treatment time steps.
#' @return A copy of Y with entries Y[i,j] for i > n_0, j > T_0 replaced with imputed values
#' @export sdid_impute
sdid_impute = function(Y, n_0, T_0, lambda.intercept = FALSE, omega.intercept = FALSE) { 
    n = nrow(Y)
    T = ncol(Y)
    for(ii in (n_0+1):n) {
        for(jj in (T_0+1):T) {
            Y[ii,jj] = sdid_predict(Y[c(1:n_0, ii), c(1:T_0, jj)], lambda.intercept=lambda.intercept, omega.intercept=omega.intercept)[1]
        }
    }
    Y
} 
