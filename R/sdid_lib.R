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
sc_weight = function(M, target, zeta = 1, intercept = FALSE, solver = c("OSQP", "ECOS", "SCS"), verbose = FALSE) {
    solver = match.arg(solver)
    if (nrow(M) != length(target)) {
        stop("invalid dimensions")
    }

    if(solver %in% c('ECOS', 'SCS')) {
        weights = CVXR::Variable(ncol(M))
        if(intercept) { 
            theintercept = CVXR::Variable(1)
            objective = zeta^2 * length(target) * sum(weights^2) + sum((M %*% weights - target - theintercept)^2)
        } else {
            objective = zeta^2 * length(target) * sum(weights^2) + sum((M %*% weights - target)^2)
        }
        constraints = list(sum(weights) == 1, weights >= 0)  
        cvx.problem = CVXR::Problem(CVXR::Minimize(objective), constraints)
        cvx.output = solve(cvx.problem, solver = solver, verbose = verbose)
        as.numeric(cvx.output$getValue(weights))     
    } else { 
        if(intercept) {
            Dmat = t(cbind(M,1)) %*% cbind(M,1) + zeta^2 * diag(rep(c(1,0),c(ncol(M), 1))) * length(target)
            dvec = t(target) %*% cbind(M,1)
            A = cbind(rbind(rep(1, ncol(M)),
                             diag(1, ncol(M))), 0)
            lb = c(1, rep(0, ncol(M)))
            ub = c(1, rep(Inf, ncol(M)))
            soln = osqp::solve_osqp(Dmat, -dvec, A, lb, ub, osqpSettings(verbose=FALSE, adaptive_rho=FALSE))
            soln$x[1:ncol(M)]

        } else {
            Dmat = t(M) %*% M + zeta^2 * diag(ncol(M)) * length(target)
            dvec = t(target) %*% M
            meq = 1
            A = rbind(rep(1, ncol(M)),
                      diag(1, ncol(M)))
            lb = c(1, rep(0, ncol(M)))
            ub = c(1, rep(Inf, ncol(M)))
            soln = osqp::solve_osqp(Dmat, -dvec, A, lb, ub, osqpSettings(verbose=FALSE, adaptive_rho=FALSE))
            soln$x
        }
    } 
}


#' use synthetic diff-in-diff and synthetic control to impute 
#' a single missing observation, element (N, T) in an N x T matrix.
#' @param Y, the observation matrix
#' @param zeta.lambda, the weight on an L2 penalty on lambda. See (6.1) in the paper. Defaults to zero.
#' @param zeta.omega,  analogous for omega. Defaults to var(Y).
#' @return a 2-vector of estimates, synthetic diff-in-diff followed by synthetic control
#' @export synthdid_impute_1
synthdid_impute_1 = function(Y, zeta.lambda=0, zeta.omega=var(as.numeric(Y)), lambda.intercept=FALSE, omega.intercept=FALSE, solver=NULL) {
    N = nrow(Y)
    T = ncol(Y)
    lambda.weight = sc_weight(Y[-N, -T, drop=FALSE], Y[-N, T], zeta = zeta.lambda, intercept = lambda.intercept, solver = solver)
    omega.weight = sc_weight(t(Y[-N, -T, drop=FALSE]), Y[N, -T], zeta = zeta.omega, intercept = omega.intercept)
    SC.transpose.est = sum(lambda.weight * Y[N, -T])
    SC.est = sum(omega.weight * Y[-N, T])
    interact.est = omega.weight %*% Y[-N, -T, drop=FALSE] %*% lambda.weight
    sdid.est = SC.est + SC.transpose.est - interact.est
    attr(sdid.est, 'sc.estimate') = SC.est
    attr(sdid.est, 'lambda') = lambda.weight
    attr(sdid.est, 'omega') = omega.weight
    return(sdid.est)
}

#' Computes synthetic diff-in-diff estimate for an average treatment effect on a treated block. 
#' See Section 4.1 of the paper. Also computes a jacknife estimate of its standard error.
#' If fast.var=TRUE, this uses the method described in Section 5 of the paper, in which lambda and omega are held fixed 
#' over jackknife replications. If fast.var=FALSE, we re-estimate lambda and omega in each jacknife replication.
#' @param Y, the observation matrix. 
#' @param N_0, the number of control units. Rows 1-N_0 of Y correspond to the control units.
#' @param T_0, the number of pre-treatment time steps. Columns 1-T_0 of Y correspond to pre-treatment time steps.
#' @param zeta.lambda, the weight on an L2 penalty on lambda. See (6.1) in the paper. Defaults to zero.
#' @param zeta.omega,  analogous for omega. Defaults to var(Y).
#' @param fast.var. Defaults to TRUE.
#' @return An average treatment effect estimate, with a standard error estimate attached as the attribute 'se'
#' @export synthdid_estimate
synthdid_estimate <- function(Y, N_0, T_0, zeta.lambda=0, zeta.omega=var(as.numeric(Y)), fast.var=TRUE, lambda.intercept=FALSE, omega.intercept=FALSE, solver=NULL){
    N = nrow(Y)
    T = ncol(Y)
    stopifnot(N > N_0, T > T_0)

    Y_00 <- Y[1:N_0,1:T_0, drop=FALSE]
    Y_10 <- Y[(N_0+1):N,1:T_0, drop=FALSE]
    Y_01 <- Y[1:N_0,(T_0+1):T, drop=FALSE]
    Y_11 <- Y[(N_0+1):N,(T_0+1):T, drop=FALSE]
    omega.weight <- sc_weight(t(Y_00), colMeans(Y_10), zeta = zeta.omega,  intercept = omega.intercept, solver = solver)
    lambda.weight <- sc_weight(Y_00, rowMeans(Y_01),   zeta = zeta.lambda, intercept = lambda.intercept, solver = solver)
    est <- synthdid_simple(omega.weight, lambda.weight, Y_00, Y_10, Y_01, Y_11)
        
    if(N == N_0 + 1) { ## if we cannot jackknife rows, return NA variance estimate
       V.hat <- NA
    } else { 
        est_jk = rep(0,N)
        if(fast.var) { 
            for (i in 1:N_0) {
                # drop the jackknifed-out row of omega and renormalize to sum to one
                # this renormalization arises from the weighted least squares interpretation of our estimator
                # it is not asymptotically necessary but helps in practice
                omega.jk  = omega.weight[-i]/sum(omega.weight[-i]) 
                est_jk[i] = synthdid_simple(omega.jk, lambda.weight, Y_00[-i,, drop=FALSE], Y_10, Y_01[-i,, drop=FALSE], Y_11)
            }
            for(i in (N_0+1):N) {
                est_jk[i] <- synthdid_simple(omega.weight, lambda.weight, Y_00, Y_10[-(i - N_0),, drop=FALSE], Y_01, Y_11[-(i - N_0),, drop=FALSE])
            }
        } else {
            for(i in 1:N) {
                est_jk[i] <- synthdid_estimate(Y[-i,, drop=FALSE], ifelse(i <= N_0, N_0 - 1, N_0), T_0, zeta.lambda=zeta.lambda, zeta.omega=zeta.omega, fast.var=T, 
                                      lambda.intercept = lambda.intercept, omega.intercept = omega.intercept)
            }
        }
        V.hat <- (N - 1) * mean((est_jk - est)^2)
    }
    attr(est, 'se')=sqrt(V.hat)
    attr(est, 'lambda') = lambda.weight
    attr(est, 'omega') = omega.weight
    return(est)
}
synthdid_simple <- function(omega.weight, lambda.weight, Y_00, Y_10, Y_01, Y_11){
    interact_est <- omega.weight %*% Y_00 %*% lambda.weight
    time_est <- sum(lambda.weight* colMeans(Y_10))
    unit_est <- sum(omega.weight*rowMeans(Y_01))
    est <- as.vector(mean(Y_11) - unit_est - time_est+ interact_est)
    return(est)
}
        

#' Returns a diff-in-diff estimate for treatment effect with a treated block. 
#' Also returns, as an attribute of the return value, a jacknife estimate of its standard error.
#' @param Y, the observation matrix. 
#' @param N_0, the number of control units. Rows 1-N_0 of Y correspond to the control units.
#' @param T_0, the number of pre-treatment time steps. Columns 1-T_0 of Y correspond to pre-treatment time steps.
#' @return An average treatment effect estimate, with a standard error estimate attached as the attribute 'se'
#' @export did_estimate
did_estimate <- function(Y, N_0, T_0){
    N = nrow(Y)
    T = ncol(Y)
    stopifnot(N > N_0, T > T_0)

    Y_00 <- Y[1:N_0,1:T_0, drop=FALSE]
    Y_10 <- Y[(N_0+1):N,1:T_0, drop=FALSE]
    Y_01 <- Y[1:N_0,(T_0+1):T, drop=FALSE]
    Y_11 <- Y[(N_0+1):N,(T_0+1):T, drop=FALSE]
        
    est <- did_simple(Y_00, Y_10, Y_01, Y_11)
    if(N == N_0 + 1) { ## if we cannot jackknife rows, return NA variance estimate
       V.hat <- NA
    } else { 
        est_jk <- rep(0, nrow=N)
        for(i in 1:N) {
            if (i <= N_0){
                Y_00_jk <- Y_00[-i,, drop=FALSE] 
                Y_10_jk <- Y_10
                Y_01_jk <- Y_01[-i,, drop=FALSE]
                Y_11_jk <- Y_11
                est_jk[i] <- did_simple(Y_00_jk, Y_10_jk, Y_01_jk, Y_11_jk)
            } else {
                Y_00_jk <- Y_00
                Y_10_jk <- Y_10[-(i - N_0),, drop=FALSE]
                Y_01_jk <- Y_01
                Y_11_jk <- Y_11[-(i - N_0),, drop=FALSE]
                est_jk[i] <- did_simple(Y_00_jk, Y_10_jk, Y_01_jk, Y_11_jk)
            }
        }
        V.hat <- (N - 1) * mean((est_jk - est)^2)
    }
    attr(est,'se')=sqrt(V.hat)
    return(est)
}
did_simple <- function(Y_00, Y_10, Y_01, Y_11){
    return(as.vector(mean(Y_11) - mean(Y_01) - mean(Y_10) + mean(Y_00)))
}


#' Imputes control potential outcomes for each element of a treated block, Y[i,j] for i > N_0, j > T_0. 
#' Each such entry is imputed separately using synthetic diff-in-diff based on control and pretreatment observations in Y.
#' Use sdid_est instead of this when estimating average treatment effects with multiple missing units,
#' as it is better to impute the mean of the missing entries than to impute each missing entry individually and average them.
#' @param Y, the observation matrix. 
#' @param N_0, the number of control units. Rows 1-N_0 of Y correspond to the control units.
#' @param T_0, the number of pre-treatment time steps. Columns 1-T_0 of Y correspond to pre-treatment time steps.
#' @param zeta.lambda, the weight on an L2 penalty on lambda. See (6.1) in the paper. Defaults to zero.
#' @param zeta.omega,  analogous for omega. Defaults to var(Y).
#' @return A copy of Y with entries Y[i,j] for i > N_0, j > T_0 replaced with imputed values
#' @export synthdid_impute
synthdid_impute = function(Y, N_0, T_0, zeta.lambda=0, zeta.omega=var(as.numeric(Y)), lambda.intercept = FALSE, omega.intercept = FALSE, solver = NULL) { 
    N = nrow(Y)
    T = ncol(Y)
    for(ii in (N_0+1):N) {
        for(jj in (T_0+1):T) {
            Y[ii,jj] = synthdid_impute_1(Y[c(1:N_0, ii), c(1:T_0, jj)], zeta.lambda=zeta.lambda, zeta.omega=zeta.omega,
                                         lambda.intercept=lambda.intercept, omega.intercept=omega.intercept, solver=solver)[1]
        }
    }
    Y
} 
