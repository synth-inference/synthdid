library(quadprog)

# Computes regularized synthetic control weights:
# gamma = argmin_g{||M g - target||_2^2 / length(target) +
#                       zeta * ||g||_2^2 : sum g = 1, g >= 0}
sc_weight = function(M, target, zeta = 1) {
    if (nrow(M) != length(target)) {
        stop("invalid dimensions")
    }
    
    # solve.QP cannot have 0 penalty for quadratic term
    if (zeta == 0) { zeta = 1e-06 }
    
    # parametrization is [weights, imbalance] 
    Dmat = diag(c(rep(zeta, ncol(M)), rep(1 / length(target), nrow(M))))
    dvec = rep(0, ncol(M) + nrow(M))
    AT = rbind(cbind(M, diag(1, nrow(M))),
               c(rep(1, ncol(M)), rep(0, nrow(M))),
               cbind(diag(1, ncol(M)), matrix(0, ncol(M), nrow(M))))
    bvec = c(target, 1, rep(0, ncol(M)))
    meq = nrow(M) + 1
    
    soln = solve.QP(Dmat, dvec, t(AT), bvec, meq = meq)
    gamma = soln$solution[1:ncol(M)]
    return(gamma)
}

# Computes regularized synthetic control weights with an intercept.
# gamma = argmin_g{||M g - target - c||_2^2 / length(target) +
#                       zeta * ||g||_2^2 : sum g = 1, g >= 0}
# This may be useful when there are global time trends; see (2.7) in the paper.
sc_weight_FE = function(M, target, zeta = 1) {
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
    bvec = c(target, 1, rep(0, ncol(M)))
    meq = nrow(M) + 1
    
    soln = solve.QP(Dmat, dvec, t(AT), bvec, meq = meq)
    gamma = soln$solution[1:ncol(M)]
    return(gamma)
}

# Makes SDID prediction for missing control outcome in the simple case where only a
# single observation (n, T) is treated. Also returns the synthetic control estimate.
sdid_predict = function(Y, zeta = 1) {
    NN = nrow(Y)
    TT = ncol(Y)
    lambda.weight = sc_weight(Y[-NN, -TT], Y[-NN, TT], zeta = zeta)
    omega.weight = sc_weight(t(Y[-NN, -T]), Y[NN, -TT], zeta = zeta)
    SC.transpose.est = sum(lambda.weight * Y[NN, -TT])
    SC.est = sum(omega.weight * Y[-NN, TT])
    interact.est = omega.weight %*% Y[-NN, -TT] %*% lambda.weight
    sdid.est = SC.est + SC.transpose.est - interact.est
    c(SDID=sdid.est, SC=SC.est)
}

# Makes a SDID estimate for treatment effect with a treated block
# standard error is estimated using the jackknife as described 
# in the footnote in Section 5 of the paper.
# If fast.var=T, we use a variant which is less computationally demanding,
# in which we do not re-estimate omega in our jackknife replications.
sdid_est <- function(Y, n_0, T_0, n, T, fast.var=T){
    	
    zeta.simple <- var(as.numeric(Y))
	Y_00 <- Y[1:n_0,1:T_0]
    Y_10 <- Y[(n_0+1):n,1:T_0]
    Y_01 <- Y[1:n_0,(T_0+1):T]
    Y_11 <- Y[(n_0+1):n,(T_0+1):T]
	omega.weight <- sc_weight(t(Y_00), colMeans(Y_10), zeta = zeta.simple)
	lambda.weight <- sc_weight(Y_00, rowMeans(Y_01),   zeta = zeta.simple)
		
	est <- sdid_simple(omega.weight, lambda.weight, Y_00, Y_10, Y_01, Y_11)
	
	est_jk <- rep(0,n-n_0)
	for(i in 1:(n-n_0)) {
		if(fast.var) {
			est_jk[i] <- sdid_simple(omega.weight, lambda.weight, Y_00, Y_10[-i,], Y_01, Y_11[-i,])
		} else {
			est_jk[i] <- sdid_est(Y[-(n_0+i),], n_0, T_0, n-1, T, fast.var=T)
		}
	}
	V.hat <- (n - n_0 - 1) * mean((est_jk - est)^2)
	attr(est, 'se')=sqrt(V.hat)
	return(est)
}
sdid_simple <- function(omega.weight, lambda.weight, Y_00, Y_10, Y_01, Y_11){
	interact_est <- omega.weight %*% Y_00 %*% lambda.weight
	time_est <- sum(lambda.weight* colMeans(Y_10))
	unit_est <- sum(omega.weight*rowMeans(Y_01))
	est <- as.vector(mean(Y_11) - unit_est - time_est+ interact_est)
	return(est)
}
		

# Makes a DID estimate for treatment effect with a treated block
# standard error is estimated using the jackknife
did_est <- function(Y, n_0, T_0, n, T){

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


