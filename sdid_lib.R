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

# Computes SDID in the simple case where only a single observation
# (n, T) is treated. Also returns the synthetic control estimate.
synth_did = function(Y, zeta = 1) {
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