# return x minimizing ||Ax - b||^2 + zeta^2 n || x ||^2                   if intercept=FALSE
#      | x minimizing min_x0 || + ||Ax + x0 - b||^2 + zeta^2 n || x ||^2  if intercept=TRUE
# here n = length(b)
simplex.least.squares =  function(A, b, zeta = 0, intercept = FALSE) {
    x = CVXR::Variable(ncol(A))
    constraints = list(sum(x) == 1, x >= 0)
    if(intercept) {
	x0 = CVXR::Variable(1)
	objective = sum((A %*% x + x0 - b)^2) + zeta^2 * length(b) * sum(x^2)
    } else {
	objective = sum((A %*% x - b)^2) + zeta^2 * length(b) * sum(x^2)
    }
    cvx.problem = CVXR::Problem(CVXR::Minimize(objective), constraints)
    cvx.output = CVXR::solve(cvx.problem, solver = 'ECOS')
    as.numeric(cvx.output$getValue(x))
}


sigma.default = function(Y, N0, T0) { sd(apply(Y[1:N0,1:T0], 1, diff)) }
epsilon = 1e-6
synthdid.reference = function(Y, N0, T0, zeta.omega=((nrow(Y)-N0)*(ncol(Y)-T0))^(1/4) * sigma.default(Y, N0, T0)) {
    N = nrow(Y); T=ncol(Y); N1 = N-N0; T1=T-T0;
    lambda = simplex.least.squares(Y[1:N0, 1:T0],    rowMeans(Y[1:N0, (T0+1):T, drop=FALSE]), zeta=epsilon*sigma.default(Y,N0,T0), intercept=TRUE)
    omega =  simplex.least.squares(t(Y[1:N0, 1:T0]), colMeans(Y[(N0+1):N, 1:T0, drop=FALSE]), zeta=zeta.omega, intercept=TRUE)
    estimate = t(c(-omega, rep(1/N1, N1))) %*% Y  %*% c(-lambda, rep(1/T1, T1))
}
sc.reference = function(Y, N0, T0, zeta.omega=1e-6 * sigma.default(Y, N0, T0)) {
    N = nrow(Y); T=ncol(Y); N1 = N-N0; T1=T-T0;
    omega =  simplex.least.squares(t(Y[1:N0, 1:T0]), colMeans(Y[(N0+1):N, 1:T0, drop=FALSE]), zeta=zeta.omega, intercept=FALSE)
    estimate = t(c(-omega, rep(1 / N1, N1))) %*% Y  %*% c(-rep(0, T0), rep(1/T1, T1))
}
did.reference = function(Y, N0, T0) {
    N = nrow(Y); T=ncol(Y); N1 = N-N0; T1=T-T0;
    estimate = t(c(-rep(1/N0, N0), rep(1 / N1, N1))) %*% Y  %*% c(-rep(1/T0, T0), rep(1 / T1, T1))
}
