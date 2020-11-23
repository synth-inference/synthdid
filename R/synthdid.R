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
  jack.se = sqrt(((n - 1)/n) * (n - 1) * var(u))

  jack.se
}

contract3 = function(X,v) {
    stopifnot(length(dim(X)) == 3, dim(X)[3] == length(v))
    out = array(0, dim=dim(X)[1:2])
    if(length(v) == 0) { return(out) }
    for(ii in 1:length(v)) {
        out = out + v[ii] * X[,,ii]
    }
    return(out)
}

# a Frank-Wolfe step for \\Ax - b||^2 + eta * ||x||^2 with x in unit simplex.
fw.step = function(A, x, b, eta, alpha=NULL) {
    Ax = A %*% x
    half.grad = t(Ax - b) %*% A + eta *  x
    i = which.min(half.grad)
    if(!is.null(alpha)) {
        x=x*(1-alpha)
        x[i] = x[i] + alpha
        return( x )
    } else {
        d.x = -x; d.x[i] = 1-x[i]
        if(all(d.x == 0)) { return(x) }
        d.err = A[,i] - Ax
        step = -t(c(half.grad)) %*% d.x / (sum(d.err^2) + eta * sum(d.x^2))
        constrained.step = min(1, max(0, step))
        return( x + constrained.step*d.x )
    }
}

# a Frank-Wolfe solver for synthetic control weights using exact line search
sc.weight.fw = function(Y, zeta, intercept=TRUE, lambda=NULL, min.decrease=1e-3, max.iter=1000) {
    T0 = ncol(Y)-1
    N0 = nrow(Y)-1
    if(is.null(lambda)) { lambda=rep(1/T0,T0) }
    if(intercept) {
        Y = apply(Y, 2, function(col) { col - mean(col) })
    }

    t=0
    vals = rep(NA, max.iter)
    A = Y[1:N0, 1:T0]
    b = Y[1:N0,T0+1]
    eta = N0*Re(zeta^2)
    while(t < max.iter && (t < 2 || vals[t-1] - vals[t] > min.decrease^2)) {
        t=t+1
        lambda.p = fw.step(A, lambda, b, eta)
        lambda = lambda.p
        err = Y[1:N0,] %*% c(lambda, -1)
        vals[t] = Re(zeta^2)*sum(lambda^2) + sum(err^2)/N0
    }
    list(lambda=lambda, vals=vals)
}

# A Frank-Wolfe + Gradient solver for lambda, omega, and beta when there are covariates
# Uses the exact line search Frank-Wolfe steps for lambda, omega and (1/t)*gradient steps for beta
# pass update.lambda=FALSE/update.omega=FALSE to fix those weights at initial values, defaulting to uniform 1/T0 and 1/N0
sc.weight.fw.covariates = function(Y, X=array(0,dim=c(dim(Y),0)), zeta.lambda = 0, zeta.omega=0,
                                   lambda.intercept=TRUE, omega.intercept=TRUE,
                                   min.decrease=1e-3, max.iter=1000,
                                   lambda=NULL, omega=NULL, beta = NULL, update.lambda=TRUE, update.omega=TRUE) {
    stopifnot(length(dim(Y))==2, length(dim(X)) == 3, all(dim(Y)==dim(X)[1:2]), all(is.finite(Y)), all(is.finite(X)))
    T0 = ncol(Y)-1
    N0 = nrow(Y)-1
    if(length(dim(X)) == 2) { dim(X) = c(dim(X),1) }
    if(is.null(lambda)) {  lambda=rep(1/T0,T0)   }
    if(is.null(omega))  {  omega=rep(1/N0,N0)    }
    if(is.null(beta))   {  beta=rep(0,dim(X)[3]) }

    update.weights = function(Y, lambda, omega) {
        Y.lambda = if(lambda.intercept) { apply(Y[1:N0,], 2, function(row) { row - mean(row) }) } else { Y[1:N0,] }
        if(update.lambda) { lambda = fw.step(Y.lambda[,1:T0], lambda, Y.lambda[,T0+1], N0*Re(zeta.lambda^2)) }
        err.lambda = Y.lambda %*% c(lambda, -1)

        Y.omega = if(omega.intercept)  { apply(t(Y[,1:T0]), 2, function(row) { row - mean(row) }) } else { t(Y[,1:T0]) }
        if(update.omega) { omega  = fw.step(Y.omega[,1:N0], omega, Y.omega[,N0+1], T0*Re(zeta.omega^2)) }
        err.omega  = Y.omega %*% c(omega,  -1)

        val = Re(zeta.omega^2) * sum(omega^2) + Re(zeta.lambda^2) * sum(lambda^2) + sum(err.omega^2) / T0 + sum(err.lambda^2) / N0
        list(val=val, lambda=lambda, omega = omega, err.lambda=err.lambda, err.omega=err.omega)
    }

    vals = rep(NA, max.iter)
    t=0
    Y.beta = Y - contract3(X,beta)
    weights = update.weights(Y.beta, lambda, omega)
    # state is kept in weights$lambda, weights$omega, beta
    while(t < max.iter && (t < 2 || vals[t-1] - vals[t] > min.decrease^2)) {
        t=t+1
        grad.beta = -if(dim(X)[3]==0) { c() } else {
            apply(X, 3, function(Xi) {
                t(weights$err.lambda) %*% Xi[1:N0,]    %*% c(weights$lambda,-1) / N0  +
                t(weights$err.omega)  %*% t(Xi[,1:T0]) %*% c(weights$omega, -1) / T0
            })
        }

        alpha = 1/t
        beta = beta - alpha * grad.beta
        Y.beta = Y - contract3(X,beta)
        weights = update.weights(Y.beta, weights$lambda, weights$omega)
        vals[t] = weights$val
    }
    list(lambda=weights$lambda, omega=weights$omega, beta=beta, vals=vals)
}

# collapse Y to an N0+1 x T0+1 vector by averaging the last N1=nrow(Y)-N0 rows and T1=ncol(Y)-T0 columns
collapsed.form = function(Y, N0, T0) {
    N = nrow(Y); T=ncol(Y)
    rbind(cbind(  Y[1:N0,1:T0, drop=FALSE],                    rowMeans(Y[1:N0,(T0+1):T, drop=FALSE])),
          cbind(  t(colMeans(Y[(N0+1):N,1:T0, drop=FALSE])),   mean(Y[(N0+1):N,(T0+1):T, drop=FALSE])))
}

# return the component-wise sum of decreasing vectors in which NA is taken to mean that the vector has stopped decreasing
# and we can use the last non-na element. Where both are NA, leave as NA.
pairwise.sum.decreasing = function(x,y) {
    na.x = is.na(x)
    na.y = is.na(y)
    x[is.na(x)] = min(x[!na.x])
    y[is.na(y)] = min(y[!na.y])
    pairwise.sum = x+y
    pairwise.sum[na.x & na.y] = NA
    pairwise.sum
}

#' Computes synthetic diff-in-diff estimate for an average treatment effect on a treated block.
#' See Section 4.1 of the paper.
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param X an optional 3-D array of time-varying covariates. Shape should be N X T X C for C covariates.
#' @param zeta.lambda Its square is weight of the ridge penalty relative to MSE. Defaults to 0.
#' @param zeta.omega Analogous for omega. Defaults to the standard deviation of first differences of Y.
#' @param lambda.intercept Binary. Use an intercept when estimating lambda.
#' @param omega.intercept Binary. Use an intercept when estimating omega.
#' @param weights a list with fields lambda and omega. If non-null weights$lambda is passed,
#'        we use them instead of estimating lambda weights. Same for weights$omega.
#' @param update.lambda If true, solve for lambda using the passed value of weights$lambda only as an initialization.
#'        If false, use it exactly as passed. Defaults to false if a non-null value of weights$lambda is passed.
#' @param update.omega  Analogous.
#' @param min.decrease Tunes a stopping criterion for our weight estimator. Stop after an iteration results in a decrease
#'		        in penalized MSE smaller than min.decrease^2.
#' @param max.iter A fallback stopping criterion for our weight estimator. Stop after this number of iterations.
#' @return An average treatment effect estimate, 'weights' and 'setup' attached as attributes.
#'         Weights contains the estimated weights lambda and omega and corresponding intercepts.
#'         If covariates X are passedas well as regression coefficients beta if X is passed
#'         Setup is a list describing the problem passed in: Y, N0, T0, X.
#' @export synthdid_estimate
#' @importFrom stats sd
synthdid_estimate <- function(Y, N0, T0, X=array(dim=c(dim(Y),0)),
                              zeta.lambda=0, zeta.omega=sd(apply(Y,1,diff)),
                              lambda.intercept=TRUE, omega.intercept=TRUE,
                              weights = list(lambda=NULL, omega=NULL, vals=NULL),
			      update.lambda=is.null(weights$lambda), update.omega = is.null(weights$omega),
			      min.decrease=1e-3*sd(apply(Y,1,diff)), max.iter=1e4) {
    stopifnot(nrow(Y) > N0, ncol(Y) > T0, length(dim(X)) %in% c(2,3), dim(X)[1:2] == dim(Y), is.list(weights),
              is.null(weights$lambda) || length(weights$lambda) == T0, is.null(weights$omega) || length(weights$omega) == N0,
	      !is.null(weights$lambda) || update.lambda, !is.null(weights$omega) || update.omega)
    if(length(dim(X)) == 2) { dim(X) = c(dim(X),1) }
    N1 = nrow(Y)-N0
    T1 = ncol(Y)-T0

    if(dim(X)[3] == 0) {
        weights$vals = NULL
        if(update.lambda) {
            Yc=collapsed.form(Y, N0, T0)
            lambda.opt = sc.weight.fw(Yc, zeta = zeta.lambda, intercept=lambda.intercept, min.decrease=min.decrease, max.iter=max.iter)
            weights$lambda = lambda.opt$lambda
            weights$vals =   lambda.opt$vals
        }
        if(update.omega) {
            Yc=collapsed.form(Y, N0, T0)
            omega.opt  = sc.weight.fw(t(Yc), zeta = zeta.omega,  intercept=omega.intercept, min.decrease=min.decrease, max.iter=max.iter)
            weights$omega = omega.opt$lambda
            if(is.null(weights$vals)) { weights$vals = omega.opt$vals }
            else { weights$vals = pairwise.sum.decreasing(weights$vals, omega.opt$vals) }
        }
    } else {
        Yc=collapsed.form(Y, N0, T0)
        Xc=apply(X, 3, function(Xi) { collapsed.form(Xi, N0, T0) })
        dim(Xc) = c(dim(Yc), dim(X)[3])
        weights = sc.weight.fw.covariates(Yc, Xc, zeta.lambda=zeta.lambda, zeta.omega=zeta.omega,
                                          lambda.intercept=lambda.intercept, omega.intercept=omega.intercept,
					  min.decrease=min.decrease, max.iter=max.iter,
                                          lambda=weights$lambda, omega=weights$omega, update.lambda=update.lambda, update.omega=update.omega)
    }

    X.beta = contract3(X, weights$beta)
    estimate = t(c(-weights$omega, rep(1/N1,N1))) %*% (Y - X.beta) %*% c(-weights$lambda, rep(1/T1, T1))

    class(estimate) = 'synthdid_estimate'
    attr(estimate, 'weights') = weights
    attr(estimate, 'setup') = list(Y=Y, X=X, N0=N0, T0=T0)
    attr(estimate, 'opts') =  list(zeta.lambda=zeta.lambda, zeta.omega=zeta.omega,
                              lambda.intercept=lambda.intercept, omega.intercept=omega.intercept)
    return(estimate)
}

#' synthdid_estimate for synthetic control estimates.
#' Takes all the same parameters, but default, passes options for synthetic control
#' with no intercept and a penalty term that defaults to the standard deviation of first differences of Y.
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param X an optional 3-D array of time-varying covariates. Shape should be N X T X C for C covariates.
#' @param zeta.lambda Its square is weight of the ridge penalty relative to MSE. Defaults to 0.
#' @param zeta.omega Analogous for omega. Defaults to the standard deviation of first differences of Y.
#' @param lambda.intercept Binary. Use an intercept when estimating lambda.
#' @param omega.intercept Binary. Use an intercept when estimating omega.
#' @param weights a list with fields lambda and omega. If non-null weights$lambda is passed,
#'        we use them instead of estimating lambda weights. Same for weights$omega.
#' @param min.decrease Tunes a stopping criterion for our weight estimator. Stop after an iteration results in a decrease
#'		        in penalized MSE smaller than min.decrease^2.
#' @param max.iter A fallback stopping criterion for our weight estimator. Stop after this number of iterations.
#' @return An average treatment effect estimate, 'weights' and 'setup' attached as attributes.
#'         Weights contains the estimated weights lambda and omega and corresponding intercepts.
#'         If covariates X are passedas well as regression coefficients beta if X is passed
#'         Setup is a list describing the problem passed in: Y, N0, T0, X.
#' @export sc_estimate
sc_estimate = function(Y, N0, T0, X=array(dim=c(dim(Y),0)),
                       zeta.lambda=0, zeta.omega=sd(apply(Y,1,diff)),
                       lambda.intercept=FALSE, omega.intercept=FALSE,
                       weights = list(lambda=rep(0,T0), omega=NULL, vals=NULL),
		       min.decrease=1e-3, max.iter=1e4) {
    synthdid_estimate(Y, N0, T0, X=X,
                       zeta.lambda=zeta.lambda, zeta.omega=zeta.omega,
                       lambda.intercept=lambda.intercept, omega.intercept=omega.intercept,
                       weights = weights, min.decrease=1e-3, max.iter=1e4)
}

#' synthdid_estimate for diff-in-diff estimates.
#' Takes all the same parameters, but default, uses constant weights lambda and omega
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param X an optional 3-D array of time-varying covariates. Shape should be N X T X C for C covariates.
#' @param zeta.lambda Its square is weight of the ridge penalty relative to MSE. Defaults to 0.
#' @param zeta.omega Analogous for omega. Defaults to the standard deviation of first differences of Y.
#' @param lambda.intercept Binary. Use an intercept when estimating lambda.
#' @param omega.intercept Binary. Use an intercept when estimating omega.
#' @param weights a list with fields lambda and omega. If non-null weights$lambda is passed,
#'        we use them instead of estimating lambda weights. Same for weights$omega.
#' @param min.decrease Tunes a stopping criterion for our weight estimator. Stop after an iteration results in a decrease
#'		        in penalized MSE smaller than min.decrease^2.
#' @param max.iter A fallback stopping criterion for our weight estimator. Stop after this number of iterations.
#' @return An average treatment effect estimate, 'weights' and 'setup' attached as attributes.
#'         Weights contains the estimated weights lambda and omega and corresponding intercepts.
#'         If covariates X are passedas well as regression coefficients beta if X is passed
#'         Setup is a list describing the problem passed in: Y, N0, T0, X.
#' @export did_estimate
did_estimate = function(Y, N0, T0, X=array(dim=c(dim(Y),0)),
                       zeta.lambda=0, zeta.omega=sd(apply(Y,1,diff)),
                       lambda.intercept=FALSE, omega.intercept=FALSE,
                       weights = list(lambda=rep(1/T0,T0), omega=rep(1/N0,N0), vals=NULL),
		       min.decrease=1e-3, max.iter=1e4) {
    synthdid_estimate(Y, N0, T0, X=X,
                       zeta.lambda=zeta.lambda, zeta.omega=zeta.omega,
                       lambda.intercept=lambda.intercept, omega.intercept=omega.intercept,
                       weights = weights, min.decrease=1e-3, max.iter=1e4)
}

#' Computes a placebo variant of our estimator using pre-treatment data only
#' @param estimate, as output by synthdid_estimate
#' @param treated.fraction, the fraction of pre-treatment data to use as a placebo treatment period
#'        Defaults to NULL, which indicates that it should be the fraction of post-treatment to pre-treatment data
#' @export synthdid_placebo
synthdid_placebo = function(estimate, treated.fraction=NULL) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    X.beta = contract3(setup$X, weights$beta)

    if(is.null(treated.fraction)) { treated.fraction = 1 - setup$T0/ncol(setup$Y) }
    placebo.T0 = floor(setup$T0*(1-treated.fraction))

    synthdid_estimate(setup$Y[,1:setup$T0], setup$N0, placebo.T0, setup$X[,1:setup$T0,],
                      zeta.lambda=opts$zeta.lambda, zeta.omega=opts$zeta.omega,
                      lambda.intercept=opts$lambda.intercept, omega.intercept=opts$omega.intercept)
}

#' Outputs the effect curve that was averaged to produce our estimate
#' @param estimate, as output by synthdid_estimate
#' @export synthdid_effect_curve
synthdid_effect_curve = function(estimate) {
    setup = attr(estimate, 'setup')
    weights = attr(estimate, 'weights')
    X.beta = contract3(setup$X, weights$beta)
    N1 = nrow(setup$Y)-setup$N0
    T1 = ncol(setup$Y)-setup$T0

    tau.sc = t(c(-weights$omega, rep(1/N1,N1))) %*% (setup$Y - X.beta)
    tau.curve = tau.sc[setup$T0+(1:T1)] - c(tau.sc[1:setup$T0] %*% weights$lambda)
    tau.curve
}

#' An estimate of the standard error of our estimator via unit-wise jackknife.
#' Should not be trusted when the number of treated units is small, and returns NA if there is only one treated unit.
#' By default, does not recompute the weights omega and lambda. Instead, weights$lambda and weights$beta will be used unchanged and
#' weights$omega will have jackknifed-out units removed then be renormalized to sum to one.
#' Pass weights=NULL to recompute weights for each jackknife replication.
#' Requires bootstrap
#' @param estimate, as output by synthdid_estimate
#' @param weights, like attr(estimate, 'weights')
#' @export synthdid_se
synthdid_se = function(estimate, weights = attr(estimate, 'weights')) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    if(setup$N0 == nrow(setup$Y)-1) { return(NA) }

    sum_normalize = function(x){ x / sum(x) }
    theta = function(ind) {
        weights.jk = weights
        if(!is.null(weights)) { weights.jk$omega = sum_normalize(weights$omega[ind[ind <= setup$N0]]) }
        estimate.jk = synthdid_estimate(setup$Y[ind,], sum(ind <= setup$N0), setup$T0, X=setup$X[ind,,],
                                        zeta.lambda=opts$zeta.lambda, zeta.omega=opts$zeta.omega,
                                        lambda.intercept=opts$lambda.intercept, omega.intercept=opts$omega.intercept,
                                        weights = weights.jk)
    }
    jackknife(1:nrow(setup$Y), theta)
}
