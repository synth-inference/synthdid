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
    squared.decrease = Inf
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
        grad.beta = if(dim(X)[3]==0) { c() } else {
            apply(X, 3, function(Xi) { 
                t(weights$err.lambda) %*% Xi[1:N0,]    %*% c(weights$lambda,-1) / N0  + 
                t(weights$err.omega)  %*% t(Xi[,1:T0]) %*% c(weights$omega, -1) / T0 
            })
        }
        
        alpha = 1/t
        beta = beta - alpha * grad.beta
        Y.beta = Y + contract3(X,beta)
        old.weights = weights
        weights = update.weights(Y.beta, weights$lambda, weights$omega)

        squared.step.length = alpha^2*sum(grad.beta^2) + sum((old.weights$lambda-weights$lambda)^2)/N0 + sum((old.weights$omega - weights$omega)^2)/T0
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
			      min.decrease=1e-3, max.iter=1e4) {
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


#' Plots treated and synthetic control trajectories and overlays a 2x2 diff-in-diff diagram of our estimator.
#' In this overlay, the treatment effect is indicated by an arrow.
#' The weights lambda defining our synthetic pre-treatment time period are plotted below.
#' If a list of estimates is passed, plots all of them. By default, does this in different facets.
#' To overlay estimates in the same facet, indicate a facet for each estimator in the argument `facet'.
#' 
#' For SC estimates (lambda=[0,0,...]), plots the trajectories and SC estimate of the effect, but no diagram.
#'
#' Requires ggplot2
#' Due to differences between ggplot and ggplotly, this will warn about an unknown aesthetic frame.
#'
#' @param estimates, a list of estimates output by synthdid_estimate. Or a single estimate.
#'          If estimates have attribute 'intercept' set (scalar in [0,1]), then plot after subtracting
#'          that fraction of the SDID adjustment for the difference between pre-treatment treated and sc curves.
#'          With intercept of almost one, this makes it easier to assess parallel-ness by making trajectories closer
#'          With intercept of one, this overlays curves, and plotting a diagram is suppressed as in the case of a SC estimate.
#' @param treated.name, the name of the treated curve that appears in the legend. Defaults to 'treated'
#' @param control.name, the name of the control curve that appears in the legend. Defaults to 'synthetic control'
#' @param force.sc, TOADD
#' @param facet, a list of the same length as estimates indicating the facet in which to plot each estimate.
#'        The values of the elements of the list are used to label the facets.
#'        If NULL, plot each estimate in a different facet. Defaults to NULL.
#' @param facet.vertical, TOADD
#' @param lambda.comparable, TRUE if the weights lambda should be plotted in such a way that the ribbons
#'        have the same mass from plot to plot, assuming the treated curve is the same. Useful for side-by-side or overlaid plots. 
#'	  Defaults to FALSE if facet is not passed, TRUE if passed.
#' @param overlay specifies a value of 'intercept' for all SDID estimates. Defaults to 0.
#'        If a vector is passed, plots at different intercept levels indicated by the 'frame' aesthetic. ggplotly will interpret this as an animation.
#' @param lambda.plot.scale determines the scale of the plot of the weights lambda. 
#' @param trajectory.linetype, the linetype of the treated and synthetic control trajectories
#' @param effect.curvature, the curvature of the arrows indicating the treatment effect. Defaults to zero. 
#'        Nonzero values help avoid overplotting when plotting multiple estimates in one facet.
#' @param trajectory.alpha determines transparency of trajectories
#' @param diagram.alpha determines transparency of diff-in-diff diagram
#' @param effect.alpha determines transparency of effect arrows
#' @param onset.alpha determines transparency of vertical lines indicating onset of treatment
#' @param alpha.multiplier, a vector of the same length as estimates, is useful for comparing multiple estimates in
#'        one facet but highlighting one or several. All plot elements associated with the estimate are displayed
#'        with alpha multiplied by the corresponding element of alpha.multiplier. Defaults to a vector of ones.
#' @export synthdid_plot
#' @import ggplot2
synthdid_plot = function(estimates, treated.name='treated', control.name='synthetic control', force.sc=FALSE, 
			 facet=NULL, facet.vertical=TRUE, lambda.comparable = !is.null(facet), overlay=0, 
			 lambda.plot.scale=3, trajectory.linetype=1, effect.curvature = 0,
			 trajectory.alpha=.4, diagram.alpha = .95, effect.alpha=.95, onset.alpha = .3, alpha.multiplier = NULL) {
    if(class(estimates) == 'synthdid_estimate') { estimates = list(estimates) } 
    if(is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }
    if(is.null(alpha.multiplier)) { alpha.multiplier = rep(1, length(estimates)) }
    treated = 1
    control = 2    
    groups = factor(c(control,treated), labels=c(control.name, treated.name))
    estimate.factors = factor(1:(length(estimates)+1), labels=c(treated.name, names(estimates)))
    facet_factors = if(is.null(facet)) {
        factor(1:length(estimates), labels=names(estimates))
    } else {
        factor(facet, levels=1:length(unique(facet)), labels=unique(facet))
    }
    grid=expand.grid(estimate=1:length(estimates), overlay=1:length(overlay))
    plot.descriptions = lapply(1:nrow(grid), function(row) {
	est = estimates[[grid$estimate[row]]]
	over = overlay[grid$overlay[row]]	
 
        setup = attr(est, 'setup')
        weights = attr(est, 'weights')
        Y = setup$Y - contract3(setup$X, weights$beta)
        N0 = setup$N0; N1 = nrow(Y)-N0
        T0 = setup$T0; T1 = ncol(Y)-T0
        
        lambda.synth = c(weights$lambda, rep(0, T1)) 
        lambda.target = c(rep(0,T0), rep(1/T1, T1))
        omega.synth  = c(weights$omega,  rep(0, N1))
        omega.target = c(rep(0,N0),  rep(1/N1, N1))

	# if we're given a synthetic control estimate, take note: we'll plot it differently
	# and if we're passed a SDID estimate with 'intercept' attribute = 1, we'll subtract the SDID intercept offset
	# and plot it like a synthetic control estimate
	is.sc = all(weights$lambda==0) || (!is.null(attr(est,'intercept')) && attr(est, 'intercept') == 1)
	if(!is.null(attr(est,'intercept'))) { over = attr(est,'intercept') }  # force overlaying curves if intercept attribute = 1
  
	intercept.offset = over * c((omega.target - omega.synth) %*% Y %*% lambda.synth)
	obs.trajectory = as.numeric(omega.target %*% Y)
        syn.trajectory = as.numeric(omega.synth %*% Y) + intercept.offset

        treated.post   = omega.target %*% Y %*% lambda.target
        treated.pre    = omega.target %*% Y %*% lambda.synth
        control.post = omega.synth %*% Y %*% lambda.target + intercept.offset
        control.pre  = omega.synth %*% Y %*% lambda.synth  + intercept.offset
        sdid.post = as.numeric(control.post + treated.pre - control.pre)

        time = as.numeric(colnames(Y))
        if(length(time) == 0) { time = 1:(T0+T1) }
        pre.time =  lambda.synth  %*% time
        post.time = lambda.target %*% time

        # construct objects on graph
        lines  = data.frame(x = rep(time,2), 
                            y = c(obs.trajectory, syn.trajectory),
                            color=rep(groups[c(treated, control)], each=length(time)))
        points   = data.frame(x = c(post.time, post.time), y=c(treated.post, sdid.post), color=groups[c(treated, control)])
        did.points = data.frame(x     =        c(pre.time,    pre.time,     post.time,     post.time),  
                                y     =        c(treated.pre, control.pre,  control.post,  treated.post),    
                                color = groups[c(treated,     control,      control,       treated)])   
        did.segments = data.frame(x    =        c(pre.time,     pre.time),      
                                  xend =        c(post.time,    post.time),      
                                  y    =        c(control.pre,  treated.pre),  
                                  yend =        c(control.post, treated.post),   
                                  color= groups[c(control, treated)])
        hallucinated.segments = data.frame(x = pre.time, xend = post.time, y = treated.pre, yend = sdid.post)
        guide.segments = data.frame(x    = c(pre.time,    post.time),    
                                    xend = c(pre.time,    post.time),    
                                    y    = c(control.pre, control.post), 
                                    yend = c(treated.pre, sdid.post))   
        arrows = data.frame(x=post.time, xend=post.time, y=sdid.post, yend=treated.post, xscale=max(time)-post.time, color=groups[control])
        
        T0s = attr(est, 'T0s')
        if(!is.null(T0s)) {
            vlines = data.frame(xintercept=time[T0s])
        } else {
            vlines = data.frame(xintercept=time[T0])
        }

        if(lambda.comparable) { 
            height = (max(c(obs.trajectory))-min(c(obs.trajectory)))/lambda.plot.scale
            bottom = min(c(obs.trajectory)) - height
            ribbons = data.frame(x=time[1:T0], ymin = rep(bottom,T0), ymax= bottom + height*lambda.synth[1:T0], color=groups[control])
        } else { 
            height = (max(c(obs.trajectory,syn.trajectory))-min(c(obs.trajectory, syn.trajectory)))/lambda.plot.scale
            bottom = min(c(obs.trajectory, syn.trajectory)) - height
            ribbons = data.frame(x=time[1:T0], ymin = rep(bottom,T0), ymax= bottom + height*lambda.synth[1:T0]/max(lambda.synth), color=groups[control])
        }
        elements = list(lines=lines, points=points, did.segments=did.segments, did.points=did.points, 
	     hallucinated.segments=hallucinated.segments, guide.segments=guide.segments,
             arrows=arrows, vlines=vlines, ribbons=ribbons)
	lapply(elements, function(x) { 
	    x$frame = over
	    x$is.sc = is.sc
	    x$estimate = estimate.factors[grid$estimate[row]+1] # offset because the treated pseudo-estimate factor is first
	    x
	})
    })
    

    one.per.facet = length(unique(facet_factors)) == length(facet_factors) 
    concatenate.field = function(field) {
        do.call(rbind, lapply(plot.descriptions, function(desc) {
                element = desc[[field]]
		estimate.factor = element$estimate[1]
                element$facet = facet_factors[as.integer(estimate.factor)-1]    # offset because the treated pseudo-estimate factor is first
		element$show = alpha.multiplier[as.integer(element$estimate)-1] # "
		element$show[element$color == groups[treated]] = 1              # show treated observations
		# if there are multiple plots per facet, color by estimator rather than by treatment/control
		# make all treated observations the same color, assuming that we're using only one treated observation per facet
		if(!one.per.facet && 'color' %in% colnames(element)) {
		    color = element$estimate
		    color[element$color == groups[treated]] = estimate.factors[1] # treated `estimate factor'
		    element$color = color
		}
		# if there are multiple plots per facet, curve treatment effect arrows so they don't lie on top of one another
		element
            }))
    }
    conc = lapply(names(plot.descriptions[[1]]), concatenate.field)
    names(conc) = names(plot.descriptions[[1]])
    no.sc = function(x) { x[!x$is.sc,] }

    p=ggplot() +
        geom_line(aes(x=x,y=y,color=color,frame=frame, alpha=trajectory.alpha*show),  data=conc$lines, linetype=trajectory.linetype) + 
        geom_point(aes(x=x,y=y,color=color,frame=frame, alpha=diagram.alpha*show), data=conc$points, shape=21) +
        geom_point(aes(x=x,y=y,color=color,frame=frame, alpha=diagram.alpha*show), data=no.sc(conc$did.points)) +
        geom_segment(aes(x=x,xend=xend,y=y,yend=yend,color=color, frame=frame, alpha=diagram.alpha*show), data=no.sc(conc$did.segments)) +
        geom_segment(aes(x=x,xend=xend,y=y,yend=yend,frame=frame, group=estimate, alpha=.5*diagram.alpha*show), data=no.sc(conc$hallucinated.segments), linetype=2,  color='black') +
        geom_segment(aes(x=x,xend=xend,y=y,yend=yend,frame=frame, group=estimate, alpha=.4*diagram.alpha*show), data=no.sc(conc$guide.segments), linetype=2, color='black') +
        geom_vline(aes(xintercept=xintercept, alpha=onset.alpha*show), data=conc$vlines, color='black') + 
        geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax, group=color, fill=color, alpha=.5*diagram.alpha*show), color='black', data = conc$ribbons, show.legend=FALSE) +
      	geom_curve(aes(x=x,xend=xend,y=y,yend=yend, alpha=effect.alpha*show),  data=conc$arrows, curvature=effect.curvature, color='black', arrow=arrow(length=unit(.2, 'cm')))
    # facet if we want multiple facets
    if(!all(conc$lines$facet == conc$lines$facet[1])) { 
        if(facet.vertical) { p = p + facet_grid(facet ~ ., scales='free_y') }
        else { p = p + facet_grid(. ~ facet) }
    }
    # if only one estimate per facet, exclude estimate-denoting linetype from legend
    if(is.null(facet)) { p = p + guides(linetype = FALSE) } 
    # use dates on x axis if provided in colnames(Y)
    p = tryCatch({
        as.Date(colnames(Y))
        p + scale_x_continuous(labels=function(time) { as.Date(time, origin='1970-01-01') })
    }, error = function(e) { p })

    p + xlab('') + ylab('') + labs(color='',fill='') + scale_alpha(guide='none') +
     theme_light() + theme(legend.direction = "horizontal", legend.position = "top")
}



#' For our estimator and a placebo, plots treated and synthetic control trajectories and overlays a 2x2 diff-in-diff diagram.
#' Requires ggplot2
#' @param estimate, as output by synthdid_estimate. 
#' @param overlay, binary, indicates whether plots should be overlaid or shown in different facets. Defaults to FALSE.
#' @param treated.fraction as in synthdid_placebo 
#' @export synthdid_placebo_plot
synthdid_placebo_plot = function(estimate, overlay=FALSE, treated.fraction=NULL) {
   estimates = list(estimate=estimate, placebo=synthdid_placebo(estimate, treated.fraction=treated.fraction))
   synthdid_plot(estimates, facet=if(overlay) { c(1,1) } else { NULL })
} 


synthdid_time_plot = function(estimate) { 
    stopifnot(class(estimate) == 'synthdid_estimate') 
    
    setup = attr(estimate, 'setup')
    weights = attr(estimate, 'weights')
    Y = setup$Y - contract3(setup$X, weights$beta)
    N0 = setup$N0; N1 = nrow(Y)-N0
    T0 = setup$T0; T1 = ncol(Y)-T0
        
    lambda.did  = c(rep(1/T0, T0),  rep(0,T1))
    lambda.synth = c(weights$lambda, rep(0, T1)) 
    lambda.target = c(rep(0,T0), rep(1/T1, T1))
    omega.synth  = c(weights$omega,  rep(0, N1))
    omega.target = c(rep(0,N0),  rep(1/N1, N1))

    points = data.frame(y=c(Y %*% lambda.target, Y %*% lambda.synth, Y %*% lambda.did),
	                x=rep(Y %*% lambda.target, 3),
	                color=rep(factor(c('target', 'synth', 'did')), each=N0+N1))
    ggplot(points) + geom_point(aes(x=x,y=y,color=color)) 
}


#' A diagnostic plot for sc.weight.fw.covariates. Plots the objective function, regularized RMSE,
#' as a function of the number of Frank-Wolfe / Gradient steps taken.
#' Requires ggplot2
#' @param estimates, a list of estimates output by synthdid_estimate. Or a single estimate.
#' @export synthdid_rmse_plot
synthdid_rmse_plot = function(estimates) { # pass an estimate or list of estimates
    if(class(estimates) == 'synthdid_estimate') { estimates = list(estimates) } 
    if(is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }
    rmse = lapply(estimates, function(est) { sqrt(attr(est, 'weights')$vals) })
    plot.data = data.frame(rmse = unlist(rmse),         
                           iteration=unlist(lapply(rmse, function(vals) { 1:length(vals) })),
                           method = unlist(mapply(function(vals, name) { rep(factor(name), length(vals)) }, rmse, names(estimates), SIMPLIFY=FALSE)))
    ggplot(plot.data[!is.na(plot.data$rmse), ]) + geom_line(aes(x=iteration, y=rmse, color=method)) + scale_y_log10() + 
        theme_light() + theme(legend.direction = "horizontal", legend.position = "top") + labs(color='')
}

#' Outputs a table of important synthetic controls and their corresponding weights omega, sorted by weight.
#' The table is truncated to exclude synthetic controls that do not matter for any estimate --- 
#' for each estimate, the truncated controls may have total weight no larger that 1-mass.
#' @param estimates, a list of estimates output by synthdid_estimate. Or a single estimate.
#' @param sort.by, the index of the estimate to sort by. Defaults to 1.
#' @param digits,  the number of digits of weight to display. Defaults to 3.
#' @param mass, which controls the length of the table. Defaults to 0.9.
#' @export synthdid_controls
synthdid_controls = function(estimates, sort.by=1, digits=3, mass=.9) { 
    if(class(estimates) == 'synthdid_estimate') { estimates = list(estimates) } 
    if(is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }
    
    omegas = do.call(cbind, lapply(estimates, function(est) { attr(est, 'weights')$omega }))
    if(is.null(dim(omegas))) { dim(omegas) = c(length(omegas), 1) }

    Y = attr(estimates[[1]], 'setup')$Y
    o = rev(order(omegas[,sort.by]))
    tab = round(omegas[o,,drop=FALSE], digits=digits)
    rownames(tab) = rownames(Y)[o]
    colnames(tab) = names(estimates)
    
    # truncate table to retain a weight sum of at least mass for each unit
    tab.len = max(apply(tab, 2, function(col) { Position(function(x){ x >= mass }, cumsum(col)) }))
    tab[1:tab.len, ]
}

## Export Methods

#' Plot a synthdid object
#' @param x The object to plot
#' @param ... Additional arguments (currently ignored).
#' @method plot synthdid_estimate
#' @export
plot.synthdid_estimate = function(x, ...) {
 synthdid_plot(x, ...)
}

#' Summarize a synthdid object
#' @param object The object to summarize
#' @param ... Additional arguments (currently ignored).
#' @method summary synthdid_estimate
#' @export
summary.synthdid_estimate = function(object, ...) {
    list(estimate = c(object), 
         se = synthdid_se(object),
         controls = synthdid_controls(object))
}

#' Print a synthdid object
#' @param x The object to print
#' @param ... Additional arguments (currently ignored).
#' @method print synthdid_estimate
#' @export
print.synthdid_estimate = function(x, ...) { cat(format(x, ...), "\n") }

#' Format a synthdid object
#' @param x The object to format
#' @param ... Additional arguments (currently ignored).
#' @method format synthdid_estimate
#' @export
format.synthdid_estimate = function(x, ...) {
    setup = attr(x, 'setup')
    weights = attr(x, 'weights')
    sprintf('synthdid x: %1.3f. Effective N0/N0 = %1.1f/%d. Effective T0/T0 = %1.1f/%d. N1,T1 = %d,%d.', 
	c(x), 1/sum(weights$omega^2), setup$N0, 1/sum(weights$lambda^2), setup$T0,
	nrow(setup$Y)-setup$N0, ncol(setup$Y) - setup$T0)
}

