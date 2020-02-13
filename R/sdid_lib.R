#' Computes regularized synthetic control weights
#  if intercept = FALSE, these weights gamma solve
#' argmin_{gamma}{||M gamma - target||_2^2 / length(target) + zeta^2 ||gamma||_2^2
#' if intercept = TRUE, these weights gamma solve
#' argmin_{gamma,c}{||M gamma + intercept - target||_2^2 / length(target) + zeta^2 ||gamma||_2^2 
#' where gamma is constrained either to the simplex or the l1 ball of radius radius
#' Allowing an intercept may be useful when there are global time trends; see (2.7) in the paper.
#' @param M
#' @param target
#' @param zeta. Defaults to 1.
#' @param intercept.  Defaults to FALSE
#' @param constraint. Defaults to simplex
#' @param radius.    Defaults to 2, used only when constraint=l1
#' @return gamma, with intercept in attribute
#' @export sc_weight
sc_weight = function(M, target, zeta = 1, intercept=FALSE, constraint = c('simplex', 'l1'), radius=2, solver = c(CVXR::installed_solvers(), 'fw'), solver.options=list()) {
    solver = match.arg(solver)
    constraint = match.arg(constraint)
    if (nrow(M) != length(target)) {
        stop("invalid dimensions")
    }

    # return constant weights for zeta=infinity, taking the limit of the solution as zeta->infinity for M finite
    if(zeta == Inf) {
        if(constraint == 'simplex') {
            return(rep(1/ncol(M), ncol(M)))
        } else {
            return(rep(0, ncol(M))); 
        }
    }

    if(solver == 'fw') {
        stopifnot(constraint == 'simplex')
        return(do.call(sc_weight.fw, c(list(M, target, zeta=zeta, intercept=intercept), solver.options)));
    }
        
    weights = CVXR::Variable(ncol(M))
    theintercept = CVXR::Variable(1)
    objective = if(!intercept) { 
        zeta^2 * length(target) * sum(weights^2) + sum((M %*% weights - target)^2)
    } else {
        zeta^2 * length(target) * sum(weights^2) + sum((M %*% weights + theintercept - target)^2)
    }
    constraints = if(constraint == 'simplex') {
        list(sum(weights) == 1, weights >= 0)
    } else {
        list(sum(abs(weights)) <= radius)
    }  
    cvx.problem = CVXR::Problem(CVXR::Minimize(objective), constraints)
    cvx.output = do.call(solve, c(list(cvx.problem, solver = solver), solver.options))
    v=as.numeric(cvx.output$getValue(weights))     
    if(intercept) { attr(v, 'intercept') = cvx.output$getValue(theintercept) }
    v
}

# frank-wolfe solver for ridge-penalized least squares over the simplex.
# essentially as described in https://arxiv.org/pdf/1911.04415.pdf
# exhibits slow convergence near optimum. Error in x relative to an exact solution tends to be about 1/max.iter
sc_weight.fw = function(M, target, zeta=1, intercept=FALSE, min.step.size=1e-3, max.iter=1000, x.init=NULL, step='linesearch') {
    lambda.reg = zeta^2*length(target)
    if(intercept) {
        target = target - mean(target)
        M = apply(M, 2, function(row) { row - mean(row) })
    } 
    
    vals = rep(NA, max.iter)
    t=0
    if(is.null(x.init)) { 
        x=rep(0,ncol(M))
        x[1]=1
    } else {
        x.init
    }

    step.size = Inf
    while(t < max.iter && step.size > min.step.size) {
        t=t+1
        Mx = M %*% x
        dif = Mx - target
        grad = lambda.reg * x + t(M) %*% dif
        i = which.min(grad)
        
        vals[t] = zeta^2*mean(x^2) + mean(dif^2) 
        if(identical(step,'linesearch')) { 
            delta = -x; delta[i] = delta[i]+1
            Mdelta = M[,i] - Mx 
            step.size = min(1, c( -t(grad) %*% delta / (lambda.reg * sum(delta^2) + sum(Mdelta^2))))
            # Check that this analytic line search is correct by doing an iterative line search
            #line.obj = function(step.size) {        
            #    x = x*(1-step.size)
            #    x[i] = x[i] + step.size
            #    lambda.reg * sum(x^2) + sum((M %*% x - target)^2)
            #}
            #approx.step.size = optimize(line.obj, c(0, 1))$minimum
            #print(c(step.size, approx.step.size))
        } else {
            step.size = 2/(t+2)
        }

        x = x*(1-step.size)
        x[i] = x[i] + step.size
    }
    attr(x,'objective') = vals
    x
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
synthdid_estimate <- function(Y, N_0, T_0, 
                              zeta.lambda=0, zeta.omega=sd(as.numeric(Y)), fast.var=TRUE, 
                              lambda.intercept=FALSE, omega.intercept=FALSE, omega.constraint = NULL, lambda.constraint = NULL, 
                              solver=NULL, solver.options=list(), standardize=FALSE,
                              lambda.solver.options=solver.options, omega.solver.options = solver.options) {
    N = nrow(Y)
    T = ncol(Y)
    stopifnot(N > N_0, T > T_0)

    Y_00 <- Y[1:N_0,1:T_0, drop=FALSE]
    Y_10 <- Y[(N_0+1):N,1:T_0, drop=FALSE]
    Y_01 <- Y[1:N_0,(T_0+1):T, drop=FALSE]
    Y_11 <- Y[(N_0+1):N,(T_0+1):T, drop=FALSE]
    
    # if standardize=TRUE, rescale when estimating omega as in Hainmueller's Synth package. 
    # Not recommended. This makes more sense in the more general framework of Hainmueller, where what we balance aren't necessarily the outcomes.
    scale = if(!standardize) { diag(T_0) } else { diag(1 / apply(rbind(Y_00,Y_10), 2, sd)) } 
        
    omega.weight <- sc_weight(t(Y_00 %*% scale), colMeans(Y_10 %*% scale), zeta = zeta.omega,  
                              intercept = omega.intercept,  constraint = omega.constraint, 
                              solver = solver, solver.options=omega.solver.options)
    lambda.weight <- sc_weight(Y_00, rowMeans(Y_01),   zeta = zeta.lambda, 
                              intercept = lambda.intercept, constraint = lambda.constraint, 
                              solver = solver, solver.options=lambda.solver.options)
    tau.curve <- synthdid_simple(omega.weight, lambda.weight, Y_00, Y_10, Y_01, Y_11)
    est <- mean(tau.curve)

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
                est_jk[i] = mean(synthdid_simple(omega.jk, lambda.weight, Y_00[-i,, drop=FALSE], Y_10, Y_01[-i,, drop=FALSE], Y_11))
            }
            for(i in (N_0+1):N) {
                est_jk[i] <- mean(synthdid_simple(omega.weight, lambda.weight, Y_00, Y_10[-(i - N_0),, drop=FALSE], Y_01, Y_11[-(i - N_0),, drop=FALSE]))
            }
        } else {
            for(i in 1:N) {
                est_jk[i] <- synthdid_estimate(Y[-i,, drop=FALSE], ifelse(i <= N_0, N_0 - 1, N_0), T_0, 
                                      zeta.lambda=zeta.lambda, zeta.omega=zeta.omega, fast.var=T, 
                                      lambda.intercept = lambda.intercept, omega.intercept = omega.intercept, 
                                      omega.constraint = omega.constraint, lambda.constraint = lambda.constraint, solver = solver, solver.options = solver.options)
            }
        }
        V.hat <- (N - 1) * mean((est_jk - est)^2)
    }
    
    class(est) = 'synthdid'
    attr(est, 'se')=sqrt(V.hat)
    attr(est, 'tau.curve') = tau.curve
    attr(est, 'lambda') = lambda.weight
    attr(est, 'omega') = omega.weight
    attr(est, 'Y') = Y
    attr(est, 'N0') = N_0
    attr(est, 'T0') = T_0
    return(est)
}

#' returns treatment effect curve tau_t for the post-treatment period
synthdid_simple <- function(omega.weight, lambda.weight, Y_00, Y_10, Y_01, Y_11){
    interact_est <- omega.weight %*% Y_00 %*% lambda.weight
    time_est <- mean(Y_10 %*% lambda.weight)
    unit_est <- omega.weight %*% Y_01
    est <- as.vector(colMeans(Y_11) - c(unit_est) - time_est + c(interact_est))
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
    
    class(est) = 'did'
    attr(est, 'se')=sqrt(V.hat)
    attr(est, 'Y') = Y
    attr(est, 'N0') = N_0
    attr(est, 'T0') = T_0

    return(est)
}
did_simple <- function(Y_00, Y_10, Y_01, Y_11){
    return(as.vector(mean(Y_11) - mean(Y_01) - mean(Y_10) + mean(Y_00)))
}

#' use synthetic diff-in-diff and synthetic control to impute 
#' a single missing observation, element (N, T) in an N x T matrix.
#' @param Y, the observation matrix
#' @param zeta.lambda, the weight on an L2 penalty on lambda. See (6.1) in the paper. Defaults to zero.
#' @param zeta.omega,  analogous for omega. Defaults to var(Y).
#' @return a 2-vector of estimates, synthetic diff-in-diff followed by synthetic control
#' @export synthdid_impute_1
synthdid_impute_1 = function(Y, zeta.lambda=0, zeta.omega=sd(as.numeric(Y)), lambda.intercept=FALSE, omega.intercept=FALSE, solver=NULL) {
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


#' Plots treated and synthetic control trajectories and overlays a 2x2 diff-in-diff diagram of our estimator.
#' In this overlay, the treatment effect is the length of the vertical dotted line. 
#' The weights lambda defining our synthetic pre-treatment time period are plotted below.
#' Requires ggplot2
#' @param est, output of synthdid_estimate
#' @param treated.name, the name of the treated curve that appears in the legend. Defaults to 'treated'
#' @param treated.name, the name of the control curve that appears in the legend. Defaults to 'synthetic control'
#' @export plot.synthdid
plot.synthdid = function(est, treated.name='treated', control.name='synthetic control', use.intercept=FALSE) { 
    library(ggplot2)
    Y = attr(est, 'Y')
    N0 = attr(est, 'N0'); N1 = nrow(Y)-N0
    T0 = attr(est, 'T0'); T1 = ncol(Y)-T0
    omega.synth  = c(attr(est, 'omega'),  rep(0, N1))
    lambda.synth = c(attr(est, 'lambda'), rep(0, T1))
    omega.target = c(rep(0,N0), rep(1/N1, N1))
    lambda.target = c(rep(0,T0), rep(1/T1, T1))
    lambda.intercept = attr(attr(est, 'lambda'), 'intercept')
    if(!use.intercept || is.null(lambda.intercept)) { lambda.intercept = 0 }

    units = factor(c(treated.name, control.name))
    treated = units[1]
    control = units[2]
    unit = rep(units, each=T0+T1)
    time = as.numeric(colnames(Y))
    linetype = 1 

    col.SC = omega.target %*% (Y %*% lambda.synth + lambda.intercept)
    row.SC = omega.synth %*% Y %*% lambda.target
    cross.term = omega.synth %*% (Y %*% lambda.synth + lambda.intercept)
    sdid = as.numeric(row.SC + col.SC - cross.term)
    post.val = omega.target %*% Y %*% lambda.target
    syn.trajectory = as.numeric(omega.synth %*% Y)
    obs.trajectory = as.numeric(omega.target %*% Y)

    pre.time = (lambda.synth / sum(lambda.synth)) %*% time # for non-simplex weights, normalize to locate the synthetic time period on the x axis
    post.time = lambda.target %*% time 
    p = ggplot() +  geom_line(aes(x=rep(time,2), y=c(obs.trajectory, syn.trajectory), color=unit), linetype=linetype, alpha=.5) +
                    geom_vline(aes(xintercept=time[T0]), color='black', linetype=1, alpha=.2) + 
                    theme_light() + xlab('') + ylab('') + labs(color='') + theme(legend.direction = "horizontal", legend.position = "top")
    if(all(lambda.synth == 0)) { # sc case
        p = p + geom_point(aes(x=post.time, y=post.val, color=treated)) + geom_point(aes(x=post.time, y=sdid,   color=control))
    } else { 
        p = p + geom_segment(aes(x=pre.time, xend=post.time, y=col.SC, yend=post.val, color=treated)) +
                geom_point(aes(x=c(pre.time, post.time), y=c(col.SC, post.val), color=treated)) +
                geom_segment(aes(x=pre.time, xend=post.time, y=cross.term, yend=row.SC, color=control)) +
                geom_point(aes(x=c(pre.time, post.time), y=c(cross.term, yend=row.SC), color=control)) +
                geom_segment(aes(x=pre.time, xend=post.time, y=col.SC, yend=sdid, color=treated), linetype=2) +
                geom_point(aes(x=c(pre.time, post.time), y=c(col.SC, yend=sdid), color=treated), shape=21) 
    }
    p = p + geom_segment(aes(x=post.time, xend=post.time, y=post.val, yend=sdid, color=treated), linetype=3)
        
    height = (max(c(obs.trajectory,syn.trajectory))-min(c(obs.trajectory, syn.trajectory)))/4
    bottom = min(c(obs.trajectory, syn.trajectory)) - height
    density.color = 'black'
    p = p + geom_ribbon(aes(x=time[1:T0], ymin = rep(bottom,T0), ymax= bottom + height*lambda.synth[1:T0]/max(lambda.synth)),
                                color=density.color, fill=density.color, alpha = .4)  
    p
}

#' Plots treated and average control trajectories and overlays a 2x2 diff-in-diff diagram of our estimator.
#' In this overlay, the treatment effect is the length of the vertical dotted line. 
#' Requires ggplot2
#' @param est, output of did_estimate
#' @export plot.did
plot.did = function(est) {
    N0 = attr(est, 'N0') 
    T0 = attr(est, 'T0') 
    attr(est, 'lambda') = rep(1/T0, T0)
    attr(est, 'omega') = rep(1/N0, N0)
    plot.synthdid(est, control.name='average control')
}

