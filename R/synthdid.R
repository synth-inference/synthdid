contract3 = function(X,v) {
    stopifnot(length(dim(X)) == 3, dim(X)[3] == length(v))
    out = array(0, dim=dim(X)[1:2])
    if(length(v) == 0) { return(out) }
    for(ii in 1:length(v)) {
        out = out + v[ii] * X[,,ii]
    }
    return(out)
}

# frank-wolfe + gradient solver
sc.weight.fw.covariates = function(Y, X=array(0,dim=c(dim(Y),0)), zeta.lambda = 0, zeta.omega=0,
                                   min.step.size=1e-3, max.iter=1000, step='linesearch', start.beta=rep(0,dim(X)[3]),
                                   lambda=NULL, omega=NULL, beta = NULL) {
    stopifnot(length(dim(Y))==2 && length(dim(X)) == 3 && all(dim(Y)==dim(X)[1:2]))
    T0 = ncol(Y)-1
    N0 = nrow(Y)-1
    if(length(dim(X)) == 2) { dim(X) = c(dim(X),1) }
    if(is.null(lambda)) {
        lambda=rep(0,T0)
        lambda[1]=1
    }    
    if(is.null(omega)) {
        omega=rep(0,N0)
        omega[1]=1
    } 
    if(is.null(beta)) {
        beta=rep(0,dim(X)[3])
    } 

    t=0
    step.size = Inf
    vals = rep(NA, max.iter)
    while(t < max.iter && step.size > min.step.size) {
        t=t+1
        X.beta = contract3(X,beta)
        Y.beta = Y - X.beta
        err.lambda = Y.beta[1:N0,] %*% c(lambda, -1) / N0
        err.omega  = t(Y.beta[,1:T0]) %*% c(omega, -1) / T0
        grad.lambda = t(err.lambda) %*% Y.beta[1:N0,1:T0]     + Re(zeta.lambda^2)  *  lambda
        grad.omega  = t(err.omega)  %*% t(Y.beta[1:N0,1:T0])  + Re(zeta.omega^2)   *  omega
        grad.beta = if(dim(X)[3]==0) { c() } else {
            -apply(X, 3, function(Xi) { 
                t(err.lambda) %*% Xi[1:N0,] %*% c(lambda,-1)   + 
                t(err.omega) %*% t(Xi[,1:T0]) %*% c(omega, -1) 
            })
        }
        i.lambda = which.min(grad.lambda)
        i.omega  = which.min(grad.omega)

        vals[t] = Re(zeta.omega^2) * sum(omega^2) + Re(zeta.lambda^2) * sum(lambda^2) + N0*sum(err.omega^2) + T0*sum(err.lambda^2)
        step.size = 2/(t+2)
        lambda = lambda*(1-step.size); lambda[i.lambda] = lambda[i.lambda] + step.size
        omega  = omega*(1-step.size);  omega[i.omega]   = omega[i.omega]   + step.size
        beta   = beta - grad.beta * step.size * (t >= start.beta)
    }
    list(lambda=lambda, omega=omega, beta=beta, vals=vals)
}

collapsed.form = function(Y, N0, T0) {
    N = nrow(Y); T=ncol(Y)
    rbind(cbind(  Y[1:N0,1:T0, drop=FALSE],                    rowMeans(Y[1:N0,(T0+1):T, drop=FALSE])),
          cbind(  t(colMeans(Y[(N0+1):N,1:T0, drop=FALSE])),   mean(Y[(N0+1):N,(T0+1):T, drop=FALSE])))
}

#' Computes synthetic diff-in-diff estimate for an average treatment effect on a treated block.
#' See Section 4.1 of the paper. 
#' @param Y, the observation matrix.
#' @param N0, the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0, the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param X, an optional 3-D array of time-varying covariates. Shape should be N X T X C for C covariates.
#' @param zeta.lambda. Its square is weight of the ridge penalty relative to MSE. Defaults to 0.
#' @param zeta.omega. Analogous for omega. Defaults to 0.
#' @return An average treatment effect estimate, 'weights' and 'setup' attached as attributes.
#'         Weights contains the estimated weights lambda and omega and corresponding intercepts.
#'         If covariates X are passedas well as regression coefficients beta if X is passed
#'         Setup is a list describing the problem passed in: Y, N0, T0, X. 
#' @export synthdid_estimate
synthdid_estimate <- function(Y, N0, T0, X=array(dim=c(dim(Y),0)),
                              zeta.lambda=0, zeta.omega=0,
                              lambda.intercept=TRUE, omega.intercept=TRUE, 
                              weights = NULL, start.beta = 0) {
    stopifnot(nrow(Y) > N0, ncol(Y) > T0, length(dim(X)) %in% c(2,3), dim(X)[1:2] == dim(Y))
    if(length(dim(X)) == 2) { dim(X) = c(dim(X),1) }

    Yc=collapsed.form(Y, N0, T0)
    Xc=apply(X, 3, function(Xi) { collapsed.form(Xi, N0, T0) })
    dim(Xc) = c(dim(Yc), dim(X)[3])
    # add intercept dummies
    if(lambda.intercept) {
        Xi = array(0,dim=dim(Yc))
        Xi[,T0+1] = 1
        Xc = array(c(Xi,Xc), dim=dim(Xc)+c(0,0,1))
    }
    if(omega.intercept) {
        Xi = array(0,dim=dim(Yc))
        Xi[N0+1,] = 1
        Xc = array(c(Xi,Xc), dim=dim(Xc)+c(0,0,1))
    }
    is.intercept = rep(FALSE, dim(Xc)[3])
    if(lambda.intercept || omega.intercept) { is.intercept[1] = TRUE }
    if(lambda.intercept && omega.intercept) { is.intercept[2] = TRUE }

    if(is.null(weights)) { 
        start.betas = rep(start.beta, dim(Xc)[3])
        start.betas[is.intercept] = 0
        weights = sc.weight.fw.covariates(Yc, Xc, zeta.lambda=zeta.lambda, zeta.omega=zeta.omega, start.beta=start.betas)
        weights.no = sc.weight.fw.covariates(Yc, Xc[,,is.intercept], zeta.lambda=zeta.lambda, zeta.omega=zeta.omega, start.beta=start.betas[is.intercept])
     }
    
    X.beta = contract3(Xc, weights$beta)
    est = t(c(weights$omega, -1)) %*% (Yc - X.beta) %*% c(weights$lambda, -1)

    # separate intercepts from regression coefficients
    if(omega.intercept)   { weights$omega0  = weights$beta[1] }
    if(lambda.intercept) { weights$lambda0 = weights$beta[sum(is.intercept)] }
    weights$beta = weights$beta[!is.intercept]

    class(est) = 'synthdid'
    attr(est, 'weights') = weights
    attr(est, 'setup') = list(Y=Y, X=X, N0=N0, T0=T0)
    return(est)
}


#' Plots treated and synthetic control trajectories and overlays a 2x2 diff-in-diff diagram of our estimator.
#' In this overlay, the treatment effect is indicated by an arrow.
#' The weights lambda defining our synthetic pre-treatment time period are plotted below.
#' Requires ggplot2
#' @param estimates, a list of estimates output by synthdid_estimate. Or a single estimate.
#' @param treated.name, the name of the treated curve that appears in the legend. Defaults to 'treated'
#' @param control.name, the name of the control curve that appears in the legend. Defaults to 'synthetic control'
#' @export synthdid_plot
synthdid_plot = function(estimates, treated.name='treated', control.name='synthetic control') {
    library(ggplot2)
    if(is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }
    
    treated = 1
    control = 2    
    groups = factor(c(treated, control), labels=c(treated.name, control.name))
    estimate.factors = factor(1:length(estimates), labels=names(estimates))

    plot.descriptions = lapply(estimates, function(est) { 
        setup = attr(est, 'setup')
        weights = attr(est, 'weights')
        Y = setup$Y - contract3(setup$X, weights$beta)
        N0 = setup$N0; N1 = nrow(Y)-N0
        T0 = setup$T0; T1 = ncol(Y)-T0

        omega.synth  = c(weights$omega,  rep(0, N1))
        lambda.synth = c(weights$lambda, rep(0, T1))
        omega.target = c(rep(0,N0),  rep(1/N1, N1))
        lambda.target = c(rep(0,T0), rep(1/T1, T1))

        treated.post   = omega.target %*% Y %*% lambda.target
        treated.pre    = omega.target %*% Y %*% lambda.synth
        control.post = omega.synth %*% Y %*% lambda.target
        control.pre  = omega.synth %*% Y %*% lambda.synth 
        sdid.post = as.numeric(control.post + treated.pre - control.pre)
        syn.trajectory = as.numeric(omega.synth %*% Y)
        obs.trajectory = as.numeric(omega.target %*% Y)
        
        time = (1-T0):T1
        pre.time =  lambda.synth  %*% time
        post.time = lambda.target %*% time

        # construct objects on graph
        lines  = data.frame(x = rep(time,2), 
                            y = c(obs.trajectory, syn.trajectory), 
                            color=rep(groups[c(treated, control)], each=length(time)))
        points = data.frame(x     =        c(pre.time,    pre.time,    post.time,     post.time,     post.time),  
                            y     =        c(treated.pre, control.pre, treated.post,  control.post,  treated.post),    
                            color = groups[c(treated,     control,     treated,       control,       treated)])   
        segments = data.frame(x    =        c(pre.time,     pre.time),      
                              xend =        c(post.time,    post.time),      
                              y    =        c(control.pre,  treated.pre),  
                              yend =        c(control.post, treated.post),   
                              color= groups[c(control, treated)])
        constructed.points   = data.frame(x = post.time, y=sdid.post)
        constructed.segments = data.frame(x = pre.time, xend = post.time, y = treated.pre, yend = sdid.post)
        faint.segments      = data.frame(x    = c(pre.time,    post.time),    
                                         xend = c(pre.time,    post.time),    
                                         y    = c(control.pre, control.post), 
                                         yend = c(treated.pre, sdid.post))    
        arrows = data.frame(x=post.time, xend = post.time, y=sdid.post, yend=treated.post)
        vlines = data.frame(xintercept=time[T0])

        height = (max(c(obs.trajectory,syn.trajectory))-min(c(obs.trajectory, syn.trajectory)))/4
        bottom = min(c(obs.trajectory, syn.trajectory)) - height
        density.color = 'black'
        ribbons = data.frame(x=time[1:T0], ymin = rep(bottom,T0), ymax= bottom + height*lambda.synth[1:T0]/max(lambda.synth))

        list(lines=lines, points=points, segments=segments, 
             constructed.points=constructed.points, constructed.segments=constructed.segments, faint.segments=faint.segments,
             arrows=arrows, vlines=vlines, ribbons=ribbons)
    })
    
   
    concatenate.field = function(descs, names, field) {
        do.call(rbind, mapply(function(desc, name) { 
                element = desc[[field]]
                element$estimate = name
                element
            }, descs, names, SIMPLIFY=FALSE))
    }
    conc = lapply(names(plot.descriptions[[1]]), function(field) { concatenate.field(plot.descriptions, estimate.factors, field) })
    names(conc) = names(plot.descriptions[[1]])
    

    p=ggplot() + facet_grid(estimate ~ ., scales='free_y') +
        geom_line(aes(x=x,y=y,color=color),  data=conc$lines) +
        geom_point(aes(x=x,y=y,color=color), data=conc$points) +
        geom_segment(aes(x=x,xend=xend,y=y,yend=yend,color=color), data=conc$segments) +
        geom_point(aes(x=x,y=y), color='black', data=conc$constructed.points, shape=21) + 
        geom_segment(aes(x=x,xend=xend,y=y,yend=yend), data=conc$constructed.segments, linetype=2) +
        geom_segment(aes(x=x,xend=xend,y=y,yend=yend), data=conc$faint.segments, alpha=.1, linetype=3) +
        geom_segment(aes(x=x,xend=xend,y=y,yend=yend), data=conc$arrows, color='black', alpha=.2, size=1, arrow=arrow(length=unit(.2, 'cm'))) +
        geom_vline(aes(xintercept=xintercept), data=conc$vlines, color='black', alpha=.4) + 
        geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax), data = conc$ribbons, color='black', fill='black', alpha=.3) + 
        xlab('') + ylab('') + labs(color='') + 
        theme_light() + theme(legend.direction = "horizontal", legend.position = "top") 
                
        p
}
# make a one-estimate variant available as plot(estimate)
plot.synthdid = synthdid_plot



#' A diagnostic plot for sc.weight.fw.covariates. Plots the objective function, regularized RMSE,
#' as a function of the number of Frank-Wolfe / Gradient steps taken.
#' Requires ggplot2
#' @param estimates, a list of estimates output by synthdid_estimate. Or a single estimate.
#' @export synthdid_rmse_plot
synthdid_rmse_plot = function(estimates) { # pass an estimate or list of estimates
    library(ggplot2)
    if(class(estimates) == 'synthdid') { estimates = list(estimates) } 
    if(is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }
    rmse = lapply(estimates, function(est) { sqrt(attr(est, 'weights')$vals) })
    plot.data = data.frame(rmse = unlist(rmse),         
                           iteration=unlist(lapply(rmse, function(vals) { 1:length(vals) })),
                           method = unlist(mapply(function(vals, name) { rep(factor(name), length(vals)) }, rmse, names(estimates), SIMPLIFY=FALSE)))
    ggplot(plot.data) + geom_line(aes(x=iteration, y=rmse, color=method)) + scale_y_log10() 
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
    if(class(estimates) == 'synthdid') { estimates = list(estimates) } 
    if(is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }
    
    omegas = do.call(cbind, lapply(estimates, function(est) { attr(est, 'weights')$omega }))
    if(length(dim(omegas))==1) { dim(omegas) = c(dim(omegas), 1) }

    Y = attr(estimates[[1]], 'setup')$Y
    o = rev(order(omegas[,sort.by]))
    tab = round(omegas[o,], digits=digits)
    rownames(tab) = rownames(Y)[o]
    colnames(tab) = names(estimates)
    
    # truncate table to retain a weight sum of at least mass for each unit
    tab.len = min(apply(tab, 2, function(col) { Position(function(x){ x >= mass }, cumsum(col)) }))
    tab[1:tab.len, ]
}
