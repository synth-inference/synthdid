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
#' @param line.width the line width.
#' @param guide.linetype (undocumented for now)
#' @param point.size (undocumented for now)
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
			 lambda.plot.scale=3, trajectory.linetype=1, effect.curvature = 0, line.width=.5, guide.linetype=2, point.size=.5,
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
        if(length(time) == 0 || !all(is.finite(time))) { time = 1:(T0+T1) }
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
        geom_line(aes(x=x,y=y,color=color,frame=frame, alpha=trajectory.alpha*show),  data=conc$lines, linetype=trajectory.linetype, size=line.width) +
        geom_point(aes(x=x,y=y,color=color,frame=frame, alpha=diagram.alpha*show), data=conc$points, shape=21, size=point.size) +
        geom_point(aes(x=x,y=y,color=color,frame=frame, alpha=diagram.alpha*show), data=no.sc(conc$did.points), size=point.size) +
        geom_segment(aes(x=x,xend=xend,y=y,yend=yend,color=color, frame=frame, alpha=diagram.alpha*show), data=no.sc(conc$did.segments), size=line.width) +
        geom_segment(aes(x=x,xend=xend,y=y,yend=yend,frame=frame, group=estimate, alpha=.6*diagram.alpha*show), data=no.sc(conc$hallucinated.segments), linetype=guide.linetype,  size=line.width, color='black') +
        geom_segment(aes(x=x,xend=xend,y=y,yend=yend,frame=frame, group=estimate, alpha=.5*diagram.alpha*show), data=no.sc(conc$guide.segments), size=line.width, linetype=guide.linetype, color='black') +
        geom_vline(aes(xintercept=xintercept, alpha=onset.alpha*show), data=conc$vlines, size=line.width, color='black') +
        geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax, group=color, fill=color, alpha=.5*diagram.alpha*show), color='black', data = conc$ribbons, size=line.width, show.legend=FALSE) +
      	geom_curve(aes(x=x,xend=xend,y=y,yend=yend, alpha=effect.alpha*show),  data=conc$arrows, curvature=effect.curvature, color='black', size=line.width, arrow=arrow(length=unit(.2, 'cm')))
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

#' Plots unit by unit difference-in-differences
#' Requires ggplot2
#' @param estimates as output by synthdid_estimate. Can be a single one or a list of them.
#' @param show.ci If TRUE, plots horizontal lines for 95\% CI as well as the point estimate. Defaults to FALSE.
#' @param negligible.threshold Unit weight threshold below which units are plotted as small, transparent xs instead of circles. Defaults to .001.
#' @param negligible.alpha Determines transparency of those xs.
#' @export synthdid_units_plot
synthdid_units_plot = function(estimates, show.ci=FALSE, negligible.threshold = .001, negligible.alpha = .3) {
    if(class(estimates) == 'synthdid_estimate') { estimates = list(estimates) }
    if(is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }
    plot.data = do.call(rbind, lapply(1:length(estimates), function(ee) {
	estimate = estimates[[ee]]
	setup = attr(estimate, 'setup')
	weights = attr(estimate, 'weights')
	Y = setup$Y - contract3(setup$X, weights$beta)
	N0 = setup$N0; N1 = nrow(Y)-N0
	T0 = setup$T0; T1 = ncol(Y)-T0

	lambda.pre = c(weights$lambda, rep(0, T1))
	lambda.post = c(rep(0,T0), rep(1/T1, T1))
	omega.control  = c(weights$omega,  rep(0, N1))
	omega.treat = c(rep(0,N0),  rep(1/N1, N1))
	difs = as.vector(t(omega.treat) %*% Y %*% (lambda.post - lambda.pre)) -  as.vector(Y[1:N0,] %*% (lambda.post - lambda.pre))
	se = if(show.ci) { synthdid_se(estimate) } else { NA }
	data.frame(y=difs, unit=rownames(Y)[1:N0], weight = omega.control[1:N0],
	           estimate = c(estimate), se = se, estimator = names(estimates)[[ee]])
    }))
    p = ggplot(plot.data) +
	geom_point(aes(x=unit,y=y,size=weight), data=plot.data[plot.data$weight >  negligible.threshold, ]) +
	geom_point(aes(x=unit,y=y,size=weight), data=plot.data[plot.data$weight <= negligible.threshold, ], alpha=negligible.alpha, shape=4,  show.legend=FALSE) +
	geom_hline(aes(yintercept=estimate), size=.75)
    if(show.ci) {
	p = p + geom_hline(aes(yintercept=estimate-1.96*se), size=.5, alpha=.5) +
		geom_hline(aes(yintercept=estimate+1.96*se), size=.5, alpha=.5)
    }
    p + facet_grid(.~estimator) + xlab('') + ylab('') + guides(shape=FALSE) +
	theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
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

#' Plot a synthdid object
#' @param x The object to plot
#' @param ... Additional arguments (currently ignored).
#' @method plot synthdid_estimate
#' @export
plot.synthdid_estimate = function(x, ...) {
 synthdid_plot(x, ...)
}
