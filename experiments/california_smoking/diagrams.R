source('setup.R')
source('../../R/synthdid.R')
library(viridis)
sdid = synthdid_estimate(Y, N0, T0)
sc = sc_estimate(Y, N0, T0)
did = did_estimate(Y, N0, T0)
sdid.i = function(i) { est = sdid; attr(est,'intercept') = i; est }


synthdid_plot(list(sdid=sdid.i(0), sc=sc, junk=sdid), facet=c(1,1,1), lambda.comparable=TRUE, 
    trajectory.linetype = 1, trajectory.alpha=.5, effect.alpha=.5, diagram.alpha=1, alpha.multiplier=c(1,.1,0), effect.curvature=-.4) + theme(legend.position='off') + scale_color_viridis_d()
ggsave('figures/smoking-parallel-diagram.pdf', width=7, height=4)
synthdid_plot(list(sdid=sdid.i(.75), sc=sc, junk=sdid), facet=c(1,1,1), lambda.comparable=TRUE, 
    trajectory.linetype = 1, trajectory.alpha=.5, effect.alpha=.5, diagram.alpha=1, alpha.multiplier=c(1,.1,0), effect.curvature=-.4) + theme(legend.position='off') + scale_color_viridis_d()
ggsave('figures/smoking-parallel-close.pdf', width=7, height=4)
synthdid_plot(list(sdid=sdid.i(1), sc=sc, junk=sdid), facet=c(1,1,1), lambda.comparable=TRUE, 
    trajectory.linetype = 1, trajectory.alpha=.5, effect.alpha=.5, diagram.alpha=1, alpha.multiplier=c(1,.1,0), effect.curvature=-.4) + theme(legend.position='off') + scale_color_viridis_d()
ggsave('figures/smoking-overlay-diagram.pdf', width=7, height=4)
synthdid_plot(list(sdid=sdid.i(1), sc=sc, junk=sdid), facet=c(1,1,1), lambda.comparable=TRUE, 
    trajectory.linetype = 1, trajectory.alpha=.5, effect.alpha=.5, diagram.alpha=1, alpha.multiplier=c(1,1,0), effect.curvature=-.4) + theme(legend.position='off') + scale_color_viridis_d()
ggsave('figures/smoking-overlay-vs-sc.pdf', width=7, height=4)


estimate = sdid
setup = attr(estimate, 'setup')
weights = attr(estimate, 'weights')
Y = setup$Y - contract3(setup$X, weights$beta)
N0 = setup$N0; N1 = nrow(Y)-N0
T0 = setup$T0; T1 = ncol(Y)-T0
    
lambda.did   = c(rep(1/T0, T0),  rep(0,T1))
lambda.synth = c(weights$lambda, rep(0, T1)) 
lambda.target = c(rep(0,T0), rep(1/T1, T1))
omega.synth  = c(weights$omega,  rep(0, N1))
omega.target = c(rep(0,N0),  rep(1/N1, N1))

points = data.frame(y=c(Y %*% lambda.target, Y %*% lambda.synth, Y %*% lambda.did),
                    state=rep(factor(rownames(Y)),3),
		    ypost = rep(Y %*% lambda.target, 3),
                    color=rep(factor(c('post', 'synth', 'did')), each=nrow(Y)))
ggplot(points) + geom_point(aes(x=state,y=y,color=color)) + theme(axis.text.x = element_text(angle = 90, hjust = 1));        ggsave('figures/time-parallel.pdf')
ggplot(points) + geom_point(aes(x=state,y=y-ypost,color=color)) + theme(axis.text.x = element_text(angle = 90, hjust = 1));  ggsave('figures/time-parallel-relative.pdf')
ggplot(points[points$state != 'California', ]) + geom_boxplot(aes(x=color, y=y-ypost));                                      ggsave('figures/time-parallel-boxplot.pdf')

Y = Y[omega.synth+omega.target > .02, ]
points = data.frame(y=c(Y %*% lambda.target, Y %*% lambda.synth, Y %*% lambda.did),
                    state=rep(factor(rownames(Y)),3),
		    ypost = rep(Y %*% lambda.target, 3),
                    color=rep(factor(c('post', 'synth', 'did')), each=nrow(Y)))
ggplot(points) + geom_point(aes(x=state,y=y,color=color)) + theme(axis.text.x = element_text(angle = 90, hjust = 1));        ggsave('figures/sc-time-parallel.pdf')
ggplot(points) + geom_point(aes(x=state,y=y-ypost,color=color)) + theme(axis.text.x = element_text(angle = 90, hjust = 1));  ggsave('figures/sc-time-parallel-relative.pdf')
ggplot(points[points$state != 'California', ]) + geom_boxplot(aes(x=color, y=y-ypost));                                      ggsave('figures/sc-time-parallel-boxplot.pdf')


source('../../R/synthdid.R')
synthdid_time_plot(sdid)
