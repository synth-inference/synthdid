source('setup.R')
source('../../R/synthdid.R')
library(viridis)
library(ggplot2)
sdid = synthdid_estimate(Y, N0, T0)
sc = sc_estimate(Y, N0, T0)
did = did_estimate(Y, N0, T0)

# wi: with intercept. Makes plots remove a fraction i in [0,1] of the SDID pre/post difference. i=0 means plot as usual, i=1 means plot overlaid (without parallelogram).
wi = function(est, i) { attr(est,'intercept') = i; est } 

### Parallel Trends Plots

synthdid_plot(list(sdid=sdid), facet.vertical=FALSE, control.name='sdid control', treated.name='california', lambda.comparable=TRUE, 
    trajectory.linetype = 1, trajectory.alpha=.8, effect.alpha=1, diagram.alpha=1, effect.curvature=-.4, onset.alpha=.7) + 
    theme(legend.position=c(.90,.90), legend.direction='vertical', legend.key=element_blank(), legend.background=element_blank())
ggsave('figures/parallel-diagram.pdf', width=7, height=4.5)

synthdid_plot(list(sdid.control=sdid, synthetic.control=sc), facet=c(1,1), treated.name='california', lambda.comparable=TRUE, 
    trajectory.linetype = 1, trajectory.alpha=.8, effect.alpha=1, diagram.alpha=1, effect.curvature=-.4, onset.alpha=.7) + 
    theme(legend.position=c(.90,.90), legend.direction='vertical', legend.key=element_blank(), legend.background=element_blank())
ggsave('figures/parallel-diagram-vs-sc.pdf', width=7, height=4.5)

synthdid_plot(list(sdid=sdid, did=did, sc=sc), facet.vertical=TRUE, control.name='control', treated.name='california', lambda.comparable=TRUE, 
    trajectory.linetype = 1, trajectory.alpha=.7, effect.alpha=.7, diagram.alpha=1, effect.curvature=-.4, onset.alpha=.7) + 
    theme(legend.position=c(.90,.97), legend.direction='vertical', legend.key=element_blank(), legend.background=element_blank())
ggsave('figures/sdid-vs-vert.pdf', width=4, height=8)

estimators = list(did, sc, sdid)
names(estimators) = c('difference in differences', 'synthetic control', 'synthetic difference in differences')
synthdid_plot(estimators, facet.vertical=FALSE, control.name='control', treated.name='california', lambda.comparable=TRUE, 
    trajectory.linetype = 1, trajectory.alpha=.7, effect.alpha=.7, diagram.alpha=1, line.width=.75, effect.curvature=-.4, onset.alpha=.7) + 
    theme(legend.position=c(.26,.07), legend.direction='horizontal', legend.key=element_blank(), legend.background=element_blank(),
    strip.background=element_blank(), strip.text.x = element_blank())
ggsave('figures/parallel-vs-horiz.pdf', width=15, height=5)

# set up the box we zoom in on in plot 5
time = as.integer(colnames(Y))
lambda = attr(sdid,'weights')$lambda
xbox.ind = c(which(lambda > .01)[1], T0+4)
xbox = time[xbox.ind] + c(-.5,.5)
ybox = range(Y[N0+1, min(xbox.ind):(max(xbox.ind))]) + c(-4,4)



estimators = function(i) { l = list(wi(sdid,i), sc, sdid); names(l)=c('sdid', 'sc', ''); l }
plot.estimators = function(ests, alpha.multiplier) {
    synthdid_plot(ests, alpha.multiplier=alpha.multiplier, facet=rep(1,length(ests)),
	trajectory.linetype = 1, trajectory.alpha=.5, effect.alpha=.5, diagram.alpha=1, effect.curvature=-.4) + 
	scale_color_viridis_d(drop=FALSE) + scale_fill_viridis_d(drop=FALSE) + scale_alpha(range=c(0,1), guide='none') 
}
plot.theme = theme(legend.position=c(.9,.85), legend.direction='vertical', legend.key=element_blank(), legend.background=element_blank())
p1 = plot.estimators(estimators(0),   alpha.multiplier=c(1,.1,0)) + plot.theme
p2 = plot.estimators(estimators(.75), alpha.multiplier=c(1,.1,0)) + plot.theme
p3 = plot.estimators(estimators(1),   alpha.multiplier=c(1,.1,0)) + plot.theme
p4 = plot.estimators(estimators(1),   alpha.multiplier=c(1, 1,0)) + plot.theme
p4.zoom = p4 + coord_cartesian(xlim=xbox, ylim=ybox) + xlab('') + ylab('') + 
    theme(axis.ticks.x= element_blank(), axis.text.x = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.position='off')
p5 = p4 + annotation_custom(ggplotGrob(p4.zoom), xmin = 1968, xmax = 1984.7,   ymin=2, ymax=95) + # manually adjusted zoom box location
	  geom_rect(aes(xmin=min(xbox), xmax=max(xbox), ymin=min(ybox), ymax=max(ybox)), color=alpha('black', .25), size=.3, fill=NA)

ggsave('figures/smoking-parallel-diagram.pdf',    plot=p1, width=7, height=4.5)
ggsave('figures/smoking-parallel-close.pdf',      plot=p2, width=7, height=4.5)
ggsave('figures/smoking-overlay-diagram.pdf',     plot=p3, width=7, height=4.5)
ggsave('figures/smoking-overlay-vs-sc.pdf',       plot=p4, width=7, height=4.5)
ggsave('figures/smoking-overlay-vs-sc-inset.pdf', plot=p5, width=7, height=4.5)

## compare using time weights vs not 
sdid.notw = synthdid_estimate(Y,N0,T0,weights=list(lambda=rep(1/T0,T0)))
plot.estimators(list(sdid=wi(sdid, 0), sdid.notw=wi(sdid.notw, 0)), alpha.multiplier=c(1,1)) + plot.theme
ggsave('figures/timeweights-vs-not.pdf', width=7, height=4)

sdid.noi = synthdid_estimate(Y,N0,T0, omega.intercept=FALSE)
sdid.notwnoi = synthdid_estimate(Y,N0,T0,weights=list(lambda=rep(1/T0,T0)), omega.intercept=FALSE)
plot.estimators(list(sdid=wi(sdid,0), sdid.noi=wi(sdid.noi, 0), sdid.notwnoi=wi(sdid.notwnoi, 0)), alpha.multiplier=c(1,1,1)) + plot.theme
ggsave('figures/timeweights-vs-not-nointercept.pdf', width=7, height=4)


### State Plots (Time Weighting vs Not)

setup = attr(sdid, 'setup')
weights = attr(sdid, 'weights')
weights.sc = attr(sc, 'weights')
weights.did = attr(did, 'weights')
Y = setup$Y - contract3(setup$X, weights$beta)
N0 = setup$N0; N1 = nrow(Y)-N0
T0 = setup$T0; T1 = ncol(Y)-T0
    
lambda.did   = c(rep(1/T0, T0),  rep(0,T1))
lambda.sdid = c(weights$lambda, rep(0, T1)) 
lambda.post = c(rep(0,T0), rep(1/T1, T1))
omega.did   = c(rep(1/N0, N0),  rep(0,N1))
omega.sdid  = c(weights$omega,  rep(0, N1))
omega.post = c(rep(0,N0),  rep(1/N1, N1))
difs.weighted =    as.vector(t(omega.post) %*% Y %*% (lambda.post - lambda.sdid)) -  as.vector(Y[1:N0,] %*% (lambda.post - lambda.sdid))
difs.unweighted =  as.vector((omega.post) %*% Y %*% (lambda.post - lambda.did))   -  as.vector(Y[1:N0,] %*% (lambda.post - lambda.did))
difs.sc         =  as.vector((omega.post) %*% Y %*% lambda.post)   -  as.vector(Y[1:N0,] %*% lambda.post)
points = data.frame(y=c(difs.weighted, difs.unweighted, difs.sc),
                    state=rep(factor(rownames(Y))[1:N0] ,3),
		    estimate = rep(c(sum(weights$omega * difs.weighted), sum(weights.did$omega * difs.unweighted), sum(weights.sc$omega * difs.sc)), each=N0),
                    estimator=rep(factor(c('sdid', 'did', 'sc')), each=N0),
		    weight = c(weights$omega, weights.did$omega, weights.sc$omega))
ggplot(points[points$state != 'California', ]) + geom_point(aes(x=state,y=y,color=estimator, size=weight)) + 
    geom_hline(aes(yintercept=estimate, color=estimator)) +
    xlab('') + ylab('') + guides(shape=FALSE) + theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave('figures/individual-differences-all.pdf', width=7, height=4)
ggplot(points[points$state != 'California' & points$estimator %in% c('sdid', 'sc'), ]) + geom_point(aes(x=state,y=y,color=estimator, size=weight)) + 
    geom_hline(aes(yintercept=estimate, color=estimator)) +
    xlab('') + ylab('') + guides(shape=FALSE) + theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave('figures/individual-differences-sc.pdf', width=7, height=4)


ggplot(points[points$state != 'California', ]) + geom_point(aes(x=state,y=y,size=weight,)) + facet_grid(.~estimator) +
    geom_hline(aes(yintercept=estimate)) + scale_size(range=c(0,5)) +
    xlab('') + ylab('') + guides(shape=FALSE) + theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
    legend.background=element_blank(), legend.title = element_blank(), legend.direction='horizontal', legend.position=c(.17,.07))
    #legend.background=element_blank(), legend.direction='horizontal', legend.position=c(.16,.07))
ggsave('figures/individual-differences-horiz.pdf', width=15, height=5)



ggplot(points[points$state != 'California', ]) + 
    geom_point(aes(x=state,y=y,size=weight), data=points[points$state != 'California' & points$weight > .001, ]) +
    geom_point(aes(x=state,y=y,size=weight), data=points[points$state != 'California' & points$weight <= .001, ], shape=1) +
    geom_hline(aes(yintercept=estimate)) + facet_grid(.~estimator) + 
    xlab('') + ylab('') + guides(shape=FALSE) + theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
    legend.background=element_blank(), legend.title = element_blank(), legend.direction='horizontal', legend.position=c(.17,.07))
    #legend.background=element_blank(), legend.direction='horizontal', legend.position=c(.16,.07))
ggsave('figures/individual-differences-horiz-nofillzero.pdf', width=15, height=5)

ggplot(points[points$state != 'California', ]) + 
    geom_point(aes(x=state,y=y,size=weight), data=points[points$state != 'California' & points$weight > .001, ]) +
    geom_point(aes(x=state,y=y,size=weight), data=points[points$state != 'California' & points$weight <= .001, ], shape=4, show.legend=FALSE) +
    geom_hline(aes(yintercept=estimate), size=.75) + facet_grid(.~estimator) + 
    xlab('') + ylab('') + guides(shape=FALSE) + theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
    legend.background=element_blank(), legend.title = element_blank(), legend.direction='horizontal', legend.position=c(.17,.07), 
    strip.background=element_blank(), strip.text.x = element_blank())
ggsave('figures/individual-differences-horiz-xzero.pdf', width=15, height=5)


