devtools::install('susanathey/MCPanel')
library(synthdid)
library(MCPanel)


## define estimators
mc_estimate = function(Y, N0, T0) {
    N1=nrow(Y)-N0
    T1=ncol(Y)-T0
    W <- outer(c(rep(0,N0),rep(1,N1)),c(rep(0,T0),rep(1,T1)))
    mc_pred <- mcnnm_cv(Y, 1-W, num_lam_L = 20)
    mc_fit  <- mc_pred$L + outer(mc_pred$u, mc_pred$v, '+')
    mc_est <- sum(W*(Y-mc_fit))/sum(W)
    mc_est
}
difp_estimate = function(Y, N0, T0) { 
    synthdid_estimate(Y, N0, T0, weights=list(lambda=rep(1/T0, T0)))  
} 	

estimators = list(sdid=synthdid_estimate, sc=sc_estimate, did=did_estimate, mc=mc_estimate, difp=difp_estimate)


## load data and define simulators
last.col = function(X) { X[, ncol(X)] }


data(CPS)
Y.logwage      = panel.matrices(CPS, treatment='min_wage', outcome='log_wage', treated.last=FALSE)$Y 
Y.hours        = panel.matrices(CPS, treatment='min_wage', outcome='hours',    treated.last=FALSE)$Y
Y.urate        = panel.matrices(CPS, treatment='min_wage', outcome='urate',    treated.last=FALSE)$Y
w.minwage      = last.col(panel.matrices(CPS, treatment='min_wage',   treated.last=FALSE)$W)
w.gunlaw       = last.col(panel.matrices(CPS, treatment='open_carry', treated.last=FALSE)$W)
w.abortion     = last.col(panel.matrices(CPS, treatment='abort_ban',  treated.last=FALSE)$W)

data(PENN)
Y.loggdp       = panel.matrices(PENN, treatment='dem', outcome='log_gdp')$Y
w.democracy    = last.col(panel.matrices(PENN, treatment='dem')$W)
w.education    = last.col(panel.matrices(PENN, treatment='educ')$W)

default=list(rank = 4, N1 = 10, T1 = 10)
cps.baseline.params  = estimate.dgp(Y.logwage, w.minwage, default$rank)

### This function creates a simulator --- a no-arg function returning a simulated dataset --- based on a modification of an estimated dgp.
##  To specify the dgp parameters, we pass in a set of baseline parameters as estimated by estimate.dgp
##  and a set of functions F, M, etc with names corresponding to those parameters which compute
##  the actual parameters for the simulator as functions of the corresponding baseline parameter
# The simulator actually uses Sigma, not ar_coef, to generate noise. 
# Updating ar_coef doesn't affect the simulations but helps us show the correct coefficients in our summary table
simulator = function(params = cps.baseline.params,
		     F=function(x){x}, M=function(x){x}, Sigma = function(x){x}, pi = function(x){x}, ar_coef = function(x){x},
		     N1=default$N1, T1=default$T1) {
    updated.params = list(F=F(params$F), M=M(params$M), Sigma=Sigma(params$Sigma), pi = pi(params$pi), ar_coef=ar_coef(params$ar_coef))
    list(params=updated.params, N1=N1, T1=T1,
	 run=function() { simulate.dgp(updated.params, N1, T1) })
}    

simulators = list(
    baseline   =  simulator(),
    # Modified outcome model
    no.corr    =  simulator(Sigma=function(Sigma) { diag(nrow(Sigma)) * norm(Sigma,'f')/sqrt(nrow(Sigma))}, ar_coef=function(coefs) { 0*coefs }),
    no.M       =  simulator(M=function(M) { 0*M }), 
    no.F       =  simulator(F=function(F) { 0*F }),
    only.noise =  simulator(M=function(M) { 0*M }, F=function(F) {0*F}),
    no.noise   =  simulator(Sigma=function(Sigma) { 0*Sigma }, ar_coef=function(coefs) { 0*coefs }),
    # Modified assignment process
    gun.law    =  simulator(estimate.dgp(Y.logwage, w.gunlaw,   default$rank)), 
    abortion   =  simulator(estimate.dgp(Y.logwage, w.abortion, default$rank)),
    random     =  simulator(pi=function(pi) { rep(.5,length(pi)) }),
    # Modified outcome variable
    hours      =  simulator(estimate.dgp(Y.hours,   w.minwage,  default$rank)),
    urate      =  simulator(estimate.dgp(Y.urate,   w.minwage,  default$rank)),
    # Assignment block size
    T1.is.one  = simulator(T1=1),
    N1.is.one  = simulator(N1=1),
    blocksize.is.one = simulator(T1=1, N1=1),
    # penn world table
    penn.democracy  = simulator(estimate.dgp(Y.loggdp,  w.democracy, default$rank)),
    penn.education  = simulator(estimate.dgp(Y.loggdp,  w.education, default$rank)),
    penn.random     = simulator(estimate.dgp(Y.loggdp,  w.education, default$rank), pi=function(pi) { rep(.5,length(pi)) }))

## run simulations to generate Tables 2 and 3
repetitions = 10
set.seed(12345)
estimate = array(NA, dim=c(repetitions, length(simulators), length(estimators)))
for(rr in 1:repetitions) { 
    for(ss in 1:length(simulators)) {
	setup = simulators[[ss]]$run()
	for(ee in 1:length(estimators)) {
	    estimate[rr,ss,ee] = estimators[[ee]](setup$Y, setup$N0, setup$T0)
	}
    }
}
rmse = apply(estimate, c(2,3), function(x) { sqrt(mean(x^2)) })
bias = apply(estimate, c(2,3), function(x) { mean(x) })

## collect info about simulator designs
sim.info = t(sapply(simulators, function(sim) { 
    data.frame(F.scale  =  round(sqrt(mean(sim$params$F^2)), 3),
               M.scale  =  round(sqrt(mean(sim$params$M^2)), 3),
	       noise.scale  =  round(sqrt(norm(sim$params$Sigma, '2')), 3),
	       ar.coefs    =  sprintf('(%.2f, %.2f)', sim$params$ar_coef[1], sim$params$ar_coef[2]))
}))

# display summary table
display.table = cbind(sim.info, as.data.frame(round(cbind(rmse, bias), 3)))
rownames(display.table) = names(simulators)
colnames(display.table) = c('F', 'M', 'Sigma', 'AR(2)', 
			    sprintf('rmse:%s', names(estimators)), 
			    sprintf('bias:%s', names(estimators)))
display.table
