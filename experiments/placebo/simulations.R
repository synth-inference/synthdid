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
sdid_reg = function(Y, N0, T0) {
    N1 = nrow(Y) - N0
    T1 = ncol(Y) - T0
    sigma = sd(apply(Y[1:N0,1:T0], 1, diff))
    synthdid_estimate(Y, N0, T0, zeta.omega = (N1*T1)^(1/4)*sigma)
}

estimators = list(sdid=synthdid_estimate, sc=sc_estimate, did=did_estimate, 
		  mc=mc_estimate, difp=difp_estimate, sdid_reg=sdid_reg)


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
Y.loggdp       = panel.matrices(PENN, treatment='dem', outcome='log_gdp', treated.last=FALSE)$Y
w.democracy    = last.col(panel.matrices(PENN, treatment='dem',  treated.last=FALSE)$W)
w.education    = last.col(panel.matrices(PENN, treatment='educ', treated.last=FALSE)$W)

default=list(rank = 4, N1 = 10, T1 = 10)
cps.baseline.params  = estimate.dgp(Y.logwage, w.minwage, default$rank)

### This function creates a simulator --- a no-arg function returning a simulated dataset --- 
### based on a modification of dgp specification as returned by estimate.dgp.
##  -To specify the dgp parameters, we pass in a set of baseline parameters as estimated by estimate.dgp
##   and a set of functions F, M, etc with names corresponding to those parameters which compute
##   the actual parameters for the simulator as functions of the corresponding baseline parameter
##  -If a numeric value of resample is passed, in each simulation resample that number of units, 
##   i.e., rows from M + F, to use in the simulation design. This is used in some rows of Table 4.
# The simulator $run uses Sigma, not ar_coef, to generate noise, so updating ar_coef doesn't affect the simulations.
# We do it here so the correct coefficients is included in $params and therefore are shown in our summary table
simulator = function(params = cps.baseline.params,
		     F=function(x){x}, M=function(x){x}, Sigma = function(x){x}, pi = function(x){x}, ar_coef = function(x){x},
		     N1=default$N1, T1=default$T1, resample=NULL) {
    updated.params = list(F=F(params$F), M=M(params$M), Sigma=Sigma(params$Sigma), pi = pi(params$pi), ar_coef=ar_coef(params$ar_coef))
    list(params=updated.params, N1=N1, T1=T1,
	 run=function() { 
	    if(!is.numeric(resample)) {
		simulate.dgp(updated.params, N1, T1) 
	    } else { 	        
		ind = sample.int(nrow(updated.params$F), size=resample, replace=TRUE)
		resampled.params=updated.params
		resampled.params$F=updated.params$F[ind,]
		resampled.params$M=updated.params$M[ind,]
		simulate.dgp(resampled.params, N1, T1) 
	    } 
	})
}    

simulators = list(
    baseline   =  simulator(),
    # Modified outcome model
    no.corr    =  simulator(Sigma=function(Sigma) { diag(nrow(Sigma)) * norm(Sigma,'f')/sqrt(nrow(Sigma))}, 
			    ar_coef=function(coefs) { 0*coefs }),
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
    # Resample rows (from Table 4)
    resample.200 = simulator(resample=200, N1=20),
    resample.400 = simulator(resample=400, N1=40),
    # penn world table (Table 3)
    penn.democracy  = simulator(estimate.dgp(Y.loggdp,  w.democracy, default$rank)),
    penn.education  = simulator(estimate.dgp(Y.loggdp,  w.education, default$rank)),
    penn.random     = simulator(estimate.dgp(Y.loggdp,  w.education, default$rank), 
				pi=function(pi) { rep(.5,length(pi)) }))

## collect summary info about simulator designs to include in tables
sim.info = do.call(rbind, lapply(simulators, function(sim) { 
    data.frame(F.scale  =  sqrt(mean(sim$params$F^2)), 
               M.scale  =  sqrt(mean(sim$params$M^2)),
	       noise.scale  =  sqrt(norm(sim$params$Sigma, '2')),
	       ar.lag1      =  sim$params$ar_coef[1],
	       ar.lag2      =  sim$params$ar_coef[2])
}))

## run simulations to estimate rmse and bias for Tables 2 and 3 and coverage for Table 4
library(abind)
library(doParallel)
library(doRNG)
cl = makeCluster(8, outfile='simulations.log')
registerDoParallel(cl, cores=8)

sim.replications  = 1000
coverage.replications = 400
set.seed(12345)
estimates = foreach(rr = 1:sim.replications, .combine=abind) %dorng% {
    library(synthdid)
    library(MCPanel)
    print(rr)

    estimate = array(NA, dim=c(length(simulators), length(estimators), 4, 1))
    for(ss in 1:length(simulators)) {
	setup = simulators[[ss]]$run()
	for(ee in 1:length(estimators)) {
	    est = estimators[[ee]](setup$Y, setup$N0, setup$T0)
	    estimate[ss,ee,1,1] = est 
	    if(names(estimators)[ee] %in% c('sdid', 'sc', 'did') && rr <= coverage.replications) {
		estimate[ss,ee,2,1] = sqrt(vcov(est, method="bootstrap")) 
		estimate[ss,ee,3,1] = sqrt(vcov(est, method="jackknife"))
		estimate[ss,ee,4,1] = sqrt(vcov(est, method="placebo")) 
	    }
	}
    }
    estimate
}
stopCluster(cl)

rmse = apply(estimates, c(1,2), function(x) { sqrt(mean(x[1,]^2)) })
bias = apply(estimates, c(1,2), function(x) { mean(x[1,]) })
coverage = lapply(2:4, function(se) {
    apply(estimates, c(1,2), function(x) { mean(abs(x[1,] / x[se,]) <= qnorm(.975))  })
})

# Display a point estimate summary table. 
# In the paper, this is split into Tables 2 and 3, with CPS sims in the former and PENN sims in the latter.
point.table = cbind(sim.info, as.data.frame(cbind(rmse, bias)))
rownames(point.table) = names(simulators)
colnames(point.table) = c('F', 'M', 'Sigma', 'ar:lag1', 'ar:lag2', 
			    sprintf('rmse:%s', names(estimators)), 
			    sprintf('bias:%s', names(estimators)))
print(point.table, digits=2)

# Display a coverage summary table. 
# A subset of these rows are in Table 4 of the paper.
coverage.table = do.call(cbind, coverage)
rownames(coverage.table) = names(simulators)
colnames(coverage.table) = c(sprintf('bootstrap:%s',  names(estimators)),
			     sprintf('jackknife:%s',  names(estimators)),
			     sprintf('placebo:%s',    names(estimators)))
print(coverage.table[,!is.na(coverage.table[1,])], digits=2)

save('estimates', file='sims.RData')
write.table(round(coverage.table,2), file='coverage_table.tab')
write.table(round(point.table,2),    file='point_table.tab')


# Plot the error distribution of the estimates. 
# The first plot, for the baseline CPS setting, is Figure 2 of the paper. 
# The rest are analogous for different dgps.
## TODO
