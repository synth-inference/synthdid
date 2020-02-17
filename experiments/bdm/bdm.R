source('../../R/synthdid.R')
library(reshape2)
library(plyr)

data = read.csv('cps_bdm.csv')
melted = melt(data, id.vars=c('state_id', 'year'))

## alex's data
Y.hours = acast(melted[melted$variable == 'uhourse',], state_id ~ year)
Y.dif.log.wage   = acast(melted[melted$variable == 'first_dif_lnwage',], state_id ~ year)
Y.prev.log.wage  = acast(melted[melted$variable == 'prev_lnwage',],      state_id ~ year)
Y.log.wage = Y.dif.log.wage + Y.prev.log.wage
Y.log.wage[,1] = Y.prev.log.wage[,2]
Y.dif.log.wage = Y.dif.log.wage[,2:ncol(Y.dif.log.wage)]

# reorder log wage data so california is in last slot
Y = Y.log.wage
T0 = which(colnames(Y) == "1998")[[1]]
N0 = 30
calif = which(rownames(Y)=='CA')
ind = c(setdiff(1:nrow(Y), calif), calif)
Y=Y[ind,]
X = Y.hours[ind,]

print(Y[nrow(Y),T0+1], digits=4) #= 6.127117 for CA 1999

sdid = synthdid_estimate(Y, N0, T0)
did = synthdid_estimate(Y, N0, T0, weights=list(omega=rep(1/N0,N0), lambda=rep(1/T0,T0)))
sc = synthdid_estimate(Y, N0, T0, omega.intercept=FALSE)
attr(sc, 'weights') = list(omega=attr(sc,'weights')$omega, lambda=rep(0,T0))

## gls
library(nlme)
reg.data = expand.grid(state = 1:nrow(Y), time = 1:ncol(Y)) 
reg.data$y = Y[t(rbind(reg.data$state, reg.data$time))]
reg.data$treated = (reg.data$state > N0) & (reg.data$time > T0)
reg.data$state.factor = as.factor(reg.data$state)
reg.data$time.factor  = as.factor(reg.data$time)
gls.fit = gls(y ~ state.factor + time.factor + treated, reg.data, corAR1(form = ~ time | state.factor))
tau.gls = gls.fit$coefficients[names(gls.fit$coefficients) == "treatedTRUE"]

estimates = list(did, sc, sdid)
names(estimates) = c('did', 'sc', 'sdid')
standard.errors = sapply(estimates, synthdid_se)
print(round(sapply(estimates, c),3))
print(round(standard.errors, 3))

summary(sdid)

## compare estimates with and without hours as covariate, placebo
sdid.hours = synthdid_estimate(Y, N0, T0, X) 
estimates = list(sdid, sdid.hours)
names(estimates) = c('log wage, no adjustment', 'log wage, adjustment for hours')
synthdid_rmse_plot(estimates)
synthdid_plot(estimates)
synthdid_placebo_plot(sdid)
synthdid_placebo_plot(sdid,overlay=TRUE)

## compare fast [holding weights fixed] vs slow [re-estimating weights] jackknife standard error estimators 
se.fast = synthdid_se(sdid)
se.jack = synthdid_se(sdid, weights=NULL)
se.fast.hours = synthdid_se(sdid.hours)
se.jack.hours = synthdid_se(sdid.hours, weights=NULL)
c(se.fast, se.jack, se.fast.hours, se.jack.hours)
