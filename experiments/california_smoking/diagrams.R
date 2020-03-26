source('../../R/synthdid.R')
data.raw = read.table("MLAB_data.txt")

STATE.NAME = c("Alabama", "Arkansas", "Colorado", "Connecticut", "Delaware",
               "Georgia", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas",
               "Kentucky", "Louisiana", "Maine", "Minnesota", "Mississippi",
               "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire",
               "New Mexico", "North Carolina", "North Dakota", "Ohio", "Oklahoma",
               "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota",
               "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "West Virginia",
               "Wisconsin", "Wyoming", "California")

STATE = data.raw[1,]
X.attr = t(data.raw[2:8,])
Y = t(data.raw[9:39,])
colnames(Y) = 1969 + 1:31 
rownames(Y) = STATE.NAME
calif = which(STATE.NAME == "California")
Y = Y[c(setdiff(1:nrow(Y), calif), calif), ]
T0 = 19
N0 = nrow(Y)-1

X.noise = rnorm(50*length(Y)); dim(X.noise)=c(dim(Y), 50)
beta.noise = rnorm(50)
estimates = list(synthdid_estimate(Y, N0, T0),
                 synthdid_estimate(Y+contract3(X.noise[,,1:2], 5*beta.noise[1:2]), N0, T0, X.noise[,,1]))
names(estimates) = c('no covariates', 'ten noise covariates')

synthdid_rmse_plot(estimates)
synthdid_plot(estimates[[3]], sc.plot=TRUE)
source('../../R/synthdid.R')

indices=which(rownames(Y) %in% c('Alabama', 'Arkansas', 'Mississippi', 'Louisiana', 'Georgia', 'Kentucky', 'Tennesee','Florida'))
estimates = list(synthdid_estimate(Y[c(indices, nrow(Y)),], length(indices), T0), synthdid_estimate(Y, N0, T0), synthdid_estimate(Y, N0, T0, omega.intercept=FALSE))
attr(estimates[[3]], 'weights')$lambda = rep(1/T0,T0)
names(estimates) = c('southeast control states', 'all control states', 'sc')
synthdid_plot(estimates, facet.vertical=FALSE)
synthdid_plot(estimates[[1]], treated.name='california', control.name='synth. california') + theme(legend.position=c(.9,.9)


Y.offset = Y[,6:ncol(Y)]
Y.offset[1:N0,] = Y[1:N0, 1:(ncol(Y)-5)]
estimates = list(synthdid_estimate(Y.offset, N0, T0-5), synthdid_estimate(Y, N0, T0))
synthdid_plot(estimates)

Y.offset = Y[,6:ncol(Y)]
Y.offset[1:N0,] = Y[1:N0, 1:(ncol(Y)-5)]
estimates = synthdid_estimate(Y.offset, N0, T0-5)
plot(estimates)

summary(estimates[[1]])
synthdid_placebo_plot(estimates[[1]])
synthdid_placebo_plot(estimates[[1]])
synthdid_placebo_plot(estimates[[1]], overlay=TRUE)

synthdid = estimates[[1]]
did = synthdid_estimate(Y, N0, T0, weights=list(lambda=rep(1/T0,T0), omega=rep(1/N0,N0)))
synthdid_plot(list(synthdid=synthdid, did=did)) 
synthdid_plot(list(synthdid=synthdid, did=did), lambda.comparable=TRUE) + facet_grid(. ~ facet)
