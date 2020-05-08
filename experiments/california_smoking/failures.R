source('setup.R')
source('../../R/synthdid.R')

indices=which(rownames(Y) %in% c('Alabama', 'Arkansas', 'Mississippi', 'Louisiana', 'Georgia', 'Kentucky', 'Tennesee','Florida'))
estimates = list(sc_estimate(Y[c(indices, nrow(Y)),], length(indices), T0), synthdid_estimate(Y, N0, T0), sc_estimate(Y, N0, T0), did_estimate(Y,N0,T0))
names(estimates) = c('southeast control states', 'all control states', 'sc', 'did')
synthdid_plot(estimates, facet.vertical=FALSE)


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
