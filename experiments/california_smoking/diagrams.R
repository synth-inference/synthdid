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
X.sqrt = sqrt(Y)
estimates = list(synthdid_estimate(Y, N0, T0),
                 synthdid_estimate(Y, N0, T0, X.noise[,,1]),
                 synthdid_estimate(Y, N0, T0, X.noise[,,1:10]),
                 synthdid_estimate(Y, N0, T0, X.sqrt))
names(estimates) = c('no covariates', 'one noise covariate', 'ten noise covariates', 'sqrt(Y) as covariate')

synthdid_rmse_plot(estimates)
synthdid_plot(estimates)

plot(estimates[[1]])
summary(estimates[[1]])
synthdid_placebo_plot(estimates[[1]])
synthdid_placebo_plot(estimates[[1]])
synthdid_placebo_plot(estimates[[1]], overlay=TRUE)

synthdid = estimates[[1]]
did = synthdid_estimate(Y, N0, T0, weights=list(lambda=rep(1/T0,T0), omega=rep(1/N0,N0)))
synthdid_plot(list(synthdid=synthdid, did=did)) 
synthdid_plot(list(synthdid=synthdid, did=did), lambda.comparable=TRUE) + facet_grid(. ~ facet)
