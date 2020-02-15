library(synthdid)
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

X.noise = rnorm(5*length(Y)); dim(X.noise)=c(dim(Y), 5)
X.sqrt = sqrt(Y)
estimates = list(synthdid_estimate(Y, N0, T0),
                 synthdid_estimate(Y, N0, T0, X.noise[,,1]),
                 synthdid_estimate(Y, N0, T0, X.noise,),
                 synthdid_estimate(Y, N0, T0, X.sqrt))
names(estimates) = c('no covariates', 'one noise covariate', 'five noise covariates', 'sqrt(Y) as covariate')

synthdid_rmse_plot(estimates)
synthdid_plot(estimates)
print(synthdid_controls(estimates[1:3], mass=.9))

