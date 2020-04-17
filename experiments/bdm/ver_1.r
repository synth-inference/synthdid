library(mvtnorm)
library(devtools)
library(xtable)
set.seed(12345)
source('../../R/synthdid.R')
source('sim_function.r')
data_lwage <- t(as.matrix(read.csv('data/lwage_cps.csv',header = FALSE)))
# Parameters

rank <- 2

# Outcome model

data_mat <- data_lwage
n <- dim(data_mat)[1]
T <- dim(data_mat)[2]

tw_st_0 <- outer(rowMeans(data_mat),rep(1,T)) +
		outer(rep(1,n),colMeans(data_mat)) - mean(data_mat)


data_mat_dm <- data_mat - tw_st_0
svd_data_mat <- svd(data_mat_dm)
factor_unit <- as.matrix(svd_data_mat$u[,1:rank]*sqrt(n))
factor_time <- as.matrix(svd_data_mat$v[,1:rank]*sqrt(T))
magnitude <- svd_data_mat$d[1:rank]/sqrt(n*T)
L_mat <- factor_unit%*%diag(magnitude, nrow = rank, ncol = rank)%*%t(factor_time)


error_mat <- data_mat_dm - L_mat
tw_st_error <- outer(rowMeans(error_mat),rep(1,T)) +
		outer(rep(1,n),colMeans(error_mat)) - mean(error_mat)
error_mat_dm <-  error_mat-tw_st_error
sigma_mat <- t(error_mat_dm)%*%error_mat_dm/n

tw_st_full <- tw_st_error +tw_st_0



noise_level <- sqrt(norm(sigma_mat, type = '2'))
signal_low_rank <- norm(L_mat, type = 'F')/sqrt(n*T)
signal_fe <- norm(tw_st_full, type = 'F')/sqrt(n*T)
lr_signal_to_noise <- signal_low_rank/noise_level

# Assignment 

W_i <- factor_unit[,1]>=quantile(factor_unit[,1],0.7)
n_0 <- n - sum(W_i)
T_0 <- floor(0.7*T)
index_order <- order(W_i)


# Simulation

B <- 100
show.results = function(results) {
    print(results$aggregates)
    synthdid_plot(results$estimates[[1]], facet.vertical=FALSE) 
}

# 1. Original design
scale_1 <- 1
scale_2 <- 1
sim_results_orig <- sim_function(L_mat,tw_st_full,index_order,sigma_mat,n_0,T_0,B,scale_1,scale_2)
show.results(sim_results_orig)

# 2. Scaled design

scale_1 <- signal_fe
scale_2 <- lr_signal_to_noise*signal_fe
sim_results_sc <- sim_function(L_mat,tw_st_full,index_order,sigma_mat,n_0,T_0,B,scale_1,scale_2)
show.results(sim_results_sc) 
show.results(sim_results_sc) + coord_cartesian(ylim=c(6,6.5))

# 3. Shuffled design

scale_1 <- 1
scale_2 <- 1
sh_index <- index_order[sample(1:n,n)]
sim_results_sh <- sim_function(L_mat,tw_st_full,sh_index,sigma_mat,n_0,T_0,B,scale_1,scale_2)
show.results(sim_results_sh)

# 4. Mixture of the two

scale_1 <- 1
scale_2 <- 1
W_base <- factor_unit[,1]>=quantile(factor_unit[,1],0.5)
W_i_mixed <- rbinom(length(W_base),1,as.numeric(W_base)/2) == 1
index_order_mixed <- order(W_i_mixed)
n_0_new <- n - sum(W_i_mixed)
sim_results_mixed <- sim_function(L_mat,tw_st_full,index_order_mixed,sigma_mat,n_0_new,T_0,B,scale_1,scale_2)
show.results(sim_results_mixed)

# 5. No low rank

scale_1 <- 0
scale_2 <- 1
sim_results_tw <- sim_function(L_mat,tw_st_full,index_order_mixed,sigma_mat,n_0_new,T_0,B,scale_1,scale_2)
show.results(sim_results_tw)

# Tex files 


xtable(sim_results_orig$aggregates, type = "latex")
xtable(sim_results_sc$aggregates, type = "latex")
xtable(sim_results_sh$aggregates, type = "latex")
xtable(sim_results_mixed$aggregates, type = "latex")
xtable(sim_results_tw$aggregates, type = "latex")



