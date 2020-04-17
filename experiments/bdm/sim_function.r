sim_function <- function(L_mat,tw_st_full,index_order,sigma_mat,n_0,T_0,B,scale_1,scale_2){

	n <- nrow(L_mat)
	T <- ncol(L_mat)
	n_1 <- n - n_0
	T_1 <- T - T_0
	
	L_mat_or_sc <- scale_1*L_mat[index_order,]
	tw_st_or <- tw_st_full[index_order,]
	
	results <- matrix(0, ncol = 3, nrow = B)

    estimates = lapply(1:B, function(b) {
		Y_b <- L_mat_or_sc + tw_st_or + scale_2*rmvnorm(n,sigma = sigma_mat)
		estimates = list(synthdid_estimate(Y_b,n_0,T_0), sc_estimate(Y_b,n_0,T_0), did_estimate(Y_b,n_0,T_0))
        names(estimates) = c('synth-did', 'sc', 'did')
        estimates 
	})
	for(b in 1:B) {
        results[b,] = sapply(estimates[[b]], c)
    }
	bias_methods <- colMeans(results)
	sd_methods <- sqrt(diag(var(results)))
	
	agg_results <- rbind(sqrt(colMeans(results^2)),bias_methods,sd_methods)
	agg_results_round <- round(agg_results,3)
	colnames(agg_results_round) <- c('sdid','did','sc')
	rownames(agg_results_round) <- c('RMSE','bias','sd')
	
	list(aggregates = agg_results_round, estimates=estimates)

}
