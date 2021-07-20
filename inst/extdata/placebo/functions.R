data_const <- function(data_list,state_id,index_full,sort_index,index_part){
	
	### This function subsets individual year-specific file from NBER

	result <- lapply(data_list, function(dta){
	
	# Focus on females from 25 to 50, 4th month of interview and not from DC

	index_fem <- as.factor(dta[,'sex']) == levels(as.factor(dta[,'sex']))[2]
	index_age <- 25 <= dta[,'age'] & dta[,'age'] <=  50
	index_minsamp <- dta[,'minsamp'] == 'MIS 4'
	index_state <-  dta[, state_id] != 'DC'
	
	## Different files have different names -- need to extract the right ones
	
	index_needed <- intersect(index_full,colnames(dta))
	
	## File with potential duplicates
	
	dta_indexed <- dta[index_fem & index_age & index_minsamp & index_state, index_needed]
	
	# Finding the duplicates using part of the information about the unit
	
	index_unique <- !duplicated(dta_indexed[,sort_index])	
	
	# organizing the data in the right way + creating the variables
	
	sd_temp <- setdiff(index_needed,index_part)
	int_temp <- intersect(index_part,index_needed)
	
	dta_final <- cbind(dta_indexed[index_unique,sd_temp],dta_indexed[index_unique,int_temp])
	colnames(dta_final) <- c(sd_temp,'emp_stat')
	dta_final[,'sex'] <- as.character(dta_final[,'sex'])
	
	return(dta_final)})
	return(result)

}



earn_const <- function(data_list,state_id){

	### This function constructs the average log earnings for units with positive earnings

	earn_res <- do.call(cbind,lapply(data_list, function(dta){
	
	ind_earn <- dta[,'earnwke'] > 0
	
	earn_out <- as.matrix(by(log(dta[ind_earn,'earnwke']),dta[ind_earn,state_id],function(data){
		return(mean(data,na.rm = TRUE))}))
	return(earn_out)}))
	
	earn_res <- earn_res[!is.na(earn_res[,1]),]
	earn_res <- earn_res[order(rownames(earn_res)),]
	return(round(earn_res,4))

}


hours_const <- function(data_list,state_id){

	### This function constructs the average working hours for units with positive earnings

	hours_res <- do.call(cbind,lapply(data_list, function(dta){
	
	ind_earn <- dta[,'earnwke'] > 0
	
	hours_out <- as.matrix(by(dta[ind_earn,'uhourse'],dta[ind_earn,state_id],function(data){
		return(mean(data,na.rm = TRUE))}))
	return(hours_out)}))
	
	hours_res <- hours_res[!is.na(hours_res[,1]),]
	hours_res <- hours_res[order(rownames(hours_res)),]
	return(round(hours_res,3))

}


urate_const <- function(data_list,state_id){

	### This function constructs the unempoyment rate for units in the labor force

	urate_res <- do.call(cbind,lapply(data_list, function(dta){
	
	ind_lf <-  !is.na(dta[,'emp_stat']) &  dta[,'emp_stat'] != 'Not In Labor Force' 
	unemp_status <- (dta[,'emp_stat'] == 'Unemployed FT') + (dta[,'emp_stat'] == 'Unemployed PT') 
	
	unemp_out <- as.matrix(by(unemp_status[ind_lf],dta[ind_lf,state_id],function(data){
		return(mean(data,na.rm = TRUE))}))
	return(unemp_out)}))
	
	urate_res <- urate_res[!is.na(urate_res[,1]),]
	urate_res <-urate_res[order(rownames(urate_res)),]
	return(round(urate_res,3))
}




