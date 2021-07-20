rm(list = ls())
library(foreign)

setwd('/Users/arkhangelsky/Dropbox/Research/double_weighted_synthetic_control/simulation_dmitry/data_collection')
source('functions.r')


# Organizing the url paths

url_place <- "https://data.nber.org/morg/annual/"


seq_1 <- 79:99
seq_2 <- 0:9
seq_3 <- 10:18

part_1 <- paste(paste('morg',seq_1,sep = ''),'dta',sep = '.')
part_2 <- paste(paste('morg0',seq_2,sep = ''),'dta',sep = '.')
part_3 <- paste(paste('morg',seq_3,sep = ''),'dta',sep = '.')

# Downloading the files 

fi_1 <- paste(url_place,part_1,sep = '')
dat_1 <- lapply(fi_1,read.dta)

fi_2 <- paste(url_place,part_2,sep = '')
dat_2 <- lapply(fi_2,read.dta)

fi_3 <- paste(url_place,part_3,sep = '')
dat_3 <- lapply(fi_3,read.dta)

data_assign <- as.matrix(read.table('state_laws.tsv')) ==1


# Subseting the data


index_full_1 <- c('hhid','hhnum','intmonth','age','sex','minsamp',
				'state','year',"lineno",'earnwke','uhourse',
				"ftpt94","ftpt89","ftpt79")
index_full_2 <- c('hhid','hhnum','intmonth','age','sex','minsamp',
				'stfips','year',"lineno",'earnwke','uhourse',
				"ftpt94","ftpt89","ftpt79")
index_part <- c("ftpt94","ftpt89","ftpt79")

sort_index_1 <- c("hhid", "hhnum", "lineno", "year", "minsamp", "intmonth", "state", "age")
sort_index_2 <- c("hhid", "hhnum", "lineno", "year", "minsamp", "intmonth", 'stfips', "age")


dat_1_sub <- data_const(dat_1,'state',index_full_1,sort_index_1,index_part)
dat_2_sub <- data_const(dat_2,'stfips',index_full_2,sort_index_2,index_part)
dat_3_sub <- data_const(dat_3,'stfips',index_full_2,sort_index_2,index_part)


# Average log wage

lwage_1 <- earn_const(dat_1_sub,'state')
lwage_2 <- earn_const(dat_2_sub,'stfips')
lwage_3 <- earn_const(dat_3_sub,'stfips')

data_mat_wage <- cbind(lwage_1,lwage_2,lwage_3)

# Average hours

hours_1 <- hours_const(dat_1_sub,'state')
hours_2 <- hours_const(dat_2_sub,'stfips')
hours_3 <- hours_const(dat_3_sub,'stfips')

# Unemployment rate

data_mat_hours <- cbind(hours_1,hours_2,hours_3)


urate_1 <- urate_const(dat_1_sub,'state')
urate_2 <- urate_const(dat_2_sub,'stfips')
urate_3 <- urate_const(dat_3_sub,'stfips')

data_mat_urate <- cbind(urate_1,urate_2,urate_3)

# Combining the data into a single file uisng assignment data

data_time_cps <- 1979:2018

#rownames(data_mat_wage) <- rownames(data_assign)

cps_matrix <- as.data.frame(
cbind(rep(rownames(data_assign), ncol(data_mat_wage)),rep(1979:2018, each = nrow(data_mat_wage)),
as.vector(data_mat_wage),as.vector(data_mat_hours),as.vector(data_mat_urate),data_assign[,1],data_assign[,2],data_assign[,3]))

names(cps_matrix) <- c('state','year','log_wage','hours','urate','min_wage','open_carry','abort_ban')
cps_matrix[1:1950,6:8] <- FALSE


write.table(cps_matrix , file = 'cps_matrix.csv',row.names = FALSE, sep=";", dec='.',quote = FALSE)





