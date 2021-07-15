### Saves the CPS and PENN data from Dropbox so that it is accessible via data(CPS) and data(PENN)
### and checks that when it is loaded and processed into matrices Y and W  as is done in the package, (via data(...) and panel.matrices(...)) 
### these matrices agree with the ones used in previous simulations.
data.dir = path.expand('~/Dropbox/double_weighted_synthetic_control/simulation_dmitry/simulation_paper/data')

## Save CPS data in data(...)-compatible format, with treatment (nominally) beginning in the last time-step
cps_data <- read.csv(sprintf('%s/../../data_for_package/cps_matrix .csv', data.dir), header = TRUE)
CPS = cps_data[,-1]
for(treatment in c('min_wage', 'open_carry', 'abort_ban')) {
    CPS[, treatment] = CPS[, treatment] & CPS$year == max(CPS$year)
}
write.table(CPS, file='../../data/CPS.csv', sep=';', row.names=FALSE, col.names=TRUE)

## Save PENN data in the same way
penn_data <- read.csv(sprintf('%s/../../data_for_package/penn_data.csv', data.dir), header = TRUE)
PENN = penn_data[penn_data$year < 2008, -1]
for(treatment in c('dem', 'educ')) {
    PENN[, treatment] = PENN[, treatment] & PENN$year == max(PENN$year)
}
write.table(PENN, file='../../data/PENN.csv', sep=';', row.names=FALSE, col.names=TRUE)

## Reload package to use this new data
devtools::install(); library(synthdid) 

## Check CPS Data as used in Dropbox sims against package version
library(foreign) 

data_mat_lwage <- as.matrix(reshape(cps_data[,2:4],direction = 'wide',timevar = 'year',idvar = 'state')[,-1])
data_mat_hours <- as.matrix(reshape(cps_data[,c(2:3,5)],direction = 'wide',timevar = 'year',idvar = 'state')[,-1])
data_mat_urate <- as.matrix(reshape(cps_data[,c(2:3,6)],direction = 'wide',timevar = 'year',idvar = 'state')[,-1])
data_assign <- cps_data[1:50,7:9]

data(CPS)
Y.logwage  = panel.matrices(CPS, treatment='min_wage', outcome='log_wage', treated.last=FALSE)$Y 
stopifnot(all(Y.logwage == data_mat_lwage))

## Check Penn Data similarly

data_mat_init <- as.matrix(reshape(penn_data[,2:4],direction = 'wide',timevar = 'year',idvar = 'country')[,-1])
T_full <- dim(data_mat_init)[2]
n <- dim(data_mat_init)[1]
data_mat <- data_mat_init[,-((T_full-9):T_full)]
assign_dem <- penn_data[1:n,5]
assign_educ <- penn_data[1:n,6]

data(PENN)
last.col = function(x) { x[,ncol(x)] }
Y.loggdp = panel.matrices(PENN, treatment='dem', outcome='log_gdp', treated.last=FALSE)$Y
w.dem    = last.col(panel.matrices(PENN, treatment='dem', outcome='log_gdp', treated.last=FALSE)$W)

# units are not in the same order and we've lost names in data_mat, so compare after sorting by first row
order1 = order(Y.loggdp[,1])
order2 = order(data_mat[,1])
stopifnot(all(Y.loggdp[order1,] == data_mat[order2,]))
stopifnot(all(w.dem[order1] == assign_dem[order2]))
