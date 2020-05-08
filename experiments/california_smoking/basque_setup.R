library(Synth)
library(reshape2)
data('basque')
Y=acast(basque[,c('gdpcap', 'regionname', 'year')], regionname ~ year, value.var='gdpcap')
treated = which(rownames(Y) == 'Basque Country (Pais Vasco)')
Y=Y[c(setdiff(1:nrow(Y), treated), treated), ]
T0=which(colnames(Y)==1967)
N0=nrow(Y)-1
Y=Y[,1:(T0+5)]



