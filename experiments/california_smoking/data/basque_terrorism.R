library(Synth)
library(reshape2)
data('basque')
Y=acast(basque[,c('gdpcap', 'regionname', 'year')], regionname ~ year, value.var='gdpcap')
treated = which(rownames(Y) == 'Basque Country (Pais Vasco)')
Y=Y[c(setdiff(1:nrow(Y), treated), treated), ]
T0=which(colnames(Y)==1967)
N0=nrow(Y)-1
Y=Y[,1:(T0+5)]


library(reshape2)
mY = melt(Y)
colnames(mY) = c('Region', 'Year', 'GDPPerCapita')
mY$treated = mY$Region == 'Basque Country (Pais Vasco)' & mY$Year >= 1967
write.csv(mY, file='basque_terrorism.csv')


