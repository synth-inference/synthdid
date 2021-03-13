## load and set up dataset
data.raw = read.table("raw/california_prop99.txt")

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
states = function(...) { which(STATE.NAME %in% c(...)) }
Y = Y[c(setdiff(1:nrow(Y), states('California')), states('California')), ]
T0 = 19
N0 = nrow(Y)-1

## write datset as csv
library(reshape2)
mY = melt(Y)
colnames(mY) = c('State', 'Year', 'PacksPerCapita')
mY$treated = as.numeric(mY$State == 'California' & mY$Year >= 1989)
write.csv(mY, file='california_prop99.csv', row.names = FALSE)

## read and compare
devtools::load_all('.')
setup = panel.matrices(read.csv('california_prop99.csv'))
if(!  all(setup$Y == Y)
   &  all(rownames(setup$Y) ==  rownames(Y))
   &  all(colnames(setup$Y) ==  colnames(Y))
   &  setup$N0 == N0
   &  setup$T0 == T0) { error('california prop 99 data matrix does not match after writing and reading csv calling make.panel') }
