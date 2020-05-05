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
states = function(...) { which(STATE.NAME %in% c(...)) }
southeast = list('Alabama', 'Arkansas', 'Mississippi', 'Louisiana', 'Georgia', 'Kentucky', 'Tennesee','Florida')
west = list('Colorado', 'Idaho', 'Montana', 'Nevada', 'New Mexico', 'North Dakota', 'South Dakota', 'Utah', 'Wyoming')
Y.southeast = Y[c(do.call(states, southeast), states('California')), ]
Y.west = Y[c(do.call(states, west), states('California')), ]
Y = Y[c(setdiff(1:nrow(Y), states('California')), states('California')), ]
T0 = 19
N0 = nrow(Y)-1

