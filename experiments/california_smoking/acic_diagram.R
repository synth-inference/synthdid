library(ggplot2)
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
X = t(data.raw[2:8,])
Y.pre = t(data.raw[9:27,])
Y.post = t(data.raw[28:39,])
pre = 1:19
post = 20:31

calif = which(STATE.NAME == "California")

# only consider 5 post-treatment periods
ss = calif
omega.weight = sc_weight(t(Y.pre[-ss,]), Y.pre[ss,], zeta = var(as.numeric(Y.pre)))
lambda.weight = sc_weight(Y.pre[-ss,], apply(Y.post[-ss,],1,mean), zeta = var(as.numeric(Y.pre)))

col.SC = sum(lambda.weight * Y.pre[ss,])
row.SC = sum(omega.weight  *  apply(Y.post[-ss,], 1, mean))
cross.term = omega.weight %*% Y.pre[-ss,] %*% lambda.weight
sdid = as.numeric(row.SC + col.SC - cross.term)
sc.trajectory = omega.weight %*% cbind(Y.pre[-ss,], Y.post[-ss,])

pre.time = sum(pre * lambda.weight) + 1969
post.time = mean(post) + 1969
calif.post = mean(Y.post[ss,])
ggplot() + geom_segment(aes(x=pre.time, xend=post.time, y=col.SC, yend=calif.post), color='red') +
           geom_point(aes(x=c(pre.time, post.time), y=c(col.SC, calif.post)), color = 'red') +
           geom_segment(aes(x=pre.time, xend=post.time, y=cross.term, yend=row.SC), color='blue') +
           geom_point(aes(x=c(pre.time, post.time), y=c(cross.term, yend=row.SC)), color='blue') +
           geom_segment(aes(x=pre.time, xend=post.time, y=col.SC, yend=row.SC + col.SC - cross.term), color='red', linetype=2) +
           geom_point(aes(x=c(pre.time, post.time), y=c(col.SC, yend=row.SC + col.SC - cross.term)), color='red', shape=21) +
           geom_line(aes(x=1969+pre, y=Y.pre[ss,]), color='red', linetype=3) +
           geom_line(aes(x=1969+1:31, y = as.numeric(sc.trajectory)), color='blue', linetype=3) + 
           geom_vline(aes(xintercept=1969+19), color='black', linetype=1, alpha=.2) +  
           xlab('') + ylab('') 


