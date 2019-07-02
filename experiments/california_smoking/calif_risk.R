rm(list = ls())
library(ggplot2)
library(lfe)
library(lmtest)
library(sandwich)
library(xtable)
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

Y.avg = apply(data.raw[9:39, -calif], 1, mean)
Y.calif = data.raw[9:39, calif]
Y.avg.sec = seq(Y.avg[1], Y.avg[31], length.out=31)
Y.calif.sec = Y.avg.sec + Y.calif[pre] - Y.avg.sec[pre]
plot.data = data.frame(y=c(Y.avg, Y.avg.sec, Y.calif, Y.calif.sec), time=rep(1960 + (9:39), 4), state=rep(c('Average', 'California'), each=31*2),
                       type = rep(rep(c('Observed', 'Secant'), each=31),2))
ggplot(plot.data) + geom_line(aes(x=time, y=y, color=state, linetype=type), size=1)

Y.avg = apply(data.raw[9:39, -calif], 1, mean)
Y.calif = data.raw[9:39, calif]
plot.data = data.frame(y=c(Y.avg, Y.avg.sec, Y.calif, Y.calif.sec), time=rep(1960 + (9:39), 4), state=rep(c('Average', 'California'), each=31*2),
                       type = rep(rep(c('Observed', 'Secant'), each=31),2))
ggplot(plot.data) + geom_line(aes(x=time, y=y, color=state, linetype=type), size=1)




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
           + xlab('') + ylab('') 

list(c(did=did, sc=sc, sdid=sdid, truth=Y.pre[ss,tt]))

pdf(fnm)
pardef = par(xpd = FALSE, mar = c(4.5, 5, 3, 3) + 0.5, cex.lab=1.4, cex.axis=1.4, cex.main=1.4, cex.sub=1.4)
y.truth = sapply(all.preds[iter,], function(vv)vv["truth"])
y.did = sapply(all.preds[iter,], function(vv)vv["did"])
y.sc = sapply(all.preds[iter,], function(vv)vv["sc"])
y.sdid = sapply(all.preds[iter,], function(vv)vv["sdid"])
plot(1969 + focal.time, y.truth, type = "l", lwd = 2,
     ylim = range(y.sc, y.did, y.sdid, y.truth),
     xlab = "year", ylab = "smoking [packs per capita]")
lines(1969 + focal.time, y.did, lwd = 2, col = 5, lty = 2)
lines(1969 + focal.time, y.sc, lwd = 2, col = 4, lty = 5)
lines(1969 + focal.time, y.sdid, lwd = 2, col = 2, lty = 4)
par(pardef)
dev.off()

row.names(state.rmse) = STATE.NAME
state.rmse.ord = state.rmse[order(STATE.NAME),]
xtab.rmse = xtable(state.rmse.ord)
print(xtab.rmse, file="statewise_rmse.tex")

state.vec = STATE.NAME[rep(1:nrow(did.err), ncol(did.err))]
state.vec = factor(state.vec)
levels(state.vec) = rev(levels(state.vec))
time.vec = 1969 + as.numeric(t(matrix(rep(focal.time, nrow(did.err)), ncol(did.err), nrow(did.err))))

pdf("calif_sdid_errors.pdf")
sdid.df = data.frame(error=as.numeric(sdid.err), state=state.vec, time=time.vec)
ggplot(data=sdid.df, aes(x=time, y=state, fill=error)) +
    geom_tile() +
    scale_fill_gradient2(low="blue", high="red", na.value="black", name="") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
dev.off()

pdf("calif_did_errors.pdf")
did.df = data.frame(error=as.numeric(did.err), state=state.vec, time=time.vec)
ggplot(data=did.df, aes(x=time, y=state, fill=error)) +
    geom_tile() +
    scale_fill_gradient2(low="blue", high="red", na.value="black", name="") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
dev.off()

pdf("calif_sc_errors.pdf")
sc.df = data.frame(error=as.numeric(sc.err), state=state.vec, time=time.vec)
ggplot(data=sc.df, aes(x=time, y=state, fill=error)) +
    geom_tile() +
    scale_fill_gradient2(low="blue", high="red", na.value="black", name="") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
dev.off()


