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

calif = which(STATE.NAME == "California")

# only consider 5 post-treatment periods
post.periods = 1
focal.time = 11:ncol(Y.pre)

Delta.pre = Y.pre[,-1] - Y.pre[,-19]

sigma.sq = var(as.numeric(Delta.pre))

all.preds = outer(1:length(STATE.NAME), focal.time, Vectorize(function(ss, tt) {
    omega.weight = sc_weight(t(Y.pre[-ss,1:(tt-1)]), Y.pre[ss,1:(tt-1)], zeta = sigma.sq)
    lambda.weight = sc_weight_FE(Y.pre[-ss,1:(tt-1)], Y.pre[-ss,tt], zeta = sigma.sq)
    col.SC = sum(lambda.weight * Y.pre[ss,1:(tt-1)])
    row.SC = sum(omega.weight * Y.pre[-ss,tt])
    cross.term = omega.weight %*% Y.pre[-ss,1:(tt-1)] %*% lambda.weight
    sdid = as.numeric(row.SC + col.SC - cross.term)
    sc = row.SC
    did = mean(Y.pre[-ss,tt]) + mean(Y.pre[ss,1:(tt-1)]) - mean(Y.pre[-ss,1:(tt-1)])
    list(c(did=did, sc=sc, sdid=sdid, truth=Y.pre[ss,tt]))
}))

did.err = apply(all.preds, 1:2, function(vv) vv[[1]]["did"] - vv[[1]]["truth"])
sc.err = apply(all.preds, 1:2, function(vv) vv[[1]]["sc"] - vv[[1]]["truth"])
sdid.err = apply(all.preds, 1:2, function(vv) vv[[1]]["sdid"] - vv[[1]]["truth"])

c(mean(did.err^2), mean(sc.err^2), mean(sdid.err^2))

state.rmse = data.frame(DID=sqrt(rowMeans(did.err^2)),
                        SC=sqrt(rowMeans(sc.err^2)),
                        SDID=sqrt(rowMeans(sdid.err^2)))

pdf("calif_sdid_vs_sc_mse.pdf")
pardef = par(xpd = FALSE, mar = c(4.5, 5, 3, 3) + 0.5, cex.lab=1.4, cex.axis=1.4, cex.main=1.4, cex.sub=1.4)
plot(state.rmse$SDID, state.rmse$SC, log = "xy", cex = 1.5,
     xlab = "synthetic diff-in-diff RMSE",
     ylab = "synthetic control RMSE")
abline(0, 1, lwd = 2, lty = 3)
points(state.rmse$SDID[calif], state.rmse$SC[calif], pch = 16, col = 4, cex = 1.5)
par(pardef)
dev.off()

for(iter in 1:length(STATE.NAME)) {
    fnm = paste0("state_plots/", STATE.NAME[iter], ".pdf")
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
}

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


