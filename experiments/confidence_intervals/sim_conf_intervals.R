library(mvtnorm)
library(synthdid)

n_0 <- 100
n_1 <- 10 
T_0 <- 120
T_1 <- 20
n <- n_0 + n_1
T <- T_0 + T_1
tau <- 1
sigma <- 2
rank <- 2
rho <- 0.7
var <- outer(1:T, 1:T, FUN=function(x, y) rho^(abs(x-y)))
W <- (1:n > n_0) %*% t(1:T > T_0)

mu.reps <- 100
noise.reps <- 25
include.slow=FALSE

results_full <- matrix(0, ncol = 3, nrow = mu.reps*noise.reps)
   
for (l in 1:mu.reps){
    print(l)
    U <- matrix(rpois(rank * n, sqrt(1:n) / sqrt(n)), n, rank)
    V <- matrix(rpois(rank * T, sqrt(1:T) / sqrt(T)), T, rank)
    mu <- U %*% t(V)
    
    for (e in 1:noise.reps){
        error <- rmvnorm(n, sigma = var, method = "chol")
        Y <- mu + tau * W  + sigma * error 

        tau.sdid = sdid_est(Y,n_0,T_0)
        t.sdid = (tau.sdid - tau)/attr(tau.sdid,'se')
        tau.did = did_est(Y,n_0,T_0)
        t.did = (tau.did - tau)/attr(tau.did,'se')
        
        if(include.slow) {
            tau.sdid.slow = sdid_est(Y,n_0,T_0, fast.var=FALSE)
            t.sdid.slow = (tau.sdid.slow - tau)/attr(tau.sdid.slow,'se')
        } else { 
            t.sdid.slow = NA 
        }

        index <- (l-1)*noise.reps + e
        results_full[index,] <- c(t.sdid, t.did, t.sdid.slow)
    }
}



pdf(sprintf('qq_sdid_cor_%d.pdf', 100*rho), width=5,height=5,paper='special') 
qqnorm(results_full[,1], pch = 1, frame = F, cex = 0.5, main = NULL)
qqline(rnorm(mu.reps*noise.reps), col = "black", lwd = 1, lty = 2)  
box(lwd=0.1)    
dev.off()

pdf(sprintf('qq_did_cor_%d.pdf', 100*rho), width=5,height=5,paper='special') 
qqnorm(results_full[,2], pch = 1, frame = F,cex = 0.5,main = NULL)
qqline(rnorm(mu.reps*noise.reps), col = "black", lwd = 1, lty = 2)
box(lwd=0.1)        
dev.off()

pdf(sprintf('qq_sdid_slow_cor_%d.pdf', 100*rho), width=5,height=5,paper='special') 
qqnorm(results_full[,3], pch = 1, frame = F,cex = 0.5,main = NULL)
qqline(rnorm(mu.reps*noise.reps), col = "black", lwd = 1, lty = 2)
box(lwd=0.1)        
dev.off()

    
print(paste('coverage', 'sdid',      round(mean(abs(results_full[,1]) < 1.96),4),
                        'sdid-slow', round(mean(abs(results_full[,3]) < 1.96),4),
                        'did',       round(mean(abs(results_full[,2]) < 1.96),4)))
    

    
   
    
    
    
  

    
    



        
        
    
