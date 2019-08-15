library(mvtnorm)
library(synthdid)
source('../../R/sdid_lib.R')

N_0 <- 100
N_1 <- 20 
T_0 <- 120
T_1 <- 5
N <- N_0 + N_1
T <- T_0 + T_1
tau <- 1
sigma <- 2
rank <- 2
rho <- 0.7
var <- outer(1:T, 1:T, FUN=function(x, y) rho^(abs(x-y)))
W <- (1:N > N_0) %*% t(1:T > T_0)

mu.reps <- 100
noise.reps <- 100
include.slow=FALSE

results_full <- matrix(0, ncol = 4, nrow = mu.reps*noise.reps)
   
for (l in 1:mu.reps){
    U <- matrix(rpois(rank * N, sqrt(1:N) / sqrt(N)), N, rank)
    V <- matrix(rpois(rank * T, sqrt(1:T) / sqrt(T)), T, rank)
    mu <- U %*% t(V)
    
    for (e in 1:noise.reps){
        print(c(l,e))
        error <- rmvnorm(N, sigma = var, method = "chol")
        Y <- mu + tau * W  + sigma * error 

        tau.sdid = synthdid_estimate(Y,N_0,T_0)
        tau.sdid.2 = synthdid_estimate(Y,N_0,T_0, zeta.omega=var(as.numeric(Y))*100)
        t.sdid = (tau.sdid - tau)/attr(tau.sdid,'se')
        t.sdid.2 = (tau.sdid.2 - tau)/attr(tau.sdid.2,'se')
        tau.did = did_estimate(Y,N_0,T_0)
        t.did = (tau.did - tau)/attr(tau.did,'se')
        
        if(include.slow) {
            tau.sdid.slow = synthdid_estimate(Y,N_0,T_0,fast.var=FALSE)
            t.sdid.slow = (tau.sdid.slow - tau)/attr(tau.sdid.slow,'se')
        } else { 
            t.sdid.slow = NA 
        }

        index <- (l-1)*noise.reps + e
        results_full[index,] <- c(t.sdid.slow, t.sdid, t.sdid.2, t.did)
    }
}

    
print(paste('coverage', 'synthdid',                        round(mean(abs(results_full[,1]) < 1.96),4),
                        'synthdid-fastvar',                round(mean(abs(results_full[,2]) < 1.96),4),
                        'synthdid-fastvar-omegapenalty',   round(mean(abs(results_full[,3]) < 1.96),4),
                        'did',                             round(mean(abs(results_full[,4]) < 1.96),4)))
    

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
    
    
  

    
    



        
        
    
