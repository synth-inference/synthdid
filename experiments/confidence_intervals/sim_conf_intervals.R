library(mvtnorm)
library(synthdid)
source('../../R/synthdid.R')

N_0 <- 100
N_1 <- 10 
T_0 <- 120
T_1 <- 1
N <- N_0 + N_1
T <- T_0 + T_1
tau <- 1
sigma <- 2
rank <- 2
rho <- 0.5
var <- outer(1:T, 1:T, FUN=function(x, y) rho^(abs(x-y)))
W <- (1:N > N_0) %*% t(1:T > T_0)

mu.reps <- 100
noise.reps <- 20
results_full <- matrix(0, ncol = 3, nrow = mu.reps*noise.reps)
   
for (l in 1:mu.reps){
    U <- matrix(rpois(rank * N, sqrt(1:N) / sqrt(N)), N, rank)
    V <- matrix(rpois(rank * T, sqrt(1:T) / sqrt(T)), T, rank)
    mu <- U %*% t(V)
    
    for (e in 1:noise.reps){
        print(c(l,e))
        error <- rmvnorm(N, sigma = var, method = "chol")
        Y <- mu + tau * W  + sigma * error 

        sdid = synthdid_estimate(Y,N_0,T_0)
        did  = synthdid_estimate(Y,N_0,T_0, weights=list(omega=rep(1/N_0,N_0), lambda=rep(1/T_0,T_0)))
        t.sdid = (c(sdid) - tau)/synthdid_se(sdid, weights=NULL)
        t.sdid.fast = (c(sdid) - tau)/synthdid_se(sdid)
        t.did = (c(did) - tau)/synthdid_se(did)
        
        index <- (l-1)*noise.reps + e
        results_full[index,] <- c(t.sdid, t.sdid.fast, t.did)
    }
}

    
print(paste('coverage', 'synthdid',                        round(mean(abs(results_full[,1]) < 1.96),4),
                        'synthdid-fastvar',                round(mean(abs(results_full[,2]) < 1.96),4),
                        'did',                             round(mean(abs(results_full[,3]) < 1.96),4)))
    

pdf(sprintf('qq_sdid_cor_%d.pdf', 100*rho), width=5,height=5,paper='special') 
qqnorm(results_full[,1], pch = 1, frame = F, cex = 0.5, main = NULL)
qqline(rnorm(mu.reps*noise.reps), col = "black", lwd = 1, lty = 2)  
box(lwd=0.1)    
dev.off()

pdf(sprintf('qq_sdid_fastvar_cor_%d.pdf', 100*rho), width=5,height=5,paper='special') 
qqnorm(results_full[,2], pch = 1, frame = F,cex = 0.5,main = NULL)
qqline(rnorm(mu.reps*noise.reps), col = "black", lwd = 1, lty = 2)
box(lwd=0.1)        
dev.off()

pdf(sprintf('qq_did_cor_%d.pdf', 100*rho), width=5,height=5,paper='special') 
qqnorm(results_full[,3], pch = 1, frame = F,cex = 0.5,main = NULL)
qqline(rnorm(mu.reps*noise.reps), col = "black", lwd = 1, lty = 2)
box(lwd=0.1)        
dev.off()
    
    
  

    
    



        
        
    
