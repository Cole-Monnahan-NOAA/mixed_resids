library(DHARMa)
library(magrittr)
library(ggplot2)
library(tidyr)

#demo of unconditional theoretical residuals skewed towards zero

# Run Simulation ===============================
#true parms
n <- 100
ng <- 5
sd.vec <- sqrt(c(2,2))
sd.vec <- sqrt(c(1,1))
eta <- 1

n.sim <- 200
pvals <- data.frame(cond.rot = rep(NA, n.sim), cond.nrot = rep(NA, n.sim),
                    uncond.rot = rep(NA, n.sim), uncond.nrot = rep(NA, n.sim)
                    )
for(ii in 1:n.sim){
  seed <- ii

  set.seed(seed)
  u <- rnorm(ng, 0, sd=sd.vec[2])
  y <- mu <- matrix(0, n, ng)
  
  for(j in 1:ng){
    for(i in 1:n){
      mu[i,j] <- exp(eta + u[j]) 
      #y[i,j] <- rgamma(1, 1/sd.vec[1]^2, scale = mu[i,j]*sd.vec[1]^2) #CV parameterization
      y[i,j] <- rpois(1, mu[i,j])
    }
  }
  
  y <- as.vector(y)
  
  #simulate data
  u.sim <- y.sim.cond <- y.sim.uncond <- matrix(0, n*ng, 1000)
  for(s in 1:1000){
    u.sim <- rnorm(ng, 0, sd=sd.vec[2])
    y.cond <- y.uncond <- mu.sim <- matrix(0, n, ng)
    for(j in 1:ng){
      for(i in 1:n){
        mu.sim[i,j] <- exp(eta + u.sim[j]) 
        #simulate data based on true u vector 
        #y.cond[i,j] <- rgamma(1, 1/sd.vec[1]^2, scale = mu[i,j]*sd.vec[1]^2) #CV parameterization
        y.cond[i,j] <- rpois(1, mu[i,j])
        #simulate data based on simulated u vector 
        #y.uncond[i,j] <- rgamma(1, 1/sd.vec[1]^2, scale = mu.sim[i,j]*sd.vec[1]^2) #CV parameterization
        y.uncond[i,j] <- rpois(1, mu.sim[i,j])
      }
    }
    y.sim.cond[,s] <- as.vector(y.cond)
    y.sim.uncond[,s] <- as.vector(y.uncond)
  }
  
  res.cond.rot <- createDHARMa(y.sim.cond, y, rotation = "estimated")
  res.uncond.rot <- createDHARMa(y.sim.uncond, y, rotation = "estimated")
  res.cond.nrot <- createDHARMa(y.sim.cond, y)
  res.uncond.nrot <- createDHARMa(y.sim.uncond, y)
  pvals$cond.rot[ii] <- ks.test(res.cond.rot$scaledResiduals, 'punif')$p.value
  pvals$cond.nrot[ii] <- ks.test(res.cond.nrot$scaledResiduals, 'punif')$p.value
  pvals$uncond.rot[ii] <- ks.test(res.uncond.rot$scaledResiduals, 'punif')$p.value
  pvals$uncond.nrot[ii] <- ks.test(res.uncond.nrot$scaledResiduals, 'punif')$p.value
  print(ii)
}
save(pvals, file = "simpleGLMM_uncond_test_Poisson.RData")
# theoretical conditional residuals are unifrom while unconditional residuals
# are skewed towards zero

pvals %>% 
  pivot_longer(., col = 1:4, names_to = "method", values_to = "pval") %>% 
  ggplot(., aes(x=pval)) + geom_histogram() + facet_wrap(~method)
0head(pvals)
