library(DHARMa)
library(magrittr)
library(ggplot2)
library(tidyr)
library(mvtnorm)

#demo of unconditional theoretical residuals skewed towards zero when distribution is positively skewed

# Run Simulation =================================================
#true parms
n <- 100
ng <- 5
sd.u <- c(0.2, 1.2)
sd.v <- c(3, 1.3)
C <- exp(-as.matrix(dist(seq(0,10,length.out = n))))
S.1 <- diag(rep(sd.v[1], n)) %*% C %*%  diag(rep(sd.v[1], n))
S.2 <- diag(rep(sd.v[2], n)) %*% C %*%  diag(rep(sd.v[2], n))

#Simulate approx normal ([1]) and skewed ([2]) gamma
alpha. <- c(100, .1) #shape parameter of gamma
# b0 <- c(2, 0.5, 2.4) #regression intercept

#Setup simulate function
sim.gamma <- function(nobs, eta, alpha){
  rgamma(nobs, alpha, scale = eta/alpha)
}

n.sim <- 100
pvals.gamma.1 <- pvals.gamma.2 <- 
  pvals.gamma.3 <- pvals.gamma.4 <-  
  data.frame(uncond.rot = rep(NA, n.sim),
             uncond.nrot = rep(NA, n.sim),
             # uncond.re.rot = rep(NA, n.sim),
             # uncond.adj = rep(NA, n.sim),
             cond.nrot = rep(NA, n.sim) )

for(ii in 1:n.sim){
  set.seed(ii); print(ii)
  u.1 <- rnorm(ng, 0, sd=sd.u[1])
  u.2 <- rnorm(ng, 0, sd=sd.u[2])
  v.1 <- as.numeric(rmvnorm(1, rep(0, ncol(C)), sigma=S.1)) 
  v.2 <- as.numeric(rmvnorm(1, rep(0, ncol(C)), sigma=S.2)) 
  
## Obs data =================================
  # Approximately normal - RE grouping term
  y.gamma.1 <- apply(as.matrix(exp(u.1)), 1, 
                     function(x) sim.gamma( n, x, alpha.[1] ))
  # Positively skewed - RE grouping term
  y.gamma.2 <- apply(as.matrix(exp(u.2)), 1, 
                     function(x) sim.gamma( n, x, alpha.[2] ))
  # Approximately normal - Banded correlation
  y.gamma.3 <- sim.gamma(n, v.1+20, alpha.[1])
  # Positively skewed - Banded correlation
  y.gamma.4 <- sim.gamma(n, exp(v.2), alpha.[2])
  
  
## Unconditional Simulation ================
  sim.gamma.1 <- replicate(1000, {
    sim.u <- rnorm(ng, 0, sd=sd.u[1])
    as.vector( apply(exp(as.matrix(sim.u)), 1, 
                   function(x) rgamma(n, alpha.[1], scale = x/alpha.[1]) 
                   ))
  })
  # sim.gamma.1.adj <- replicate(1000, {
  #   sim.u <- rnorm(ng, 0, sd=sd.u[1])
  #   sim.y <- apply(exp(as.matrix(sim.u)), 1, 
  #                  function(x) rgamma(n, alpha.[1], scale = x/alpha.[1]) )
  #   as.vector(apply(sim.y, 2, function(x) (x-mean(x))/sd(x)))
  # })
  sim.gamma.2 <- replicate(1000, {
    sim.u <- rnorm(ng, 0, sd=sd.u[2])
    as.vector( apply(exp(as.matrix(sim.u)), 1, 
                     function(x) rgamma(n, alpha.[2], scale = x/alpha.[2]) 
    ))
  })
  # sim.gamma.2.adj <- replicate(1000, {
  #   sim.u <- rnorm(ng, 0, sd=sd.u[2])
  #   sim.y <- apply(exp(as.matrix(sim.u)), 1, 
  #                  function(x) rgamma(n, alpha.[2], scale = x/alpha.[2]) )
  #   as.vector(apply(sim.y, 2, function(x) (x-mean(x))/sd(x)))
  # })
  # sim.norm.1 <- replicate(1000, {
  #   rep(rnorm(ng, 0, sd=sd.u), each = n)
  # })
  # sim.gamma.1 <- replicate(1000, {
  #   sim.u <- rep(rnorm(ng, 0, sd=sd.u), each = n)
  #   sim.gamma(n*ng, exp(sim.u), alpha.[1])
  # })
  # sim.gamma.1 <- apply(exp(sim.norm.1), 2, 
  #                      function(x) sim.gamma(n*ng, x, alpha.[1]))
  # sim.gamma.1 <- replicate(1000, {
  #   sim.u <- rnorm(ng, 0, sd=sd.u)
  #   sim.y <- apply(exp(as.matrix(sim.u)), 1, function(x) rgamma(n, alpha.[1], scale = x/alpha.[1]) )
  #   as.vector(apply(sim.y, 2, function(x) (x-mean(x))/sd(x)))
  #   })
  # sim.norm.2 <- replicate(1000, {
  #   rep(rnorm(ng, 0, sd=sd.u), each = n)
  # })
  # sim.gamma.2 <- apply(exp(sim.norm.2), 2, 
  #                      function(x) sim.gamma(n*ng, x, alpha.[2]))
  # 
  sim.norm.3 <- replicate(1000, {
    as.numeric(rmvnorm(1, rep(0, ncol(C)), sigma=S.1))# + b0[3]
  })
  sim.gamma.3 <- apply(sim.norm.3+20, 2, 
                       function(x) sim.gamma(n, x, alpha.[1]))
  sim.norm.4 <- replicate(1000, {
    as.numeric(rmvnorm(1, rep(0, ncol(C)), sigma=S.2))# + b0[3]
  })
  sim.gamma.4 <- apply(exp(sim.norm.4), 2, 
                       function(x) sim.gamma(n, x, alpha.[2]))
  
  
  # Calculate p-values ===================
  ## Rotated - estimated using covariance of simulated data
  pvals.gamma.1$uncond.rot[ii] <- 
    ks.test(
      createDHARMa(sim.gamma.1, as.vector(y.gamma.1),
                   rotation = "estimated")$scaledResiduals, 
      'punif')$p.value
  pvals.gamma.2$uncond.rot[ii] <-
    ks.test(
      createDHARMa(sim.gamma.2, as.vector(y.gamma.2),
                   rotation = "estimated")$scaledResiduals,
      'punif')$p.value
  
  pvals.gamma.3$uncond.rot[ii] <-
    ks.test(
      createDHARMa(sim.gamma.3, y.gamma.3,
                   rotation = "estimated")$scaledResiduals,
      'punif')$p.value
  pvals.gamma.4$uncond.rot[ii] <-
    ks.test(
      createDHARMa(sim.gamma.4, y.gamma.4,
                   rotation = "estimated")$scaledResiduals,
      'punif')$p.value
  
  ## Rotated - estimated using covariance of simulated RE
  # covar.1 <- Matrix::nearPD(cov(t(sim.norm.1)))$mat
  # pvals.gamma.1$uncond.re.rot[ii] <-
  #   ks.test(
  #     createDHARMa(sim.gamma.1, y.gamma.1,
  #                  rotation = as.matrix(covar.1))$scaledResiduals,
  #     'punif')$p.value
  # covar.2 <- Matrix::nearPD(cov(t(sim.norm.2)))$mat
  # pvals.gamma.2$uncond.re.rot[ii] <-
  #   ks.test(
  #     createDHARMa(sim.gamma.2, y.gamma.2,
  #                  rotation = as.matrix(covar.2))$scaledResiduals,
  #     'punif')$p.value
  # covar.3 <- Matrix::nearPD(cov(t(sim.norm.3)))$mat
  # pvals.gamma.3$uncond.re.rot[ii] <- 
  #   ks.test(
  #     createDHARMa(sim.gamma.3, y.gamma.3,
  #                  rotation = as.matrix(covar.3))$scaledResiduals,
  #     'punif')$p.value
  # 
  # ## Adjusted, Not rotated
  # pvals.gamma.1$uncond.adj[ii] <- 
  #   ks.test(
  #     createDHARMa(sim.gamma.1.adj,
  #                  apply(y.gamma.1, 2, 
  #                        function(x) (x-mean(x))/sd(x)))$scaledResiduals, 
  #     'punif')$p.value
  # pvals.gamma.2$uncond.adj[ii] <- 
  #   ks.test(
  #     createDHARMa(sim.gamma.2.adj,
  #                  apply(y.gamma.2, 2, 
  #                        function(x) (x-mean(x))/sd(x)))$scaledResiduals, 
  #     'punif')$p.value
  
  
  ## Not Rotated
  pvals.gamma.1$uncond.nrot[ii] <- 
    ks.test(
      createDHARMa(sim.gamma.1, y.gamma.1)$scaledResiduals, 
      'punif')$p.value
  pvals.gamma.2$uncond.nrot[ii] <-
    ks.test(
      createDHARMa(sim.gamma.2, y.gamma.2)$scaledResiduals,
      'punif')$p.value
  pvals.gamma.3$uncond.nrot[ii] <-
    ks.test(
      createDHARMa(sim.gamma.3, y.gamma.3)$scaledResiduals,
      'punif')$p.value
  pvals.gamma.4$uncond.nrot[ii] <-
    ks.test(
      createDHARMa(sim.gamma.4, y.gamma.4)$scaledResiduals,
      'punif')$p.value
  
  
# Conditional Simulation ==============================
  sim.gamma.1 <- replicate(1000, {
    sim.gamma( n*ng, exp( rep(u.1, each = n)), alpha.[1] )
  })
  sim.gamma.2 <- replicate(1000, {
    sim.gamma( n*ng, exp( rep(u.2, each = n)), alpha.[2] )
  })
  sim.gamma.3 <- replicate(1000,{
    sim.gamma(n, v.1+20, alpha.[1])
  })
  sim.gamma.4 <- replicate(1000,{
    sim.gamma(n, exp(v.2), alpha.[2])
  })

  # Calculate p-values ===================
  
  pvals.gamma.1$cond.nrot[ii] <- 
    ks.test(
      createDHARMa(sim.gamma.1, y.gamma.1)$scaledResiduals, 
      'punif')$p.value
  pvals.gamma.2$cond.nrot[ii] <-
    ks.test(
      createDHARMa(sim.gamma.2, y.gamma.2)$scaledResiduals,
      'punif')$p.value
  pvals.gamma.3$cond.nrot[ii] <-
    ks.test(
      createDHARMa(sim.gamma.3, y.gamma.3)$scaledResiduals,
      'punif')$p.value
  pvals.gamma.4$cond.nrot[ii] <-
    ks.test(
      createDHARMa(sim.gamma.4, y.gamma.4)$scaledResiduals,
      'punif')$p.value
}

save(pvals.gamma.1, file = "code/unconditional_test/pvals_gamma_1.RData")
save(pvals.gamma.2, file = "code/unconditional_test/pvals_gamma_2.RData")
save(pvals.gamma.3, file = "code/unconditional_test/pvals_gamma_3.RData")
save(pvals.gamma.4, file = "code/unconditional_test/pvals_gamma_4.RData")

# Plot results ==========================================
load("code/unconditional_test/pvals_gamma_1.RData")
load("code/unconditional_test/pvals_gamma_2.RData")
load("code/unconditional_test/pvals_gamma_3.RData")
load("code/unconditional_test/pvals_gamma_4.RData")

pvals.gamma.1 <- pivot_longer(pvals.gamma.1, 1:3, names_to = "method",
                              values_to = "pvalue")
pvals.gamma.2 <- pivot_longer(pvals.gamma.2, 1:3, names_to = "method",
                              values_to = "pvalue")
pvals.gamma.3 <- pivot_longer(pvals.gamma.3, 1:3, names_to = "method",
                              values_to = "pvalue")
pvals.gamma.4 <- pivot_longer(pvals.gamma.4, 1:3, names_to = "method",
                              values_to = "pvalue")

pvals.gamma.1$case <- "Approx. Normal Gamma w/ Grouping RE"
pvals.gamma.2$case <- "Skewed Gamma w/ Grouping RE"
pvals.gamma.3$case <- "Approx. Normal Gamma w/ Banded Correlation"
pvals.gamma.4$case <- "Skewed Gamma w/ Banded Correlation"

pval.results <- rbind(
  pvals.gamma.1, pvals.gamma.2, pvals.gamma.3, pvals.gamma.4
)

method.labs <- c(
  uncond.rot = "Unconditional, Rotated", 
  uncond.nrot = "Unconditional, Not Rotated",
  # uncond.re.rot = "Unconditional, RE Rotated",
  # uncond.adj = "Unconditional, Adjusted",
  cond.nrot = "Conditional, Not Rotated")

pval.results %>%
  ggplot(., aes(x=pvalue)) + geom_histogram(bins = 50) + 
  facet_grid(case~method, 
             labeller = labeller(method = method.labs, 
                                 .default = label_wrap_gen(17))) +
  scale_x_continuous(
    breaks=c(0,0.5, 1),
    labels = scales::number_format(accuracy = 0.1)) +
  theme_bw()
