library(DHARMa)
library(magrittr)
library(ggplot2)
library(tidyr)

#demo of unconditional theoretical residuals skewed towards zero when distribution is positively skewed

# Run Simulation =================================================
#true parms
n <- 100
ng <- 5
sd.u <- 1
beta <- c(2.3, 0.4)

#Simulate approx normal ([1]) and skewed ([2]) gamma
y.var <- c(1, 4.5)

n.sim <- 100
pvals.gamma.1 <- pvals.gamma.2 <-
  data.frame(uncond.rot = rep(NA, n.sim),
             uncond.nrot = rep(NA, n.sim),
             uncond.re.rot = rep(NA, n.sim),
             cond.nrot = rep(NA, n.sim) )

for(ii in 1:n.sim){
  set.seed(ii); print(ii)
  u <- rnorm(ng, 0, sd=sd.u)
  
## Obs data =================================
  # Approximately normal
  y.gamma.1 <- as.matrix(u + beta[1]) %>% exp() %>%
    apply(., 1, function(x) rgamma(n, x^2/y.var[1], scale = y.var[1]/x)) %>% 
    as.vector() 
  # Positively skewed
  y.gamma.2 <- as.matrix(u + beta[2]) %>% exp() %>%
    apply(., 1, function(x) rgamma(n, x^2/y.var[2], scale = y.var[2]/x)) %>% 
    as.vector() 
 
## Unconditional Simulation ================
  sim.norm.1 <- replicate(1000, {
    rep(rnorm(ng, beta[1], sd=sd.u), each = n)
  })
  sim.gamma.1 <- sim.norm.1 %>%
      exp() %>%
      apply(., 2, function(x) rgamma(n*ng, x^2/y.var[1], scale = y.var[1]/x))
  
  sim.norm.2 <- replicate(1000, {
    rep(rnorm(ng, beta[2], sd=sd.u), each = n)
  })
  sim.gamma.2 <- sim.norm.2 %>%
    exp() %>%
    apply(., 2, function(x) rgamma(n*ng, x^2/y.var[2], scale = y.var[2]/x))
  
  # Calculate p-values ===================
  ## Rotated - estimated using covariance of simulated data
  pvals.gamma.1$uncond.rot[ii] <- 
    ks.test(
      createDHARMa(sim.gamma.1, y.gamma.1,
                   rotation = "estimated")$scaledResiduals, 
      'punif')$p.value
  pvals.gamma.2$uncond.rot[ii] <-
    ks.test(
      createDHARMa(sim.gamma.2, y.gamma.2,
                   rotation = "estimated")$scaledResiduals,
      'punif')$p.value
  
  ## Rotated - estimated using covariance of simulated RE
  covar.1 <- Matrix::nearPD(cov(t(sim.norm.1)))$mat
  pvals.gamma.1$uncond.re.rot[ii] <- 
    ks.test(
      createDHARMa(sim.gamma.1, y.gamma.1,
                   rotation = as.matrix(covar.1))$scaledResiduals,
      'punif')$p.value
  covar.2 <- Matrix::nearPD(cov(t(sim.norm.2)))$mat
  pvals.gamma.2$uncond.re.rot[ii] <- 
    ks.test(
      createDHARMa(sim.gamma.2, y.gamma.2,
                   rotation = as.matrix(covar.2))$scaledResiduals,
      'punif')$p.value
  
  
  ## Not Rotated
  pvals.gamma.1$uncond.nrot[ii] <- 
    ks.test(
      createDHARMa(sim.gamma.1, y.gamma.1)$scaledResiduals, 
      'punif')$p.value
  pvals.gamma.2$uncond.nrot[ii] <-
    ks.test(
      createDHARMa(sim.gamma.2, y.gamma.2)$scaledResiduals,
      'punif')$p.value
  
  
  
# Conditional Simulation ==============================
  sim.gamma.1 <- replicate(1000, {
    as.matrix(u + beta[1])  %>% exp() %>%
      apply(., 1, function(x) rgamma(n, x^2/y.var[1], scale = y.var[1]/x)) %>%
      as.vector()
  })
  sim.gamma.2 <- replicate(1000, {
    as.matrix(u + beta[2])%>% exp() %>%
      apply(., 1, function(x) rgamma(n, x^2/y.var[2], scale = y.var[2]/x)) %>%
      as.vector()
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
}

# Rotation correctly adjusts unconditional residuals when distribution approx. normal
pvals.gamma.1 %>% 
  pivot_longer(., col = 1:4, names_to = "method", values_to = "pval") %>% 
  ggplot(., aes(x=pval)) + geom_histogram() + facet_wrap(~method)
# Rotation does not correctly adjust unconditional residuals when distribution skewed
pvals.gamma.2 %>% 
  pivot_longer(., col = 1:4, names_to = "method", values_to = "pval") %>% 
  ggplot(., aes(x=pval)) + geom_histogram() + facet_wrap(~method)

save(pvals.gamma.1, file = "code/unconditional_test/pvals_gamma_1.RData")
save(pvals.gamma.2, file = "code/unconditional_test/pvals_gamma_2.RData")
