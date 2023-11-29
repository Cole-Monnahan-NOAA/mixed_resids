# Reproduce analysis for the simpleGLMM model 
# with a uniformally distributed covariate

library(goftest)
library(DHARMa)
library(TMB)
library(tmbstan)
source("R/startup.R")
dyn.load(dynlib('src/simpleGLMM'))

simulate.simpleGLMM <- function(seed = NULL, b0 = 4, 
                                sig.y = 0.5,
                                sig.u = 2,
                                ngroups=5, nobs=100){

  #reproduce simulation steps and simulation seeds
  set.seed(seed)
  X <- cbind(rep(1, nobs), runif(nobs, -.5, .5))
  
  set.seed(seed)
  u <- rnorm(ngroups, 0, sig.u)
  
  set.seed(seed*2)
  y <- matrix(0, nobs, ngroups)
  for(i in 1:nobs){
    for(j in 1:ngroups){
      y[,j] <- rnorm(nobs, b0 + u[j], sig.y)
    }
  }

  Dat <- list(y = as.vector(y), X = X,
              group = rep(1:ngroups, each = nobs) - 1, #C++ starts indexing at 0
              obs = rep(1:nobs, ngroups) - 1, #C++ starts indexing at 0
              family = 0, link = 2, #family = gaussian(link = identity)
              sim_re = 0)
  Par <- list(beta = c(0,0), ln_sig_y = 0, ln_sig_u = 0, ln_sig_v = numeric(0),
              u=rep(0, ngroups), v = rep(0, nobs))
  return(list(Data = Dat, Par = Par, u = u, y = y, sdvec = c(sig.u, sig.y)))
}

#Single model run
n.gr=5; n.obs=100
out <- simulate.simpleGLMM(seed = 1, ngroups=n.gr, nobs=n.obs)
obj <- MakeADFun(out$Data, out$Par, 
                 random = "u", 
                 map = list(v = rep(factor(NA), n.obs)), #fix v at 0
                 DLL = "simpleGLMM", silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
report <- obj$report(obj$env$last.par.best)

#ecdf residuals from DHARMa
#Conditional
dharma.cond <- createDHARMa(replicate(1000, obj$simulate()$y),
                            observedResponse = out$Data$y,
                            fittedPredictedResponse = report$fpr)
ks.test(dharma.cond$scaledResiduals, "punif")

#Unconditional
#turn on RE simulation in TMB model 
#(corresponds to running line 90 in simpleGLMM.cpp)
obj$env$data$sim_re <- 1
dharma.uncond <- createDHARMa(replicate(1000, obj$simulate()$y),
                            observedResponse = out$Data$y,
                            fittedPredictedResponse = report$fpr)
ks.test(dharma.uncond$scaledResiduals, "punif")

#Unconditional with rotation
dharma.uncond.rot <- createDHARMa(replicate(1000, obj$simulate()$y),
                              observedResponse = out$Data$y,
                              fittedPredictedResponse = report$fpr,
                              rotation = "estimated")
ks.test(dharma.uncond.rot$scaledResiduals, "punif")


## Now quick simulation testing of pvalues
nsim <- 5000
parcond.pval <- parcond2.pval <- uncond.pval <- cond.pval <- fg.pval <- posttruth.pval <- postmle.pval <- rep(0,nsim)
for(i in 1:nsim){
  out <- simulate.simpleGLMM(i, ngroups=n.groups, nobs=n.obs)
  obj <- MakeADFun(out$Data, out$Par, random = "u", DLL = "simpleGLMM", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  rep <- obj$report()
  osa <- oneStepPredict(obj,observation.name = 'y',
                        method = 'fullGaussian')
  fg.pval[i] <- ad.test(osa$residual, null = 'pnorm', estimated = TRUE)$p.value[[1]]
  ## Conditional DHARMa
  obj$env$data$sim_re <- 0 #turn off RE simulation
  dharma.cond <- createDHARMa(replicate(1000, obj$simulate()$y),
                              observedResponse=out$Data$y,
                              fittedPredictedResponse = rep$mu)
  resids.cond <- residuals(dharma.cond, quantileFunction = qnorm,
                           outlierValues = c(-7,7))
  cond.pval[i] <- ad.test(resids.cond, null = 'pnorm', estimated = TRUE)$p.value[[1]]
  ## Unconditional DHARMa
  obj$env$data$sim_re <- 1 #turn on RE simulation
  dharma.uncond <- createDHARMa(replicate(1000, obj$simulate()$y),
                                observedResponse=out$Data$y,
                                fittedPredictedResponse = rep$mu)
  resids.uncond <- residuals(dharma.uncond, quantileFunction = qnorm,
                             outlierValues = c(-7,7))
  uncond.pval[i] <- ad.test(resids.uncond, null = 'pnorm', estimated = TRUE)$p.value[[1]]
  ## Partial conditioned in two ways. First is simulating from
  ## full joint posterior (i.e. w/ fixed effects)
  sdr <-  sdreport(obj, getJointPrecision=TRUE)
  dat <- out$Dat
  obj$env$data$simRE <- 0 # turn off RE simulation
  joint.mle <- obj$env$last.par.best
  tmp <- replicate(1000, {
        newpar <- rmvnorm_prec(mu=joint.mle, prec=sdr$jointPrecision)
        obj$simulate(par=newpar)[['y']]
  })
  dharma <- createDHARMa(tmp, dat$y, fittedPredictedResponse=rep$mu)
  resids <- residuals(dharma, quantileFunction = qnorm, outlierValues = c(-7,7))
  parcond.pval[i] <- goftest::ad.test(resids,'pnorm', estimated = TRUE)$p.value
  ## To drop the fixed effects I Think I need to rebuild the obj
  ## with them fixed
  map <- list(b0=factor(NA), ln_sig_u=factor(NA), ln_sig_y=factor(NA))
  pars <- c(as.list(opt$par)); pars$u <- rep(0, len=n.groups)
  objmle <- MakeADFun(out$Data, parameters=pars, map = map, silent=TRUE)
  objmle$env$data$simRE <- 0 # turn off RE simulation
  optmle <- nlminb(objmle$par, objmle$fn, objmle$gr)
  sdrmle <- sdreport(objmle)
  joint.mle <- objmle$env$last.par.best
  tmp <- replicate(1000, {
    newpar <- mvtnorm::rmvnorm(1, mean=joint.mle, sigma=sdrmle$cov.fixed)
    objmle$simulate(par=newpar)[['y']]
  })
  dharma <- createDHARMa(tmp, dat$y, fittedPredictedResponse=rep$mu)
  resids <- residuals(dharma, quantileFunction = qnorm, outlierValues = c(-7,7))
  parcond2.pval[i] <- goftest::ad.test(resids,'pnorm', estimated = TRUE)$p.value
  ## Use posterior draws via tmbstan, after fixing the FE at
  ## either MLE or truth
  ## start with the MLE one
  fitmle <- tmbstan(objmle, chains=1, warmup=500, iter=501, refresh=-1)
  postmle <- as.numeric(as.matrix(fitmle))
  Fx <- pnorm(q=out$Data$y,
              mean=pars$b0+ postmle[out$Data$group+1],
              sd=exp(pars$ln_sig_y))
  residual <- qnorm(Fx)
  postmle.pval[i] <- ad.test(residual, null='pnorm', estimated=TRUE)$p.value[1]
  ## now the truth
  pars <- list(b0=4, ln_sig_u=log(sqrt(out$sig2[2])), ln_sig_y=log(sqrt(out$sig2[1])), u=rep(0, n.groups))
  objtruth <- MakeADFun(out$Data, parameters=pars, map = map)
  fittruth <- tmbstan(objtruth, chains=1, warmup=500, iter=501, refresh=-1)
  posttruth <- as.numeric(as.matrix(fittruth))
  Fx <- pnorm(q=out$Data$y,
              mean=pars$b0+ posttruth[out$Data$group+1],
              sd=exp(pars$ln_sig_y))
  residual <- qnorm(Fx)
  posttruth.pval[i] <- ad.test(residual, null='pnorm', estimated=TRUE)$p.value[1]
}

pvals <-
  tibble(fg.pval, parcond.pval, parcond2.pval, cond.pval, uncond.pval, posttruth.pval,
                postmle.pval) %>%
  pivot_longer(cols=everything(), names_to='method',
                values_to='pvalue') %>%
                mutate(method=gsub('.pval', '', method))
g <- ggplot(pvals, aes(pvalue)) + facet_wrap('method') +
  geom_histogram(bins=20)
ggsave('plots/simpleGLMM_test.png', g, width=7, height=5)


## Explore the simulated data in the parcond2 vs posterior ways
fit <- tmbstan(objmle, iter=1000)
post <- as.data.frame(fit)[,1:10]
tmp <- replicate(2000, {
  newpar <- mvtnorm::rmvnorm(1, mean=joint.mle, sigma=sdrmle$cov.fixed)
  #objmle$simulate(par=newpar)[['y']]
})[1,,] %>% t %>% as.data.frame
names(post) <- names(tmp)
xx <- bind_rows(data.frame(type='post', post),
                data.frame(type='hessian', tmp))
pairs(xx[,-1], col=ifelse(xx$type=='post', 1,2))
## I'm confused why these are the same, but above the residuals
## aren't right for parcond2
