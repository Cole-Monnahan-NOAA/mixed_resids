## estimate and validate a random walk model with and without
## drift. based on code example provided by uffe h?gsbro thygesen
## and kasper kristensen, 2016, in the tmb package examples:
## randomwalkvalidation.r.

## modified starting 11/2020 by cole

## randomwalk is conditional, randomwalk2 is unconditional
run.randomwalk.iter <- function(ii){
  library(TMB)
  library(DHARMa)
  library(INLA)
  library(dplyr)
  library(tidyr)
  library(R.utils)
  dyn.load(dynlib("models/randomwalk"))

  ## simulate data with these parameters
  message(ii, ": Simulating data...")
  set.seed(ii)
  mu <- 0.75
  sigma <- 1
  s <- 1
  huge <- 1e3
  ## simulate random track
  nt <- 50
  X <- c(0,cumsum(rnorm(nt-1,mean=mu,sd=sigma)))
  ## simulate random measurements
  Y <- X + rnorm(nt,sd=s)
  data <- list(y=Y,huge=huge)
  parameters <- list(x=X, mu=0, logsigma=log(sigma), logs=log(s))

  ## H0: mu=0 so no drift, underspecified model
  message(ii, ": Optimizing two competing models...")
  obj0 <- MakeADFun(data, parameters, random=c("x"),
                    dll="randomwalk", map=list(mu=factor(NA)))
  trash <- obj0$env$beSilent()
  opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
  opt0 <- add_aic(opt0, n=length(Y))
  sdr0 <- sdreport(obj0, getJointPrecision=TRUE)
  rep0 <- obj0$report(obj0$env$last.par.best)
  ## H1: mu estimated, correctly specified model
  obj1 <- MakeADFun(data, parameters, random=c("x"), dll="randomwalk")
  trash <- obj1$env$beSilent()
  opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
  opt1 <- add_aic(opt1, n=length(Y))
  sdr1 <- sdreport(obj1, getJointPrecision=TRUE)
  rep1 <- obj1$report(obj1$env$last.par.best)
  opt0$xpred <- obj0$report()$xpred
  opt1$xpred <- obj1$report()$xpred

  ## Save MLEs to test for properties. These are the true pars as
  ## parameterized in the TMB model
  truepars <- c(log(sigma), log(s))
  mles <- rbind(
    data.frame(version='m0', rep=ii, mle=opt0$par,
               par=names(obj0$par), true=truepars),
    data.frame(version='m1', rep=ii, mle=opt1$par,
               par=names(obj1$par), true=c(mu,truepars)))
  dir.create('results/randomwalk_mles', showWarnings=FALSE)
  saveRDS(mles, file=paste0('results/randomwalk_mles/mles_', ii, '.RDS'))


  message(ii, ": Calculating residuals..")
  osa0 <- calculate.osa(obj0, methods=c('fg', 'osg', 'cdf'), observation.name='y')
  osa1 <- calculate.osa(obj1, methods=c('fg', 'osg', 'cdf'), observation.name='y')

  ## DHARMa resids, both conditional and unconditional
  ## hack to get this to evaluate in a function
  expr <- expression(obj$simulate()$x2)
  sim0_cond <-
    calculate.dharma(obj0, expr, obs=Y, fpr=rep0$xpred)
  obj0$env$data$simRE <- 1 #turn on RE simulation
  sim0_uncond <-
    calculate.dharma(obj0, expr, obs=Y, fpr=rep0$xpred)
  sim1_cond <-
    calculate.dharma(obj1, expr, obs=Y, fpr=rep1$xpred)
  obj1$env$data$simRE <- 1 #turn on RE simulation
  sim1_uncond <-
    calculate.dharma(obj1, expr, obs=Y, fpr=rep1$xpred)

  ## Try adding residuals from the joint precisions matrix
  sim0_parcond <- calculate.jp(obj0, sdr0, opt0, Y, 'y', fpr=rep0$xpred)
  sim1_parcond <- calculate.jp(obj1, sdr1, opt1, Y, 'y', fpr=rep1$xpred)

  ## Combine together in tidy format for analysis and plotting later
  r0 <- data.frame(model='randomwalk', replicate=ii, y=Y,
                   ypred=opt0$xpred, version='m0',
                   osa.cdf = osa0$cdf, osa.gen = osa0$gen,
                   osa.fg=osa0$fg, osa.osg=osa0$osg,
                   sim_cond=sim0_cond$resids,
                   sim_uncond=sim0_uncond$resids,
                   sim_parcond=sim0_parcond$resids,
                   maxgrad=max(abs(obj0$gr(opt0$par))),
                   AIC=opt0$AIC, AICc=opt0$AICc)
  r1 <- data.frame(model='randomwalk', replicate=ii, y=Y,
                   ypred=opt1$xpred, version='m1',
                   osa.cdf = osa1$cdf, osa.gen = osa1$gen,
                   osa.fg=osa1$fg, osa.osg=osa1$osg,
                   sim_cond=sim1_cond$resids, sim_uncond=sim1_uncond$resids,
                   sim_parcond=sim1_parcond$resids,
                   maxgrad=max(abs(obj1$gr(opt1$par))),
                   AIC=opt1$AIC, AICc=opt1$AICc)
  resids <- rbind(r0, r1)

  ## Calculate p-values. Dharma and JPdone already above
  osa.pvals0 <- calc.osa.pvals(osa0)
  osa.pvals1 <- calc.osa.pvals(osa1)

  pvals0 <- make.pval.df(osa.pvals0, sim0_cond, sim0_uncond, sim0_parcond)
  pvals0$version <- 'm0'
  pvals1 <- make.pval.df(osa.pvals1, sim1_cond, sim1_uncond, sim1_parcond)
  pvals1$version <- 'm1'
  pvals <- rbind(pvals0, pvals1)
  pvals$replicate <- ii; pvals$model <- 'randomwalk'

  ## save to file in case it crashes can recover what did run
  dir.create('results/randomwalk_pvals', showWarnings=FALSE)
  dir.create('results/randomwalk_resids', showWarnings=FALSE)
  saveRDS(pvals, file=paste0('results/randomwalk_pvals/pvals_', ii, '.RDS'))
  saveRDS(resids, file=paste0('results/randomwalk_resids/resids_', ii, '.RDS'))
  return(invisible(list(pvals=pvals, resids=resids)))
}


## Clean up the old runs
unlink('results/randomwalk_resids', TRUE)
unlink('results/randomwalk_pvals', TRUE)
unlink('results/randomwalk_mles', TRUE)

message("Preparing workspace to run ", Nreps, " iterations in parallel...")
TMB::compile("models/randomwalk.cpp") # modified for simulation
sfInit( parallel=TRUE, cpus=cpus )
## sfExport('run.randomwalk.iter', 'sim.omega', 'cMatern', 'sim.data',
##          'rmvnorm_prec', 'add_aic')
sfExportAll()

message("Starting parallel runs...")
results <- sfLapply(1:Nreps, function(ii) run.randomwalk.iter(ii))

## ## Read results back in from file
## fs <- list.files('results/randomwalk_pvals/', full.names=TRUE)
## ## Sometimes they fail to run for some unkonwn reason so try
## ## rerunning those ones once
## if(length(fs)<Nreps){
##   message("Rerunning some failed runs...")
##   bad <- which.failed(Nreps)
##   results <- sfLapply(bad, function(ii) run.randomwalk.iter(ii))
##   fs <- list.files('results/randomwalk_pvals/', full.names=TRUE)
## }
## bad <- which.failed(Nreps)
## if(length(bad)>0) warning(length(bad), " runs failed")

message("Randomwalk: processing and saving final results...")
## Read results back in from file
fs <- list.files('results/randomwalk_pvals/', full.names=TRUE)
pvals <- lapply(fs, readRDS) %>% bind_rows %>% filter(!is.na(pvalue))
saveRDS(pvals, file='results/randomwalk_pvals.RDS')
## Read in residuals
fs <- list.files('results/randomwalk_resids/', full.names=TRUE)
resids <- lapply(fs, readRDS) %>% bind_rows
saveRDS(resids, file='results/randomwalk_resids.RDS')
fs <- list.files('results/randomwalk_mles/', full.names=TRUE)
mles <- lapply(fs, readRDS) %>% bind_rows
saveRDS(mles, file='results/randomwalk_mles.RDS')

