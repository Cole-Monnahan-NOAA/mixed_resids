
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
  mu <- .75
  sigma <- 1
  tau <- 1
  huge <- 1e3
  ## simulate random track
  nt <- 100
  X <- rnorm(nt,mean=0,sd=tau)
  Ypred <- rep(NA, nt)
  Ypred[1] <- X[1]
  for(t in 2:nt) Ypred[t] <- Ypred[t-1]+X[t]+mu
  ## simulate random measurements
  Y <- rnorm(nt, Ypred, sd=sigma)
  data <- list(y=Y, simRE=0, huge=huge)
  parameters <- list(x=X, mu=0, logsigma=log(sigma), logtau=log(tau))

  message(ii, ": Optimizing two competing models...")
  ## H0: mu estimated, correctly specified model
  obj0 <- MakeADFun(data, parameters, random=c("x"), DLL="randomwalk")
  trash <- obj0$env$beSilent()
  opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
  opt0 <- add_aic(opt0, n=length(Y))
  sdr0 <- sdreport(obj0, getJointPrecision=TRUE)
  rep0 <- obj0$report(obj0$env$last.par.best)
  ## H1: mu=0 so no drift, underspecified model. Note opposite
  ## order from original script.
  obj1 <- MakeADFun(data, parameters, random=c("x"),
                    DLL="randomwalk", map=list(mu=factor(NA)))
  trash <- obj1$env$beSilent()
  ## obj1$fn()
  ## plot(1:nt, Ypred, type='l')
  ## points(1:nt, Y, col=2)
  ## lines(1:nt, obj1$report()$ypred, col=3)
  ## plot(1:nt, obj1$report()$x)
  opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
  opt1 <- add_aic(opt1, n=length(Y))
  sdr1 <- sdreport(obj1, getJointPrecision=TRUE)
  rep1 <- obj1$report(obj1$env$last.par.best)
  opt0$xpred <- obj0$report()$xpred
  opt1$xpred <- obj1$report()$xpred

  ## Save MLEs to test for properties. These are the true pars as
  ## parameterized in the TMB model
  truepars <- c(log(sigma), log(tau))
  mles <- rbind(
    data.frame(version='m0', rep=ii, mle=opt1$par,
               par=names(obj1$par), true=truepars),
    data.frame(version='m0', rep=ii, mle=opt0$par,
               par=names(obj0$par), true=c(mu,truepars)))
  dir.create('results/randomwalk_mles', showWarnings=FALSE)
  saveRDS(mles, file=paste0('results/randomwalk_mles/mles_', ii, '.RDS'))


  message(ii, ": Calculating residuals..")
  osa0 <- calculate.osa(obj0, methods=c('gen','fg', 'osg', 'cdf'), observation.name='y')
  osa1 <- calculate.osa(obj1, methods=c('gen','fg', 'osg', 'cdf'), observation.name='y')

  ## DHARMa resids, both conditional and unconditional
  ## hack to get this to evaluate in a function
  expr <- expression(obj$simulate()$y)
  sim0_cond <-
    calculate.dharma(obj0, expr, obs=Y, fpr=rep0$ypred)
  obj0$env$data$simRE <- 1 #turn on RE simulation
  sim0_uncond <-
    calculate.dharma(obj0, expr, obs=Y, fpr=rep0$ypred)
  sim1_cond <-
    calculate.dharma(obj1, expr, obs=Y, fpr=rep1$ypred)
  obj1$env$data$simRE <- 1 #turn on RE simulation
  sim1_uncond <-
    calculate.dharma(obj1, expr, obs=Y, fpr=rep1$ypred)

  ## Try adding residuals from the joint precisions matrix
  sim0_parcond <- calculate.jp(obj0, sdr0, opt0, Y, 'y', fpr=rep0$ypred)
  sim1_parcond <- calculate.jp(obj1, sdr1, opt1, Y, 'y', fpr=rep1$ypred)
  ## Combine together in tidy format for analysis and plotting later
  r0 <- data.frame(model='randomwalk', replicate=ii, y=Y,
                   ypred=rep0$ypred, version='m0',
                   osa.cdf = osa0$cdf, osa.gen = osa0$gen,
                   osa.fg=osa0$fg, osa.osg=osa0$osg,
                   sim_cond=sim0_cond$resids,
                   sim_uncond=sim0_uncond$resids,
                   sim_parcond=sim0_parcond$resids,
                   maxgrad=max(abs(obj0$gr(opt0$par))),
                   AIC=opt0$AIC, AICc=opt0$AICc)
  r1 <- data.frame(model='randomwalk', replicate=ii, y=Y,
                   ypred=rep1$ypred, version='m1',
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
