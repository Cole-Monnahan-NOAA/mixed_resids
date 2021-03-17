simulate.simpleGLMM <- function(seed, n.j=3, n.i=10){
  ## n.j <- 3 #number of subjects
  ## n.i <- 10 #number of observations
  b0 <- 4
  sig2.y <- .1 #obs variance
  sig2.u <- 5 # between group variance
  set.seed(seed)
  u <- rnorm(n.j, mean=0, sd=sqrt(sig2.u))
  y <- matrix(0, n.i, n.j)
  for(j in 1:n.j){
    y[,j] <- rnorm(n.i, b0 + u[j], sqrt(sig2.y))
  }
  Dat <- data.frame(y = as.vector(y), group = rep(1:n.j, each = n.i))
  Data <- list(y = Dat[,1], group = Dat[,2]-1, sim_re = 0)
  Par <- list(b0=0, ln_sig_u=0, ln_sig_y=0, u=rep(0, n.j))
  return(list(Data=Data, Par=Par))
}



run.simpleGLMM.iter <- function(ii){
  library(TMB)
  library(DHARMa)
  library(INLA)
  library(dplyr)
  library(tidyr)
  library(R.utils)
  dyn.load(dynlib("models/simpleGLMM"))

  ## simulate data with these parameters
  message(ii, ": Simulating data...")
  out <- simulate.simpleGLMM(ii, 5, 8)
  dat0 <- out$Data
  ## Add an outlier to each group
  dat1 <- dat0
  ind <- which(!duplicated(dat0$group))
  set.seed(ii)
  dat1$y[ind] <- dat1$y[ind]+sample(c(-2,2), size=length(ind), replace=TRUE)

  message(ii, ": Optimizing two competing models...")
  ## H0: b0 fixed at zero, underspecified model
  obj0 <- MakeADFun(dat0, out$Par, random = 'u', DLL = 'simpleGLMM')
  trash <- obj0$env$beSilent()
  opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
  opt0 <- add_aic(opt0, n=length(dat0$y))
  sdr0 <- sdreport(obj0, getJointPrecision=TRUE)
  rep0 <- obj0$report(obj0$env$last.par.best)
  ## H1: mu estimated, correctly specified model
  obj1 <- MakeADFun(dat1, out$Par, random = 'u', DLL = 'simpleGLMM')
  trash <- obj1$env$beSilent()
  opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
  opt1 <- add_aic(opt1, n=length(dat1$y))
  sdr1 <- sdreport(obj1, getJointPrecision=TRUE)
  rep1 <- obj1$report(obj1$env$last.par.best)

  ## Save MLEs to test for properties. These are the true pars as
  ## parameterized in the TMB model
  truepars <- c(4, log(sqrt(3)), log(sqrt(.1)))
  mles <- rbind(
    data.frame(version='m0', rep=ii, mle=opt0$par,
               par=names(obj0$par), true=truepars),
    data.frame(version='m1', rep=ii, mle=opt1$par,
               par=names(obj1$par), true=truepars))
  dir.create('results/simpleGLMM_mles', showWarnings=FALSE)
  saveRDS(mles, file=paste0('results/simpleGLMM_mles/mles_', ii, '.RDS'))


  message(ii, ": Calculating residuals..")
  osa0 <- calculate.osa(obj0, methods=c('gen','fg', 'osg', 'cdf'), observation.name='y')
  osa1 <- calculate.osa(obj1, methods=c('gen','fg', 'osg', 'cdf'), observation.name='y')

  ## DHARMa resids, both conditional and unconditional
  ## hack to get this to evaluate in a function
  expr <- expression(obj$simulate()$y)
  sim0_cond <-
    calculate.dharma(obj0, expr, obs=dat0$y, fpr=rep0$mu)
  obj0$env$data$sim_re <- 1 #turn on RE simulation
  sim0_uncond <-
    calculate.dharma(obj0, expr, obs=dat0$y, fpr=rep0$mu)
  sim1_cond <-
    calculate.dharma(obj1, expr, obs=dat1$y, fpr=rep1$mu)
  obj1$env$data$sim_re <- 1 #turn on RE simulation
  sim1_uncond <-
    calculate.dharma(obj1, expr, obs=dat1$y, fpr=rep1$mu)

  ## Try adding residuals from the joint precisions matrix
  sim0_parcond <- calculate.jp(obj0, sdr0, opt0, dat0$y, 'y', fpr=rep0$mu)
  sim1_parcond <- calculate.jp(obj1, sdr1, opt1, dat1$y, 'y', fpr=rep1$mu)

  ## Combine together in tidy format for analysis and plotting later
  r0 <- data.frame(model='simpleGLMM', replicate=ii, y=dat0$y,
                   ypred=rep0$mu, version='m0',
                   osa.cdf = osa0$cdf, osa.gen = osa0$gen,
                   osa.fg=osa0$fg, osa.osg=osa0$osg,
                   sim_cond=sim0_cond$resids,
                   sim_uncond=sim0_uncond$resids,
                   sim_parcond=sim0_parcond$resids,
                   maxgrad=max(abs(obj0$gr(opt0$par))),
                   AIC=opt0$AIC, AICc=opt0$AICc)
  r1 <- data.frame(model='simpleGLMM', replicate=ii, y=dat0$y,
                   ypred=rep1$mu, version='m1',
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
  pvals$replicate <- ii; pvals$model <- 'simpleGLMM'

  ## save to file in case it crashes can recover what did run
  dir.create('results/simpleGLMM_pvals', showWarnings=FALSE)
  dir.create('results/simpleGLMM_resids', showWarnings=FALSE)
  saveRDS(pvals, file=paste0('results/simpleGLMM_pvals/pvals_', ii, '.RDS'))
  saveRDS(resids, file=paste0('results/simpleGLMM_resids/resids_', ii, '.RDS'))
  if(ii==1){
    message("Making plots for replicate 1...")
    library(ggplot2)
    resids.long <- cbind(resids,x=1:length(dat0$y), group=c(dat0$group, dat1$group)) %>%
      pivot_longer(c('osa.cdf', 'osa.gen', 'osa.fg', 'osa.osg',
                     'sim_cond', 'sim_uncond', 'sim_parcond')) %>%
      filter(!is.na(value))
    theme_set(theme_bw())
    ## Plot of data
    png('plots/simpleGLMM_simdata.png', width=7, height=4, units='in', res=200)
    par(mfrow=c(1,2))
    boxplot(y~group, dat0, ylim=range(dat1$y))
    boxplot(y~group, dat1, ylim=range(dat1$y))
    dev.off()
    ## plot of resids
    g <- ggplot(resids.long, aes(x, y=value, color=factor(group))) +
      geom_point() + facet_grid(version~name) +
      labs(x='Order', y='Residual')
    ggsave('plots/simpleGLMM_resids_by_group.png', g, width=9, height=6)
    g <- GGally::ggpairs(resids, columns=6:12, mapping=aes(color=version), title='Random Walk')
    ggsave('plots/simpleGLMM_resids_pairs.png', g, width=7, height=5)
    ## Plot of  DHARMa simulated data look like
    ff <- function(x, v, re) data.frame(x=1:length(dat0$y), x$sims, group=dat0$group, version=v, method=re)
    g <- rbind(ff(sim0_cond, 'm0', 'cond'),
               ff(sim0_parcond, 'm0', 'parcond'),
               ff(sim0_uncond, 'm0', 'uncond'),
               ff(sim1_cond, 'm1', 'cond'),
               ff(sim1_parcond, 'm1', 'parcond'),
               ff(sim1_uncond, 'm1', 'uncond')) %>%
      pivot_longer(cols=c(-x, -group, -version, -method), names_prefix="X",
                   names_to='replicate', values_to='y') %>%
      mutate(replicate=as.numeric(replicate)) %>%
      ggplot(aes(group, y, fill=factor(group))) +
      geom_violin() +
      ## geom_jitter(alpha=.5, pch='.', width=.3, height=0) +
      facet_grid(version~method)
    g <- g+geom_jitter(alpha=.5, width=.2, data=rbind(data.frame(dat0), data.frame(dat1)))
    ggsave('plots/simpleGLMM_simdata.png', g, width=9, height=6)
  }
  return(invisible(list(pvals=pvals, resids=resids)))
}
