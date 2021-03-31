## functions for simulating data
cMatern <- function(H, Nu, Kap) {
  ifelse(H > 0, besselK(H*Kap, Nu) * (H*Kap)^Nu, 1) / gamma(Nu) * 2^(1-Nu)
}

# Simulate spatial field
sim.omega <- function(Range, sig2, Dmat, Nu = 1, method, mesh){
  Kappa <- sqrt(8)/Range
  Phi <- Range
  Tau <- sqrt(1/(4*pi*Kappa^2*sig2))
  n <- dim(Dmat)[1]

  #Simulate random field and obs
  if(method == 'TMB.matern'){
  #  dyn.load(dynlib('spatial'))
    dat <- list(y = rep(0,n), X = matrix(1, n,1),
                dd = Dmat, nu = Nu, v_i = (1:n)-1, simRE = 1,
                family = 000, link = 2, reStruct = 00)
    dat$spde <- INLA::inla.spde2.matern(mesh)$param.inla[c('M0', 'M1', 'M2')]
    obj <- TMB::MakeADFun(data =  dat,
                     parameters = list(beta = 0, theta = 0, log_tau = log(Tau),
                                       log_kappa = log(Kappa),
                                       omega = rep(0,n)),
                     random = 'omega',
                     DLL = 'spatial')
    sim <- obj$simulate()
    omega <- sim$omega
 #   dyn.unload(dynlib('spatial'))
  }
  if(method == 'TMB.spde'){
 #   dyn.load(dynlib('spatial'))
    dat <- list(y = rep(0,n), X = matrix(1, n,1),
                dd = Dmat, nu = 1, v_i = mesh$idx$loc-1, simRE = 1,
                family = 000, link = 2, reStruct = 10)
    dat$spde <- INLA::inla.spde2.matern(mesh)$param.inla[c('M0', 'M1', 'M2')]
    obj <- TMB::MakeADFun(data =  dat,
                     parameters = list(beta = 0, theta = 0, log_tau = log(Tau),
                                       log_kappa = log(Kappa),
                                       omega = rep(0,mesh$n)),
                     random = 'omega',
                     DLL = 'spatial')
    sim <- obj$simulate()
    omega <- sim$omega
 #   dyn.unload(dynlib('spatial'))
  }
  return(omega)
}


# Simulate data
sim.data <- function(X, Beta, omega, parm, fam, link, Loc){
  Xbeta <- X %*% Beta
  eta <- as.numeric(Xbeta + omega)
  N <- nrow(X)
  if(link == 'identity'){
    mu <- eta
  }
  if(link == 'log'){
    mu <- exp(eta)
  }
  if(fam == 'Gamma'){
    Y <- rgamma(N,1/parm[1]^2,scale=mu*parm[1]^2)
  }
  if(fam == 'Poisson'){
    Y <- rpois(N, mu)
  }
  if(fam == 'Tweedie'){
    Y <- tweedie::rtweedie(N, mean, parm[1], parm[2])
  }
  if(fam == 'Normal'){
    ## assuming parm[1] is SD of sampling process
    Y <- rnorm(N, mu, parm[1])
  }
  return(Y)
}



## Wrapper function to run a single simulation iteration. Called
## in spatial.R using parallel hence the extra stuff
run.spatial.iter <- function(ii){
  library(TMB)
  library(DHARMa)
  library(INLA)
  library(dplyr)
  library(tidyr)
  library(R.utils)
  library(goftest)
  dyn.load(TMB::dynlib("models/spatial")) ## simulate data with these parameters

  message(ii, ": Simulating data...")
  set.seed(ii)
  n <- 100
  sp.var <- 0.5
  sd.obs <- .5
  Range <- 20
  ## Simulate spatial random effects
  Loc <- matrix(runif(n*2,0,100),ncol=2)
  dmat <- as.matrix(dist(Loc))
  mesh <- try(
    withTimeout( INLA::inla.mesh.2d(Loc, max.edge = c(Range/3, Range), offset = c(2, Range*.75)),
                 timeout = 30, onTimeout = 'silent' ))
  if(is.character(mesh)){
    system("Taskkill /IM fmesher.exe /F")
    warning("mesh failed in rep=", ii)
    return(NULL)
  }
  Omega <- sim.omega(Range,sp.var,dmat,method="TMB.spde",mesh=mesh)
  ## simulate random measurements
  ## True beta: an interecept and single covariate
  Beta <- c(1,2)
  X1 <- rep(1, nrow(Loc))               # intercept
  X2 <- rnorm(n=nrow(Loc), 0, 1)        # random
  X <- as.matrix(cbind(X1, X2))
  y0 <- sim.data(X=X, Beta=Beta, omega=Omega[mesh$idx$loc],
                parm=sd.obs, fam='Normal', link='identity')
  ## Add outliers to 5 random values to create second data set
  noutlier <- 5
  y1 <- y0
  ind <- sample(1:n, size=noutlier)
  y1[ind] <- y0[ind]+rnorm(noutlier, 0, sd.obs*4)
  par <- list(beta = 0*Beta, theta = 0, log_tau = 0, log_kappa = 0,
              omega = rep(0,mesh$n))
  dat <- list(y=y0, X=X,
              dd=dmat, nu=1,
              v_i=mesh$idx$loc-1,
              simRE=0, family=000, link=2, reStruct=10)
  dat$spde <- INLA::inla.spde2.matern(mesh)$param.inla[c('M0', 'M1', 'M2')]

  message(ii, ": Optimizing two competing models...")
  ## H0: Space w/ normal; properly specified
  obj0 <- TMB::MakeADFun(dat, par, random=c('omega'), dll="spatial")
  trash <- obj0$env$beSilent()
  opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
  opt0 <- add_aic(opt0, n=length(dat$y))
  sdr0 <- sdreport(obj0, getJointPrecision=TRUE)
  rep0 <- obj0$report(obj0$env$last.par.best)
  ## H1: Space w/ normal and unmodeled outliers (underspecified)
  dat$y <- y1
  obj1 <- TMB::MakeADFun(dat, par, random=c("omega"), dll="spatial")
  trash <- obj1$env$beSilent()
  opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
  opt1 <- add_aic(opt1, n=length(dat$y))
  sdr1 <- sdreport(obj1, getJointPrecision=TRUE)
  rep1 <- obj1$report(obj1$env$last.par.best)
  opt0$ypred <- obj0$report()$mu
  opt1$ypred <- obj1$report()$mu
  ## Save MLEs to test for properties. These are the true pars as
  ## parameterized in the TMB model
  logkappa_true <- log(sqrt(8)/Range)
  logtau_true <- log(sqrt(1/(4*pi*(sqrt(8/Range)^2*sp.var))))
  truepars <- c(Beta, log(sd.obs), logtau_true, logkappa_true)
  mles <- rbind(
    data.frame(version='m0', rep=ii, mle=opt0$par,
               par=names(obj0$par), true=truepars),
    data.frame(version='m1', rep=ii, mle=opt1$par,
               par=names(obj1$par), true=truepars))
  dir.create('results/spatial_mles', showWarnings=FALSE)
  saveRDS(mles, file=paste0('results/spatial_mles/mles_', ii, '.RDS'))

  message(ii, ": Calculating residuals..")
  osa0 <- calculate.osa(obj0, methods=c('gen', 'fg', 'osg', 'cdf'), observation.name='y')
  osa1 <- calculate.osa(obj1, methods=c('gen', 'fg', 'osg', 'cdf'), observation.name='y')

  ## DHARMa resids, both conditional and unconditional
  ## hack to get this to evaluate in a function
  expr <- expression(obj$simulate()$y)
  sim0_cond <-
    calculate.dharma(obj0, expr, obs=y0, fpr=rep0$Xbeta)
  obj0$env$data$simRE <- 1 #turn on RE simulation
  sim0_uncond <-
    calculate.dharma(obj0, expr, obs=y0, fpr=rep0$Xbeta)
  sim1_cond <-
    calculate.dharma(obj1, expr, obs=y1, fpr=rep1$Xbeta)
  obj1$env$data$simRE <- 1 #turn on RE simulation
  sim1_uncond <-
    calculate.dharma(obj1, expr, obs=y1, fpr=rep1$Xbeta)

  ## Try adding residuals from the joint precisions matrix
  sim0_parcond <- calculate.jp(obj0, sdr0, opt0, y0, 'y', fpr=rep0$Xbeta)
  sim1_parcond <- calculate.jp(obj1, sdr1, opt1, y1, 'y', fpr=rep1$Xbeta)

  ## Combine together in tidy format for analysis and plotting later
  r0 <- data.frame(model='spatial', replicate=ii, y0=y0, y1=y1,
                   ypred=opt0$ypred, x=Loc[,1], y=Loc[,2], version='m0',
                   osa.cdf = osa0$cdf, osa.gen = osa0$gen,
                   osa.fg=osa0$fg, osa.osg=osa0$osg,
                   sim_cond=sim0_cond$resids,
                   sim_uncond=sim0_uncond$resids,
                   sim_parcond=sim0_parcond$resids,
                   maxgrad=max(abs(obj0$gr(opt0$par))),
                   AIC=opt0$AIC, AICc=opt0$AICc)
  r1 <- data.frame(model='spatial', replicate=ii, y0=y0, y1=y1,
                   ypred=opt1$ypred, x=Loc[,1], y=Loc[,2], version='m1',
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
  ## Extra ones for spatial model. Only test for positive
  ## correlation
  sac0_cond <- sac0_uncond  <- sac0_parcond <-   sac1_cond <- sac1_uncond  <- sac1_parcond <- NA
  if(!is.na(sim0_cond$resids[1])) sac0_cond <- testSpatialAutocorrelation(sim0_cond$resids, x=Loc[,1], y=Loc[,2], plot=FALSE, alternative='greater', )$p.value
  if(!is.na(sim0_uncond$resids[1]))  sac0_uncond <- testSpatialAutocorrelation(sim0_uncond$resids, x=Loc[,1], y=Loc[,2], plot=FALSE, alternative='greater', )$p.value
  if(!is.na(sim0_parcond$resids[1]))  sac0_parcond <- testSpatialAutocorrelation(sim0_parcond$resids, x=Loc[,1], y=Loc[,2], plot=FALSE, alternative='greater', )$p.value
  if(!is.na(sim1_cond$resids[1]))  sac1_cond <- testSpatialAutocorrelation(sim1_cond$resids, x=Loc[,1], y=Loc[,2], plot=FALSE, alternative='greater', )$p.value
  if(!is.na(sim1_uncond$resids[1]))  sac1_uncond <- testSpatialAutocorrelation(sim1_uncond$resids, x=Loc[,1], y=Loc[,2], plot=FALSE, alternative='greater', )$p.value
  if(!is.na(sim1_parcond$resids[1]))  sac1_parcond <- testSpatialAutocorrelation(sim1_parcond$resids, x=Loc[,1], y=Loc[,2], plot=FALSE, alternative='greater', )$p.value

  ## calculate Moran's I by hand for osa
  w <- 1/dmat;  diag(w) <- 0
  sac0 <- lapply(osa0, function(x) calc.sac(x, w))
  sac1 <- lapply(osa1, function(x) calc.sac(x, w))

  pvals0 <- make.pval.df(osa.pvals0, sim0_cond, sim0_uncond, sim0_parcond)
  ## Add on the SAC ones
  pvals0 <- rbind( pvals0,
                  data.frame(method='osa.fg', test='sac', pvalue=sac0$fg),
                  data.frame(method='osa.osg', test='sac', pvalue=sac0$osg),
                  data.frame(method='osa.cdf', test='sac', pvalue=sac0$cdf),
                  data.frame(method='osa.gen', test='sac', pvalue=sac0$gen),
                  data.frame(method='cond', test='sac', pvalue=sac0_cond),
                  data.frame(method='uncond', test='sac', pvalue=sac0_uncond),
                  data.frame(method='parcond', test='sac', pvalue=sac0_parcond))
  pvals0$version <- 'm0'
  pvals1 <- make.pval.df(osa.pvals1, sim1_cond, sim1_uncond, sim1_parcond)
  ## Add on the SAC ones
  pvals1 <- rbind( pvals1,
                  data.frame(method='osa.fg', test='sac', pvalue=sac1$fg),
                  data.frame(method='osa.osg', test='sac', pvalue=sac1$osg),
                  data.frame(method='osa.cdf', test='sac', pvalue=sac1$cdf),
                  data.frame(method='osa.gen', test='sac', pvalue=sac1$gen),
                  data.frame(method='cond', test='sac', pvalue=sac1_cond),
                  data.frame(method='uncond', test='sac', pvalue=sac1_uncond),
                  data.frame(method='parcond', test='sac', pvalue=sac1_parcond))
  pvals1$version <- 'm1'
  pvals <- rbind(pvals0, pvals1)
  pvals$replicate <- ii; pvals$model <- 'spatial'

  ## Exploratory plots for first replicate
  if(ii==1){
    message("Making plots for replicate 1...")
    library(ggplot2)
    resids.long <- resids %>%
      pivot_longer(c('osa.cdf', 'osa.gen', 'osa.fg', 'osa.osg',
                     'sim_cond', 'sim_uncond', 'sim_parcond')) %>%
      filter(!is.na(value))
    theme_set(theme_bw())
    ## Plot of data
    png('plots/spatial_simdata.png', width=7, height=4, units='in', res=200)
    par(mfrow=c(1,3))
    hist(y0, xlim=range(c(y0,y1)))
    hist(y1, xlim=range(c(y0,y1)))
    plot(y0, y1)
    dev.off()
    g <- rbind(data.frame(x=Loc[,1], y=Loc[,2], z=y0, version='m0'),
               data.frame(x=Loc[,1], y=Loc[,2], z=y1, version='m1')) %>%
      ggplot(aes(x,y, size=abs(z), color=z>0)) + geom_point(alpha=.5) + facet_wrap('version')
    ggsave('plots/spatial_data_example.png', g, width=7, height=5)
    ## plot of resids
    g <- ggplot(resids.long, aes(x, y, size=abs(value), color=value<0)) +
      geom_point(alpha=.5) + facet_grid(version~name)
    ggsave('plots/spatial_resids_by_space.png', g, width=9, height=6)
    g <- GGally::ggpairs(resids, columns=9:15, mapping=aes(color=version), title='Random Walk')
    ggsave('plots/spatial_resids_pairs.png', g, width=7, height=5)
    ## ## Plot of  DHARMa simulated data look like
    ## ff <- function(x, v, re) data.frame(x=Loc[,1], y=Loc[,2], version=v, method=re, x$simulatedResponse[,1:4])
    ## g <- rbind(ff(sim0_cond$resids, 'm0', 'cond'),
    ##            ff(sim0_parcond$resids, 'm0', 'parcond'),
    ##            ff(sim0_uncond$resids, 'm0', 'uncond'),
    ##            ff(sim1_cond$resids, 'm1', 'cond'),
    ##            ff(sim1_parcond$resids, 'm1', 'parcond'),
    ##            ff(sim1_uncond$resids, 'm1', 'uncond')) %>%
    ##   pivot_longer(cols=c(-x,-y, -version, -method), names_prefix="X",
    ##                names_to='replicate', values_to='z') %>%
    ##   mutate(replicate=as.numeric(replicate))  %>%
    ##   ggplot(aes(x, y, size=log(z))) + geom_point(alpha=.2) +
    ##   facet_grid(version+method~replicate)
    ## ggsave('plots/spatial_simdata.png', g, width=9, height=9)
  }
  ## save to file in case it crashes can recover what did run
  dir.create('results/spatial_pvals', showWarnings=FALSE)
  dir.create('results/spatial_resids', showWarnings=FALSE)
  saveRDS(pvals, file=paste0('results/spatial_pvals/pvals_', ii, '.RDS'))
  saveRDS(resids, file=paste0('results/spatial_resids/resids_', ii, '.RDS'))
  return(invisible(list(pvals=pvals, resids=resids)))
}

