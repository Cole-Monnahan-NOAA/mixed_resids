library(MASS) # rmvnorm
library(TMB)
library(INLA)
library(DHARMa)
## library(VAST)
library(ggplot2)
library(dplyr)
library(tidyr)
## Has ggpairs for comparing residual types
library(GGally)
library(snowfall)
library(R.utils)
library(Matrix)

## Some global settings
ggwidth <- 7
ggheight <- 5
theme_set(theme_bw())

message("Loading functions...")
## Function to simulate parameters from the joint precisions
## matrix (fixed + random effects). Modified from
## FishStatsUtils::simulate_data
rmvnorm_prec <- function(mu, prec ) {
  ##set.seed( random_seed )
  z = matrix(rnorm(length(mu)), ncol=1)
  L = Matrix::Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.vector(z)
  return(mu + z)
}

## Quick fn to check for failed runs by looking at results output
## that doesn't exist
which.failed <- function(Nreps){
  success <- gsub('results/spatial_pvals/pvals_|.RDS', "", x=fs) %>%
    as.numeric()
  fail <- which(! 1:Nreps %in% success)
  fail
}

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
                                       log_kappa = log(Kappa), log_zeta=0,
                                       omega = rep(0,n), u=rep(0,n)),
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
                                       log_kappa = log(Kappa), log_zeta=0,
                                       omega = rep(0,mesh$n), u=rep(0,n)),
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
  dyn.load(TMB::dynlib("models/spatial"))
  ## simulate data with these parameters
  message(ii, ": Simulating data...")
  set.seed(ii)
  n <- 100
  sp.var <- 0.5
  CV <- 1
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
  ## True beta, an interecept and single covariate
  Beta <- c(1,2)
  X1 <- rep(1, nrow(Loc))               # intercept
  X2 <- rnorm(n=nrow(Loc), 0, 1)        # random
  X <- as.matrix(cbind(X1, X2))
  y0 <- sim.data(X=X, Beta=Beta, omega=Omega[mesh$idx$loc],
                parm=CV, fam='Gamma', link='log')
  ## Add overdispersion in form of lognormal
  u <- rlnorm(n=n, meanlog=0, sdlog=1)
  y <- y0*u
  ## par(mfrow=c(1,2))
  ## hist(log(y0), xlim=range(log(y)));
  ## hist(log(y), xlim=range(log(y)))
  par <- list(beta = 0*Beta, theta = 0, log_tau = 0, log_kappa = 0,
              log_zeta=0, omega = rep(0,mesh$n), u=rep(0,n))
  dat <- list(y=y, X=X,
              dd=dmat, nu=1,
              v_i=mesh$idx$loc-1,
              simRE=0, family=100, link=0, reStruct=10)
  dat$spde <- INLA::inla.spde2.matern(mesh)$param.inla[c('M0', 'M1', 'M2')]

  message(ii, ": Optimizing two competing models...")
  ## H0: Space but no overdispersion (underspecified)
  map <- list(log_zeta=factor(NA), u=factor(NA*par$u))
  obj0 <- TMB::MakeADFun(dat, par, random=c('omega', 'u'),
                         dll="spatial", map=map)
  trash <- obj0$env$beSilent()
  opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
  opt0 <- add_aic(opt0, n=length(dat$y))
  sdr0 <- sdreport(obj0, getJointPrecision=TRUE)
  rep0 <- obj0$report(obj0$env$last.par.best)
  ## estimate states and parameters under h1: spatial variance
  map <- list()
  obj1 <- TMB::MakeADFun(dat, par, random=c("omega", 'u'),
                         dll="spatial", map=map)
  trash <- obj1$env$beSilent()
  opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
  opt1 <- add_aic(opt1, n=length(dat$y))
  sdr1 <- sdreport(obj1, getJointPrecision=TRUE)
  rep1 <- obj1$report(obj1$env$last.par.best)
  opt0$ypred <- obj0$report()$mu
  opt1$ypred <- obj1$report()$mu
  opt0$joint.mle <- obj0$env$last.par.best
  opt1$joint.mle <- obj1$env$last.par.best

  message(ii, ": Calculating residuals..")
  ## OSA residuals

  ## Full Gaussian method
  # osa0.fg <- tryCatch(
  #   oneStepPredict(obj0, observation.name="y",
  #                  method="fullGaussian", trace=FALSE)$residual,
  #   error=function(e) 'error')
  # osa1.fg <- tryCatch(
  #   oneStepPredict(obj1, observation.name="y",
  #                  method="fullGaussian", trace=FALSE)$residual,
  #   error=function(e) 'error')
  #
  # if(is.character(osa0.fg) | is.character(osa1.fg)){
  #   warning("OSA Full Gaussian failed in rep=", ii)
  #   return(NULL)
  # }

  ## one step Gaussian method
  # osa0.osg <- tryCatch(
  #   oneStepPredict(obj0, observation.name="y",
  #                  data.term.indicator='keep' ,
  #                  method="oneStepGaussian", trace=FALSE)$residual,
  #   error=function(e) 'error')
  # osa1.osg <- tryCatch(
  #   oneStepPredict(obj1, observation.name="y",
  #                  data.term.indicator='keep' ,
  #                  method="oneStepGaussian", trace=FALSE)$residual,
  #   error=function(e) 'error')
  # if(is.character(osa0.osg) | is.character(osa1.osg)){
  #   warning("OSA one Step Gaussian failed in rep=", ii)
  #   return(NULL)
  # }

  ## CDF method
  osa0.cdf <- tryCatch(
    oneStepPredict(obj0, observation.name="y",
                   data.term.indicator='keep' ,
                   method="cdf", trace=FALSE)$residual,
      error=function(e) 'error')
  osa1.cdf <- tryCatch(
    oneStepPredict(obj1, observation.name="y",
                   data.term.indicator='keep' ,
                   method="cdf", trace=FALSE)$residual,
    error=function(e) 'error')
  if(is.character(osa0.cdf) | is.character(osa1.cdf)){
      warning("OSA CDF failed in rep=", ii)
      return(NULL)
    }
  if(any(is.infinite(osa0.cdf) | any(is.infinite(osa1.cdf)))){
    warning("OSA failed some residuals are infinite in rep=",ii)
    return(NULL)
  }
  ## one step Generic method
  osa0.gen <- tryCatch(
  oneStepPredict(obj0, observation.name="y",
                 data.term.indicator='keep' , range = c(0,Inf),
                 method="oneStepGeneric", trace=FALSE)$residual,
    error=function(e) 'error')
  osa1.gen <- tryCatch(
    oneStepPredict(obj1, observation.name="y",
                   data.term.indicator='keep' , range = c(0,Inf),
                   method="oneStepGeneric", trace=FALSE)$residual,
    error=function(e) 'error')
  if(is.character(osa0.gen) | is.character(osa1.gen)){
    warning("OSA Generic failed in rep=", ii)
    return(NULL)
  }

### DHARMa resids, both conditional and unconditional
  tmp <- replicate(1000, {obj0$simulate()$y})
  dharma0_cond <- createDHARMa(tmp, y, fittedPredictedResponse = rep0$Xbeta)
  obj0$env$data$simRE <- 1 #turn on RE simulation
  tmp <- replicate(1000, {obj0$simulate()$y})
  dharma0_uncond <- createDHARMa(tmp, y, fittedPredictedResponse = rep0$Xbeta)
  tmp <- replicate(1000, {obj1$simulate()$y})
  dharma1_cond <- createDHARMa(tmp, y, fittedPredictedResponse = rep1$Xbeta)
  obj1$env$data$simRE <- 1 #turn on RE simulation
  tmp <- replicate(1000, {obj1$simulate()$y})
  dharma1_uncond <- createDHARMa(tmp, y, fittedPredictedResponse = rep1$Xbeta)
  ## warning("don't know the right way to calculate DHARMa resids")
  sim0_cond <- residuals(dharma0_cond, quantileFunction = qnorm, outlierValues = c(-7,7))
  sim0_uncond <- residuals(dharma0_uncond, quantileFunction = qnorm, outlierValues = c(-7,7))
  sim1_cond <- residuals(dharma1_cond, quantileFunction = qnorm, outlierValues = c(-7,7))
  sim1_uncond <- residuals(dharma1_uncond, quantileFunction = qnorm, outlierValues = c(-7,7))
  ## Try adding residuals from the joint precisions matrix
  test0 <- tryCatch(Matrix::Cholesky(sdr0$jointPrecision, super=TRUE),
                    error=function(e) 'error')
  test1 <- tryCatch(Matrix::Cholesky(sdr1$jointPrecision, super=TRUE),
                    error=function(e) 'error')
  if(is.character(test0) | is.character(test1)){
    warning("Joint-Precision approach failed b/c Chol factor failed for rep=",ii)
    return(NULL)
  }
  jp.sim0 <- function(){
    newpar <- rmvnorm_prec(mu=opt0$joint.mle, prec=sdr0$jointPrecision)
    obj0$env$data$simRE <- 0 # turn off RE simulation
    obj0$simulate(par=newpar)$y
  }
  jp.sim1 <- function(){
    newpar <- rmvnorm_prec(mu=opt1$joint.mle, prec=sdr1$jointPrecision)
    obj1$env$data$simRE <- 0 # turn off RE simulation
    obj1$simulate(par=newpar)$y
  }
  tmp0 <- replicate(1000, {jp.sim0()})
  tmp1 <- replicate(1000, {jp.sim1()})
  dharma0_parcond <- createDHARMa(tmp0, y, fittedPredictedResponse=NULL)
  dharma1_parcond <- createDHARMa(tmp1, y, fittedPredictedResponse=NULL)
  sim0_parcond <- residuals(dharma1_parcond, quantileFunction = qnorm, outlierValues = c(-7,7))
  sim1_parcond <- residuals(dharma1_parcond, quantileFunction = qnorm, outlierValues = c(-7,7))

  ## Combine together in tidy format for analysis and plotting later
  d0 <- data.frame(model='spatial', replicate=ii, ytrue=dat$y,
                   ypred=opt0$ypred,
                   x=Loc[,1], y=Loc[,2], version='m0',
                   #osa=osa0, sim_cond=sim0_cond,
                   #osa.fg = osa0.fg, osa.osg = osa0.osg,
                   osa.cdf = osa0.cdf, osa.gen = osa0.gen,
                   sim_cond=sim0_cond,
                   sim_uncond=sim0_uncond,
                   sim_parcond=sim0_parcond,
                   maxgrad=max(abs(obj0$gr(opt0$par))),
                   AIC=opt0$AIC, AICc=opt0$AICc)
  d1 <- data.frame(model='spatial', replicate=ii, ytrue=dat$y,
                   ypred=opt1$ypred,
                   x=Loc[,1], y=Loc[,2], version='m1',
                   #osa.fg = osa1.fg, osa.osg = osa1.osg,
                   osa.cdf = osa1.cdf, osa.gen = osa1.gen,
                   sim_cond=sim1_cond,
                   sim_uncond=sim1_uncond,
                   sim_parcond=sim0_parcond,
                   maxgrad=max(abs(obj1$gr(opt1$par))),
                   AIC=opt1$AIC, AICc=opt1$AICc)
  resids <- rbind(d0, d1)

  ## Extract p-values calculated by DHARMa
  ##
  ## Note: Type binomial for continuous, if integer be careful. Not
  ## sure if we want two-sided for dispersion? Using defaults for
  ## now.
  ## AMH: change to alternative = 'greater' when testing for overdispersion in positive only distributions
  ## AMH: Add significance tests
  disp0_uncond <- testDispersion(dharma0_uncond, alternative = 'greater', plot=FALSE)
  outlier0_uncond <- testOutliers(dharma0_uncond, alternative = 'greater',
                                  margin = 'upper', type='binomial', plot=FALSE)
  pval0_uncond <- suppressWarnings(ks.test(dharma0_uncond$scaledResiduals,'punif')$p.value)
  disp1_uncond <- testDispersion(dharma1_uncond, alternative = 'greater', plot=FALSE)
  outlier1_uncond <- testOutliers(dharma1_uncond, alternative = 'greater',
                                  margin = 'upper', type='binomial', plot=FALSE)
  pval1_uncond <- suppressWarnings(ks.test(dharma1_uncond$scaledResiduals,'punif')$p.value)
  disp0_cond <- testDispersion(dharma0_cond, alternative = 'greater', plot=FALSE)
  outlier0_cond <- testOutliers(dharma0_cond, alternative = 'greater',
                                margin = 'upper', type='binomial', plot=FALSE)
  sac0_cond <- testSpatialAutocorrelation(dharma0_cond, x=Loc[,1], y=Loc[,2], plot=FALSE, alternative='greater', ) #only test for positive correlation
  pval0_cond <- suppressWarnings(ks.test(dharma0_cond$scaledResiduals,'punif')$p.value)
  disp1_cond <- testDispersion(dharma1_cond, alternative = 'greater', plot=FALSE)
  outlier1_cond <- testOutliers(dharma1_cond, alternative = 'greater',
                                margin = 'upper', type='binomial', plot=FALSE)
  sac1_cond <- testSpatialAutocorrelation(dharma1_cond, x=Loc[,1], y=Loc[,2], plot=FALSE, alternative='greater') #only test for positive correlation
  pval1_cond <- suppressWarnings(ks.test(dharma1_cond$scaledResiduals,'punif')$p.value)
  ## The joint precision resids
  disp0_parcond <- testDispersion(dharma0_parcond, alternative = 'greater', plot=FALSE)
  outlier0_parcond <- testOutliers(dharma0_parcond, alternative = 'greater',
                                margin = 'upper', type='binomial', plot=FALSE)
  sac0_parcond <- testSpatialAutocorrelation(dharma0_parcond, x=Loc[,1], y=Loc[,2], plot=FALSE, alternative='greater') #only test for positive correlation
  pval0_parcond <- suppressWarnings(ks.test(dharma0_parcond$scaledResiduals,'punif')$p.value)
  disp1_parcond <- testDispersion(dharma1_parcond, alternative = 'greater', plot=FALSE)
  outlier1_parcond <- testOutliers(dharma1_parcond, alternative = 'greater',
                                margin = 'upper', type='binomial', plot=FALSE)
  sac1_parcond <- testSpatialAutocorrelation(dharma1_parcond, x=Loc[,1], y=Loc[,2], plot=FALSE, alternative='greater') #only test for positive correlation
  pval1_parcond <- suppressWarnings(ks.test(dharma1_parcond$scaledResiduals,'punif')$p.value)
  ## osa
  ## osa
  # pval0_osa.fg <- suppressWarnings(ks.test(osa0.fg,'pnorm')$p.value)
  # pval0_osa.osg <- suppressWarnings(ks.test(osa0.osg,'pnorm')$p.value)
  pval0_osa.cdf <- suppressWarnings(ks.test(osa0.cdf,'pnorm')$p.value)
  pval0_osa.gen <- suppressWarnings(ks.test(osa0.gen,'pnorm')$p.value)
  #calculate Moran's I by hand for osa
  w <- 1/dmat
  diag(w) <- 0
  # sac0_osa.fg <- ape::Moran.I(osa0.fg, w, alternative = 'greater') #only test for positive correlation
  # sac0_osa.osg <- ape::Moran.I(osa0.osg, w, alternative = 'greater') #only test for positive correlation
  sac0_osa.cdf <- ape::Moran.I(osa0.cdf, w, alternative = 'greater') #only test for positive correlation
  sac0_osa.gen <- ape::Moran.I(osa0.gen, w, alternative = 'greater') #only test for positive correlation
  # pval1_osa.fg <- suppressWarnings(ks.test(osa1.fg,'pnorm')$p.value)
  # pval1_osa.osg <- suppressWarnings(ks.test(osa1.osg,'pnorm')$p.value)
  pval1_osa.cdf <- suppressWarnings(ks.test(osa1.cdf,'pnorm')$p.value)
  pval1_osa.gen <- suppressWarnings(ks.test(osa1.gen,'pnorm')$p.value)
  # sac1_osa.fg <- ape::Moran.I(osa1.fg, w, alternative = 'greater') #only test for positive correlation
  # sac1_osa.osg <- ape::Moran.I(osa1.osg, w, alternative = 'greater') #only test for positive correlation
  sac1_osa.cdf <- ape::Moran.I(osa1.cdf, w, alternative = 'greater') #only test for positive correlation
  sac1_osa.gen <- ape::Moran.I(osa1.gen, w, alternative = 'greater') #only test for positive correlation
  pvals <- rbind(
    data.frame(version='m0', RE='parcond', test='outlier', pvalue=outlier0_parcond$p.value),
    data.frame(version='m0', RE='cond', test='outlier', pvalue=outlier0_cond$p.value),
    data.frame(version='m0', RE='uncond', test='outlier', pvalue=outlier0_uncond$p.value),
    data.frame(version='m0', RE='osa', test='outlier', pvalue=NA),
    data.frame(version='m0', RE='parcond', test='disp', pvalue=disp0_parcond$p.value),
    data.frame(version='m0', RE='cond', test='disp', pvalue=disp0_cond$p.value),
    data.frame(version='m0', RE='uncond', test='disp', pvalue=disp0_uncond$p.value),
    data.frame(version='m0', RE='osa', test='disp', pvalue=NA),
    data.frame(version='m0', RE='parcond', test='sac', pvalue=sac0_parcond$p.value),
    data.frame(version='m0', RE='cond', test='sac', pvalue=sac0_cond$p.value),
    data.frame(version='m0', RE='uncond', test='sac', pvalue=NA),
    # data.frame(version='m0', RE='osa.fg', test='sac', pvalue=sac0_osa.fg$p.value),
    # data.frame(version='m0', RE='osa.osg', test='sac', pvalue=sac0_osa.osg$p.value),
    data.frame(version='m0', RE='osa.cdf', test='sac', pvalue=sac0_osa.cdf$p.value),
    data.frame(version='m0', RE='osa.gen', test='sac', pvalue=sac0_osa.gen$p.value),
    data.frame(version='m0', RE='parcond', test='GOF', pvalue=pval0_parcond),
    data.frame(version='m0', RE='cond', test='GOF', pvalue=pval0_cond),
    data.frame(version='m0', RE='uncond', test='GOF', pvalue=pval0_uncond),
    # data.frame(version='m0', RE='osa.fg', test='GOF', pvalue=pval0_osa.fg),
    # data.frame(version='m0', RE='osa.osg', test='GOF', pvalue=pval0_osa.osg),
    data.frame(version='m0', RE='osa.cdf', test='GOF', pvalue=pval0_osa.cdf),
    data.frame(version='m0', RE='osa.gen', test='GOF', pvalue=pval0_osa.gen),
    data.frame(version='m1', RE='parcond', test='outlier', pvalue=outlier1_parcond$p.value),
    data.frame(version='m1', RE='cond', test='outlier', pvalue=outlier1_cond$p.value),
    data.frame(version='m1', RE='uncond', test='outlier', pvalue=outlier1_uncond$p.value),
    data.frame(version='m1', RE='osa', test='outlier', pvalue=NA),
    data.frame(version='m1', RE='parcond', test='disp', pvalue=disp1_parcond$p.value),
    data.frame(version='m1', RE='cond', test='disp', pvalue=disp1_cond$p.value),
    data.frame(version='m1', RE='uncond', test='disp', pvalue=disp1_uncond$p.value),
    data.frame(version='m1', RE='osa', test='disp', pvalue=NA),
    data.frame(version='m1', RE='parcond', test='sac', pvalue=sac1_parcond$p.value),
    data.frame(version='m1', RE='cond', test='sac', pvalue=sac1_cond$p.value),
    data.frame(version='m1', RE='uncond', test='sac', pvalue=NA),
    # data.frame(version='m1', RE='osa.fg', test='sac', pvalue=sac1_osa.fg$p.value),
    # data.frame(version='m1', RE='osa.osg', test='sac', pvalue=sac1_osa.osg$p.value),
    data.frame(version='m1', RE='osa.cdf', test='sac', pvalue=sac1_osa.cdf$p.value),
    data.frame(version='m1', RE='osa.gen', test='sac', pvalue=sac1_osa.gen$p.value),
    data.frame(version='m1', RE='parcond', test='GOF', pvalue=pval1_parcond),
    data.frame(version='m1', RE='cond', test='GOF', pvalue=pval1_cond),
    data.frame(version='m1', RE='uncond', test='GOF', pvalue=pval1_uncond),
    # data.frame(version='m1', RE='osa.fg', test='GOF', pvalue=pval1_osa.fg),
    # data.frame(version='m1', RE='osa.osg', test='GOF', pvalue=pval1_osa.osg),
    data.frame(version='m1', RE='osa.cdf', test='GOF', pvalue=pval1_osa.cdf),
    data.frame(version='m1', RE='osa.gen', test='GOF', pvalue=pval1_osa.gen))
  pvals$replicate <- ii; pvals$model <- 'spatial'

  ## Exploratory plots for first replicate
  if(ii==1){
    message("Making plots for replicate 1...")
    library(ggplot2)
    resids.long <- resids %>%
      pivot_longer(c('osa.cdf', 'osa.gen', 'sim_cond', 'sim_uncond', 'sim_parcond'))
    theme_set(theme_bw())
    ## Plot of data
    g <- data.frame(x=Loc[,1], y=Loc[,2], z=y) %>%
      ggplot(aes(x,y, size=z)) + geom_point(alpha=.5)
    ggsave('plots/spatial_data_example.png', g, width=7, height=5)
    ## plot of resids
    g <- ggplot(resids.long, aes(x, y, size=abs(value), color=value<0)) +
      geom_point(alpha=.5) + facet_grid(version~name)
    ggsave('plots/spatial_resids_by_space.png', g, width=9, height=6)
    g <- GGally::ggpairs(resids, columns=8:11, mapping=aes(color=version), title='Random Walk')
    ggsave('plots/spatial_resids_pairs.png', g, width=7, height=5)
    ## Plot of  DHARMa simulated data look like
    ff <- function(x, v, re) data.frame(x=Loc[,1], y=Loc[,2], version=v, RE=re, x$simulatedResponse[,1:4])
    g <- rbind(ff(dharma0_cond, 'm0', 'cond'),
               ff(dharma0_parcond, 'm0', 'parcond'),
               ff(dharma0_uncond, 'm0', 'uncond'),
               ff(dharma1_cond, 'm1', 'cond'),
               ff(dharma1_parcond, 'm1', 'parcond'),
               ff(dharma1_uncond, 'm1', 'uncond')) %>%
      pivot_longer(cols=c(-x,-y, -version, -RE), names_prefix="X",
                   names_to='replicate', values_to='z') %>%
      mutate(replicate=as.numeric(replicate))  %>%
      ggplot(aes(x, y, size=log(z))) + geom_point(alpha=.2) +
      facet_grid(version+RE~replicate)
    ggsave('plots/spatial_simdata.png', g, width=9, height=9)
  }
  ## save to file in case it crashes can recover what did run
  dir.create('results/spatial_pvals', showWarnings=FALSE)
  dir.create('results/spatial_resids', showWarnings=FALSE)
  saveRDS(pvals, file=paste0('results/spatial_pvals/pvals_', ii, '.RDS'))
  saveRDS(resids, file=paste0('results/spatial_resids/resids_', ii, '.RDS'))
  return(invisible(pvals))
}

add_aic <- function(opt,n){
  opt$AIC <- TMBhelper::TMBAIC(opt, n=Inf)
  opt$AICc <- TMBhelper::TMBAIC(opt, n=n)
  opt$BIC <- TMBhelper::TMBAIC(opt, p=log(n))
  opt
}
