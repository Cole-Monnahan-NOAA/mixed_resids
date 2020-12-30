## estimate and validate a spatial model with and without
## drift. based on code example provided by uffe h?gsbro thygesen
## and kasper kristensen, 2016, in the tmb package examples:
## randomwalkvalidation.r.

## modified starting 11/2020 by cole
## modified to a spatial model by Andrea 12/2020

## Wrapper function to use parallel. Runs a single simulation
## iteration. Called below
run.iter <- function(ii){
  library(TMB)
  library(DHARMa)
  library(INLA)
  library(dplyr)
  library(tidyr)
  dyn.load(TMB::dynlib("models/spatial"))
  ## simulate data with these parameters
  message(ii, ": Simulating data...")
  set.seed(ii)
  n <- 50
  sp.var <- 0.5
  Beta <- 1
  CV <- 0.5
  Range <- 60
  ## Simulate spatial random effects
  Loc <- matrix(runif(n*2,0,100),ncol=2)
  dmat <- as.matrix(dist(Loc))
  mesh <- INLA::inla.mesh.2d(Loc, max.edge = c(Range, Range/3), offset = c(2, Range*.75))
  Omega <- sim.omega(Range,sp.var,dmat,method="TMB.spde",mesh=mesh)
  ## simulate random measurements
  y <- sim.data(X = matrix(1, nrow(Loc),1), Beta=Beta, omega = Omega[mesh$idx$loc],
                parm = CV, fam = 'Gamma', link = 'log')
  dat <- list(y = y, X = matrix(1, n,1),
              dd = dmat, nu = 1, v_i = mesh$idx$loc-1,
              simRE = 0, family = 100, link = 0, reStruct = 10)
  dat$spde <- INLA::inla.spde2.matern(mesh)$param.inla[c('M0', 'M1', 'M2')]
  par <- list(beta = 0, theta = 0, log_tau = 0, log_kappa = 0,
              omega = rep(0,mesh$n))

  ## estimate states and parameters under h0: Cov(omega)~Exp
  message(ii, ": Optimizing two competing models...")
  dummy.mesh <- INLA::inla.mesh.create(matrix(runif(4), ncol=2))
  dat <- list(y = y, X = matrix(1, n,1),
              dd = dmat, nu = 0.5, v_i = (1:n)-1,
              simRE = 0, family = 100, link = 0, reStruct = 00)
  dat$spde <- INLA::inla.spde2.matern(dummy.mesh)$param.inla[c('M0', 'M1', 'M2')]
  par = list(beta = 0, theta = 0, log_tau = 0, log_kappa = 0,
             omega = rep(0,n))
  obj0 <- TMB::MakeADFun(dat, par, random=c("omega"), dll="spatial")
  trash <- obj0$env$beSilent()
  opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
  sdr0 <- sdreport(obj0)
  estOmega0 <- summary(sdr0,"random")
  ## estimate states and parameters under h1: Cov(Omega)~Matern
  dat$nu = 1
  obj1 <- TMB::MakeADFun(dat, par, random=c("omega"), dll="spatial")
  trash <- obj1$env$beSilent()
  opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
  sdr1 <- sdreport(obj1)
  estOmega1 <- summary(sdr1,"random")

  message(ii, ": Calculating residuals..")
  ## OSA residuals
  osa0 <- tryCatch(
    oneStepPredict(obj0, observation.name="y",
                   data.term.indicator='keep' ,
                   method="cdf", trace=FALSE)$residual,
    error=function(e) 'error')
  osa1 <- tryCatch(
    oneStepPredict(obj1, observation.name="y",
                   data.term.indicator='keep' ,
                   method="cdf", trace=FALSE)$residual,
    error=function(e) 'error')
  if(is.character(osa0) | is.character(osa1)){
    warning("OSA failed in rep=", ii)
    next
  }

### DHARMa resids, both conditional and unconditional
  tmp <- replicate(1000, {obj0$simulate()$y})
  dharma0_cond <- createDHARMa(tmp, y, fittedPredictedResponse = rep(opt0$par['beta'],n))
  obj0$env$data$simRE <- 1 #turn on RE simulation
  tmp <- replicate(1000, {obj0$simulate()$y})
  dharma0_uncond <- createDHARMa(tmp, y, fittedPredictedResponse = rep(opt0$par['beta'],n))
  tmp <- replicate(1000, {obj1$simulate()$y})
  dharma1_cond <- createDHARMa(tmp, y, fittedPredictedResponse = rep(opt1$par['beta'], n))
  obj1$env$data$simRE <- 1 #turn on RE simulation
  tmp <- replicate(1000, {obj1$simulate()$y})
  dharma1_uncond <- createDHARMa(tmp, y, fittedPredictedResponse = rep(opt1$par['beta'], n))
  ## warning("don't know the right way to calculate DHARMa resids")
  sim0_cond <- residuals(dharma0_cond, quantileFunction = qnorm, outlierValues = c(-7,7))
  sim0_uncond <- residuals(dharma0_uncond, quantileFunction = qnorm, outlierValues = c(-7,7))
  sim1_cond <- residuals(dharma1_cond, quantileFunction = qnorm, outlierValues = c(-7,7))
  sim1_uncond <- residuals(dharma1_uncond, quantileFunction = qnorm, outlierValues = c(-7,7))

### Combine together in tidy format for analysis and plotting
  d0 <- data.frame(x=Loc[,1], y=Loc[,2], version='m0',# pearson=pearson0,
                   osa=osa0, sim_cond=sim0_cond,
                   sim_uncond=sim0_uncond)
  d1 <- data.frame(x=Loc[,1], y=Loc[,2], version='m1', #pearson=pearson1,
                   osa=osa1, sim_cond=sim1_cond,
                   sim_uncond=sim1_uncond)
  resids <- rbind(d0, d1)
  resids.long <- resids %>% pivot_longer(-c(x, y, version))

### Extract p-values calculated by DHARMa
  ## Note: Type binomial for continuous, if integer be careful. Not
  ## sure if we want two-sided for dispersion? Using defaults for
  ## now.
  ## AMH: change to alternative = 'greater' when testing for overdispersion in positive only distributions
                                        #AMH: Add significance tests
  disp0_uncond <- testDispersion(dharma0_uncond, plot=FALSE)
  outlier0_uncond <- testOutliers(dharma0_uncond, type='binomial', plot=FALSE)
  pval0_uncond <- suppressWarnings(ks.test(dharma0_uncond$scaledResiduals,'punif')$p.value)
  disp1_uncond <- testDispersion(dharma1_uncond, plot=FALSE)
  outlier1_uncond <- testOutliers(dharma1_uncond, type='binomial', plot=FALSE)
  pval1_uncond <- suppressWarnings(ks.test(dharma1_uncond$scaledResiduals,'punif')$p.value)
  disp0_cond <- testDispersion(dharma0_cond, plot=FALSE)
  outlier0_cond <- testOutliers(dharma0_cond, type='binomial', plot=FALSE)
  pval0_cond <- suppressWarnings(ks.test(dharma0_cond$scaledResiduals,'punif')$p.value)
  disp1_cond <- testDispersion(dharma1_cond, plot=FALSE)
  outlier1_cond <- testOutliers(dharma1_cond, type='binomial', plot=FALSE)
  pval1_cond <- suppressWarnings(ks.test(dharma1_cond$scaledResiduals,'punif')$p.value)
                                        #osa
  pval0_osa <- suppressWarnings(ks.test(osa0,'pnorm')$p.value)
  pval1_osa <- suppressWarnings(ks.test(osa1,'pnorm')$p.value)

  pvals <- rbind(
    data.frame(version='m0', RE='cond', test='outlier', pvalue=outlier0_cond$p.value),
    data.frame(version='m0', RE='uncond', test='outlier', pvalue=outlier0_uncond$p.value),
    data.frame(version='m0', RE='osa', test='outlier', pvalue=NA),
    data.frame(version='m0', RE='cond', test='disp', pvalue=disp0_cond$p.value),
    data.frame(version='m0', RE='uncond', test='disp', pvalue=disp0_uncond$p.value),
    data.frame(version='m0', RE='osa', test='disp', pvalue=NA),
    data.frame(version='m0', RE='cond', test='GOF', pvalue=pval0_cond),
    data.frame(version='m0', RE='uncond', test='GOF', pvalue=pval0_uncond),
    data.frame(version='m0', RE='osa', test='GOF', pvalue=pval0_osa),
    data.frame(version='m1', RE='cond', test='outlier', pvalue=outlier1_cond$p.value),
    data.frame(version='m1', RE='uncond', test='outlier', pvalue=outlier1_uncond$p.value),
    data.frame(version='m1', RE='osa', test='outlier', pvalue=NA),
    data.frame(version='m1', RE='cond', test='disp', pvalue=disp1_cond$p.value),
    data.frame(version='m1', RE='uncond', test='disp', pvalue=disp1_uncond$p.value),
    data.frame(version='m1', RE='osa', test='disp', pvalue=NA),
    data.frame(version='m1', RE='cond', test='GOF', pvalue=pval1_cond),
    data.frame(version='m1', RE='uncond', test='GOF', pvalue=pval1_uncond),
    data.frame(version='m1', RE='osa', test='GOF', pvalue=pval1_osa))
  pvals$replicate <- ii

  ## Exploratory plots for first replicate
  if(ii==1){
    library(ggplot2)
    g <- ggplot(resids.long, aes(x, y, size=abs(value), color=value<0)) +
      geom_point(alpha=.5) + facet_wrap('version', ncol=1)
    labs(x='Time Step', y='Residual Value', title='Random Walk')
    ggsave('plots/spatial_resids_by_space.png', g, width=7, height=8)
    g <- GGally::ggpairs(resids, columns=4:6, mapping=aes(color=version), title='Random Walk')
    ggsave('plots/spatial_resids_pairs.png', g, width=7, height=5)

    ## What do the simulated data look like?
    ff <- function(x, v, re) data.frame(x=Loc[,1], y=Loc[,2], version=v, RE=re, x$simulatedResponse[,1:4])
    g <- rbind(ff(dharma0_cond, 'm0', 'cond'),
               ff(dharma0_uncond, 'm0', 'uncond'),
               ff(dharma1_cond, 'm1', 'cond'),
               ff(dharma1_uncond, 'm1', 'uncond')) %>%
      pivot_longer(cols=c(-x,-y, -version, -RE), names_prefix="X",
                   names_to='replicate', values_to='z') %>%
      mutate(replicate=as.numeric(replicate))  %>%
      ggplot(aes(x, y, size=z)) + geom_point(alpha=.2) +
      facet_grid(version+RE~replicate)
    ggsave('plots/spatial_simdata.png', g, width=9, height=9)
  }
  saveRDS(pvals, file=paste0('results/spatial_pvals/pvals_', ii, '.RDS'))
  return(pvals)
}

library(snowfall)
TMB::compile("models/spatial.cpp") # modified for simulation
sfInit( parallel=TRUE, cpus=cpus )
sfExport('run.iter', 'sim.omega', 'cMatern', 'sim.data')
results <- sfLapply(1:Nreps, function(ii) run.iter(ii))

## Read results back in from file
fs <- list.files('results/spatial_pvals', full.names=TRUE)
results <- lapply(fs, readRDS) %>% do.call(rbind, .) %>%
  filter(!is.na(pvalue))
## Did any fail to run? Try to rerun
bad <- which(!1:Nreps %in% results$replicate)
results <- sfLapply(bad, function(ii) run.iter(ii))

## Read results back in from file
fs <- list.files('results/spatial_pvals', full.names=TRUE)
results <- lapply(fs, readRDS) %>% do.call(rbind, .) %>%
  filter(!is.na(pvalue))

g <- ggplot(filter(results, test=='outlier') , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_outlier.png', g, width=5, height=5)
g <- ggplot(filter(results, test=='disp') , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_disp.png', g, width=5, height=5)
g <- ggplot(filter(results, test=='GOF') , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_GOF.png', g, width=5, height=5)


