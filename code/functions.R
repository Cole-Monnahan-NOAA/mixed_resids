message("Loading global functions...")

make.pval.df <- function(osa, sim_cond, sim_uncond, sim_parcond){
  pvals <- rbind(
    data.frame(method='parcond', test='outlier', pvalue=sim_parcond$outlier),
    data.frame(method='cond', test='outlier', pvalue=sim_cond$outlier),
    data.frame(method='uncond', test='outlier', pvalue=sim_uncond$outlier),
    data.frame(method='parcond', test='disp', pvalue=sim_parcond$disp),
    data.frame(method='cond', test='disp', pvalue=sim_cond$disp),
    data.frame(method='uncond', test='disp', pvalue=sim_uncond$disp),
    data.frame(method='parcond', test='GOF', pvalue=sim_parcond$pval),
    data.frame(method='cond', test='GOF', pvalue=sim_cond$pval),
    data.frame(method='uncond', test='GOF', pvalue=sim_uncond$pval),
    data.frame(method='osa.fg', test='GOF', pvalue=osa$fg),
    data.frame(method='osa.osg', test='GOF', pvalue=osa$osg),
    data.frame(method='osa.cdf', test='GOF', pvalue=osa$cdf),
    data.frame(method='osa.gen', test='GOF', pvalue=osa$gen))
  return(pvals)
}


calc.sac <- function(x, w){
  y <- NA
  if(is.numeric(x)){
    ## only test for positive correlationa
    y <- ape::Moran.I(x, w, alternative = 'greater')$p.value
  }
  return(y)
}


calc.osa.pvals.ks <- function(osa){
  fg <- osg <- cdf <- gen <- NA
  if(is.numeric(osa$fg)) fg <- suppressWarnings(ks.test(osa$fg,'pnorm')$p.value)
  if(is.numeric(osa$osg)) osg <- suppressWarnings(ks.test(osa$osg,'pnorm')$p.value)
  if(is.numeric(osa$cdf)) cdf <- suppressWarnings(ks.test(osa$cdf,'pnorm')$p.value)
  if(is.numeric(osa$gen)) gen <- suppressWarnings(ks.test(osa$gen,'pnorm')$p.value)
  return(list(fg=fg, osg=osg, cdf=cdf, gen=gen))
}

calc.osa.pvals <- function(osa){
  fg <- osg <- cdf <- gen <- NA
  if(is.numeric(osa$fg))
    fg <- goftest::ad.test(osa$fg,'pnorm', estimated = TRUE)$p.value
  if(is.numeric(osa$osg)) osg <- goftest::ad.test(osa$osg,'pnorm', estimated = TRUE)$p.value
  if(is.numeric(osa$cdf)) cdf <- goftest::ad.test(osa$cdf,'pnorm', estimated = TRUE)$p.value
  if(is.numeric(osa$gen)) gen <- goftest::ad.test(osa$gen,'pnorm', estimated = TRUE)$p.value
  return(list(fg=fg, osg=osg, cdf=cdf, gen=gen))
}

calc.dharma.pvals.ks <- function(dharma, alternative = c("two.sided", "greater",
                                                      "less")){
  ## Extract p-values calculated by DHARMa
  ##
  ## Note: Type binomial for continuous, if integer be careful. Not
  ## sure if we want two-sided for dispersion? Using defaults for
  ## now.
  ## AMH: change to alternative = 'greater' when testing for overdispersion in positive only distributions
  ## AMH: Add significance tests
  alternative <- match.arg(alternative)
  disp <- testDispersion(dharma, alternative, plot=FALSE)
  outlier <- testOutliers(dharma, alternative,
                          margin = 'upper', type='binomial', plot=FALSE)
  pval <-
    suppressWarnings(ks.test(dharma$scaledResiduals,'punif')$p.value)
  return(list(disp=disp, outlier=outlier, pval=pval))
}

calc.dharma.pvals <- function(dharma, alternative = c("two.sided", "greater",
                                                      "less")){
  ## Extract p-values calculated by DHARMa
  ##
  ## Note: Type binomial for continuous, if integer be careful. Not
  ## sure if we want two-sided for dispersion? Using defaults for
  ## now.
  ## AMH: change to alternative = 'greater' when testing for overdispersion in positive only distributions
  ## AMH: Add significance tests
  alternative <- match.arg(alternative)
  disp <- testDispersion(dharma, alternative, plot=FALSE)
  outlier <- testOutliers(dharma, alternative,
                          margin = 'upper', type='binomial', plot=FALSE)
  resids <- residuals(dharma, quantileFunction = qnorm, outlierValues = c(-7,7))
  pval <- goftest::ad.test(resids,'pnorm', estimated = TRUE)$p.value
  return(list(disp=disp, outlier=outlier, pval=pval))
}



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


add_aic <- function(opt,n){
  opt$AIC <- TMBhelper::TMBAIC(opt, n=Inf)
  opt$AICc <- TMBhelper::TMBAIC(opt, n=n)
  opt$BIC <- TMBhelper::TMBAIC(opt, p=log(n))
  opt
}

calculate.jp <- function(obj, sdr, opt, obs, data.name, fpr, N=1000, random = TRUE,
                         alternative = c("two.sided", "greater","less")){
  alternative = match.arg(alternative)
  joint.mle <- obj$env$last.par.best
  if(random){
  test <- tryCatch(Matrix::Cholesky(sdr$jointPrecision, super=TRUE),
                   error=function(e) 'error')
  if(is.character(test)){
    warning("Joint-Precision approach failed b/c Chol factor failed")
    return(list(sims=NA, resids=NA, disp=NA, outlier=NA, pval=NA))
  }
  jp.sim <- function(){
    newpar <- rmvnorm_prec(mu=joint.mle, prec=sdr$jointPrecision)
    obj$env$data$simRE <- 0 # turn off RE simulation
    obj$simulate(par=newpar)[[data.name]]
  }
  ## newpars <- replicate(1000, {rmvnorm_prec(mu=joint.mle, prec=sdr$jointPrecision)})
  ## pairs(t(newpars))
  newpar <- rmvnorm_prec(mu=joint.mle, prec=sdr$jointPrecision)
  } else {
    jp.sim <- function(){
      newpar <- mvtnorm::rmvnorm(1, sdr$par.fixed, sdr$cov.fixed)
      obj$env$data$simRE <- 0 # turn off RE simulation
      obj$simulate(par=newpar)[[data.name]]
    }
  }
  tmp <- replicate(N, {jp.sim()})
  if(any(is.nan(tmp))){
    warning("NaN values in JP simulated data")
    return(list(sims=NA, resids=NA, disp=NA, outlier=NA, pval=NA))
  }
  dharma <- createDHARMa(tmp, obs, fittedPredictedResponse=fpr)
  resids <- residuals(dharma, quantileFunction = qnorm, outlierValues = c(-7,7))
  disp <- testDispersion(dharma, alternative = alternative, plot=FALSE)
  outlier <- testOutliers(dharma, alternative = alternative,
                          margin = 'upper', type='binomial', plot=FALSE)
  pval <- goftest::ad.test(resids,'pnorm', estimated = TRUE)$p.value
  return(list(sims=tmp, resids=resids, disp=disp$p.value, outlier=outlier$p.value, pval=pval))
}



#AMH: repetitive with calc.dharma.pvals...can we condense?
calculate.dharma <- function(obj, expr, N=1000, obs, fpr,
                             alternative = c("two.sided", "greater","less")){
  alternative = match.arg(alternative)
  tmp <- replicate(N, eval(expr))
  dharma <- createDHARMa(tmp, obs, fittedPredictedResponse = fpr)
  resids <- residuals(dharma, quantileFunction = qnorm,
                      outlierValues = c(-7,7))

  ## Extract p-values calculated by DHARMa
  ##
  ## Note: Type binomial for continuous, if integer be careful. Not
  ## sure if we want two-sided for dispersion? Using defaults for
  ## now.
  ## AMH: change to alternative = 'greater' when testing for overdispersion in positive only distributions
  ## AMH: Add significance tests
  disp <- testDispersion(dharma, alternative = alternative, plot=FALSE)
  outlier <- testOutliers(dharma, alternative = alternative,
                          margin = 'upper', type='binomial', plot=FALSE)
  pval <- goftest::ad.test(resids,'pnorm', estimated = TRUE)$p.value
  return(list(sims=tmp, resids=resids, disp=disp$p.value, outlier=outlier$p.value, pval=pval))
}

calculate.osa <- function(obj, methods, observation.name,
                          data.term.indicator='keep',
                          Range = c(-Inf,Inf)){
  ## OSA residuals
  fg <- osg <- cdf <- gen <- NA

  if('fg' %in% methods){
    fg <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     method="fullGaussian", trace=FALSE)$residual,
      error=function(e) 'error')
    if(is.character(fg)){
      warning("OSA Full Gaussian failed")
      fg <- NA
    }
  }
  ## one step Gaussian method
  if('osg' %in% methods){
    osg <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     data.term.indicator='keep' ,
                     method="oneStepGaussian", trace=FALSE)$residual,
      error=function(e) 'error')
    if(is.character(osg)){
      warning("OSA one Step Gaussian failed")
      osg <- NA
    }
  }
  ## cdf method
  if('cdf' %in% methods){
    cdf <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     data.term.indicator='keep' ,
                     method="cdf", trace=FALSE)$residual,
      error=function(e) 'error')
    if(is.character(cdf) | any(!is.finite(cdf))){
      warning("OSA cdf failed")
      cdf <- NA
    }
  }
  ## one step Generic method
  if('gen' %in% methods){
    gen <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     data.term.indicator='keep' ,
                     range = Range,
                     ##! range = c(0,Inf) only when obs>0 ,
                     method="oneStepGeneric", trace=FALSE)$residual,
      error=function(e) 'error')
    if(is.character(gen) | (!is.character(gen) & any(!is.finite(gen)))){
      warning("OSA Generic failed")
      gen <- NA
    }
  }
  return(list(gen=gen, fg=fg, osg=osg, cdf=cdf))
}


