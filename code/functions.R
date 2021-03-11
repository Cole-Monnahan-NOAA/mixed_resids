message("Loading global functions...")

calc.osa.pvals <- function(osa){
  fg <- osg <- cdf <- gen <- NA
  if(!is.null(osa$fg)) fg <- suppressWarnings(ks.test(osa.fg,'pnorm')$p.value)
  if(!is.null(osa$osg)) osg <- suppressWarnings(ks.test(osa.osg,'pnorm')$p.value)
  if(!is.null(osa$cdf)) cdf <- suppressWarnings(ks.test(osa.cdf,'pnorm')$p.value)
  if(!is.null(osa$gen)) gen <- suppressWarnings(ks.test(osa.gen,'pnorm')$p.value)
  return(list(fg=fg, osg=osg, cdf=cdf, gen=gen))
}

calc.dharma.pvals <- function(dharma){
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
  pval0_uncond <-
    suppressWarnings(ks.test(dharma0_uncond$scaledResiduals,'punif')$p.value)
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

calculate.jp <- function(sdr, opt, N=1000){
  joint.mle <- obj$env$last.par.best
  test <- tryCatch(Matrix::Cholesky(sdr$jointPrecision, super=TRUE),
                   error=function(e) 'error')
  if(is.character(test)){
    warning("Joint-Precision approach failed b/c Chol factor failed")
    return(NULL)
  }
  jp.sim <- function(){
    newpar <- rmvnorm_prec(mu=joint.mle, prec=sdr$jointPrecision)
    obj$env$data$simRE <-  # turn off RE simulation
      obj$simulate(par=newpar)$y
  }
  tmp <- replicate(N, {jp.sim()})
  dharma_parcond <- createDHARMa(tmp, y, fittedPredictedResponse=NULL)
  resids <- residuals(dharma1_parcond, quantileFunction = qnorm, outlierValues = c(-7,7))
  return(resids)
}

calculate.dharma <- function(obj, expr, N=1000, fpr){
  tmp <- replicate(N, {expr})
  dharma <- createDHARMa(tmp, y, fittedPredictedResponse = fpr)
  resids <- residuals(dharma, quantileFunction = qnorm,
                      outlierValues = c(-7,7))
  return(resids)
}

calculate.osa <- function(obj, methods, ii, observation.name,
                          data.term.indicator='keep'){

  ## OSA residuals
  fg <- osg <- cdf <- gen <- NULL

  if('fg' %in% methods){
    fg <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     method="fullGaussian", trace=FALSE)$residual,
      error=function(e) 'error')
    if(is.character(fg)){
      warning("OSA Full Gaussian failed in rep=", ii)
      fg <- NULL
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
      warning("OSA one Step Gaussian failed in rep=", ii)
      osg <- NULL
    }
  }
  ## cdf method
  if('cdf' %in% methods){
    cdf <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     data.term.indicator='keep' ,
                     method="cdf", trace=FALSE)$residual,
      error=function(e) 'error')
    if(is.character(cdf) | any(is.infinite(cdf))){
      warning("OSA cdf failed in rep=", ii)
      cdf <- NULL
    }
  }
  ## one step Generic method
  if('gen' %in% methods){
    gen <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     data.term.indicator='keep' , range = c(0,Inf),
                     method="oneStepGeneric", trace=FALSE)$residual,
      error=function(e) 'error')
    if(is.character(gen)){
      warning("OSA Generic failed in rep=", ii)
      gen <- NULL
    }
  }
  return(list(gen=gen, fg=fg, osg=osg, cdf=cdf))
}


