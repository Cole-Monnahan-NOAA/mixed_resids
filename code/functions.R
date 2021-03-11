
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
  FG <- OSG <- CDF <- GEN <- NULL

  if('FG' %in% methods){
    FG <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     method="fullGaussian", trace=FALSE)$residual,
      error=function(e) 'error')
    if(is.character(FG)){
      warning("OSA Full Gaussian failed in rep=", ii)
      FG <- NULL
    }
  }
  ## one step Gaussian method
  if('OSG' %in% methods){
    OSG <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     data.term.indicator='keep' ,
                     method="oneStepGaussian", trace=FALSE)$residual,
      error=function(e) 'error')
    if(is.character(OSG)){
      warning("OSA one Step Gaussian failed in rep=", ii)
      OSG <- NULL
    }
  }
  ## CDF method
  if('CDF' %in% methods){
    CDF <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     data.term.indicator='keep' ,
                     method="cdf", trace=FALSE)$residual,
      error=function(e) 'error')
    if(is.character(CDF) | any(is.infinite(CDF))){
      warning("OSA CDF failed in rep=", ii)
      CDF <- NULL
    }
  }
  ## one step Generic method
  if('GEN' %in% methods){
    GEN <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     data.term.indicator='keep' , range = c(0,Inf),
                     method="oneStepGeneric", trace=FALSE)$residual,
      error=function(e) 'error')
    if(is.character(GEN)){
      warning("OSA Generic failed in rep=", ii)
      GEN <- NULL
    }
  }
  return(list(GEN=GEN, FG=FG, OSG=OSG, CDF=CDF))
}

}
