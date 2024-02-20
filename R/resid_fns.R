calculate.osa <- function(obj, methods, observation.name,
                          data.term.indicator='keep',
                          Range = c(-Inf,Inf), Discrete = NULL,
                          Subset = NULL){
  ## OSA residuals
  fg <- osg <- cdf <- gen <- pears <- NA
  runtime.fg <- runtime.osg <- runtime.cdf <- runtime.gen <- NA
  if('fg' %in% methods){
    t0 <- Sys.time()
    fg <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     method="fullGaussian", trace=FALSE,
                     discrete = Discrete,
                     subset = Subset)$residual,
      error=function(e) 'error')
    runtime.fg <- as.numeric(Sys.time()-t0, 'secs')
    if(is.character(fg)){
      warning("OSA Full Gaussian failed")
      fg <- NA; runtime.fg <- NA
    }
  }
  ## one step Gaussian method
  if('osg' %in% methods){
    t0 <- Sys.time()
    osg <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     data.term.indicator='keep' ,
                     method="oneStepGaussian", trace=FALSE,
                     discrete = Discrete,
                     subset = Subset)$residual,
      error=function(e) 'error')
    runtime.osg <- as.numeric(Sys.time()-t0, 'secs')
    if(is.character(osg)){
      warning("OSA one Step Gaussian failed")
      osg <- NA; runtime.osg <- NA
    }
  }
  ## cdf method
  if('cdf' %in% methods){
    t0 <- Sys.time()
    cdf <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     data.term.indicator='keep' ,
                     method="cdf", trace=FALSE,
                     discrete = Discrete,
                     subset = Subset)$residual,
      error=function(e) 'error')
    runtime.cdf <- as.numeric(Sys.time()-t0, 'secs')
    if(is.character(cdf)){# | any(!is.finite(cdf))){
      warning("OSA cdf failed")
      cdf <- NA; runtime.cdf <- NA
    }
  }
  ## one step Generic method
  if('gen' %in% methods){
    t0 <- Sys.time()
    gen <- tryCatch(
      oneStepPredict(obj, observation.name=observation.name,
                     data.term.indicator='keep' ,
                     range = Range,
                     method="oneStepGeneric", trace=FALSE,
                     discrete = Discrete,
                     subset = Subset)$residual,
      error=function(e) 'error')
    runtime.gen <- as.numeric(Sys.time()-t0, 'secs')
    if(is.character(gen) | (!is.character(gen) & any(!is.finite(gen)))){
      warning("OSA Generic failed")
      gen <- NA; runtime.gen <- NA
    }
  }
  #Calculate Pearson's residuals 
  pears.df <- length(obj$env$data$y) - length(obj$par)
  if('pears' %in% methods){
    report <- obj$report()
    if(Discrete == TRUE){
      if(obj$env$data$family == 200){#Poisson model
        pears <- (obj$env$data$y - report$exp_val)/sqrt(report$exp_val)
      }
    } else {
      sig <- if(is.null(report$sig)) report$sig_y else report$sig
      pears <- (obj$env$data$y - report$exp_val)/sig
    }
  }
  return(list(gen=gen, fg=fg, osg=osg, cdf=cdf, 
              pears = pears, pears.df = pears.df,
              runtime.gen=runtime.gen, runtime.fg=runtime.fg,
              runtime.osg=runtime.osg, runtime.cdf=runtime.cdf))
}

calculate.dharma <- function(obj, expr, N=1000, obs, idx, fpr, int.resp, rot){
  #alternative <- match.arg(alternative)
  t0 <- Sys.time()
  tmp <- replicate(N, eval(expr)[idx])
  dharma <- createDHARMa(tmp, obs, fittedPredictedResponse = fpr,
                         integerResponse = int.resp,
                         rotation = rot)
  resids <- residuals(dharma, quantileFunction = qnorm,
                      outlierValues = c(-7,7))
  runtime <- as.numeric(Sys.time()-t0, 'secs')
  ## Extract p-values calculated by DHARMa
  ##
  ## Note: Type binomial for continuous, if integer be careful. Not
  ## sure if we want two-sided for dispersion? Using defaults for
  ## now.
  ## AMH: change to alternative = 'greater' when testing for overdispersion in positive only distributions
  ## AMH: Add significance tests
  # disp <- testDispersion(dharma, alternative = alt, plot=FALSE)
  # outlier <- testOutliers(dharma, alternative = alt,
  #                         margin = 'upper', type='binomial', plot=FALSE)
  # pval.ks <-
  #   suppressWarnings(ks.test(dharma$scaledResiduals,'punif')$p.value)
  # pval.ad <- goftest::ad.test(resids,'pnorm', estimated = TRUE)$p.value
  return(list(sims=tmp, resids=resids, out=dharma, runtime=runtime))#, disp=disp$p.value,
              # outlier=outlier$p.value, pval.ks=pval.ks,
              # pval.ad=pval.ad))
}

calc.quantile <- function(fam, report, y){
  out <- NA
  if(fam == 'Gaussian'){
    Fx <- pnorm(q = y, mean = report$exp_val, sd = report$sig_y)
    out <-  qnorm(squeeze(Fx))
  }
  if(fam == 'Lognormal'){
    Fx <- pnorm(q = log(y), mean = report$exp_val, sd = report$sig_y)
    out <-  qnorm(squeeze(Fx))
  }
  if(fam == 'Gamma'){
    Fx <- pgamma(q = y, shape = 1/report$sig_y^2, scale = report$exp_val*report$sig_y^2)
    out <-  qnorm(squeeze(Fx))
  }
  if(fam == 'Poisson'){
    Fx <- ppois(y, report$exp_val)
    px <- dpois(y, report$exp_val)
    u <- runif(length(Fx))
    out <- qnorm(squeeze(Fx) - u * px)
  }
  if(fam == 'NB'){
    Fx <- pnbinom(y, size = report$size, mu = report$exp_val)
    px <- dnbinom(y, size = report$size, mu = report$exp_val)
    u <- runif(length(Fx))
    out <- qnorm(squeeze(Fx) - u * px)
  }  
  if(fam == "Tweedie"){
    grp.idx <- group + 1
    p <- report$power
    phi <- report$sig_y
    Fx <- ptweedie(q = y, mu = report$exp_val, phi = phi, power = p)
    zero.idx <- which(y == 0)
    #as implemented in statmod::qres.tweedie
    Fx[zero.idx] <- runif(length(zero.idx), 0, Fx[zero.idx])
    out <- qnorm(squeeze(Fx))
  }
  if(fam == "Delta_Gamma"){
    zero.idx <- which(y == 0)
    pos.idx <- which(y > 0)
    sig <- tmp$sig_y
    Fx <- rep(0,length(y))
    Fx[pos.idx] <- pgamma(y[pos.idx], shape = 1/report$sig_y^2, 
                          scale = report$exp_val[pos.idx]*report$sig_y^2)
    Fx[zero.idx] <- runif(length(zero.idx), 0, 1)
    out <- qnorm(squeeze(Fx))
  }
  return(out)
}

#duplicated functions, rewrite each test (eg. dispersion, spatial outlier) separately rather than use this code
# calc.dharma.pvals <-
#   function(dharma, alternative = c("two.sided", "greater", "less")){
#     ## Extract p-values calculated by DHARMa
#     ##
#     ## Note: Type binomial for continuous, if integer be careful. Not
#     ## sure if we want two-sided for dispersion? Using defaults for
#     ## now.
#     ## AMH: change to alternative = 'greater' when testing for overdispersion in positive only distributions
#     ## AMH: Add significance tests
#     alternative <- match.arg(alternative)
#     disp <- testDispersion(dharma, alternative, plot=FALSE)
#     outlier <- testOutliers(dharma, alternative,
#                             margin = 'upper', type='binomial', plot=FALSE)
#     resids <- residuals(dharma, quantileFunction = qnorm, outlierValues = c(-7,7))
#     pval.ks <-
#       suppressWarnings(ks.test(dharma$scaledResiduals,'punif')$p.value)
#     pval.ad <- goftest::ad.test(resids,'pnorm', estimated = TRUE)$p.value
#     return(list(disp=disp, outlier=outlier, pval.ks=pval.ks, pval.ad=pval.ad))
#   }
#

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


calculate.jp <- function(obj, sdr, opt, obs, data.name, fpr, N=1000, random = TRUE){
  t0 <- Sys.time()
  joint.mle <- obj$env$last.par.best
  if(random){
    test <- tryCatch(Matrix::Cholesky(sdr$jointPrecision, super=TRUE),
                     error=function(e) 'error')
    if(is.character(test)){
      warning("Joint-Precision approach failed b/c Chol factor failed")
      return(list(sims=NA, runtime=NA, resids=NA, out = NA))#disp=NA, outlier=NA,
                  #pval.ks=NA, pval.ad=NA))
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
    return(list(sims=NA, runtime=NA, resids=NA, out = NA))#disp=NA, outlier=NA,
    #pval.ks=NA, pval.ad=NA))
  }
  dharma <- createDHARMa(tmp, obs, fittedPredictedResponse=fpr)
  resids <- residuals(dharma, quantileFunction = qnorm, outlierValues = c(-7,7))
  runtime <- as.numeric(Sys.time()-t0, 'secs')
  # disp <- testDispersion(dharma, alternative = alternative, plot=FALSE)
  # outlier <- testOutliers(dharma, alternative = alternative,
  #                         margin = 'upper', type='binomial', plot=FALSE)
  # pval.ks <-
  #   suppressWarnings(ks.test(dharma$scaledResiduals,'punif')$p.value)
  # pval.ad <- goftest::ad.test(resids,'pnorm', estimated = TRUE)$p.value
  return(list(sims=tmp, runtime=runtime, resids=resids, out = dharma#, disp=disp$p.value,
              #outlier=outlier$p.value,
              #pval.ks=pval.ks, pval.ad=pval.ad))
  ))
}


calc.sac <- function(res.type, dat, res.obj, version){
  
  if(res.type == 'osa'){
    res.names <- c('cdf', 'gen', 'fg', 'mcmc', 'osg',
                   'pears')
  }
  if(res.type == 'sim'){
    res.names <- c('cond', 'uncond', 'cond_nrot', 'uncond_nrot')
  }
  df <- data.frame(res.type = character(), method = character(), model = character(),
                   test = character(), version = character(), pvalue = numeric())
  
  dmat.obs <- as.matrix(dist(dat$loc, upper = TRUE))
  wt.obs<- 1/dmat.obs; diag(wt.obs) <- 0

 # y <- NA
  
  for(m in 1:length(res.obj)){
    nms <- names(res.obj)[m]
    if(res.type == "osa"){
      x <- res.obj[[m]]
    }
    if(res.type == "sim"){
      x <- res.obj[[m]]$out$scaledResiduals
    }
    if (nms %in% res.names) {
        wt <- wt.obs
      if(is.numeric(x)){
        ## only test for positive correlationa
        y <- ape::Moran.I(x, wt, alternative = 'greater')$p.value
        df <- rbind(df,data.frame(res.type= res.type, 
                         method = names(res.obj)[m], 
                         model='spatial', 
                         test='SAC', 
                         version = version,
                         pvalue = y))
      } 
    }
  }
  if(nrow(df) == 0){
    df <- data.frame(res.type = res.type, method = NA, model = NA,
                     test = 'SAC', version = version, pvalue = NA)
  }

   return(df)
  
}

calc.auto <- function(res.type, res.obj, version){
  
  if(res.type == 'osa'){
    res.names <- c('cdf', 'gen', 'fg', 'mcmc', 'osg',
                   'pears')
  }
  if(res.type == 'sim'){
    res.names <- c('cond', 'uncond', 'cond_nrot', 'uncond_nrot')
  }
  df <- data.frame(res.type = character(), method = character(), model = character(),
                   test = character(), version = character(), pvalue = numeric())
  
  
  for(m in 1:length(res.obj)){
    nms <- names(res.obj)[m]
    if(res.type == "osa"){
      x <- res.obj[[m]]
    }
    if(res.type == "sim"){
      x <- res.obj[[m]]$out$scaledResiduals
    }
    if (nms %in% res.names) {
      if(is.numeric(x)){
        p <- lmtest::dwtest(x[is.finite(x)] ~ 1)$p.value
        df <- rbind(df,data.frame(res.type= res.type, 
                                  method = names(res.obj)[m], 
                                  model='randomwalk', 
                                  test='Auto', 
                                  version = version,
                                  pvalue = p))
      } 
    }
  }
  if(nrow(df) == 0){
    df <- data.frame(type = type, method = NA, model = NA,
                     test = 'Auto', version = version, pvalue = NA)
  }
  
  return(df)
  
}

calc.cat <- function(res.type, group, res.obj, version){
  
  if(res.type == 'osa'){
    res.names <- c('cdf', 'gen', 'fg', 'mcmc', 'osg',
                   'pears')
  }
  if(res.type == 'sim'){
    res.names <- c('cond', 'uncond', 'cond_nrot', 'uncond_nrot')
  }
  df <- data.frame(res.type = character(), method = character(), model = character(),
                   test = character(), version = character(), pvalue = numeric())
  
  for(m in 1:length(res.obj)){
    nms <- names(res.obj)[m]
    if(res.type == "osa"){
      x <- res.obj[[m]]
    }
    if(res.type == "sim"){
      x <- res.obj[[m]]$out$scaledResiduals
    }
    if (nms %in% res.names) {
      if(is.numeric(x)){
        ## only test for positive correlation
        idx <- which(is.finite(x))
        y <- leveneTest(x[idx] ~ factor(group[idx]))$`Pr(>F)`[1]
        df <- rbind(df,data.frame(res.type = res.type, 
                                  method = names(res.obj)[m], 
                                  model = 'simpleGLMM', 
                                  test = 'Cat', 
                                  version = version,
                                  pvalue = y))
      } 
    }
  }
  if(nrow(df) == 0){
    df <- data.frame(type = type, method = NA, model = NA,
                     test = 'Cat', version = version, pvalue = NA)
  }
  
  return(df)
  
}



calc.pvals <- function(res.type, method, mod, res.obj, version, fam, doTrue){
  df <- data.frame(res.type = character(), method = character(), model = character(),
                   test = character(), version = character(), pvalue = numeric())
  if(!is.null(method)){
    if(res.type == 'osa'){
      for(m in 1:length(method)){
        if(is.numeric(res.obj[[method[m]]])){
          #outlier, disp, GOF.ad, GOF.ks !outlier test not available yet for osa
          #remove NA values from residual vector
          res.real <- res.obj[[method[m]]][!is.na(res.obj[[method[m]]])]
          if(doTrue){
            ad <- goftest::ad.test(res.real,'pnorm', estimated = FALSE)$p.value #assume mean=0,sd=1?
          } else {
            ad <- goftest::ad.test(res.real,'pnorm', estimated = TRUE)$p.value
          }
          ks <- suppressWarnings(ks.test(res.obj[[method[m]]],'pnorm')$p.value)
          df <- rbind(df, data.frame(res.type='osa', method=method[m], model=mod, test='GOF.ad', version = version, pvalue = ad))
          df <- rbind(df, data.frame(res.type='osa', method=method[m], model=mod, test='GOF.ks', version = version, pvalue = ks))
        }
      }
      if(!is.null(fam)){
        if(all(fam == 'Poisson' & !is.na(res.obj$pears))){
          disp <- 1 - pchisq(sum(res.obj$pears^2), res.obj$pears.df)
          df <- rbind(df, data.frame(res.type='osa', method='pears', model=mod, test='disp',
                                     version = version, pvalue = disp))
        }
      }
  
    }
    if(res.type == 'sim'){
      alt <- 'two.sided'
      marg <- 'both'
      if(!is.null(fam)){
        if(fam == 'Poisson' | fam == 'Gamma' | fam == 'Lognormal' | fam == 'NB'){
          alt <- 'greater'
          marg <- 'upper'
        }
      }
  
      for(m in 1:length(method)){
        if( all( !is.na(res.obj[[method[m]]]) ) ) {
  
          if(!is.null(fam)){
            if(fam == 'Poisson'){
              disp <- testDispersion(res.obj[[method[m]]]$out, alternative = alt, plot=FALSE)$p.value
  
              #outlier test type = 'bootstrap' when discrete
              outlier <- testOutliers(res.obj[[method[m]]]$out, alternative = alt,
                                      margin = marg, type='bootstrap', plot=FALSE)$p.value
  
              df <- rbind(df, data.frame(res.type='sim', method=method[m], model=mod, test='disp',version = version, pvalue = disp))
              df <- rbind(df, data.frame(res.type='sim', method=method[m], model=mod, test='outlier',version = version, pvalue = outlier))
            } else {
              #outlier test type = 'binomial' only appropriate for continuous distributions
              outlier <- testOutliers(res.obj[[method[m]]]$out, alternative = alt,
                                      margin = marg, type='binomial', plot=FALSE)$p.value
              df <- rbind(df, data.frame(res.type='sim', method=method[m], model=mod, test='outlier',version = version, pvalue = outlier))
            }
          } else {
            #outlier test type = 'binomial' only appropriate for continuous distributions
            outlier <- testOutliers(res.obj[[method[m]]]$out, alternative = alt,
                                    margin = marg, type='binomial', plot=FALSE)$p.value
            df <- rbind(df,  data.frame(res.type='sim',method=method[m], model=mod, test='outlier',version = version, pvalue = outlier))
          }
          
          if(doTrue){
            #use squeeze [0,1] -> (0,1) for ad.test
            ad <- goftest::ad.test(squeeze(res.obj[[method[m]]]$out$scaledResiduals),'punif')$p.value 
          } else {
            ad <- goftest::ad.test(squeeze(res.obj[[method[m]]]$out$scaledResiduals),'punif', estimated = TRUE)$p.value
          }
          ks <- suppressWarnings(ks.test(res.obj[[method[m]]]$out$scaledResiduals,'punif')$p.value)
  
          df <- rbind(df, data.frame(res.type='sim', method=method[m], model=mod, test='GOF.ad', version = version, pvalue = ad))
          df <- rbind(df, data.frame(res.type='sim', method=method[m], model=mod, test='GOF.ks', version = version, pvalue = ks))
        }
      }
    }
  } 
  if(is.null(method) | nrow(df) == 0){
    df <- data.frame(res.type = res.type, method = NA, model = NA,
                     test = NA, version = version, pvalue = NA)
  }
  return(df)
}


