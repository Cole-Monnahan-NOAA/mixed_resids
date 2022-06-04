setup_trueparms <- function(mod, misp){
  true.comp <- list()
  if(mod=='linmod'){
    theta <- c(4,-5)
    sd.vec <- 1
    true.comp[[1]] <- true.comp[[2]] <-
      list(beta_1 = theta[1], beta_2 = theta[2],
                      ln_sig = log(sd.vec))
    sp.parm <- 0
    fam <- NULL
    link <- NULL
    if(misp == 'overdispersion'){
      #sd.vec <- c(sd.vec, sd.vec*log(4))
      sd.vec <- c(sd.vec, 1)
    }
    if(misp == 'misscov'){
      true.comp[[2]] <- list(beta = theta[1],
                          ln_sig = log(sd.vec[1]))
    }
  }
  if(mod=='randomwalk'){
    theta <- 2
    sd.vec <- c(1,1)
    sp.parm <- 0
    fam <- NULL
    link <- NULL
    true.comp[[1]] <- list(mu = theta, ln_sig = log(sd.vec[1]),
           ln_tau = log(sd.vec[2]))
    if(misp == 'mu0'){
      true.comp[[2]] <- list(ln_sig = log(sd.vec[1]),
                             ln_tau = log(sd.vec[2]))
    }
    if(misp == 'outliers'){
      true.comp[[2]] <- true.comp[[1]]
    }
  }
  if(mod=='simpleGLMM'){
    #parms when fam = 'Gaussian';link='identity'
    theta <- 1
    sd.vec <- sqrt(c(2,10))
    fam <- 'Tweedie'
    pow <- 1.5
    link <- 'log'

    #parms when fam = 'Poisson'; link = 'log'
    # theta <- 1.5
    # sd.vec <- sqrt(c(.5,2))
    # fam <- 'Poisson'
    # link <- 'log'

    sp.parm <- 0
    if(misp=='misscov'){
       #parms when fam = 'Gaussian';link='identity
       #theta <- c(4,-5)
       theta <- c(4,-.4)
       true.comp[[1]] <- list(beta_1 = theta[1], beta_2 = theta[2],
                              ln_sig_y = log(sd.vec[1]),
                              ln_sig_u = log(sd.vec[2]))
       true.comp[[2]] <- list(beta = theta[1],
                              ln_sig_y = log(sd.vec[1]),
                              ln_sig_u = log(sd.vec[2]))
    }
    if(misp=='overdispersion') {
      sd.vec <- c(sd.vec, 1)
      true.comp[[1]] <- list(beta = theta,
                             ln_sig_y = log(sd.vec[1]),
                             ln_sig_u = log(sd.vec[2]),
                             ln_sig_v = log(sd.vec[3]))
      true.comp[[2]] <- list(beta = theta[1],
                             ln_sig_y = log(sd.vec[1]),
                             ln_sig_u = log(sd.vec[2]))
    }
    if(misp == 'outliers'){
      true.comp[[1]] <- true.comp[[2]] <-
        list(beta = theta[1],
             ln_sig_y = log(sd.vec[1]),
             ln_sig_u = log(sd.vec[2]))
    }
    if(fam == 'Poisson'){
      true.comp[[1]]$ln_sig_y <- NULL
      true.comp[[2]]$ln_sig_y <- NULL
    }
    if(fam == 'Tweedie'){
      sd.vec[3] <- pow
      if(misp == 'deltagamma'){
        #convert sd of tweedie to sd of gamma, tweedie: var(y) = mu^pow*phi
        sd.y <- sqrt(theta^pow*sd.vec[1]) 
        true.comp[[1]] <- list(beta = theta,
                               ln_sig_y = log(sd.vec[1]),
                               ln_sig_u = log(sd.vec[2]))
        true.comp[[2]] <- list(beta = theta,
                               ln_sig_y = log(sd.y/theta), #cv
                               ln_sig_u = log(sd.vec[2]))
      }
      true.comp[[1]]$ln_sig_y[2] <- log((pow-1)/(2-pow))
      #lambda definition based on Compound Poisson-Gamma relationship to Tweedie
      #Prob(Y==0) = exp(-lambda) based on Poisson 
      lambda <- (exp(theta)^(2 - pow))/((2-pow)*sd.vec[1])
      pz <- exp(-lambda)
      #logit transform of the probability of zero
      true.comp[[2]]$ln_sig_y[2] <- log(pz/(1-pz))
    }
  }
  if(mod=='spatial'){
    theta=0.5
    sd.vec <- c(.5,sqrt(2))
    sp.parm <- 50
    fam <- 'Poisson'
    link <- 'log'

    true.comp[[1]] <- true.comp[[2]] <-
      list(beta = theta, theta = log(sd.vec[1]),
           ln_tau = log(1/(2*sqrt(pi)*sqrt(8)/sp.parm*sd.vec[2])), #1/(2*sqrt(pi)*kappa*sp.sd))
           ln_kappa = log(sqrt(8)/sp.parm))

    if(misp=='misscov'){
      theta <- c(1,2)
      true.comp[[1]] <-
        list(beta_1 = theta[1], beta_2 = theta[2], theta = log(sd.vec[1]),
             ln_tau = log(1/(2*sqrt(pi)*sqrt(8)/sp.parm*sd.vec[2])), #1/(2*sqrt(pi)*kappa*sp.sd))
             ln_kappa = log(sqrt(8)/sp.parm))
      true.comp[[2]] <-
        list(beta = theta[1], theta = log(sd.vec[1]),
             ln_tau = log(1/(2*sqrt(pi)*sqrt(8)/sp.parm*sd.vec[2])), #1/(2*sqrt(pi)*kappa*sp.sd))
             ln_kappa = log(sqrt(8)/sp.parm))
    }
    if(misp=='overdispersion'){
      sd.vec <- c(sd.vec, sd.vec[2]*.5) #overdispersion variance = 0.5*spatial variance
      true.comp[[1]]$ln_sig_v <- log(sd.vec[3])
    }
    if(misp == 'dropRE'){
      true.comp[[2]] <- list(beta = theta, theta = log(sd.vec[1]))
    }
    if(misp == 'mispomega'){
      sd.vec[2] <- sd.vec[2]/3
    }
    if(fam == "Poisson"){
      true.comp[[1]]$theta <- NULL
      true.comp[[2]]$theta <- NULL
    }
  }
  true.pars <- list(theta=theta, sd.vec=sd.vec, sp.parm=sp.parm,
                    fam=fam, link=link, true.comp=true.comp)
  return(true.pars)
}

run_iter <- function(ii, n=100, ng=0, mod, cov.mod = 'norm', misp, do.true = FALSE, savefiles=TRUE){
  library(TMB)
  library(DHARMa)
  library(INLA)
  library(dplyr)
  library(tidyr)
  library(R.utils)
  library(goftest)
  library(tweedie)

  if(mod == 'linmod'){
    Random <- FALSE
  } else {
    Random <- TRUE
  }

  setupTMB(mod)
  true.parms <- setup_trueparms(mod,misp)

  ## simulate data with these parameters
  message(ii, ": Simulating data...")
  sim.dat <- simdat(n, ng, mod, cov.mod, true.parms, misp, ii)
  if(is.null(sim.dat)){
    set.seed(ii)
    while(is.null(sim.dat)){
      new.seed = round(runif(1, 5000, 8000))
      sim.dat <-  simdat(n, ng, mod, cov.mod, true.parms, misp, new.seed)
    }
  }

  init.dat <- mkTMBdat(sim.dat, true.parms, mod, misp)
  init.par <- mkTMBpar(true.parms, sim.dat, mod, misp, do.true)
  init.random <- mkTMBrandom(mod, misp, do.true)
  init.map <- mkTMBmap(init.par, mod, misp, true.parms$fam, do.true)
  mod.out <- osa.out <- dharma.out <- list(h0 = NULL, h1 = NULL)
  pvals <- data.frame(id = character(), type = character(), method = character(),
                      model = character(), test = character(), 
                      version = character(), pvalue = numeric())
  mles <-  r <- out <- list()

  for(h in 1:2){

    message(ii, ": Optimizing two competing models...")

    id <- paste0(mod, '_', misp, '_', do.true, '_h', h-1, '_', ii)

    init.obj <- list(data = init.dat[[h]], parameters = init.par[[h]], map = init.map[[h]], random = init.random[[h]], DLL = mod)
    
    if(is.null(init.random[[h]])) Random <- FALSE
    
    mod.out[[h]] <- try(
      fit_tmb(
        obj.args = init.obj,
        control = list(run.model = !do.true, do.sdreport = TRUE)
        )
    )
    if(class(mod.out[[h]])!='try-error'){
      if(!do.true){
        ## if estimating, return MLE values
        tmp1 <- true.parms$true.comp[[h]]
        tmp2 <- mod.out[[h]]$opt$par
        stopifnot(length(tmp2)>0)
        ##if('fam' %in% names(tmp1)){
        tmp1$fam <- NULL
        tmp1$link <- NULL
        ##}
        tmp1 <- unlist(tmp1)
        ## super hacky way to get unique names when there are vectors
        names(tmp2) <- unlist(sapply(unique(names(tmp2)), function(x) {
          y <- names(tmp2)[names(tmp2)==x]
          if(length(y)>1) paste(y, 1:length(y), sep="_") else y
        }))
        ## Save the true values and estimated ones to file
        mles[[h]] <- data.frame(h=h-1, par=names(tmp2),
                                mle=as.numeric(tmp2),
                                true=as.numeric(tmp1),
                                bias = tmp2-tmp1)
        mles[[h]] <- cbind(mles[[h]], id=id, replicate=ii,
                           do.true=do.true,  model=mod, misp=misp,
                           n=n)
      } else {
        ## otherwise just NULL b/c nothing estimated
        mles[[h]] <- NULL
      }

      message(ii, ": Calculating residuals..")
      disc <- FALSE; ran <- c(-Inf,Inf); rot <- NULL
      if(!is.null(true.parms$fam)){
        if(true.parms$fam == 'Poisson'){
          disc <- TRUE
          ran <- c(0,Inf)
        }
        if(true.parms$fam == 'Gamma') ran <- c(0,Inf)
      }
      if(mod == 'randomwalk' | mod == 'spatial') rot <- "estimated"
      
      osa.out[[h]] <- calculate.osa(mod.out[[h]]$obj, methods=osa.methods, observation.name='y', Discrete = disc, Range = ran)

      expr <- expression(obj$simulate()$y)
     
      if('cond' %in% dharma.methods){
        dharma.out[[h]]$cond <- 
          calculate.dharma(mod.out[[h]]$obj, expr, obs=sim.dat[[h]],
                           fpr=mod.out[[h]]$report$fpr, 
                           int.resp = disc, rot = NULL)
      } else {
        dharma.out[[h]]$cond <- list()
        dharma.out[[h]]$cond$resids <- NA
        dharma.out[[h]]$cond$runtime <- NA
      }
      
      if('uncond_nrot' %in% dharma.methods){
        mod.out[[h]]$obj$env$data$sim_re <- 1 #turn on RE simulation
        dharma.out[[h]]$uncond_nrot <- 
          calculate.dharma(mod.out[[h]]$obj, expr, obs=sim.dat[[h]], 
                           fpr=mod.out[[h]]$report$fpr, 
                           int.resp = disc, rot = NULL)
      } else {
        dharma.out[[h]]$uncond_nrot <- list()
        dharma.out[[h]]$uncond_nrot$resids <- NA
        dharma.out[[h]]$uncond_nrot$runtime <- NA
      }

      if('uncond' %in% dharma.methods){
        mod.out[[h]]$obj$env$data$sim_re <- 1 #turn on RE simulation
        dharma.out[[h]]$uncond <- 
          calculate.dharma(mod.out[[h]]$obj, expr, obs=sim.dat[[h]], 
                           fpr=mod.out[[h]]$report$fpr, 
                           int.resp = disc, rot = rot)
      } else {
        dharma.out[[h]]$uncond <- list()
        dharma.out[[h]]$uncond$resids <- NA
        dharma.out[[h]]$uncond$runtime <- NA
      }
      
      if('re_uncond' %in% dharma.methods & mod != 'linmod' & Random ){
        mod.out[[h]]$obj$env$data$sim_re <- 1 #turn on RE simulation
        if(mod == 'spatial'){
          expr <- expression(obj$simulate()$omega)  
          obs <- mod.out[[h]]$obj$env$parList()$omega
        } else {
          #center simulation and obs
          expr <- expression(obj$simulate()$u )
          obs <- mod.out[[h]]$obj$env$parList()$u
        }
        if(mod == 'randomwalk'){
          fpr <- mod.out[[h]]$report$ypred
        } else {
          fpr <- rep(0, length(obs))
        }
        dharma.out[[h]]$re_uncond <- 
          calculate.dharma(mod.out[[h]]$obj, expr, obs=obs, 
                           fpr = fpr, 
                           int.resp = disc, rot = rot)
      } else {
        dharma.out[[h]]$re_uncond <- list()
        dharma.out[[h]]$re_uncond$resids <- NA
        dharma.out[[h]]$re_uncond$runtime <- NA
      }
      
      if('re_uncond_nrot' %in% dharma.methods & mod != 'linmod' & Random ){
        mod.out[[h]]$obj$env$data$sim_re <- 1 #turn on RE simulation
        if(mod == 'spatial'){
          expr <- expression(obj$simulate()$omega)  
          obs <- mod.out[[h]]$obj$env$parList()$omega
          fpr <- rep(0, length(obs))
        } else {
          #center simulation and obs
          expr <- expression(obj$simulate()$u )
          obs <- mod.out[[h]]$obj$env$parList()$u
          if(mod == 'randomwalk'){
            fpr <- mod.out[[h]]$report$ypred
          } else {
            fpr <- rep(0, length(obs))
          }
        }
        dharma.out[[h]]$re_uncond_nrot <- 
          calculate.dharma(mod.out[[h]]$obj, expr, obs=obs, 
                           fpr = fpr, 
                           int.resp = disc, rot = NULL)
      } else {
        dharma.out[[h]]$re_uncond_nrot <- list()
        dharma.out[[h]]$re_uncond_nrot$resids <- NA
        dharma.out[[h]]$re_uncond_nrot$runtime <- NA
      }

      ## only makes sense to run when MLE is estimated and there
      ## are RE - turning on when do.true == TRUE (AMH, 6/3/2022)
      if('mcmc' %in% osa.methods & Random ){
        t0 <- Sys.time()
        ## Build up TMB obj again
        if(do.true)  FE <- mod.out[[h]]$obj$par # true FE
        if(!do.true) FE <- mod.out[[h]]$opt$par # estimated FE
        ## make into list; https://stackoverflow.com/questions/46251725/convert-named-vector-to-list-in-r/46251794
        FE <- split(unname(FE),names(FE))
        MLE <- modifyList(init.par[[h]], FE) #
        ## Get FE and map them off
        xx <- names(MLE)[-which(names(MLE) %in% init.random[[h]])]
        map <- lapply(names(FE), function(x) factor(FE[[x]]*NA))
        names(map) <- names(FE)
        map <- c(map, init.map[[h]])
        ## Rebuild with original FE mapped off and RE as FE
        objmle <- MakeADFun(data=init.dat[[h]], parameters=MLE,
                            map=map, DLL=mod.out[[h]]$obj$env$DLL)
        fitmle <- tmbstan::tmbstan(objmle, chains=1, warmup=300, iter=301, seed=ii, refresh=-1)
        postmle <- as.numeric(as.matrix(fitmle)) ## single sample
        postmle <- postmle[-length(postmle)] # drop lp__ value
        ## Calculate residuals given the sample
        tmp <- objmle$report(postmle)
        if(is.null(true.parms$fam)){
          sig <- if(is.null(tmp$sig)) tmp$sig_y else tmp$sig
          Fx <- pnorm(q=init.dat[[h]]$y, mean= tmp$exp_val, sd=sig)
          osa.out[[h]]$mcmc <-  qnorm(Fx)
        } else {
          if(true.parms$fam == 'Poisson'){
            Fx <- ppois(init.dat[[h]]$y, tmp$exp_val)
            px <- dpois(init.dat[[h]]$y, tmp$exp_val)
            u <- runif(length(Fx))
            osa.out[[h]]$mcmc <- qnorm(Fx - u * px)
          }
          if(true.parms$fam == 'Gamma'){
            sig <- if(is.null(tmp$sig)) tmp$sig_y else tmp$sig
            Fx <- pgamma(init.dat[[h]]$y, shape = 1/sig^2, scale = tmp$exp_val*sig^2)
            osa.out[[h]]$mcmc <-  qnorm(Fx)
          }
          if(true.parms$fam == 'Gaussian'){
            sig <- if(is.null(tmp$sig)) tmp$sig_y else tmp$sig
            Fx <- pnorm(q=init.dat[[h]]$y, mean= tmp$exp_val, sd=sig)
            osa.out[[h]]$mcmc <-  qnorm(Fx)
          }
          if(init.dat[[h]]$family == 400){
            grp.idx <- init.dat[[h]]$group + 1
            res <- NULL
            p <- tmp$power
            phi <- tmp$sig_y
            #calculate by group as ptweedie relies on pop min/max
            for(j in 1:ng){
              idx <- which(grp.idx == j)
              y <- init.dat[[h]]$y[idx]
              mu <- tmp$exp_val[idx]
              Fx <- ptweedie(q = y, mu = mu, 
                             phi = phi, power = p)
              zero.idx <- which(y == 0)
              #as implemented in statmod::qres.tweedie
              Fx[zero.idx] <- runif(length(zero.idx), 0, Fx[zero.idx])
              res <- c(res, qnorm(Fx))
            }
            osa.out[[h]]$mcmc <- res
          }
          if(init.dat[[h]]$family == 500){
            zero.idx <- which(init.dat[[h]]$y == 0)
            pos.idx <- which(init.dat[[h]]$y > 0)
            sig <- tmp$sig_y
            Fx <- rep(0,length(init.dat[[h]]$y))
            Fx[pos.idx] <- pgamma(init.dat[[h]]$y[pos.idx], shape = 1/sig^2, 
                                     scale = tmp$exp_val[pos.idx]*sig^2)
            Fx[zero.idx] <- runif(length(zero.idx), 0, 1)
            osa.out[[h]]$mcmc <- qnorm(Fx)
          }
        }
        osa.out[[h]]$runtime.mcmc <- as.numeric(Sys.time()-t0, 'secs')
        
        if('re_mcmc' %in% osa.methods){
          if(mod == "randomwalk"){
            sig <- tmp$tau
            mu <- tmp$ypred
            Fx <- pnorm(q = postmle, mean = mu, sd = sig) #needs rotation?
            osa.out[[h]]$re_mcmc <- qnorm(Fx)
          } 
          if(mod == "simpleGLMM"){
            u.idx <- grep('u', names(objmle$par))
            sig <- tmp$sig_u
            mu <- rep(0,length(u.idx))
            Fx <- pnorm(q = postmle[u.idx], mean = mu, sd = sig) 
            osa.out[[h]]$re_mcmc <- qnorm(Fx)
          }
          if(mod == "spatial"){
            omega.idx <- grep('omega', names(objmle$par))
            mu <- rep(0, length(postmle[omega.idx]))
            #rotation using spatial covariance matrix
            Sig <- solve(tmp$Q * exp(2 * FE$ln_tau))
            r1 <- as.vector(solve(t(chol(Sig)), postmle[omega.idx]))
            Fx <- pnorm(q = r1, mean = mu, sd = 1) 
            osa.out[[h]]$re_mcmc <- qnorm(Fx)
          }
          if(mod == "spatial"){#filter on omegas associated with obs
            omega.idx <- omega.idx[init.dat[[h]]$mesh_i+1]
            r2 <- r1[omega.idx]
            mu <- rep(0,length(omega.idx))
            Fx <- pnorm(q = r2, mean = mu, sd = 1)
            osa.out[[h]]$re_obs_mcmc <- qnorm(Fx)
            if(h==1) osa.methods <- c(osa.methods, 're_obs_mcmc')
          } else {
            osa.out[[h]]$re_obs_mcmc <- NA
          }
          
          if(misp == 'overdispersion' & mod != 'linmod' & h == 1){
            mu <- rep(0, n)
            sig <- tmp$sig_v
            Fx <- pnorm(q = postmle[grep('v', names(objmle$par))], 
                        mean = mu, sd = sig)
            osa.out[[h]]$re_mcmc_v <- qnorm(Fx)
            osa.methods <- c(osa.methods, 're_mcmc_v')
          } else {
            osa.out[[h]]$re_mcmc_v <- NA
          }
        }
        
      } else {
        osa.out[[h]]$mcmc <- osa.out[[h]]$runtime.mcmc <- NA
        osa.out[[h]]$re_mcmc <- NA
        osa.out[[h]]$re_obs_mcmc <- NA
      }
    
      AIC <- ifelse(do.true, NA, mod.out[[h]]$aic$AIC) #!doesn't work if doTrue == TRUE
      AICc <- ifelse(do.true, NA, mod.out[[h]]$aic$AICc) #!doesn't work if doTrue == TRUE
      maxgrad <- ifelse(do.true,NA, max(abs(mod.out[[h]]$obj$gr(mod.out[[h]]$opt$par))) ) #!doesn't work if doTrue == TRUE
      converge <- NA
      if(!do.true){
        if(mod.out[[h]]$opt$convergence == 1 | mod.out[[h]]$sdr$pdHess == FALSE ){
          converge <- 1
        } else {
          converge <- 0
        }
      }
      r[[h]] <- data.frame(id=id, model=mod, misp = misp, version=names(mod.out)[h],
                           replicate=ii, y=sim.dat[[h]],
                           ypred=mod.out[[h]]$report$exp_val,
                           osa.cdf = osa.out[[h]]$cdf, 
                           osa.gen = osa.out[[h]]$gen,
                           osa.fg=osa.out[[h]]$fg, 
                           osa.osg=osa.out[[h]]$osg,
                           osa.mcmc=osa.out[[h]]$mcmc, 
                           pears=osa.out[[h]]$pears,
                          # osa.re_mcmc=osa.out[[h]]$re_mcmc, #difficult to include b/c diff dimension from data
                           sim_cond= dharma.out[[h]]$cond$resids,
                           sim_uncond=dharma.out[[h]]$uncond$resids)#,
                           #sim_re_uncond=dharma.out[[h]]$re_uncond$resids) #difficult to include b/c diff dimension from data
      out[[h]] <- data.frame(id=id, model=mod, misp = misp, version=names(mod.out)[h],
                             replicate = ii,
                             runtime_cond=dharma.out[[h]]$cond$runtime,
                             runtime_uncond=dharma.out[[h]]$uncond$runtime,
                             # runtime_parcond=dharma.out[[h]]$parcond$runtime,
                             runtime.cdf=osa.out[[h]]$runtime.cdf,
                             runtime.fg=osa.out[[h]]$runtime.fg,
                             runtime.osg=osa.out[[h]]$runtime.osg,
                             runtime.gen=osa.out[[h]]$runtime.gen,
                             runtime.mcmc=osa.out[[h]]$runtime.mcmc,
                             maxgrad=maxgrad, converge=converge,
                             AIC=AIC, AICc=AICc)

      pvals <- rbind(pvals, cbind(id,calc.pvals( type = 'osa', method = osa.methods, mod = mod,
                                        res.obj = osa.out[[h]], version = names(mod.out)[h],
                                        fam = true.parms$fam, do.true )))
      pvals <- rbind(pvals, cbind(id,calc.pvals( type = 'sim', method = dharma.methods, mod = mod,
                                        res.obj = dharma.out[[h]], version = names(mod.out)[h],
                                        fam = true.parms$fam, do.true )))
      if(mod == 'spatial'){
        if(!is.null(osa.methods)){
          sac.pvals <- calc.sac( type = 'osa', 
                                 dat = sim.dat, 
                                 res.obj = osa.out[[h]],
                                 version = names(mod.out)[h])
          pvals <- rbind(pvals, cbind(id, sac.pvals))
        }    
      }
          
         
    }
  }
  resids <- rbind(r[[1]], r[[2]])
  stats <- rbind(out[[1]], out[[2]])
  pvals$replicate <- ii
  pvals$misp <- misp
  ## Tack this on for plotting later
  resids$do.true <- do.true
  pvals$do.true <- do.true
  stats$do.true <- do.true
  mles <- do.call(rbind, mles)
  if(savefiles){
    if(do.true) mod <- paste0(mod, "_true")
    dir.create(paste0('results/', mod, '_', misp, '_pvals'), showWarnings=FALSE)
    dir.create(paste0('results/', mod, '_', misp, '_resids'), showWarnings=FALSE)
    saveRDS(pvals, file=paste0('results/', mod, '_', misp, '_pvals/pvals_', ii, '.RDS'))
    saveRDS(resids, file=paste0('results/', mod, '_', misp, '_resids/resids_', ii, '.RDS'))
    if(length(mles)>0){
      dir.create(paste0('results/', mod, '_', misp, '_mles'), showWarnings=FALSE)
      saveRDS(mles, file=paste0('results/', mod, '_', misp, '_mles/mles_', ii, '.RDS'))
    }
    dir.create(paste0('results/', mod, '_', misp, '_stats'), showWarnings=FALSE)
    saveRDS(stats, file=paste0('results/', mod, '_', misp, '_stats/stats_', ii, '.RDS'))
  }


  #
  # if(ii==1 & savefiles){
  #   message("Making plots for replicate 1...")
  #   library(ggplot2)
    ## Need to generalize plots
    # if(mod == 'simpleGLMM'){
    #
    #   resids.long <- cbind(resids,x=1:length(dat0$y), group=c(dat0$group, dat1$group)) %>%
    #     pivot_longer(c('osa.cdf', 'osa.gen', 'osa.fg', 'osa.osg',
    #                    'sim_cond', 'sim_uncond', 'sim_parcond')) %>%
    #     filter(!is.na(value))
    #   theme_set(theme_bw())
    ## Plot of data
    # png(paste0('plots/',mod,'_simdata.png'), width=7, height=4, units='in', res=200)
    # par(mfrow=c(1,2))
    # boxplot(y~group, dat0, ylim=range(dat1$y))
    # boxplot(y~group, dat1, ylim=range(dat1$y))
    # dev.off()
    # } else {
    #   resids.long <- resids %>%
    #     pivot_longer(c('osa.cdf', 'osa.gen', 'osa.fg', 'osa.osg',
    #                    'sim_cond')) %>%
    #     filter(!is.na(value))
    # }

#
#
#     ## plot of resids
#     g <- ggplot(resids.long, aes(x, y=value, color=factor(group))) +
#       geom_point() + facet_grid(version~name) +
#       labs(x='Order', y='Residual')
#     ggsave('plots/simpleGLMM_resids_by_group.png', g, width=9, height=6)
#     g <- GGally::ggpairs(resids, columns=6:12, mapping=aes(color=version), title='Random Walk')
#     ggsave('plots/simpleGLMM_resids_pairs.png', g, width=7, height=5)
#     ## Plot of  DHARMa simulated data look like
#     ff <- function(x, v, re) data.frame(x=1:length(dat0$y), x$sims, group=dat0$group, version=v, method=re)
#     g <- rbind(ff(sim0_cond, 'm0', 'cond'),
#                ff(sim0_parcond, 'm0', 'parcond'),
#                ff(sim0_uncond, 'm0', 'uncond'),
#                ff(sim1_cond, 'm1', 'cond'),
#                ff(sim1_parcond, 'm1', 'parcond'),
#                ff(sim1_uncond, 'm1', 'uncond')) %>%
#       pivot_longer(cols=c(-x, -group, -version, -method), names_prefix="X",
#                    names_to='replicate', values_to='y') %>%
#       mutate(replicate=as.numeric(replicate)) %>%
#       ggplot(aes(group, y, fill=factor(group))) +
#       geom_violin() +
#       ## geom_jitter(alpha=.5, pch='.', width=.3, height=0) +
#       facet_grid(version~method)
#     g <- g+geom_jitter(alpha=.5, width=.2, data=rbind(data.frame(dat0), data.frame(dat1)))
#     ggsave('plots/simpleGLMM_simdata.png', g, width=9, height=6)
  return(invisible(list(pvals=pvals, resids=resids, mles=mles, stats=stats)))
}

mkTMBdat <- function(Dat, Pars, Mod, Misp){
  if(Mod == 'linmod'){
    dat0 <- list(y = Dat$y0, X = Dat$x)
    dat1 <- list(y = Dat$y1, X = Dat$x)
  }
  if(Mod == 'randomwalk'){
    dat0 <- list(y = Dat$y0, sim_re = 0)
    dat1 <- list(y = Dat$y1, sim_re = 0)
  }
  if(Mod == 'simpleGLMM'){
    ng <- length(Dat$random$u)
    ni <- length(Dat$y0)/ng
    dat0 <- list(y=Dat$y0, X = Dat$x,
                 group = rep(0:(ng-1), each = ni),
                 obs = rep(0:(ni-1), ng),
                 family = fam_enum(Pars$fam),
                 link = link_enum(Pars$link),
                 sim_re = 0)
    dat1 <- dat0
    dat1$y = Dat$y1
    if(Misp=="deltagamma"){
      dat1$family = fam_enum("Delta_Gamma")
    }
  }
  if(Mod == 'spatial'){
    loc <- Dat$loc
    mesh <- Dat$mesh
    dd <- as.matrix(dist(loc))
    dat0 <- list(y = Dat$y0, X = Dat$x, dd = dd, nu = 1,
                 mesh_i = mesh$idx$loc-1, sim_re = 0,
                 family = fam_enum(Pars$fam),
                 link = link_enum(Pars$link), reStruct = 10)
    dat0$spde <-  INLA::inla.spde2.matern(mesh)$param.inla[c('M0', 'M1', 'M2')]
    dat1 <- dat0
    dat1$y <- Dat$y1
  }
  if(Misp=='misscov'){
    dat1$X <- as.matrix(dat1$X[,1])
  }
  if(Misp=='dropRE'){
    dat1$reStruct = 20 #do not fit spatial likelihood
  }
  out = list(h0 = dat0, h1 = dat1)
  return(out)
}

mkTMBpar <- function(Pars, Dat, Mod, Misp, doTrue){

  if(Mod == 'linmod'){
    if(doTrue){
      par0 <-  list(beta = Pars$theta, ln_sig = log(Pars$sd.vec[1]))
    } else {
      par0 <- list(beta = c(0,0), ln_sig = 0)
    }
    par1 <- par0

    if(Misp == 'misscov'){
      par1$beta <- par1$beta[1]
    }
  }

  if(Mod == 'randomwalk'){
    if(doTrue){
      par0 <-  list(u = Dat$random$u, mu = Pars$theta, ln_sig = log(Pars$sd.vec[1]), ln_tau = log(Pars$sd.vec[2]))
    } else {
      par0 <- list(u = rep(1, length(Dat$random$u)), mu = 0, ln_sig = 0, ln_tau = 0)
    }
    par1 <- par0
    if(Misp == 'mu0'){
      par1$mu = 0
    }
  }

  if(Mod == 'simpleGLMM'){
    if(doTrue){
      par0 <- list(beta = Pars$theta, 
                   ln_sig_y = log(Pars$sd.vec[1]),
                   ln_sig_u = log(Pars$sd.vec[2]), ln_sig_v = numeric(0),
                   u = Dat$random$u, v = rep(0,length(Dat$y0)/length(Dat$random$u)))
      if(Pars$fam == 'Tweedie'){
        par0$ln_sig_y <- c(par0$ln_sig_y, 
                           log((Pars$sd.vec[3]-1)/(2-Pars$sd.vec[3])))
      }
    } else {
      par0 <- list(beta = 0, ln_sig_y = 0, ln_sig_u = 0, ln_sig_v = numeric(0),
                   u = rep(0, length(Dat$random$u)),
                   v = rep(0, length(Dat$y0)/length(Dat$random$u)))
      if(Pars$fam == 'Tweedie'){
        par0$ln_sig_y <- c(0,0)
      }
    }
    if(Pars$fam == "Poisson") par0$ln_sig_y = numeric(0)
    par1 <- par0
    if(Misp == 'overdispersion'){
      if(doTrue){
        par0$ln_sig_v <- log(Pars$sd.vec[3])
        par0$v <- Dat$random$v
      } else {
        par0$ln_sig_v <- 0
      }
    }
    if(Misp == 'misscov'){
      par1$beta <- par1$beta[1]
      if(!doTrue) par0$beta <- c(0,0)
    }
  }

  if(Mod == 'spatial'){
    if(doTrue){
      omega.true <- rep(0, Dat$mesh$n)
      omega.true[Dat$mesh$idx$loc] <- as.vector(Dat$random$omega)
      par0 <- list(beta = Pars$theta, theta = log(Pars$sd.vec[1]),
                   ln_tau = log(1/(2*sqrt(pi)*sqrt(8)/Pars$sp.parm*sd.vec[2])),
                   ln_kappa = log(sqrt(8)/Pars$sp.parm),
                   ln_sig_v = numeric(0),
                   omega = omega.true,
                   v = rep(0,length(Dat$y0)))
    } else {
      par0 <- list(beta = 0, theta = 0, ln_tau = 0, ln_kappa = 0,
                   ln_sig_v = numeric(0),
                   omega = rep(0, Dat$mesh$n),
                   v = rep(0,length(Dat$y0)))
      if(Misp == 'misscov'){
        par0$beta <- c(0,0)
      }
    }
    if(Pars$fam == "Poisson") par0$theta = 0
    par1 <- par0
    if(Misp == 'overdispersion'){
      if(doTrue){
        par0$ln_sig_v <- log(Pars$sd.vec[3])
        par0$v <- Dat$random$v
      } else {
        par0$ln_sig_v <- 0
      }
    }
    if(Misp == 'misscov'){
      par1$beta <- par1$beta[1]
    }
    if(Misp == 'dropRE'){
      par1$ln_tau = 0
      par1$ln_kappa = 0
    }
  }
  out = list(h0 = par0, h1 = par1)
  return(out)
}

mkTMBrandom <- function(Mod, Misp, doTrue){
  if(Mod == 'linmod'){
    Random.h0 = NULL
    Random.h1 = NULL
  }
  if(Mod == 'randomwalk'){
    Random.h0 <- 'u'
    Random.h1 <- 'u'
  }
  if(Mod == 'simpleGLMM'){
    Random.h1 <- 'u'
    if(Misp == 'overdispersion'){
      Random.h0 <- c('u', 'v')
    } else {
      Random.h0 <- 'u'
    }
  }
  if(Mod == 'spatial'){
    Random.h1 <- 'omega'
    if(Misp == 'overdispersion'){
      Random.h0 <- c('omega', 'v')
    } else {
      Random.h0 <- 'omega'
    }
    if(Misp == 'dropRE'){
      Random.h1 <- c()
    }
  }
  out <- list(h0 = Random.h0, h1 = Random.h1)
  return(out)
}

mkTMBmap <- function(Pars, Mod, Misp, Fam, doTrue){
  if(Mod == 'linmod'){
    map.h0 <- list()
    map.h1 <- list()
  }
  if(Mod == 'randomwalk'){
    map.h0 <- list()
    if(Misp == 'mu0'){
      map.h1 <- list(mu = factor(NA))
    } else {
      map.h1 <- list()
    }
  }
  if(Mod == 'simpleGLMM' | Mod == 'spatial'){
    map.h1 <- list(v = rep(factor(NA), length(Pars$h1$v)))
    if(Misp == 'overdispersion'){
      map.h0 <- list()
    } else {
      map.h0 <- list(v = rep(factor(NA), length(Pars$h0$v)))
    }
    if(Fam == "Poisson"){
      map.h0$theta <- factor(NA)
      map.h1$theta <- factor(NA)
    }
    if(Misp == 'dropRE'){
      map.h1$ln_kappa = factor(NA)
      map.h1$ln_tau = factor(NA)
      map.h1$omega <- rep(factor(NA), length(Pars$h1$omega))
    }
  }
  out <- list(h0 = map.h0, h1 = map.h1)
  return(out)
}


fam_enum <- function(fam){
  if(fam == 'Gaussian') out <- 000
  if(fam == 'Gamma') out <- 100
  if(fam == 'Poisson') out <- 200
  if(fam == 'lognormal') out <- 300
  if(fam == 'Tweedie') out <- 400
  if(fam == 'Delta_Gamma') out <- 500
  return(out)
}

link_enum <- function(link){
  if(link == 'log') out <- 0
  if(link == 'logit') out <- 1
  if(link == 'identity') out <- 2
  if(link == 'probit') out <- 3
  return(out)
}




