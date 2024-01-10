setup_trueparms <- function(mod, misp, fam, link, type){
  true.comp <- list()

  if(mod == 'linmod'){
    true.pars <- setup_linmod(mod, misp, fam, link)
  }

  if(mod == 'randomwalk'){
    true.pars <- setup_randomwalk(mod, misp, fam, link, type)
  }

  if(mod == 'simpleGLMM'){
    true.pars <- setup_simpleGLMM(mod, misp, fam, link, type)
  }

  if(mod == 'spatial'){
    true.pars <- setup_spatial(mod, misp, fam, link, type)
  } 
  return(true.pars)
}

setup_linmod <- function(mod, misp, fam, link){
  beta <- c(4,-5)
  sd.vec <- 1
  true.comp <- list()
  true.comp[[1]] <- true.comp[[2]] <-
    list(beta_1 = beta[1], beta_2 = beta[2],
      ln_sig = log(sd.vec))
  
  if(misp == 'overdispersion'){
    sd.vec <- c(sd.vec, 1)
  }
  if(misp == 'misscov'){
    true.comp[[2]] <- list(beta = beta[1],
                           ln_sig = log(sd.vec[1]))
  }
  true.parms <- list(beta=beta, sd.vec=sd.vec, 
    fam=fam, link=link, true.comp=true.comp)
  return(true.parms)
}

setup_randomwalk <- function(mod, misp, fam, link, type){
  if(type == "LMM"){
    beta <- 6
    drift <- 2
    sd.vec <- c(1,1)
  }
  if(type == "GLMM"){
    beta <- 0.1
    drift <- 0.05
    sd.vec <- c(0.5, 0.05)
  }
  
  true.comp <- list()
  true.comp[[1]] <- list(mu = beta, ln_sig_y = log(sd.vec[1]),
          ln_sig_u = log(sd.vec[2]))
  for(i in 1:length(misp)){
    if(misp[i] == 'mu0'){
      true.comp[[i+1]] <- list(ln_sig_y = log(sd.vec[1]),
                               ln_sig_u = log(sd.vec[2]))
    } 
    if(misp[i] == 'missre'){
      true.comp[[i+1]] <- list(mu = beta, ln_sig_y = log(sd.vec[1]))
    }
    if(misp[i] == 'gamma-lognorm' | misp[i] == 'normal-lognorm'){
      true.comp[[i+1]] <- true.comp[[1]]
    }
  }

  true.parms <- list(beta=beta, sd.vec=sd.vec, drift = drift,
    fam=fam, link=link, true.comp=true.comp)

  return(true.parms)
}

setup_simpleGLMM <- function(mod, misp, fam, link, type){
  theta <- NA
  true.comp <- list()
  if(type == "LMM"){
    #currently only misp = misscovunif or misscovnorm is implemented
    beta <- c(4,-8)
    sd.vec <- c(0.5, 2)
    true.comp[[1]] <- list(beta_1 = beta[1], beta_2 = beta[2],
                           ln_sig_y = log(sd.vec[1]),
                           ln_sig_u = log(sd.vec[2]))
   
  }

  if(type == "GLMM"){
    sd.vec <- rep(NA, 2)

    if(fam == "NB"){
      beta <- log(2)
      size <- 1
      sd.vec[2] <- 1 #sig_u
      theta <- size

      true.comp[[1]] <- list(beta = beta,
                             theta = log(size),
                             ln_sig_u = log(sd.vec[2]))

    }

    if(fam == "Poisson"){
      beta <- 0.5
      sd.vec[2] <- sqrt(0.5) #sig_u

       true.comp[[1]] <- list(beta = beta[1],
                              ln_sig_u = log(sd.vec[2]))
    }

    if(fam == "Tweedie"){
      beta <- 1.5
      pow <- 1.3
      sd.vec <- c(1.7,1) #phi, sig_u
      theta <- pow

      true.comp[[1]] <- list(beta = beta,
                             ln_sig_y = log(sd.vec[1]),
                             theta = log((pow-1)/(2-pow)),
                             ln_sig_u = log(sd.vec[2]))
    }
  }

  for(i in 1:length(misp)){
  true.comp[[i+1]] <- unlist(true.comp[1])
    if(misp[i] == 'missunifcov' | misp[i] == 'missnormcov'){
      true.comp[[i+1]]$beta_1 <- beta[1]
      true.comp[[i+1]]$beta_2 <- NULL
      names(true.comp[[i+1]])[1] = "beta"
    }
    if(misp[i] == 'missre'){
      true.comp[[i+1]]$ln_sig_u <- NULL
    }
    if(misp[i] == 'nb-pois'){
      true.comp[[i+1]]$theta <- NULL
    }
  }
  true.parms <- list(beta = beta, sd.vec=sd.vec, theta = theta,
    fam=fam, link=link, true.comp=true.comp)

  return(true.parms)
}


setup_spatial<- function(mod, misp, fam, link, type){
  sp.parm <- 50
  true.comp <- list()
  if(type == "LMM"){
    beta <- 20
    sd.vec <- c(1, 1)

    true.comp[[1]] <- list(beta = beta, theta = log(sd.vec[1]),
           ln_tau = log(1/(2*sqrt(pi)*sqrt(8)/sp.parm*sd.vec[2])), #1/(2*sqrt(pi)*kappa*sp.sd))
           ln_kappa = log(sqrt(8)/sp.parm))
  }

  if(type == "GLMM"){
    if(fam == "Poisson"){
      beta <- 0.5
      sd.vec <- c(NA, sqrt(0.25))  
    }
    true.comp[[1]] <- list(beta = beta, 
           ln_tau = log(1/(2*sqrt(pi)*sqrt(8)/sp.parm*sd.vec[2])), #1/(2*sqrt(pi)*kappa*sp.sd))
           ln_kappa = log(sqrt(8)/sp.parm))
  }

  for(i in 1:length(misp)){
    true.comp[[i+1]] <- true.comp[1]
    if(misp[i] == "missre"){
      true.comp[[i+1]]$ln_tau <- NULL
      true.comp[[i+1]]$ln_kappa <- NULL
    }
  }

  
  true.parms <- list(beta=beta, sd.vec=sd.vec, sp.parm=sp.parm,
                    fam=fam, link=link, true.comp=true.comp)
  return(true.parms)
}

run_iter <- function(ii, n=100, ng=0, mod, cov.mod = NULL, misp, type,
                     family, link, do.true = FALSE, savefiles=TRUE){
  library(TMB)
  library(DHARMa)
  library(fmesher)
  library(dplyr)
  library(tidyr)
  library(R.utils)
  library(goftest)
  library(tweedie)
  
  if(do.true){
    mod.name <- paste0(mod, "_true")
  } else {
    mod.name <- mod
  }
  res.name <- paste0('results/', mod.name, '_', type)

  if(mod == 'linmod'){
    Random <- FALSE
  } else {
    Random <- TRUE
  }
  
  # if(misp == "missnormcov"){
  #   cov.mod <- "norm"
  #   misp <- "misscov"
  # }
  # if(misp == "missunifcov"){
  #   cov.mod <- "unif"
  #   misp <- "misscov"
  # }

  setupTMB(mod)
  true.parms <- setup_trueparms(mod, misp, family, link, type)

  ## simulate data with these parameters
  message(ii, ": Simulating data...")
  sim.dat <- simdat(n, ng, mod, cov.mod, type, true.parms, misp, ii)
  if(is.null(sim.dat)){
    set.seed(ii)
    while(is.null(sim.dat)){
      new.seed = round(runif(1, 5000, 8000))
      sim.dat <-  simdat(n, ng, mod, cov.mod, type, true.parms, misp, new.seed)
    }
  }

  init.dat <- mkTMBdat(sim.dat, true.parms, mod, misp, type)
  init.par <- mkTMBpar(true.parms, sim.dat, mod, misp, type, do.true)
  init.random <- mkTMBrandom(mod, misp)
  init.map <- mkTMBmap(init.par, mod, misp, type)
  mod.out <- osa.out <- dharma.out <- list()
  pvals <- data.frame(id = character(), type = character(), misp = character(),
                      res.type = character(), method = character(),
                      model = character(), test = character(), 
                      version = character(), pvalue = numeric())
  mles <-  r <- out <- list()

  for(h in 1:(length(misp)+1)){
    message(ii, ": Optimizing  models...")
    if(h == 1){
      misp.name <- "correct"
      id <- paste0(mod, '_', do.true, '_', type, '_correct_h0_', ii)
    } else {
      id <- paste0(mod, '_', do.true, '_', type, '_', misp[h-1], '_h', h-1, '_', ii)
      misp.name <- misp[h-1]
    }

    # Set the distribution based on the true or mis-specified likelihood
    mod.fam <- "Gaussian"
    if(!is.null(true.parms$fam)){
      mod.fam <- true.parms$fam
    }
    if(h > 1){
      if(misp[h-1] == "nb-pois"){
        mod.fam <- "Poisson"
      }
      if(misp[h-1] == "gamma-lognorm" | misp[h-1] == "normal-lognorm"){
        mod.fam <- "Lognormal"
      }
      if(misp[h-1] == "norm-gamma"){
        mod.fam <- "Gamma"
      }
    }

    if(h == 1 | mod == "linmod"){
      init.obj <- list(data = init.dat[[h]], parameters = init.par[[h]], 
                      map = init.map[[h]], random = init.random[[h]], DLL = mod)
      if(is.null(init.random[[h]])) Random <- FALSE
    } else {
      init.obj <- list(data = init.dat$h1[[h-1]], parameters = init.par$h1[[h-1]], 
                      map = init.map$h1[[h-1]], random = init.random$h1[[h-1]], DLL = mod)
      if(is.null(init.random$h1[[h-1]])) Random <- FALSE
    }
    
    if(Random) init.obj$hessian <- TRUE
    
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
                           do.true=do.true,  model=mod, misp=misp.name,
                           n=n)
      } else {
        ## otherwise just NULL b/c nothing estimated
        mles[[h]] <- NULL
      }

      message(ii, ": Calculating residuals..")
      disc <- FALSE; ran <- c(-Inf,Inf); rot <- NULL
      if(mod.fam == 'Poisson' | mod.fam == 'NB' | mod.fam == "Tweedie"){
          disc <- TRUE
          ran <- c(0,Inf)
      }
      if(mod.fam == 'Gamma' | mod.fam == 'Lognormal'){
        ran <- c(0,Inf)
      } 
        
      if(mod == 'randomwalk' | mod == 'spatial' | mod == 'simpleGLMM') rot <- "estimated"
      
      osa.out[[h]] <- calculate.osa(mod.out[[h]]$obj, methods=osa.methods, 
                                    observation.name='y', Discrete = disc, 
                                    Range = ran)

      expr <- expression(obj$simulate()$y)
      
      if(mod == 'simpleGLMM'){
        n_ <- n*ng
      } else {
        n_ <- n
      }

      if(h == 1 | mod == "linmod"){
        init.obs <- sim.dat[[h]]
      } else {
        init.obs <- sim.dat$y1[[h-1]]
      }
      dharma.out[[h]] <- list()
      if('cond' %in% dharma.methods){
        dharma.out[[h]]$cond <- 
          calculate.dharma(mod.out[[h]]$obj, expr, obs = init.obs,
                           idx = 1:n_, fpr=mod.out[[h]]$report$fpr, 
                           int.resp = disc, rot = rot)
      } else {
        dharma.out[[h]]$cond <- list()
        dharma.out[[h]]$cond$resids <- NA
        dharma.out[[h]]$cond$runtime <- NA
      }
      
      if('cond_nrot' %in% dharma.methods){
        dharma.out[[h]]$cond_nrot <- 
          calculate.dharma(mod.out[[h]]$obj, expr, obs = init.obs,
                           idx = 1:n_, fpr=mod.out[[h]]$report$fpr, 
                           int.resp = disc, rot = NULL)
      } else {
        dharma.out[[h]]$cond_nrot <- list()
        dharma.out[[h]]$cond_nrot$resids <- NA
        dharma.out[[h]]$cond_nrot$runtime <- NA
      }
      
      if('uncond_nrot' %in% dharma.methods | 
         'uncond' %in% dharma.methods){
        mod.out[[h]]$obj$env$data$sim_re <- 1 #turn on RE simulation
        #retape, do not reset parameters to initial values
        mod.out[[h]]$obj$retape(set.defaults = FALSE) 
      }
      
      if('uncond_nrot' %in% dharma.methods){
        dharma.out[[h]]$uncond_nrot <- 
          calculate.dharma(mod.out[[h]]$obj, expr, obs = init.obs, 
                           idx = 1:n_, fpr=mod.out[[h]]$report$fpr, 
                           int.resp = disc, rot = NULL)
      } else {
        dharma.out[[h]]$uncond_nrot <- list()
        dharma.out[[h]]$uncond_nrot$resids <- NA
        dharma.out[[h]]$uncond_nrot$runtime <- NA
      }

      if('uncond' %in% dharma.methods){
        dharma.out[[h]]$uncond <- 
          calculate.dharma(mod.out[[h]]$obj, expr, obs = init.obs, 
                           idx = 1:n_, fpr=mod.out[[h]]$report$fpr, 
                           int.resp = disc, rot = rot)
      } else {
        dharma.out[[h]]$uncond <- list()
        dharma.out[[h]]$uncond$resids <- NA
        dharma.out[[h]]$uncond$runtime <- NA
      }
      
      ## only makes sense to run when MLE is estimated and there
      ## are RE - turning on when do.true == TRUE (AMH, 6/3/2022)
      if('mcmc' %in% osa.methods){
        t0 <- Sys.time()
        if(!(misp.name == "missre")){
          ## Build up TMB obj again
          if(do.true)  FE <- mod.out[[h]]$obj$par # true FE
          if(!do.true) FE <- mod.out[[h]]$opt$par # estimated FE
          ## make into list; https://stackoverflow.com/questions/46251725/convert-named-vector-to-list-in-r/46251794
          FE <- split(unname(FE),names(FE))
          if(h == 1){
            MLE <- modifyList(init.par[[h]], FE) #
            ## Get FE and map them off
            xx <- names(MLE)[-which(names(MLE) %in% init.random[[h]])]
            map <- lapply(names(FE), function(x) factor(FE[[x]]*NA))
            names(map) <- names(FE)
            map <- c(map, init.map[[h]])
            ## Rebuild with original FE mapped off and RE as FE
            objmle <- MakeADFun(data=init.dat[[h]], parameters=MLE,
                                map=map, DLL=mod.out[[h]]$obj$env$DLL)
          } else {
            MLE <- modifyList(init.par$h1[[h-1]], FE) #
            ## Get FE and map them off
            xx <- names(MLE)[-which(names(MLE) %in% init.random$h1[[h-1]])]
            map <- lapply(names(FE), function(x) factor(FE[[x]]*NA))
            names(map) <- names(FE)
            map <- c(map, init.map$h1[[h-1]])
            objmle <- MakeADFun(data=init.dat$h1[[h-1]], parameters=MLE,
                                map=map, DLL=mod.out[[h]]$obj$env$DLL)
          }
          
          fitmle <- tmbstan::tmbstan(objmle, chains=1, warmup=300, iter=301, seed=ii, refresh=-1)
          postmle <- as.numeric(as.matrix(fitmle)) ## single sample
          postmle <- postmle[-length(postmle)] # drop lp__ value
          ## Calculate residuals given the sample
          tmp <- objmle$report(postmle)
        }
        #if no random effect in model mispecification, return the quantile residual
        if(misp.name == "missre"){ 
          tmp <- mod.out[[h]]$report
        }
        # Calculate quantile residual
        if(h == 1){
          osa.out[[h]]$mcmc <- calc.quantile(mod.fam, tmp, init.dat[[h]]$y)
        } else {
          osa.out[[h]]$mcmc <- calc.quantile(mod.fam, tmp, init.dat$h1[[h-1]]$y)
        }
        osa.out[[h]]$runtime.mcmc <- as.numeric(Sys.time()-t0, 'secs')
      
      } else {
        osa.out[[h]]$mcmc <- osa.out[[h]]$runtime.mcmc <- NA
      }
    
      AIC <- ifelse(do.true, NA, mod.out[[h]]$aic$AIC) #!doesn't work if doTrue == TRUE
      AICc <- ifelse(do.true, NA, mod.out[[h]]$aic$AICc) #!doesn't work if doTrue == TRUE
      maxgrad <- ifelse(do.true, NA, 
                        max(abs(mod.out[[h]]$obj$gr(mod.out[[h]]$opt$par))) 
                        ) #!doesn't work if doTrue == TRUE
      convergestatus <- ifelse(do.true, NA, 
                               mod.out[[h]]$opt$convergence)
      convergehessian <- ifelse(do.true, NA, 
                                mod.out[[h]]$sdr$pdHess)
      
      r[[h]] <- data.frame(id=id, model=mod, type = type, misp = misp.name, 
                           version = paste0("h", h-1),
                           replicate=ii, y=init.obs,
                           ypred=mod.out[[h]]$report$exp_val,
                           osa.cdf = osa.out[[h]]$cdf, 
                           osa.gen = osa.out[[h]]$gen,
                           osa.fg=osa.out[[h]]$fg, 
                           osa.osg=osa.out[[h]]$osg,
                           osa.mcmc=osa.out[[h]]$mcmc, 
                           pears=osa.out[[h]]$pears,
                           sim_cond_rot= dharma.out[[h]]$cond$resids,
                           sim_uncond_rot=dharma.out[[h]]$uncond$resids,
                           sim_cond_nrot= dharma.out[[h]]$cond_nrot$resids,
                           sim_uncond_nrot=dharma.out[[h]]$uncond_nrot$resids) 
     
      out[[h]] <- data.frame(id=id, model=mod, type = type, misp = misp.name, 
                             version = paste0("h", h-1), 
                             replicate = ii, 
                             converge.maxgrad = maxgrad, 
                             converge.status = convergestatus,
                             converge.hessian = convergehessian,
                             AIC = AIC, 
                             AICc = AICc)
      out[[h]][names(osa.out[[h]])[grep("runtime", names(osa.out[[h]]))]] <- 
        sapply(grep("runtime", names(osa.out[[h]])), function(x) osa.out[[h]][[x]] )
      out[[h]][paste0("runtime.", names(dharma.out[[1]])[!grepl("re", names(dharma.out[[1]]))])] <- 
        sapply(names(dharma.out[[1]])[!grepl("re", names(dharma.out[[1]]))], function(x) dharma.out[[h]][[x]]$runtime)
        
    
      pvals <- rbind(pvals, cbind(id = id, type = type, misp = misp.name,
                                  calc.pvals( res.type = 'osa', method = osa.methods, mod = mod,
                                              res.obj = osa.out[[h]], version = paste0("h", h-1),
                                              fam = true.parms$fam, do.true )))
      pvals <- rbind(pvals, cbind(id = id, type = type, misp = misp.name, 
                                  calc.pvals( res.type = 'sim', method = dharma.methods, mod = mod,
                                              res.obj = dharma.out[[h]], version = paste0("h", h-1),
                                              fam = true.parms$fam, do.true )))
      if(mod == 'spatial'){
        if(!is.null(osa.methods)){
          sac.pvals <- calc.sac( res.type = 'osa', 
                                 dat = sim.dat, 
                                 res.obj = osa.out[[h]],
                                 version = paste0("h", h-1))
          pvals <- rbind(pvals, cbind(id = id, type = type, 
                                      misp = misp.name,  sac.pvals))
        } 
        if(!is.null(dharma.methods)){
          sac.pvals <- calc.sac( res.type = 'sim', 
                                 dat = sim.dat, 
                                 res.obj = dharma.out[[h]],
                                 version = paste0("h", h-1))
          pvals <- rbind(pvals, cbind(id = id, type = type, 
                                      misp = misp.name,  sac.pvals))
        } 
      }
          
         
    }
  }
  resids <- rbind(r[[1]], r[[2]])
  stats <- rbind(out[[1]], out[[2]])
  pvals$replicate <- ii
  ## Tack this on for plotting later
  resids$do.true <- do.true
  pvals$do.true <- do.true
  stats$do.true <- do.true
  mles <- do.call(rbind, mles)
  if(savefiles){
    dir.create(paste0(res.name, '_pvals'), showWarnings=FALSE)
    dir.create(paste0(res.name, '_resids'), showWarnings=FALSE)
    saveRDS(pvals, file=paste0(res.name, '_pvals/pvals_', ii, '.RDS'))
    saveRDS(resids, file=paste0(res.name, '_resids/resids_', ii, '.RDS'))
    if(length(mles)>0){
      dir.create(paste0(res.name, '_mles'), showWarnings=FALSE)
      saveRDS(mles, file=paste0(res.name, '_mles/mles_', ii, '.RDS'))
    }
    dir.create(paste0(res.name, '_stats'), showWarnings=FALSE)
    saveRDS(stats, file=paste0(res.name, '_stats/stats_', ii, '.RDS'))
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

mkTMBdat <- function(Dat, Pars, Mod, Misp, Type){
  if(Mod == 'linmod'){
    dat0 <- list(y = Dat$y0, X = as.matrix(Dat$x))
    dat1 <- list(y = Dat$y1, X = as.matrix(Dat$x))
  }
  if(Mod == 'randomwalk'){
    dat0 <- list(y = Dat$y0, mod = 0, sim_re = 0) #default: type = LMM
    if(Type == "GLMM"){
      dat0$mod <- 2
    }
    dat1 <- list()
    if(length(Misp) != length(Dat$y1)){
      "Stop, length of mispecifications does not match list of misspecified data"
    }
    for(m in 1:length(Misp)){
      dat1[[m]] <- dat0
      dat1[[m]]$y <- Dat$y1[[m]]
      if(Misp[m] == 'normal-lognormal' | Misp[m] == 'gamma-lognorm'){
        dat1[[m]]$mod <- 1 # Lognormal
      }
    }
    
  }
  if(Mod == 'simpleGLMM'){
    ng <- length(Dat$random$u0)
    ni <- length(Dat$y0)/ng
    dat0 <- list(y=Dat$y0, X = as.matrix(Dat$x),
                 group = rep(0:(ng-1), each = ni),
                 obs = rep(0:(ni-1), ng),
                 family = fam_enum(Pars$fam),
                 link = link_enum(Pars$link),
                 sim_re = 0)
    dat1 <- list()
    if(length(Misp) != length(Dat$y1)){
      "Stop, length of mispecifications does not match list of misspecified data"
    }
    for(m in 1:length(Misp)){
      dat1[[m]] <- dat0
      dat1[[m]]$y <- Dat$y1[[m]]
      if(Misp[m] == "missunifcov" | Misp[m] == "missnormcov"){
        dat1[[m]]$X <- as.matrix(Dat$x[,1])
      }
      if(Misp[m] == "nb-pois"){
        dat1[[m]]$family <- fam_enum("Poisson")
      }
    }
  }
  if(Mod == 'spatial'){
    loc <- Dat$loc
    mesh <- Dat$mesh
    dd <- as.matrix(dist(loc))
    dat0 <- list(y = Dat$y0, X = as.matrix(Dat$x), dd = dd, nu = 1,
                 mesh_i = mesh$idx$loc-1, sim_re = 0,
                 family = fam_enum(Pars$fam),
                 link = link_enum(Pars$link), reStruct = 10)
    spde <- fmesher::fm_fem(mesh, order = 2)
    dat0$spde <-  list(M0 = spde$c0,
                       M1 = spde$g1,
                       M2 = spde$g2)
    dat1 <- list()
    if(length(Misp) != length(Dat$y1)){
      "Stop, length of mispecifications does not match list of misspecified data"
    }
    for(m in 1:length(Misp)){
      dat1[[m]] <- dat0
      
      dat1[[m]]$y <- Dat$y1[[m]]
      if(Misp[m] == "missre" ){
        dat1[[m]]$reStruct <- 20 #do not fit a spatial likelihood
      }
      if(Misp[m] == "norm-gamma"){
        dat1[[m]]$family <- fam_enum("Gamma")
      }
    }
  }
  
  out = list(h0 = dat0, h1 = dat1)
  return(out)
}

mkTMBpar <- function(Pars, Dat, Mod, Misp, Type, doTrue){
  if(Mod == 'linmod'){
    out <- mkTMBpar_linmod(Pars, Dat, Mod, Misp, Type, doTrue)
  }
  if(Mod == 'randomwalk'){
    out <- mkTMBpar_randomwalk(Pars, Dat, Mod, Misp, Type, doTrue)
  }
  if(Mod == 'simpleGLMM'){
    out <- mkTMBpar_simpleGLMM(Pars, Dat, Mod, Misp, Type, doTrue)
  }
  if(Mod == 'spatial'){
    out <- mkTMBpar_spatial(Pars, Dat, Mod, Misp, Type, doTrue)
  }
  return(out)
}

mkTMBpar_linmod <- function(Pars, Dat, Mod, Misp, Type, doTrue){
  if(doTrue){
    par0 <-  list(beta = Pars$beta, ln_sig_y = log(Pars$sd.vec[1]))
  } else {
    par0 <- list(beta = c(0,0), ln_sig_y = 0)
  }
  par1 <- par0

  if(Misp == 'misscov'){
    par1$beta <- par1$beta[1]
  }
  
  out = list(h0 = par0, h1 = par1)
  return(out)
}

mkTMBpar_randomwalk <- function(Pars, Dat, Mod, Misp, Type, doTrue){
  if(doTrue){
    par0 <-  list(beta = Pars$beta, mu = Pars$drift, ln_sig_y = log(Pars$sd.vec[1]), 
                  ln_sig_u = log(Pars$sd.vec[2]), u = Dat$random$u0)
  } else {
    par0 <-  list(beta = 0, mu = 0, ln_sig_y = 0, ln_sig_u = 0, 
                  u = rep(1, length(Dat$random$u0)))
  }
  par1 <- list()
  for(m in 1:length(Misp)){
    par1[[m]] <- par0
    if(Misp[m] == "missre"){
      par1[[m]]$ln_sig_u = numeric(0)
      par1[[m]]$u = rep(0, length(Dat$random$u0))
    }
    if(Misp[m] == "mu0"){
      par1[[m]]$mu = 0
      if(doTrue){
        par1[[m]]$u = Dat$random$u1[[m]]
      }
    }
  }

  out = list(h0 = par0, h1 = par1)
  return(out)
}

mkTMBpar_simpleGLMM <- function(Pars, Dat, Mod, Misp, Type, doTrue){
  if(doTrue){
    par0 <- list(beta = Pars$beta, 
                 ln_sig_y = log(Pars$sd.vec[1]),
                 theta = 0,
                 ln_sig_u = log(Pars$sd.vec[2]), 
                 u = Dat$random$u0)
    if(Pars$fam == 'Tweedie'){
      par0$theta <- log((Pars$theta-1)/(2-Pars$theta))
    }
    if(Pars$fam == 'NB'){
      par0$ln_sig_y = 0
      par0$theta <- log(Pars$theta)
    }
  } else {
    par0 <- list(beta = rep(0,length(Pars$beta)), ln_sig_y = 0, 
                 theta = 0, ln_sig_u = 0,
                 u = rep(0, length(Dat$random$u0)))
    if(Pars$fam == 'Tweedie' | Pars$fam == 'NB'){
      par0$theta <- 0
    }
  }

  par1 <- list()
  for(m in 1:length(Misp)){
    par1[[m]] <- par0
    if(Misp[m] == "missunifcov" | Misp[m] == "missnormcov"){
      par1[[m]]$beta <- par1[[m]]$beta[1]
    }
    if(Misp[m] == "missre"){
      par1[[m]]$ln_sig_u <- numeric(0)
      par1[[m]]$u <- rep(0, length(Dat$random$u0))
    }
    if(Misp[m] == "nb-pois"){
      par1[[m]]$theta <- 0
    }
    if(doTrue & Misp[m] == "mispre"){
      par1[[m]]$u <- Dat$random$u1[[m]]
    }
  }
  
  out = list(h0 = par0, h1 = par1)
  return(out)
}

mkTMBpar_spatial <- function(Pars, Dat, Mod, Misp, Type, doTrue){
  if(doTrue){
    omega.true <- rep(0, Dat$mesh$n)
    omega.true[Dat$mesh$idx$loc] <- as.vector(Dat$random$omega0)
    par0 <- list(beta = Pars$beta, theta = log(Pars$sd.vec[1]),
                   ln_tau = log(1/(2*sqrt(pi)*sqrt(8)/Pars$sp.parm*Pars$sd.vec[2])),
                   ln_kappa = log(sqrt(8)/Pars$sp.parm),
                   omega = omega.true)
  } else {
    par0 <- list(beta = 0, theta = 0, ln_tau = 0, ln_kappa = 0,
                   omega = rep(0, Dat$mesh$n))
     
  }
  if(Pars$fam == "Poisson") par0$theta = 0
  
  par1 <- list()
  for(m in 1:length(Misp)){
    par1[[m]] <- par0
    if(Misp[m] == "missre"){
      par1[[m]]$ln_tau = 0
      par1[[m]]$ln_kappa = 0
      par1[[m]]$omega = rep(0, Dat$mesh$n)
    }
  }

  out = list(h0 = par0, h1 = par1)
  return(out)
}

mkTMBrandom <- function(Mod, Misp){
  if(Mod == 'linmod'){
    Random.h0 = NULL
    Random.h1 = NULL
  }
  if(Mod == 'randomwalk' | Mod == "simpleGLMM"){
    Random.h0 <- 'u'
    Random.h1 <- list()
    for(m in 1:length(Misp)){
      Random.h1[[m]] <- 'u'
    }
  }
  if(Mod == 'spatial'){
    Random.h0 <- 'omega'
    Random.h1 <- list()
    for(m in 1:length(Misp)){
      Random.h1[[m]] <- 'omega'
    }
  }
  if("missre" %in% Misp){
    Random.h1[which(Misp == "missre")] <- list(NULL)
  }
  out <- list(h0 = Random.h0, h1 = Random.h1)
  return(out)
}

mkTMBmap <- function(Pars, Mod, Misp, Type){
  map.h0 <- list()
  map.h1 <- list()
  if(Mod == 'randomwalk'){
    for(m in 1:length(Misp)){
      map.h1[[m]] <- list()
      if(Misp[m] == "missre"){
        map.h1[[m]]$u = rep(factor(NA), length(Pars$h0$u))
      }
      if(Misp[m] == "mu0"){
        map.h1[[m]]$mu = factor(NA)
      }
    }
  }
  if(Mod == 'simpleGLMM'){
    if(Type == "LMM"){
      map.h0$theta = factor(NA)
    }
    if(Type == "GLMM"){
      map.h0$ln_sig_y = factor(NA)
    }
    for(m in 1:length(Misp)){
      map.h1[[m]] <- list()
      if(Type == "LMM"){
        map.h1[[m]]$theta = factor(NA)
      }
      if(Type == "GLMM"){
        map.h1[[m]]$ln_sig_y = factor(NA)
        if(Misp[[m]] == "nb-pois"){
          map.h1[[m]]$theta = factor(NA)
        }
      }
      
      if(Misp[m] == "missre"){
        map.h1[[m]]$u = rep(factor(NA), length(Pars$h0$u))
      }
    }
  } 
  if(Mod == 'spatial'){
    if(Type == "GLMM"){
      map.h0 <- list(theta = factor(NA))
    } 
    for(m in 1:length(Misp)){
      map.h1[[m]] <- map.h0
      if(Misp[m] == "missre"){
        map.h1[[m]]$ln_tau = factor(NA)
        map.h1[[m]]$ln_kappa = factor(NA)
        map.h1[[m]]$omega = rep(factor(NA), length(Pars$h0$omega))
      }
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
  if(fam == 'NB') out <- 600
  return(out)
}

link_enum <- function(link){
  if(link == 'log') out <- 0
  if(link == 'logit') out <- 1
  if(link == 'identity') out <- 2
  if(link == 'probit') out <- 3
  return(out)
}




