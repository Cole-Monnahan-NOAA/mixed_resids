setup_trueparms <- function(mod, misp){
  if(mod=='linmod'){
    theta <- c(4,-5)
    sd.vec <- 1
    sp.parm <- 0
    fam <- NULL
    link <- NULL
    if(misp == 'overdispersion'){
      #sd.vec <- c(sd.vec, sd.vec*log(4))
      sd.vec <- c(sd.vec, 1)
    }
  }
  if(mod=='randomwalk'){
    theta <- .75
    sd.vec <- c(1,1)
    sp.parm <- 0
    fam <- NULL
    link <- NULL
  }
  if(mod=='simpleGLMM'){
    #parms when fam = 'Gaussian';link='identity'
    theta <- 4
    sd.vec <- sqrt(c(.5,10))
    fam <- 'Gaussian'
    link <- 'identity'

    #parms when fam = 'Poisson'; link = 'log'
    # theta <- 1.5
    # sd.vec <- sqrt(c(.5,2))
    # fam <- 'Poisson'
    # link <- 'log'

    sp.parm <- 0
    if(misp=='misscov'){
      #parms when fam = 'Gaussian';link='identity
       theta <- c(4,-5)
     #  theta <- c(1.5,-2)
    }
    if(misp=='overdispersion') {
      sd.vec <- c(sd.vec, 1)
    }
  }
  if(mod=='spatial'){
    theta=2
    sd.vec <- c(.5,sqrt(0.5))
    sp.parm <- 20
    fam <- 'Poisson'
    link <- 'log'
    if(misp=='misscov'){
      theta <- c(1,2)
    }
    if(misp=='overdispersion'){
      sd.vec <- c(sd.vec, sd.vec[1]*2)
    }
  }
  true.pars <- list(theta=theta, sd.vec=sd.vec, sp.parm=sp.parm,
                    fam=fam, link=link)
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
  pvals <- data.frame(type = character(), method = character(),
                      test = character(), version = character(),
                      pvalue = numeric())
  mles <-  r <- out <- list()

  for(h in 1:2){

    message(ii, ": Optimizing two competing models...")
    init.obj <- list(data = init.dat[[h]], parameters = init.par[[h]], map = init.map[[h]], random = init.random[[h]], DLL = mod)
    mod.out[[h]] <- fit_tmb(obj.args = init.obj, control = list(run.model = !do.true, do.sdreport = TRUE))
    if(!do.true){
      ## if estimating, return MLE values
      tmp1 <- true.parms;
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
      mles[[h]] <- rbind(data.frame(h=h-1, type='true', par=names(tmp1), value=as.numeric(tmp1)),
                         data.frame(h=h-1, type='mle', par=names(tmp2), value=as.numeric(tmp2)))
      mles[[h]] <- cbind(mles[[h]], replicate=ii,
                         do.true=do.true,  model=mod, misp=misp)
    } else {
      ## otherwise just NULL b/c nothing estimated
      mles[[h]] <- NULL
    }

    message(ii, ": Calculating residuals..")
    disc <- FALSE; ran <- c(-Inf,Inf)
    if(!is.null(true.parms$fam)){
      if(true.parms$fam == 'Poisson'){
        disc <- TRUE
        ran <- c(0,Inf)
      }
      if(true.parms$fam == 'Gamma') ran <- c(0,Inf)
    }
    osa.out[[h]] <- calculate.osa(mod.out[[h]]$obj, methods=osa.methods, observation.name='y', Discrete = disc, Range = ran)

    expr <- expression(obj$simulate()$y)
    if('cond' %in% dharma.methods){
      dharma.out[[h]]$cond <- calculate.dharma(mod.out[[h]]$obj, expr, obs=sim.dat[[h]], fpr=mod.out[[h]]$report$fpr)
    }

    if('uncond' %in% dharma.methods){
      mod.out[[h]]$obj$env$data$sim_re <- 1 #turn on RE simulation
      dharma.out[[h]]$uncond <- calculate.dharma(mod.out[[h]]$obj, expr, obs=sim.dat[[h]], fpr=mod.out[[h]]$report$fpr)
    }

    ## only makes sense to run when MLE is estimated and there
    ## are no RE
    if('mcmc' %in% osa.methods & !do.true & Random ){
      t0 <- Sys.time()
      ## Build up TMB obj again
      FE <- mod.out[[h]]$opt$par # estimated FE
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
      if(mod=='spatial'){
        Fx <- ppois(init.dat[[h]]$y, tmp$exp_val)
        px <- dpois(init.dat[[h]]$y, tmp$exp_val)
        u <- runif(length(Fx))
        osa.out[[h]]$mcmc <- qnorm(Fx - u * px)
      } else {
        sig <- if(is.null(tmp$sig)) tmp$sig_y else tmp$sig
        Fx <- pnorm(q=init.dat[[h]]$y, mean= tmp$exp_val, sd=sig)
        osa.out[[h]]$mcmc <-  qnorm(Fx)
      }
      osa.out[[h]]$runtime.mcmc <- as.numeric(Sys.time()-t0, 'secs')
    } else {
      osa.out[[h]]$mcmc <- osa.out[[h]]$runtime.mcmc <- NA
    }

    AIC <- ifelse(do.true, NA, mod.out[[h]]$aic$AIC) #!doesn't work if doTrue == TRUE
    AICc <- ifelse(do.true, NA, mod.out[[h]]$aic$AICc) #!doesn't work if doTrue == TRUE
    maxgrad <- ifelse(do.true,NA, max(abs(mod.out[[h]]$obj$gr(mod.out[[h]]$opt$par))) ) #!doesn't work if doTrue == TRUE
    converge <- NA
    if(!do.true){
      if(mod.out[[h]]$opt$convergence == 0 & mod.out[[h]]$sdr$pdHess == TRUE){
        converge <- 0
      } else {
        converge <- 1
      }
    }
    r[[h]] <- data.frame(model=mod, replicate=ii, y=sim.dat[[h]],
                         ypred=mod.out[[h]]$report$exp_val, version=names(mod.out)[h],
                         osa.cdf = osa.out[[h]]$cdf, osa.gen = osa.out[[h]]$gen,
                         osa.fg=osa.out[[h]]$fg, osa.osg=osa.out[[h]]$osg,
                         osa.mcmc=osa.out[[h]]$mcmc,
                         sim_cond= dharma.out[[h]]$cond$resids,
                         sim_uncond=dharma.out[[h]]$uncond$resids)#,
                         # sim_parcond=dharma.out[[h]]$parcond$resids)
    out[[h]] <- data.frame(model=mod, replicate = ii, version=names(mod.out)[h],
                           runtime_cond=dharma.out[[h]]$cond$runtime,
                         runtime_uncond=dharma.out[[h]]$uncond$runtime,
                        # runtime_parcond=dharma.out[[h]]$parcond$runtime,
                         runtime.cdf=osa.out[[h]]$runtime.cdf,
                         runtime.fg=osa.out[[h]]$runtime.fg,
                         runtime.osg=osa.out[[h]]$runtime.osg,
                         runtime.gen=osa.out[[h]]$runtime.gen,
                         runtime.mcmc=osa.out[[h]]$runtime.mcmc,
                         maxgrad=maxgrad, converge=converge,AIC=AIC, AICc=AICc)

    pvals <- rbind(pvals, calc.pvals( type = 'osa', method = osa.methods, mod = mod,
                                      res.obj = osa.out[[h]], version = names(mod.out)[h],
                                      fam = true.parms$fam, do.true ))
    pvals <- rbind(pvals, calc.pvals( type = 'sim', method = dharma.methods, mod = mod,
                                      res.obj = dharma.out[[h]], version = names(mod.out)[h],
                                      fam = true.parms$fam, do.true ))
    if(mod == 'spatial'){
      dmat <- as.matrix(dist(sim.dat$loc, upper = TRUE))
      wt <- 1/dmat;  diag(wt) <- 0
      for(m in 1:length(osa.methods)){
        pvals <- rbind(pvals, data.frame(type='osa', method=osa.methods[m], model=mod, test='SAC', version = names(mod.out)[h],
                       pvalue = calc.sac(osa.out[[h]][[osa.methods[m]]], wt)))
      }
      for(m in 1:length(osa.methods)){
        pvals <- rbind(pvals, data.frame(type='sim', method=dharma.methods[m], model=mod, test='SAC', version = names(mod.out)[h],
                       pvalue = calc.sac(dharma.out[[h]][[dharma.methods[m]]]$resids, wt)))
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
  return(invisible(list(pvals=pvals, resids=resids, mles=mles)))
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
  }

  if(Mod == 'simpleGLMM'){
    if(doTrue){
      par0 <- list(beta = Pars$theta, ln_sig_y = log(Pars$sd.vec[1]),
                   ln_sig_u = log(Pars$sd.vec[2]), ln_sig_v = numeric(0),
                   u = Dat$random$u, v = rep(0,length(Dat$y0)/length(Dat$random$u)))
    } else {
      par0 <- list(beta = 0, ln_sig_y = 0, ln_sig_u = 0, ln_sig_v = numeric(0),
                   u = rep(0, length(Dat$random$u)),
                   v = rep(0, length(Dat$y0)/length(Dat$random$u)))
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
      par0 <- list(beta = Pars$theta, theta = log(Pars$sd.vec[1]),
                   ln_tau = log(1/(2*sqrt(pi)*sqrt(8)/Pars$sp.parm*sd.vec[2])),
                   ln_kappa = log(sqrt(8)/Pars$sp.parm),
                   ln_sig_v = numeric(0),
                   omega = Dat$random$omega,
                   v = rep(0,length(Dat$y0)))
    } else {
      par0 <- list(beta = 0, theta = 0, ln_tau = 0, ln_kappa = 0,
                   ln_sig_v = numeric(0),
                   omega = rep(0, length(Dat$random$omega)),
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
  if(fam == 'Binomial') out <- 500
  return(out)
}

link_enum <- function(link){
  if(link == 'log') out <- 0
  if(link == 'logit') out <- 1
  if(link == 'identity') out <- 2
  if(link == 'probit') out <- 3
  return(out)
}




