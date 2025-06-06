#' Function used to simulate data for each case study
#'
#' @param n sample size
#' @param ng group size, specific to simpleGLMM case
#' @param mod string defining case study, options are:
#' * linmod: linear model, no random effects
#' * randomwalk: AR1 process with drift term
#' * simpleGLMM: normally distrbuted data with global mean and group-specific random effect
#' * spatial: spatial random effect model
#' @param cov.mod string specifying if covariates are normal or uniform
#' @param theta beta coefficients
#' @param sd.vec vector of standard deviations for obs and random effects
#' @param sp.parm spatial parameters, specific to spatial case
#' @param sp.fam observation family, specific to spatial case
#' @param sp.link observation link, specific to spatial case
#' @param misp model misspecification, options are:
#' * missunifcov: missing uniform coviariate
#' * missnormcov: missing normal covariate
#' * missre: missing random effect
#' * mu0: missing drift term, specific to linmod
#' * mispre: mispecified random effect, e.g. exp(re)
#' * nb-pois: misspecified distribution (nb data fit to poisson)
#' * overdispersion: misspecified distribution (lognormal data fit to normal)
#' * lognorm-gamma: misspecified distribution (lognormal data fit to gamma)
#' * hsk: heteroskadacity introduced to simulated data
#' * mispzip: zero-inflated data fit to non-zero inflated data
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
simdat <- function(n, ng=0, mod, cov.mod = NULL, type = NULL,
                   trueparms = list(beta, sd.vec, sp.parm, fam, link),
                   misp, seed) {
  for(i in 1:length(misp)){
  if(!(misp[i] %in% c('misscov', 'missunifcov', 'missnormcov', 'missre', 
                   'mu0', 'mispre', 'nb-pois', 'overdispersion', 
                   'gamma-normal', 'hsk', 'pois-zip', 'aniso',
                   'normal-gamma', 'ln-error', 'identity-log'))){
    stop('incorrect mis-specification name')
  } }
  if(!(mod %in% c('linmod', 'pois_glm', 'randomwalk', 
                  'simpleGLMM', 'spatial', 'phylo'))) stop('incorrect model name')

  #Simulate Covariate
  if(!is.null(cov.mod)){
    set.seed(seed)
    if(cov.mod == 'norm') X <- cbind(rep(1,n), rnorm(n))
    if(cov.mod == 'unif') X <- cbind(rep(1,n), runif(n, -.5,.5))
  } else {
    X <- matrix(1, n, 1)
  }

  if(mod == 'linmod'){
    dat.out <- simdat.linmod(n, mod, cov.mod, type,
                   trueparms, misp, seed, X)
  }
  if(mod == 'pois_glm'){
    dat.out <- simdat.glm(n, mod, cov.mod, type,
                             trueparms, misp, seed)
  }
  if(mod == 'randomwalk'){
    dat.out <- simdat.randomwalk(n, mod, cov.mod, type,
                   trueparms, misp, seed)
  }
  if(mod == 'simpleGLMM'){
    dat.out <- simdat.simpleGLMM(n, ng, mod, cov.mod, type,
                   trueparms, misp, seed, X)
  }
  if(mod == 'spatial'){
    dat.out <- simdat.spatial(n, mod, type,
                   trueparms, misp, seed, X)
  }
  if(mod == 'phylo'){
    dat.out <- simdat.phylo(n, mod, type,
                            trueparms, misp, seed, X)
  }

  return(dat.out)
}

simdat.linmod <- function(n, mod, cov.mod, type=NULL, 
  trueparms, misp, seed, X){
    if(misp != 'overdispersion') stop("Misspecification not available for linmod")
    list2env(trueparms, envir = environment(simdat.linmod))
    mu <- X %*% beta
    #Correctly specified data
    set.seed(seed*2)
    eps <- rnorm(n, 0, sd.vec[1])
    y0 <- eps + mu
    # Mis-specified data
    set.seed(seed*3)
    v <- rnorm(n,0,sd.vec[2])  #sd.vec[2] = sd.vec[1]*log(4)
    y1 <- y0*exp(v)

    dat.out <- list(y0 = y0, y1 = y1, x = X)
    return(dat.out)
}

simdat.glm <- function(n, mod, cov.mod, type = NULL,
                       trueparms, misp, seed){
  if(misp != 'pois-zip') stop("Misspecification not available for linmod")
  list2env(trueparms, envir = environment(simdat.glm))
  set.seed(seed)
  y0 <- rpois(n, exp(beta))
  set.seed(seed)
  y1 <- rbinom(n, 1, 0.7)*y0
  dat.out <- list(y0 = y0, y1 = y1)
  return(dat.out)
}

simdat.randomwalk <- function(n, mod, cov.mod, type, 
  trueparms, misp, seed){
  for(i in 1:length(misp)){
    if(!(misp[i] %in% c('missre', 'hsk', 'gamma-normal', 'normal-gamma',
                        'mu0', 'mispre'))) {
      stop("Misspecification not available for random walk")
    }
  }
  list2env(trueparms, envir = environment(simdat.randomwalk))
  ## Simulate random measurements
  set.seed(seed)
  u0 <- c(0, cumsum( rnorm(n-1, drift, sd.vec[2]) ) ) 
  set.seed(seed*10)
  y0 <- sim_y(Eta = u0, omega = rep(0, n), parm = sd.vec[1], 
              fam = fam, link = link)
  
  u1 <- list()
  y1 <- list()
  for(i in 1:length(misp)){
    y1[[i]] <- y0
    u1[[i]] <- u0
    
    if(misp[i] == "mispre"){
      #simulate non-normal random effect
      if(type == "LMM"){
        set.seed(seed)
        u.misp <- cumsum( exp(rnorm(n, drift, sd.vec[2])) )
      }
      if(type == "GLMM"){
        set.seed(seed)
       # u.misp <- cumsum( rgamma(n, 0.5, 20) )
        u.misp <- cumsum( rgamma(n, 0.5, 10) )
      }
      set.seed(seed)
      y.misp <- sim_y(Eta = u.misp, omega = rep(0, n), parm = sd.vec[1], 
                      fam = fam, link = link)
      
      y1[[i]] <- y.misp
      u1[[i]] <- u.misp
    }
    
    if(misp[i] == "hsk"){
      #simulate data with heterosckadicity
      set.seed(123)
      var.hsk <- 5 + c(rep(30, n/4), rep(-4.5, n/4), rep(30, n/4), rep(-4.5, n/4))
      y.misp <- sim_y(Eta = u0, omega = rep(0, n), parm = sqrt(var.hsk), 
                  fam = fam, link = link)
      y1[[i]] <- y.misp
    }

    
    if(misp[i] == "mu0"){
      set.seed(seed)
      u1[[i]] <- c(0, cumsum( rnorm(n-1, 0, sd.vec[2]) ) ) 
    }
  }
  
  
  
  

  
 
    
  random <- list(u0 = u0, u1 = u1)

  dat.out <- list(y0 = y0, y1 = y1, random = random)
  return(dat.out)
}
 

simdat.simpleGLMM <- function(n, ng, mod, cov.mod, type, 
  trueparms, misp, seed, X){
  for(i in 1:length(misp)){
    if(!(misp[i] %in% c('missre', 'nb-pois', 'mispre', 
      'missunifcov', 'misscovnorm'))) {
      stop("Misspecification not available for simpleGLMM")
    }
  }
  list2env(trueparms, envir = environment(simdat.simpleGLMM))
  mu <- X %*% beta 
  set.seed(seed)
  u0 <- rnorm(ng, 0, sd=sd.vec[2])
  y0 <- eta <- matrix(0, n, ng)

  inits <- NA
  if(type == "LMM"){
    inits <- sd.vec[1]
  }
  if(type == "GLMM"){
    if(fam == "NB"){
      inits <- theta
    }
    if(fam == "Tweedie"){
      inits <- c(sd.vec[1], theta)
    }
  }

  set.seed(seed*2)
  for(j in 1:ng){
    for(i in 1:n){
      eta[i,j] <- mu[i] + u0[j] 
      y0[i,j] <- sim_y(as.matrix(eta[i,j]), 0, parm=inits, fam=fam, link=link)
    }
  }
  y0 <- as.vector(y0)
  u1 <- list()
  y1 <- list()

  for(m in 1:length(misp)){
    if(misp[m] == "mispre"){
      set.seed(seed)
      u1[[m]] <- rgamma(ng, 1, 1)
      ymisp <- eta <- matrix(0, n, ng)
      set.seed(seed*2)
      for(j in 1:ng){
        for(i in 1:n){
          eta[i,j] <- mu[i] + u1[[m]][j] 
          ymisp[i,j] <- sim_y(as.matrix(eta[i,j]), 0, parm=inits, 
                                fam=fam, link=link)
        }
      }
     
      y1[[m]] <- as.vector(ymisp)
    } else {
      u1[[m]] <- u0
      y1[[m]] <- y0
    }
  }
  random <- list(u0=u0, u1=u1)
  dat.out <- list(y0 = y0, y1 = y1, random = random, x = X)
  return(dat.out)
}

simdat.spatial <- function(n, mod, type, trueparms, 
                           misp, seed, X){
  for(i in 1:length(misp)){
    if(!(misp[i] %in% c('missre', 'pois-zip', 'mispre', 'aniso', 'ln-error',
      'normal-gamma'))) {
      stop(paste0("Misspecification", misp[i]," not available for spatial"))
    }
  }
  list2env(trueparms, envir = environment(simdat.spatial))
  mu <- X %*% beta 
  set.seed(seed)
  #simulate Omega with at least a 500x500 spatial "grid"
  if(n < 500){
    nspobs <- 500
  } else {
    nspobs <- n
  }
  #simulate spatial
  Loc <- matrix(runif(nspobs*2,0,100),ncol=2)
  dmat <-  as.matrix(dist(Loc))
  Sigma <- sd.vec[2]^2 * cMatern(dmat,1,sqrt(8)/sp.parm)
  set.seed(seed)
  Omega <- t(mvtnorm::rmvnorm(1, rep(0,nrow(Sigma)), 
                              sigma = Sigma, method = 'chol'))
  samp.idx <- sample(1:nspobs, n)
  loc <- Loc[samp.idx,]
  omega0 <- Omega[samp.idx]
  #generate mesh using lower resolution location matrix
  mesh <- try(
    R.utils::withTimeout( 
      fmesher::fm_mesh_2d(loc, max.edge = c(sp.parm/3,sp.parm), 
                          offset = c(sp.parm/20,sp.parm), min.angle = 26),
      timeout = 30, onTimeout = 'silent' ))
    
  if("aniso" %in% misp){
    Loc.aniso <- Loc
    rad <- 45*pi/180
    R <- rbind(c(cos(rad), sin(rad)), c(-sin(rad), cos(rad)))
    S <- rbind(c(sqrt(15), 0), c(0, 1/sqrt(15)))
    for(s in 1:nrow(Loc)){
      Loc.aniso[s,] <- Loc[s,]%*% solve(S) %*% solve(R)
    }
    loc.aniso <- Loc.aniso[samp.idx,]
    mesh.aniso <- try(
      R.utils::withTimeout( 
        fmesher::fm_mesh_2d(loc.aniso, max.edge = c(sp.parm/3,sp.parm), 
                            offset = c(sp.parm/10,sp.parm), min.angle = 26),
        timeout = 30, onTimeout = 'silent' ))
  }
  
  if(is.character(mesh)){
    system("Taskkill /IM fmesher.exe /F")
    # warning("mesh failed in rep=", ii)
    return(NULL)
  }
  
  set.seed(seed)
  y0 <- sim_y(Eta = mu, omega=omega0,
              parm=sd.vec[1], fam=fam, link=link)
  omega1 <- list()
  y1 <- list()

  for(m in 1:length(misp)){
    omega1[[m]] <- omega0
    y1[[m]] <- y0
    if(misp[m] == "mispre"){
      set.seed(seed)
      y1[[m]] <- sim_y(Eta = mu, omega=exp(omega0),
              parm=sd.vec[1], fam=fam, link=link)
     
      omega1[[m]] <- exp(omega0)
    } 
    if(misp[[m]] == "pois-zip"){
      set.seed(seed)
      y1[[m]] <- rbinom(n, 1, 0.7) * y1[[m]]
    }
    if(misp[[m]] == "ln-error"){
      #simulate data with lognormal error
      set.seed(123)
      y1[[m]] <- y0 * exp(rnorm(n))
    }
  
  }

  random <- list(omega0 = omega0, omega1 = omega1)
  dat.out <- list(y0 = y0, y1 = y1, random = random, x = X)
  dat.out$loc = loc
  dat.out$mesh = mesh
  if("aniso" %in% misp){
    dat.out$mesh.aniso = mesh.aniso
  }
  return(dat.out)
}

simdat.phylo <- function(n, mod, type, trueparms, 
                         misp, seed, X){
  for(i in 1:length(misp)){
    if(!(misp[i] %in% c('missre', 'nb-pois', 'mispre', 'identity-log', 'misscov'))) {
      stop(paste0("Misspecification", misp[i]," not available for phlyo"))
    }
  }
  list2env(trueparms, envir = environment(simdat.phylo))
  
  mu <- X %*% beta
  set.seed(seed)
  tree = rtree(n)
  set.seed(seed*2)
  u0 <- rTrait(n = 1, phy = tree, model="BM",
               parameters = list(ancestral.state = 0, sigma2 = sd.vec[2]^2))
  set.seed(seed*3)
  if(type == "LMM"){
    parm=sd.vec[1]
  }
  if(type == "GLMM"){
   #parm <- size
    parm <- sd.vec[1]
  }
  y0 <-  sim_y(Eta = mu, omega=u0,
               parm=parm, fam=fam, link=link)
  y1 <- u1 <- list()
  for(m in 1:length(misp)){
    u1[[m]] <- u0
    y1[[m]] <- y0
    
    if(misp[m] == "identity-log"){
      y1[[m]] <- sim_y(Eta = mu, omega=u0,
                  parm=sd.vec[1], fam=fam, link="log")
    }
    if(misp[m] == "mispre"){ 
      set.seed(seed*2)
      if(type == "LMM"){
        u1[[m]] <- exp(u0)
      }
      if(type == "GLMM"){       
        # u1[[m]] <- rTrait(n = 1, phy = tree, model="OU",
        #                                            parameters = list(trend = 2,
        #                                                              optimal.value = 2,
        #                                                              alpha = -.1,
        #                                                              sigma2 = sd.vec[2]^2))
        u1[[m]] <- rTrait(n = 1, phy = tree, model="EB",
                          parameters = list(
                            ancestral.state = 1,
                            rate = -2,
                            sigma2 = sd.vec[2]^2))
        # u1[[m]] <- rTrait(n = 1, phy = tree, model="lambda",
        #                   parameters = list(lambda = 0.5,
        #                                     sigma2 = sd.vec[2]^2))
      }
      y1[[m]] <- sim_y(Eta = mu, omega=u1[[m]],
                       parm=parm, fam=fam, link=link)
    }
  }
  
  random <- list(u0=u0, u1=u1)
  dat.out <- list(y0 = y0, y1 = y1, random = random, x = X, tree = tree)
}

## functions for simulating spatial data
cMatern <- function(H, Nu, Kap) {
  ifelse(H > 0, besselK(H*Kap, Nu) * (H*Kap)^Nu, 1) / gamma(Nu) * 2^(1-Nu)
}

# Simulate spatial field
sim_omega <- function(Range, sig2, Dmat, Nu = 1, method, mesh=NULL){
  Kappa <- sqrt(8)/Range
  Phi <- Range
  Tau <- sqrt(1/(4*pi*Kappa^2*sig2))
  n <- dim(Dmat)[1]

  #Simulate random field and obs
  if(method == 'R.matern'){
    Sigma <- sig2 * cMatern(Dmat,1,Kappa)
    omega <- t(mvtnorm::rmvnorm(1, rep(0,nrow(Sigma)), 
                                sigma = Sigma, method = 'chol'))
  }
  
  if(method == 'TMB.matern'){
    #  dyn.load(dynlib('spatial'))
    dat <- list(y = rep(0,n), X = matrix(1, n,1),
                dd = Dmat, nu = Nu, mesh_i = (1:n)-1, sim_re = 1,
                family = 000, link = 2, reStruct = 00)
    dat$spde <- INLA::inla.spde2.matern(mesh)$param.inla[c('M0', 'M1', 'M2')]
    obj <- TMB::MakeADFun(data =  dat,
                          parameters = list(beta = 0, theta = 0, ln_tau = log(Tau),
                                            ln_kappa = log(Kappa), ln_sig_v = 0,
                                            omega = rep(0,n), v = rep(0,n)),
                          random = 'omega',
                          DLL = 'spatial')
    sim <- obj$simulate()
    omega <- sim$omega
    #   dyn.unload(dynlib('spatial'))
  }
  if(method == 'TMB.spde'){
    #   dyn.load(dynlib('spatial'))
    dat <- list(y = rep(0,n), X = matrix(1, n,1),
                dd = Dmat, nu = 1, mesh_i = mesh$idx$loc-1, sim_re = 1,
                family = 000, link = 2, reStruct = 10)
    dat$spde <- INLA::inla.spde2.matern(mesh)$param.inla[c('M0', 'M1', 'M2')]
    obj <- TMB::MakeADFun(data =  dat,
                          parameters = list(beta = 0, theta = 0, ln_tau = log(Tau),
                                            ln_kappa = log(Kappa), ln_sig_v = 0,
                                            omega = rep(0,mesh$n), v = rep(0,n)),
                          random = 'omega',
                          DLL = 'spatial')
    sim <- obj$simulate()
    omega <- sim$omega
    #   dyn.unload(dynlib('spatial'))
  }
  return(omega)
}


# Simulate data
sim_y <- function(Eta, omega, parm, fam, link){
  Eta <- Eta + omega
  if(is.vector(Eta)) N <- length(Eta)
  if(is.matrix(Eta)) N <- nrow(Eta)
  if(link == 'identity'){
    mu <- Eta
  }
  if(link == 'log'){
    mu <- exp(Eta)
  }
  if(fam == 'Gamma'){
    #parm = CV
    Y <- rgamma(N,1/parm[1]^2,scale=mu*parm[1]^2)
  }
  if(fam == 'Poisson'){
    Y <- rpois(N, mu)
  }
  if(fam == 'Tweedie'){
    #parm = (phi, power)
    Y <- tweedie::rtweedie(N, mu = mu, phi = parm[1], power = parm[2])
  }
  if(fam == 'Gaussian'){
    ## assuming parm[1] is SD of sampling process
    Y <- rnorm(N, mu, parm)
  }
  if(fam == 'NB'){
    #parm = size
    Y <- rnbinom(N, size = parm[1], mu = mu)
  }
  return(Y)
}
