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
#' * outliers: randomly added outlier terms
#' * misscov: missing coviariate
#' * overdispersion: observation level random error added to observations
#' * mu0: missing drift term, specific to linmod
#' * mispomega: non-normal spatial effect, specific to spatial
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
simdat <- function(n, ng=0, mod, cov.mod = 'norm',
                   trueparms = list(theta, sd.vec, sp.parm, fam, link),
                   misp, seed){
  if(!(misp %in% c('outliers', 'misscov', 'overdispersion', 
                   'mu0', 'normal', 'mispomega', 'dropRE', 
                   'deltagamma', 'aniso'))){
    stop('incorrect mis-specification name')
  } 
  if(!(mod %in% c('linmod', 'randomwalk', 'simpleGLMM', 'spatial'))) stop('incorrect model name')

  list2env(trueparms, envir = environment(simdat))

  #Simulate Covariate
  set.seed(seed)
  N <- n
  if(mod == 'simpleGLMM') N <- n*ng
  if(misp == 'misscov' | mod == 'linmod'){
    if(cov.mod == 'norm') X <- cbind(rep(1,n), rnorm(n))
    if(cov.mod == 'unif') X <- cbind(rep(1,n), runif(n, 5,10))
  } else {
    X <- matrix(1, n, 1)
  }
  mu <- X %*% theta

  if(mod == 'linmod'){
    set.seed(seed*2)
    eps <- rnorm(N, 0, sd.vec[1])
    y0 <- eps + mu
    v <- NULL
    if(misp == 'mu0' | misp == 'mispomega') stop("Misspecification not available for linmod")
    if(misp == 'overdispersion'){
      set.seed(seed*3)
      v <- rnorm(N,0,sd.vec[2])  #sd.vec[2] = sd.vec[1]*log(4)
      y1 <- y0*exp(v)
    }
    if(misp == 'misscov'){
      y1 <- y0
    }
    random <- list(v=v)
  }

  if(mod == 'randomwalk'){
    if(misp == 'overdispersion' | misp == 'misscov'| misp == 'mispomega'){
      stop("Misspecification not available for random walk")
    }
    set.seed(seed)
    ## Simulate random measurements
    u <- c(0,cumsum(rnorm(N-1,mean=mu,sd=sd.vec[2])))
    if(misp == "normal"){
      y0 <- u + rlnorm(N,sd=sd.vec[1])
    } else {
      y0 <- u + rnorm(N,sd=sd.vec[1])
    }
    
    y1 <- y0
    
    random <- list(u=u)
  }

  if(mod == 'simpleGLMM'){ 
    if(misp == 'mu0'| misp == 'mispomega') stop("Misspecification not available for simpleGLMM")
    set.seed(seed)
    u <- rnorm(ng, 0, sd=sd.vec[2])
    y0 <- eta <- matrix(0, n, ng)
    v <- rep(0, n)
    if(misp == 'overdispersion'){
      set.seed(seed*3)
      v <- rnorm(n,0,sd.vec[3]) #sd.vec[3] = 1
    }
    set.seed(seed*2)
    inits <- sd.vec[1]
    for(j in 1:ng){
      for(i in 1:n){
        eta[i,j] <- mu[i] + u[j] + v[i]
        # will not work if Tweedie and misp = overdispersion
        if(fam == 'Tweedie') inits <- c(inits, sd.vec[3])
        y0[i,j] <- sim_y(as.matrix(eta[i,j]), 0, parm=inits, fam=fam, link=link)
      }
    }
    y0 <- as.vector(y0)
    y1 <- y0
    random <- list(u=u,v=v)
  }

  if(mod == 'spatial'){
    ## Simulate spatial random effects
    set.seed(seed)
    Loc <- matrix(runif(1000*2,0,100),ncol=2)
    dmat <-  as.matrix(dist(Loc))
    Sigma <- sd.vec[2]^2 * cMatern(dmat,1,sqrt(8)/sp.parm)
    Omega <- t(mvtnorm::rmvnorm(1, rep(0,nrow(Sigma)), 
                                sigma = Sigma, method = 'chol'))
    samp.idx <- sample(1:1000, N)
    loc <- Loc[samp.idx,]
    omega <- Omega[samp.idx]
    
    mesh <- try(
      R.utils::withTimeout( 
        INLA::inla.mesh.2d(loc, max.edge = c(sp.parm/3,sp.parm), 
                           offset = c(sp.parm/10,sp.parm*2), min.angle = 26),
        timeout = 30, onTimeout = 'silent' ))
    
    if(misp == "aniso"){
      loc.aniso <- loc
      rad <- 45*pi/180
      R <- rbind(c(cos(rad), sin(rad)), c(-sin(rad), cos(rad)))
      S <- rbind(c(sqrt(10), 0), c(0, 1/sqrt(10)))
      for(s in 1:nrow(loc)){
        loc.aniso[s,] <- loc[s,]%*% solve(S) %*% solve(R)
      }
      mesh.aniso <- try(
        R.utils::withTimeout( 
          INLA::inla.mesh.2d(loc.aniso, max.edge = c(sp.parm/3,sp.parm), 
                             offset = c(sp.parm/10,sp.parm*2), min.angle = 26),
          timeout = 30, onTimeout = 'silent' ))
    }
   
    if(is.character(mesh)){
      system("Taskkill /IM fmesher.exe /F")
     # warning("mesh failed in rep=", ii)
      return(NULL)
    }
    # Omega <- sim_omega(sp.parm,sd.vec[2]^2,dmat,method="R.matern",mesh=mesh)
    # if(length(Omega)>n){
    #   omega <- Omega[mesh$idx$loc]
    # } else {
    #   omega <- Omega
    # }
    v <- rep(0, n)
    if(misp == 'overdispersion'){
      set.seed(seed)
      v <- rnorm(N,0,sd.vec[3]) #sd.vec[3] = sd.vec[1]*2
      mu <- mu + v
    }
    set.seed(seed)
    y0 <- sim_y(Eta = mu, omega=omega,
                  parm=sd.vec[1], fam=fam, link=link)


    if(misp == 'mu0') stop("Misspecification not available for spatial")
    if(misp == 'misscov' | misp == 'overdispersion' | 
       misp == 'dropRE' | misp == "aniso"){
      y1 <- y0
    }
    if(misp=='mispomega'){
      set.seed(seed)
      y1 <- sim_y(Eta = mu, omega=exp(omega),
                  parm=sd.vec, fam=fam, link=link)
    }
    random <- list(omega=omega, v=v)
  }

  if(misp == 'outliers'){
    set.seed(seed)
    noutlier <- 5
    y1 <- y0
    ind <- sample(1:N, size=noutlier)
    y1[ind] <- y0[ind]+rnorm(noutlier, 0, sd.vec[1]*4)
  }

  dat.out <- list(y0=y0, y1=y1, x=X, random=random)
  if(mod == 'spatial'){
    dat.out$loc = loc
    dat.out$mesh = mesh
    if(misp == "aniso"){
      dat.out$mesh.aniso = mesh.aniso
    }
  }
  return(dat.out)
}


## functions for simulating spatial data
cMatern <- function(H, Nu, Kap) {
  ifelse(H > 0, besselK(H*Kap, Nu) * (H*Kap)^Nu, 1) / gamma(Nu) * 2^(1-Nu)
}

# Simulate spatial field
sim_omega <- function(Range, sig2, Dmat, Nu = 1, method, mesh){
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
  N <- nrow(Eta)
  if(link == 'identity'){
    mu <- Eta
  }
  if(link == 'log'){
    mu <- exp(Eta)
  }
  if(fam == 'Gamma'){
    Y <- rgamma(N,1/parm[1]^2,scale=mu*parm[1]^2)
  }
  if(fam == 'Poisson'){
    Y <- rpois(N, mu)
  }
  if(fam == 'Tweedie'){
    Y <- tweedie::rtweedie(N, mu = mu, phi = parm[1], power = parm[2])
  }
  if(fam == 'Gaussian'){
    ## assuming parm[1] is SD of sampling process
    Y <- rnorm(N, mu, parm[1])
  }
  return(Y)
}
