library(MASS) # rmvnorm
library(TMB)
library(INLA)
library(DHARMa)
library(VAST)
library(ggplot2)
library(dplyr)
library(tidyr)
## Has ggpairs for comparing residual types
library(GGally)

## Some global settings
ggwidth <- 7
ggheight <- 5
theme_set(theme_bw())

## functions for simulating data
cMatern <- function(H, Nu, Kap) {
  ifelse(H > 0, besselK(H*Kap, Nu) * (H*Kap)^Nu, 1) / gamma(Nu) * 2^(1-Nu)
}

# Simulate spatial field
sim.omega <- function(Range, sig2, Dmat, Nu = 1, method, mesh){
  Kappa <- sqrt(8)/Range
  Phi <- Range
  Tau <- sqrt(1/(4*pi*Kappa^2*sig2))
  n <- dim(Dmat)[1]
  
  #Simulate random field and obs
  if(method == 'TMB.matern'){
  #  dyn.load(dynlib('spatial'))
    dat <- list(y = rep(0,n), X = matrix(1, n,1),
                dd = Dmat, nu = Nu, v_i = (1:n)-1, simRE = 1, 
                family = 000, link = 2, reStruct = 00)
    dat$spde <- inla.spde2.matern(mesh)$param.inla[c('M0', 'M1', 'M2')]
    obj <- MakeADFun(data =  dat,
                     parameters = list(beta = 0, theta = 0, log_tau = log(Tau),  
                                       log_kappa = log(Kappa),
                                       omega = rep(0,n)),
                     random = 'omega',
                     DLL = 'spatial')
    sim <- obj$simulate()
    omega <- sim$omega
 #   dyn.unload(dynlib('spatial'))
  }
  if(method == 'TMB.spde'){
 #   dyn.load(dynlib('spatial'))
    dat <- list(y = rep(0,n), X = matrix(1, n,1),
                dd = Dmat, nu = 1, v_i = mesh$idx$loc-1, simRE = 1, 
                family = 000, link = 2, reStruct = 10)
    dat$spde <- inla.spde2.matern(mesh)$param.inla[c('M0', 'M1', 'M2')]
    obj <- MakeADFun(data =  dat,
                     parameters = list(beta = 0, theta = 0, log_tau = log(Tau),  
                                       log_kappa = log(Kappa),
                                       omega = rep(0,mesh$n)),
                     random = 'omega',
                     DLL = 'spatial')
    sim <- obj$simulate()
    omega <- sim$omega
 #   dyn.unload(dynlib('spatial'))
  }
  return(omega)
}


# Simulate data
sim.data <- function(X, Beta, omega, parm, fam, link, Loc){
  Xbeta <- X%*%matrix(Beta, nrow=1)
  eta <- Xbeta + omega
  N <- nrow(X)
  if(link == 'identity'){
    mu <- eta
  }
  if(link == 'log'){
    mu <- exp(eta)
  }
  if(fam == 'Gamma'){
    Y <- rgamma(N,1/parm[1]^2,scale=mu*parm[1]^2)
  }
  if(fam == 'Poisson'){
    Y <- rpois(N, mu)
  }
  if(fam == 'Tweedie'){
    Y <- tweedie::rtweedie(N, mean, parm[1], parm[2])
  }
  return(Y)
}

