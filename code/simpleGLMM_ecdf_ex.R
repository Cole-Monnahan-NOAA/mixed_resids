# Reproduce analysis for the simpleGLMM model 
# with a uniformally distributed covariate

library(DHARMa)
library(TMB)
compile('src/simpleGLMM.cpp')
dyn.load(dynlib('src/simpleGLMM'))

simulate.simpleGLMM <- function(seed = NULL, b0 = 4, 
                                sig.y = 0.5,
                                sig.u = 2,
                                ngroups=5, nobs=100){

  #reproduce simulation steps and simulation seeds
  set.seed(seed)
  X <- cbind(rep(1, nobs), runif(nobs, -.5, .5))
  
  set.seed(seed)
  u <- rnorm(ngroups, 0, sig.u)
  
  set.seed(seed*2)
  y <- matrix(0, nobs, ngroups)
  for(i in 1:nobs){
    for(j in 1:ngroups){
      y[,j] <- rnorm(nobs, b0 + u[j], sig.y)
    }
  }

  Dat <- list(y = as.vector(y), X = X,
              group = rep(1:ngroups, each = nobs) - 1, #C++ starts indexing at 0
              obs = rep(1:nobs, ngroups) - 1, #C++ starts indexing at 0
              family = 0, link = 2, #family = gaussian(link = identity)
              sim_re = 0)
  Par <- list(beta = c(0,0), ln_sig_y = 0, ln_sig_u = 0, ln_sig_v = numeric(0),
              u=rep(0, ngroups), v = rep(0, nobs))
  return(list(Data = Dat, Par = Par, u = u, y = y, sdvec = c(sig.u, sig.y)))
}

#Single model run
n.gr=5; n.obs=100
out <- simulate.simpleGLMM(seed = 1, ngroups=n.gr, nobs=n.obs)
obj <- MakeADFun(out$Data, out$Par, 
                 random = "u", 
                 map = list(v = rep(factor(NA), n.obs)), #fix v at 0
                 DLL = "simpleGLMM", silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
report <- obj$report(obj$env$last.par.best)

#ecdf residuals from DHARMa
#Conditional
ecdf.cond <- createDHARMa(replicate(1000, obj$simulate()$y),
                          observedResponse = out$Data$y,
                          fittedPredictedResponse = report$fpr,
                          integerResponse = FALSE)
ks.test(ecdf.cond$scaledResiduals, "punif")

#Unconditional
#turn on RE simulation in TMB model 
#(corresponds to running line 90 in simpleGLMM.cpp)
obj$env$data$sim_re <- 1
ecdf.uncond <- createDHARMa(replicate(1000, obj$simulate()$y),
                            observedResponse = out$Data$y,
                            fittedPredictedResponse = report$fpr,
                            integerResponse = FALSE)
ks.test(ecdf.uncond$scaledResiduals, "punif")

#Unconditional with rotation
ecdf.uncond.rot <- createDHARMa(replicate(1000, obj$simulate()$y),
                              observedResponse = out$Data$y,
                              fittedPredictedResponse = report$fpr,
                              integerResponse = FALSE,
                              rotation = "estimated")
ks.test(ecdf.uncond.rot$scaledResiduals, "punif")

#Confirm RE simulation flag is working
#Estimated RE
obj$env$parList()$u

#Turn off RE simulation
obj$env$data$sim_re <- 0
replicate(3, obj$simulate()$u)

#Turn on RE simulation
obj$env$data$sim_re <- 1
replicate(3, obj$simulate()$u)
