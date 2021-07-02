simulate.simpleGLMM <- function(seed = NULL, b0 = 4, sig2 = c(.5,10), ngroups=3, nobs=10){
  sig2.y <- sig2[1] # within group variance
  sig2.u <- sig2[2] # between group variance
  set.seed(seed)
  u <- rnorm(ngroups, 0, sqrt(sig2.u))
  y <- matrix(0, nobs, ngroups)
  for(j in 1:ngroups){
    y[,j] <- rnorm(nobs, b0 + u[j], sqrt(sig2.y))
  }

  Dat <- data.frame(y = as.vector(y), group = rep(1:ngroups, each = nobs))
  Data <- list(y = Dat[,1], group = Dat[,2]-1, sim_re = 0)
  Par <- list(b0=0, ln_sig_u=0, ln_sig_y=0, u=rep(0, ngroups))
  return(list(Data=Data, Par=Par, u = u))
}

library(TMB)
library(goftest)
compile('models/simpleGLMM.cpp')
dyn.load(dynlib('models/simpleGLMM'))

#check consistency
out <- simulate.simpleGLMM(seed = 1)
out$Data$sim_re <- 1
obj <- MakeADFun(out$Data, out$Par, random = "u", DLL = "simpleGLMM", silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
print(checkConsistency(obj))

nsim <- 500
fg.pval <- rep(0,nsim)
set.seed(123)
for(i in 1:nsim){
  out <- simulate.simpleGLMM()
  obj <- MakeADFun(out$Data, out$Par, random = "u", DLL = "simpleGLMM", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  osa <- oneStepPredict(obj,observation.name = 'y',
                        data.term.indicator = 'keep',
                        method = 'fullGaussian')
  fg.pval[i] <- ad.test(osa$residual, null = 'pnorm', estimated = TRUE)$p.value[[1]]
}
hist(fg.pval, xlim = c(0,1),breaks= 500)
