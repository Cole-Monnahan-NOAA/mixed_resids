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
source("code/startup.R")

#check consistency
out <- simulate.simpleGLMM(seed = 123)
out$Data$sim_re <- 1
obj <- MakeADFun(out$Data, out$Par, random = "u", DLL = "simpleGLMM", silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
print(checkConsistency(obj))

nsim <- 500
uncond.pval <- cond.pval <- fg.pval <- rep(0,nsim)
set.seed(123)
for(i in 1:nsim){
  out <- simulate.simpleGLMM(i)
  obj <- MakeADFun(out$Data, out$Par, random = "u", DLL = "simpleGLMM", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  rep <- obj$report()
  osa <- oneStepPredict(obj,observation.name = 'y',
                        method = 'fullGaussian')
  fg.pval[i] <- ad.test(osa$residual, null = 'pnorm', estimated = TRUE)$p.value[[1]]
  ## Conditional DHARMa
  obj$env$data$sim_re <- 0 #turn off RE simulation
  dharma.cond <- createDHARMa(replicate(1000, obj$simulate()$y),
                              observedResponse=out$Data$y,
                              fittedPredictedResponse = rep$mu)
  resids.cond <- residuals(dharma.cond, quantileFunction = qnorm,
                           outlierValues = c(-7,7))
  cond.pval[i] <- ad.test(resids.cond, null = 'pnorm', estimated = TRUE)$p.value[[1]]
  ## Unconditional DHARMa
  obj$env$data$sim_re <- 1 #turn on RE simulation
  dharma.uncond <- createDHARMa(replicate(1000, obj$simulate()$y),
                                observedResponse=out$Data$y,
                                fittedPredictedResponse = rep$mu)
  resids.uncond <- residuals(dharma.uncond, quantileFunction = qnorm,
                             outlierValues = c(-7,7))
  uncond.pval[i] <- ad.test(resids.uncond, null = 'pnorm', estimated = TRUE)$p.value[[1]]
}

par(mfrow=c(3,1))
hist(fg.pval, xlim = c(0,1))
hist(cond.pval, xlim = c(0,1))
hist(uncond.pval, xlim = c(0,1))


## Try Kasper's simulation residuals using tmbstan
library(tmbstan)
out <- simulate.simpleGLMM(1, ngroups=10, nobs=5)
obj <- MakeADFun(out$Data, out$Par, random = "u", DLL = "simpleGLMM", silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Remake Obj with fixed effects at their truth
FIXED <- factor(NA)
map <- list(b0=FIXED, ln_sig_u=FIXED, ln_sig_y=FIXED)
pars <- list(b0=4, ln_sig_u=log(10), ln_sig_y=log(.5), u=rep(0, 10))
obj2 <- MakeADFun(out$Data, parameters=pars, map = map)
fit2 <- tmbstan(obj2)
u.post2 <- as.matrix(fit2)[,1:10]

mean <- pars$b0+u.post2[,out$Data$group+1]
Fx <- t(apply(mean, 1, function(x) pnorm(out$Data$y, x, sd=exp(pars$ln_sig_y))))
residual <- qnorm(Fx)
qqnorm(residual); qqline(residual)
boxplot(residual)
## Check correlations of resids
pairs(residual, horInd=1:8, verInd=1:8)

