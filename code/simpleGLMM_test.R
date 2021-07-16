library(TMB)
library(goftest)
library(DHARMa)
compile('models/simpleGLMM.cpp')
dyn.load(dynlib('models/simpleGLMM'))
source("code/startup.R")

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
  return(list(Data=Data, Par=Par, u = u, y=y, sig2=sig2))
}

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
dev.off()

## Try Kasper's simulation residuals using tmbstan
library(tmbstan)

#Laplace approx can be biased when n.obs small
n.groups <- 5
n.obs <- 50
out <- simulate.simpleGLMM(1, ngroups=n.groups, nobs=n.obs)
out$Data$sim_re <- 1
obj <- MakeADFun(out$Data, out$Par, random = "u", DLL = "simpleGLMM", silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr
chk <- checkConsistency(obj)
chk
s <- summary(chk)
s$marginal$p.value

report <- obj$report()
res.fg <- oneStepPredict(obj, observation.name = 'y', method = 'fullGaussian')
qqnorm(res.fg$residual); qqline(res.fg$residual)


## Remake Obj with fixed effects at their truth
FIXED <- factor(NA)
map <- list(b0=FIXED, ln_sig_u=FIXED, ln_sig_y=FIXED)
pars <- list(b0=4, ln_sig_u=log(sqrt(out$sig2[2])), ln_sig_y=log(sqrt(out$sig2[1])), u=rep(0, n.groups))
obj2 <- MakeADFun(out$Data, parameters=pars, map = map)
fit2 <- tmbstan(obj2)
u.post2 <- as.matrix(fit2)[,1:n.groups]
#posteriors centered around true values of u with some bias
hist(u.post2[,1], breaks = 100); abline(v=out$u[1],col='red')
hist(u.post2[,2], breaks = 100); abline(v=out$u[2],col='red')
hist(u.post2[,3], breaks = 100); abline(v=out$u[3],col='red')
hist(u.post2[,4], breaks = 100); abline(v=out$u[4],col='red')
hist(u.post2[,5], breaks = 100); abline(v=out$u[5],col='red')
pairs(u.post2)

#mean. <- pars$b0+u.post2[,out$Data$group+1]
#Fx <- t(apply(mean., 1, function(x) pnorm(out$Data$y, x, sd=exp(pars$ln_sig_y))))
mean. <- pars$b0+apply(u.post2,2,median)[out$Data$group+1]
Fx <- pnorm(out$Data$y, mean., sd=exp(pars$ln_sig_y))
residual <- qnorm(Fx)
qqnorm(residual); qqline(residual)
ad.test(residual, 'pnorm')
#boxplot(residual)
plot(out$Data$y, residual)

## Check correlations of resids
# pairs(residual, horInd=1:8, verInd=1:8)

## unconditional DHARMa:
u <- t(replicate(4000, {rnorm(n.groups,0,sd=exp(pars$ln_sig_u))}))
#simulations centered around 0
hist(u[,1], breaks = 100); abline(v=out$u[1],col='red')
hist(u[,2], breaks = 100); abline(v=out$u[2],col='red')
hist(u[,3], breaks = 100); abline(v=out$u[3],col='red')
hist(u[,4], breaks = 100); abline(v=out$u[4],col='red')
hist(u[,5], breaks = 100); abline(v=out$u[5],col='red')
pairs(u)
mean. <- pars$b0 + u[,out$Data$group+1]
y.sim <- sapply(1:4000, function(x) rnorm(n.groups*n.obs, mean.[x,], sd=exp(pars$ln_sig_y)))
res.dharma <- createDHARMa(y.sim,out$Data$y, fittedPredictedResponse = rep(pars$b0,n.groups*n.obs))
qqnorm(res.dharma$scaledResiduals); qqline(res.dharma$scaledResiduals)
plot( out$Data$y, res.dharma$scaledResiduals)

#use DHARMa mean in pnorm method
mean. <- pars$b0+apply(u,2,median)[out$Data$group+1]
Fx <- pnorm(out$Data$y, mean., sd=exp(pars$ln_sig_y))
residual <- qnorm(Fx)
qqnorm(residual[is.finite(residual)]); qqline(residual[is.finite(residual)])
plot(out$Data$y, residual)

#use stan mean in DHARMa method
mean. <-  pars$b0 + u.post2[,out$Data$group+1]
y.sim <- sapply(1:4000, function(x) rnorm(n.groups*n.obs, mean.[x,], sd=exp(pars$ln_sig_y)))
res.dharma <- createDHARMa(y.sim,out$Data$y, fittedPredictedResponse = rep(pars$b0,n.groups*n.obs))
qqnorm(res.dharma$scaledResiduals); qqline(res.dharma$scaledResiduals)
plot(out$Data$y, res.dharma$scaledResiduals)


## conditional DHARMa:
mean. <- pars$b0 + out$u[out$Data$group+1]
y.sim <- sapply(1:4000, function(x) rnorm(n.groups*n.obs, mean., sd=exp(pars$ln_sig_y)))
res.dharma <- createDHARMa(y.sim,out$Data$y, fittedPredictedResponse = rep(pars$b0,n.groups*n.obs))
qqnorm(res.dharma$scaledResiduals); qqline(res.dharma$scaledResiduals)
plot(out$Data$y, res.dharma$scaledResiduals )
