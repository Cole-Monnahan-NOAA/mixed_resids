library(TMB)
library(corrplot)
## This has a function to return the covar matrix of predicted
## data for the FG method stripped down
source('extractOSA.R')
## simulate a simple LMM

set.seed(101)
logsigma <- 1
logtau <- -1
ngroups <- 4
nobs <- 3
means <- rnorm(ngroups, 0, exp(logsigma))
group <- rep(1:ngroups, each=nobs)
obs <- rnorm(ngroups*nobs, mean=means[group],
             sd=exp(logtau))
boxplot(obs~group)
data <- list(obs=obs, group=group)
pars <- list(beta=rep(0,ngroups), logsigma=1)

## A simple linear model (no random effects)
compile('../src/demo_lm.cpp')
dyn.load('../src/demo_lm')
obj.lm <- MakeADFun(data, pars, DLL='demo_lm')
opt.lm <- with(obj.lm, nlminb(par, fn, gr))
cor.lm <- cov2cor(extractOSAcovar(obj.lm, 'obs'))
osa.lm <- oneStepPredict(obj.lm, method='fullGaussian', observation.name='obs')$residual
pea.lm <- (obs-obj.lm$report()$mu)/exp(opt.lm$par[ngroups+1])
osa.lm-pea.lm # identical

## Repeat with random group effects
compile('../src/demo_lmm.cpp')
dyn.load('../src/demo_lmm')
pars$logtau <- 1
obj.lmm <- MakeADFun(data, pars, random='beta', DLL='demo_lmm')
opt.lmm <- with(obj.lmm, nlminb(par, fn, gr))
cor.lmm <- cov2cor(extractOSAcovar(obj.lmm, 'obs'))
osa.lmm <- oneStepPredict(obj.lmm, method='fullGaussian', observation.name='obs')$residual
pea.lmm <- (obs-obj.lmm$report()$mu)/exp(opt.lmm$par[2])
osa.lmm-pea.lmm # not identical


## very similar fits and data predictions
boxplot(obs~group)
points(1:ngroups, obj.lm$report()$beta, pch=16, col=2)
points(.25+1:ngroups, obj.lmm$report()$beta, pch=16, col=3)

## but very different correlations
par(mfrow=c(1,2))
corrplot(cor.lm, diag=TRUE, type='lower')
corrplot(cor.lmm, diag=TRUE, type='lower')

out <- data.frame(obs=obs, osa.lm, pea.lm, osa.lmm, pea.lmm)
