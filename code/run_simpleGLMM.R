library(TMB)
library(magrittr)
library(dplyr)
library(DHARMa)
compile('models/simpleGLMM.cpp')
dyn.load(dynlib('models/simpleGLMM'))

n.j <- 3 #number of subjects
n.i <- 10 #number of observations
b0 <- 4
sig2.y <- .1 #obs variance
sig2.u <- 3 # between group variance
set.seed(123)
u <- c(-2.75, -.3, 3.25)
y <- matrix(0, n.i, n.j)

for(j in 1:n.j){
  y[,j] <- rnorm(n.i, b0 + u[j], sqrt(sig2.y))
}

Dat <- data.frame(y = as.vector(y), group = rep(1:3, each = n.i))

Data <- list(y = Dat[,1], group = Dat[,2]-1, sim_re = 0)
Par <- list(b0 = 0, ln_sig_u = 0,# ln_sig_v = c(0,0,0),
            ln_sig_y = 0, u = rep(0, n.j))#, v = rep(0, nrow(Dat)))
obj <- MakeADFun(Data, Par, random = 'u', DLL = 'simpleGLMM')
opt <- nlminb(obj$par, obj$fn, obj$gr)
b0; opt$par[1]
c(sd(u), sqrt(sig2.y)); exp(opt$par[2:3])

osa.fg <- oneStepPredict(obj, observation.name = 'y', data.term.indicator = 'keep',
                          method = 'fullGaussian')
osa.osg <- oneStepPredict(obj, observation.name = 'y', data.term.indicator = 'keep', 
                           method = "oneStepGaussian")
osa.cdf <- oneStepPredict(obj, observation.name = 'y', data.term.indicator = 'keep',
                           method = 'cdf')
osa.gen <- oneStepPredict(obj, observation.name = 'y', data.term.indicator = 'keep',
                           method = 'oneStepGeneric')
sim.c <- replicate(1000, {obj$simulate()$y})
obj$env$data$sim_re <- 1
sim.u <- replicate(1000, {obj$simulate()$y})
res.c <- createDHARMa(sim.c, y)
res.u <- createDHARMa(sim.u, y)

ks.test(osa.fg$residual, 'pnorm')$p.value
ks.test(osa.osg$residual, 'pnorm')$p.value
ks.test(osa.cdf$residual, 'pnorm')$p.value
ks.test(osa.gen$residual, 'pnorm')$p.value
ks.test(res.c$scaledResiduals, 'punif')$p.value
ks.test(res.u$scaledResiduals, 'punif')$p.value
