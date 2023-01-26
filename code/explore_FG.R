rm(list=ls())
## https://github.com/kaskr/adcomp/blob/master/TMB/inst/examples/simple.cpp
#runExample('simple')

source("R/startup.R")

mod <- "randomwalk"
misp <- "mu0"
setupTMB(mod)
true.parms <- setup_trueparms(mod, misp, fam = "Gaussian", link = "identity")
sim.dat <- simdat(n=100, ng=0, mod, 
                  cov.mod="norm", true.parms, 
                  misp, seed = 1)
init.dat <- mkTMBdat(sim.dat, true.parms, mod, misp)
init.par <- mkTMBpar(true.parms, sim.dat, mod, misp, FALSE) #use FALSE b/c estimating model
init.random <- mkTMBrandom(mod, misp, FALSE)
init.map <- mkTMBmap(init.par, mod, misp, true.parms$fam, FALSE)
h <- 1 #fit correct model
obj <- MakeADFun(data = init.dat[[h]], parameters = init.par[[h]], 
                 map = init.map[[h]], random = init.random[[h]], DLL = mod)
opt <- nlminb(obj$par, obj$fn, obj$gr)

observation.name <- 'y'
obs <- as.vector(obj$env$data[[observation.name]])
## Rebuild arg list from the original fit
args0 <- args <- as.list(obj$env)[intersect(names(formals(MakeADFun)), ls(obj$env))]
## Update par list with the MLEs
args$parameters <- obj$env$parList(par = obj$env$last.par.best)
## Turn off estimation of fixed effects
names.random <- unique(names(obj$env$par[obj$env$random]))
names.fixed <- setdiff(names(args$parameters), names.random)
args$map <- lapply(args$parameters[names.fixed], function(x) factor(x * NA))
## Move observations from data into parameter list, delcare as
## random effects, turn on estimation of it via map
args$random <- c(names.random, observation.name)
args$parameters[observation.name] <- args$data[observation.name]
args$data[observation.name] <- NULL
## Build new object
newobj2 <- do.call("MakeADFun", args)
## mode before optimization the RE
mode0 <- newobj2$env$last.par
newobj2$fn()
## mode after
mode <- newobj2$env$last.par
## Get joint hessian of data and random effects. The inversion
## includes the real random effects and this is what is
## "unconditional" about it.

#rotate then subset
Sigma <- solve(as.matrix(newobj2$env$spHess(mode, random = TRUE)))
L <- t(chol(Sigma))
i <- which(names(newobj2$env$par[newobj2$env$random]) == observation.name)
## Do the rotation of just the data.
#L <- t(chol(Sigma[i,i]))
#res <- obs - mode[i]
res <- mode0 - mode
r1 <- as.vector(solve(L, res))[i] #subset on obs here
#pred <- data.frame(residual = as.vector(solve(L, res)))

##TMB way (subset then rotate):
h <- newobj2$env$spHess(mode, random = TRUE)
i <- which(names(newobj2$env$par[newobj2$env$random]) == observation.name)
Sigma2 <- solve(as.matrix(GMRFmarginal(h, i)))
res <- obs - mode[i]
L <- t(chol(Sigma2))
r2 <- as.vector(solve(L, res))

plot(r1, r2);abline(0,1)
sum(r1-r2)


## explore what it's doing a bit
names(args0$data)
names(args$data)
str(args0$parameters)
str(args$parameters)
str(args0$map)
str(args$map)
str(args0$random)                       # u only
str(args$random)

modedf <- data.frame(par=names(mode0), mode0, mode)
u.ind <- which(modedf$par == 'u')
x.ind <- which(modedf$par == 'y')
ggplot(modedf, aes(mode0, mode)) + facet_wrap('par', scales='free') + geom_point()

corr <- cov2cor(Sigma)
par <- paste0(modedf$par, c(1:length(x.ind), 1:length(u.ind)))
dimnames(corr) <- list(par,par)
n <- 10
ind <- c(head(x.ind,n),head(u.ind,n))
library(corrplot)
corrplot(corr[ind,ind], type='upper')
corrplot(corr[x.ind[1:10], x.ind[1:10]], type='upper', diag=FALSE)

## Try doing it conditionally. Is this right?
# i <- which(names(newobj2$env$par[newobj2$env$random]) == observation.name)
# Sigma2 <- solve(as.matrix(newobj2$env$spHess(mode, random = TRUE))[i,i])
corr2 <- cov2cor(Sigma2)
corrplot(corr2[1:10, 1:10], type='upper', diag=FALSE)

#visualize with INLA
library(INLA)
#both matrices show correlation
image(Sigma) #joint correlation of RE and data
image(Sigma2) #correlation of just data

GMRFmarginal <- function(Q, i, ...) {
  ind <- 1:nrow(Q)
  i1 <- (ind)[i]
  i0 <- setdiff(ind, i1)
  if (length(i0) == 0)
    return(Q)
  Q0 <- as(Q[i0, i0, drop = FALSE], "symmetricMatrix")
  L0 <- Cholesky(Q0, ...)
  ans <- Q[i1, i1, drop = FALSE] - Q[i1, i0, drop = FALSE] %*%
    solve(Q0, Q[i0, i1, drop = FALSE])
  ans
}

Sigma <- solve(as.matrix(GMRFmarginal(h, i)))
