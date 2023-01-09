rm(list=ls())
## https://github.com/kaskr/adcomp/blob/master/TMB/inst/examples/simple.cpp
runExample('simple')
observation.name <- 'x'
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
Sigma <- solve(as.matrix(newobj2$env$spHess(mode, random = TRUE)))
i <- which(names(newobj2$env$par[newobj2$env$random]) == observation.name)
res <- obs - mode[i]
## Do the rotation of just the data.
L <- t(chol(Sigma[i,i]))
pred <- data.frame(residual = as.vector(solve(L, res)))

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
x.ind <- which(modedf$par == 'x')
ggplot(modedf, aes(mode0, mode)) + facet_wrap('par', scales='free') + geom_point()

corr <- cov2cor(Sigma)
par <- paste0(modedf$par, c(1:length(x.ind), 1:length(u.ind)))
dimnames(corr) <- list(par,par)
n <- 10
ind <- c(head(x.ind,n),head(u.ind,n))
library(corrplot)
corrplot(corr[ind,ind], type='upper')
corrplot(corr[x.ind[1:10], x.ind[1:10]], type='upper', diag=FALSE)

## Try doing it conditionally
i <- which(names(newobj2$env$par[newobj2$env$random]) == observation.name)
Sigma2 <- solve(as.matrix(newobj2$env$spHess(mode, random = TRUE))[i,i])
corr2 <- cov2cor(Sigma2)
corrplot(corr2[1:10, 1:10], type='upper')
