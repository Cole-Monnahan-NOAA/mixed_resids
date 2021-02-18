
## Pull out pieces of OSA calcs for FullGaussian and recreate
## minimally to explore.
rm(list=ls())
source('code/run_simpleGLMM.R')
## Rebuild obj just in case
obj <- MakeADFun(Data, Par, random = 'u', DLL = 'simpleGLMM')
opt <- nlminb(obj$par, obj$fn, obj$gr)


pred.FG <- oneStepPredict(obj, 'y', method='fullGaussian')
pred.OSG <- oneStepPredict(obj, 'y', data.term.indicator='keep', method='fullGaussian')
all.equal(pred.FG, pred.OSG)
## So yes they give exactly teh same residuals in this case

### ------------------------------------------------------------
### First try visualizing the Full Guassian method
## Manualyl specify the input arguments
observation.name = 'y'
data.term.indicator = NULL#'keep'
method=c("oneStepGaussianOffMode", "fullGaussian",
         "oneStepGeneric", "oneStepGaussian", "cdf")[2]
subset <- NULL# = 1:100
conditional = NULL
discrete = NULL
## discreteSupport = NULL
range = c(-Inf, Inf)
seed = 123
parallel = FALSE
trace = TRUE
reverse = (method == "oneStepGaussianOffMode")
obs <- as.vector(obj$env$data[[observation.name]])
subset <- 1:length(obs)
unconditional <- setdiff(1:length(obs), union(subset, conditional))
args <- as.list(obj$env)[intersect(names(formals(MakeADFun)), ls(obj$env))]
if (length(obj$env$random)) {
  args$parameters <- obj$env$parList(par = obj$env$last.par.best)
} else {
  args$parameters <- obj$env$parList(obj$env$last.par.best)
}
names.random <- unique(names(obj$env$par[obj$env$random]))
names.all <- names(args$parameters)
fix <- setdiff(names.all, names.random)
map <- lapply(args$parameters[fix], function(x) factor(x * NA))
ran.in.map <- names.random[names.random %in% names(args$map)]
if (length(ran.in.map)) map <- c(map, args$map[ran.in.map])
args$map <- map
args$random <- names.random
args$regexp <- FALSE
args$parameters[observation.name] <- args$data[observation.name]
args$data[observation.name] <- NULL
if (!is.null(data.term.indicator)) {
  one <- rep(1, length(obs))
  zero <- rep(0, length(obs))
  if (method == "cdf") {
    args$parameters[[data.term.indicator]] <- cbind(one,
                                                    zero, zero)
  }
  else {
    args$parameters[[data.term.indicator]] <- cbind(one)
  }
}
if (length(unconditional) > 0) {
  if (is.null(data.term.indicator))
    stop("Failed to disable some data terms (because 'data.term.indicator' missing)")
  args$parameters[[data.term.indicator]][unconditional, 1] <- 0
}
if (length(conditional) > 0) {
  if (is.null(data.term.indicator))
    stop("Failed to enable some data terms (because 'data.term.indicator' missing)")
  args$parameters[[data.term.indicator]][conditional, 1] <- 1
}
makeFac <- function(x) {
  fac <- as.matrix(x)
  fac[] <- 1:length(x)
  fac[conditional, ] <- NA
  fac[unconditional, ] <- NA
  fac[subset, ] <- 1:(length(subset) * ncol(fac))
  factor(fac)
}
map <- list()
map[[observation.name]] <- makeFac(obs)
if (!is.null(data.term.indicator)) {
  map[[data.term.indicator]] <- makeFac(args$parameters[[data.term.indicator]])
}
args$map <- c(args$map, map)
args$silent <- TRUE
newobj <- do.call("MakeADFun", args)
nm <- names(newobj$par)
obs.pointer <- which(nm == observation.name)
if (method == "cdf") {
  tmp <- matrix(which(nm == data.term.indicator), ncol = 3)
  data.term.pointer <- tmp[, 1]
  lower.cdf.pointer <- tmp[, 2]
  upper.cdf.pointer <- tmp[, 3]
} else {
  data.term.pointer <- which(nm == data.term.indicator)
  lower.cdf.pointer <- NULL
  upper.cdf.pointer <- NULL
}
observation <- local({
  obs.local <- newobj$par
  i <- 1:length(subset)
  function(k, y = NULL, lower.cdf = FALSE, upper.cdf = FALSE) {
    obs.local[data.term.pointer[k < i]] <- 0
    if (!is.null(y))
      obs.local[obs.pointer[k]] <- y
    if (lower.cdf | upper.cdf) {
      obs.local[data.term.pointer[k]] <- 0
      if (lower.cdf)
        obs.local[lower.cdf.pointer[k]] <- 1
      if (upper.cdf)
        obs.local[upper.cdf.pointer[k]] <- 1
    }
    obs.local
  }
})
tracefun <- function(k) if (trace) print(k)
applyMethod <- function(oneStepMethod) {
  ord <- seq_along(subset)
  if (reverse)
    ord <- rev(ord)
  pred <- do.call("rbind", lapply(ord, oneStepMethod))
  pred <- as.data.frame(pred)[ord, ]
  pred$Fx <- 1/(1 + exp(pred$nlcdf.lower - pred$nlcdf.upper))
  pred$px <- 1/(exp(-pred$nlcdf.lower + pred$nll) + exp(-pred$nlcdf.upper +
                                                        pred$nll))
  if (discrete) {
    if (!is.null(seed)) {
      Random.seed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- Random.seed)
      set.seed(seed)
    }
    U <- runif(nrow(pred))
  }
  else {
    U <- 0
  }
  pred$residual <- qnorm(pred$Fx - U * pred$px)
  pred
}
## Done preparing object
##if (method == "fullGaussian") {
args2 <- args
args2$random <- c(args2$random, observation.name)
fix <- data.term.indicator
args2$map[fix] <- lapply(args2$map[fix], function(x) factor(NA * unclass(x)))
newobj2 <- do.call("MakeADFun", args2)
newobj2$fn()
mode <- newobj2$env$last.par
GMRFmarginal <- function(Q, i, ...) {
  ind <- 1:nrow(Q)
  i1 <- (ind)[i]
  i0 <- setdiff(ind, i1)
  if (length(i0) == 0) return(Q)
  Q0 <- as(Q[i0, i0, drop = FALSE], "symmetricMatrix")
  L0 <- Matrix::Cholesky(Q0, ...)
  ans <- Q[i1, i1, drop = FALSE] - Q[i1, i0, drop = FALSE] %*%
    Matrix::solve(Q0, Q[i0, i1, drop = FALSE])
  ans
}
h <- newobj2$env$spHess(mode, random = TRUE)
i <- which(names(newobj2$env$par[newobj2$env$random]) ==
           observation.name)
Sigma <- solve(as.matrix(GMRFmarginal(h, i)))
res <- obs[subset] - mode[i]
L <- t(chol(Sigma))
pred <- data.frame(residual = as.vector(solve(L, res)))
##}
## Hack way to generate random draws and plot the real data vs
## them.
library(mvtnorm)
nsim <- 1000
ysub <- 5
x <- rmvnorm(n=nsim, mean=mode[i], sigma=Sigma)
x <- rbind(x, obs[subset], mode[i])
ind <- c(1,2,3, 11,12,13, 21,22,23)
labs <- paste0('y[',ind,']\n', round(pred[ind,1],2))
png('plots/glmm_fullGaussian_mvn.png', width=7, height=7,
    units='in', res=500)
pairs(x[,ind], upper.panel=NULL, labels=labs,
      col=c(rep(rgb(0,0,0,.1),nsim), rgb(1,0,0), rgb(0,1,0)),
      pch=16, gap=0)
dev.off()
## Now look at Cholesky rotated version
xseq <- seq(-5,5, len=1000); yseq <- dnorm(xseq)
png('plots/glmm_fullGaussian_znorm.png', width=9, height=5,
    units='in', res=500)
par(mfrow=c(1,2), mgp=c(1.5,.5,0), tck=-.02, mar=c(3,3,.5,.5))
plot(xseq, yseq, type='l', xlab='Residual', ylab='Density')
mycol <- rep(c(1,2,3), each=10)
text(pred[,1], y=dnorm(pred[,1]), labels=paste0('y[', 1:30,']'),
     cex=.7, col=mycol)
boxplot(pred[,1]~Dat[,2], xlab='Group', ylab='Residual')
dev.off()



### ------------------------------------------------------------
## Repeat with oneStepGauss
observation.name = 'y'
data.term.indicator = 'keep'
method=c("oneStepGaussianOffMode", "fullGaussian",
         "oneStepGeneric", "oneStepGaussian", "cdf")[4]
subset <- NULL# = 1:100
conditional = NULL
discrete = NULL
## discreteSupport = NULL
range = c(-Inf, Inf)
seed = 123
parallel = FALSE
trace = TRUE
reverse = (method == "oneStepGaussianOffMode")
obs <- as.vector(obj$env$data[[observation.name]])
subset <- 1:length(obs)
unconditional <- setdiff(1:length(obs), union(subset, conditional))
args <- as.list(obj$env)[intersect(names(formals(MakeADFun)), ls(obj$env))]
if (length(obj$env$random)) {
  args$parameters <- obj$env$parList(par = obj$env$last.par.best)
} else {
  args$parameters <- obj$env$parList(obj$env$last.par.best)
}
names.random <- unique(names(obj$env$par[obj$env$random]))
names.all <- names(args$parameters)
fix <- setdiff(names.all, names.random)
map <- lapply(args$parameters[fix], function(x) factor(x * NA))
ran.in.map <- names.random[names.random %in% names(args$map)]
if (length(ran.in.map)) map <- c(map, args$map[ran.in.map])
args$map <- map
args$random <- names.random
args$regexp <- FALSE
args$parameters[observation.name] <- args$data[observation.name]
args$data[observation.name] <- NULL
if (!is.null(data.term.indicator)) {
  one <- rep(1, length(obs))
  zero <- rep(0, length(obs))
  if (method == "cdf") {
    args$parameters[[data.term.indicator]] <- cbind(one,
                                                    zero, zero)
  }
  else {
    args$parameters[[data.term.indicator]] <- cbind(one)
  }
}
if (length(unconditional) > 0) {
  if (is.null(data.term.indicator))
    stop("Failed to disable some data terms (because 'data.term.indicator' missing)")
  args$parameters[[data.term.indicator]][unconditional, 1] <- 0
}
if (length(conditional) > 0) {
  if (is.null(data.term.indicator))
    stop("Failed to enable some data terms (because 'data.term.indicator' missing)")
  args$parameters[[data.term.indicator]][conditional, 1] <- 1
}
makeFac <- function(x) {
  fac <- as.matrix(x)
  fac[] <- 1:length(x)
  fac[conditional, ] <- NA
  fac[unconditional, ] <- NA
  fac[subset, ] <- 1:(length(subset) * ncol(fac))
  factor(fac)
}
map <- list()
map[[observation.name]] <- makeFac(obs)
if (!is.null(data.term.indicator)) {
  map[[data.term.indicator]] <- makeFac(args$parameters[[data.term.indicator]])
}
args$map <- c(args$map, map)
args$silent <- TRUE
newobj <- do.call("MakeADFun", args)
nm <- names(newobj$par)
obs.pointer <- which(nm == observation.name)
if (method == "cdf") {
  tmp <- matrix(which(nm == data.term.indicator), ncol = 3)
  data.term.pointer <- tmp[, 1]
  lower.cdf.pointer <- tmp[, 2]
  upper.cdf.pointer <- tmp[, 3]
} else {
  data.term.pointer <- which(nm == data.term.indicator)
  lower.cdf.pointer <- NULL
  upper.cdf.pointer <- NULL
}
observation <- local({
  obs.local <- newobj$par
  i <- 1:length(subset)
  function(k, y = NULL, lower.cdf = FALSE, upper.cdf = FALSE) {
    obs.local[data.term.pointer[k < i]] <- 0
    if (!is.null(y))
      obs.local[obs.pointer[k]] <- y
    if (lower.cdf | upper.cdf) {
      obs.local[data.term.pointer[k]] <- 0
      if (lower.cdf)
        obs.local[lower.cdf.pointer[k]] <- 1
      if (upper.cdf)
        obs.local[upper.cdf.pointer[k]] <- 1
    }
    obs.local
  }
})
tracefun <- function(k) if (trace) print(k)
applyMethod <- function(oneStepMethod) {
  ord <- seq_along(subset)
  if (reverse)
    ord <- rev(ord)
  pred <- do.call("rbind", lapply(ord, oneStepMethod))
  pred <- as.data.frame(pred)[ord, ]
  pred$Fx <- 1/(1 + exp(pred$nlcdf.lower - pred$nlcdf.upper))
  pred$px <- 1/(exp(-pred$nlcdf.lower + pred$nll) + exp(-pred$nlcdf.upper + pred$nll))
  if (discrete) {
    if (!is.null(seed)) {
      Random.seed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- Random.seed)
      set.seed(seed)
    }
    U <- runif(nrow(pred))
  }
  else {
    U <- 0
  }
  pred$residual <- qnorm(pred$Fx - U * pred$px)
  pred
}

## Done preparing object
##  if (method == "oneStepGaussian") {
p <- newobj$par
newobj$fn(p)
oneStepGaussian <- function(k) {
  tracefun(k)
  index <- subset[k]
  f <- function(y) {
    newobj$fn(observation(k, y))
  }
  g <- function(y) {
    newobj$gr(observation(k, y))[obs.pointer[k]]
  }
  cbind(obs, observation(k,obs[index]) %>% as.numeric)
  opt <- nlminb(obs[index], f, g)
  H <- optimHess(opt$par, f, g)
  c(observation = obs[index], mean = opt$par, sd = sqrt(1/H))
}
ord <- seq_along(subset)
if (reverse) ord <- rev(ord)
pred <- do.call("rbind", lapply(ord, oneStepGaussian))
pred <- as.data.frame(pred)[ord, ]
pred$residual <- (pred$observation - pred$mean)/pred$sd
## }

png('plots/glmm_oneStepGauss_marginals.png', width=7, height=7,
    units='in', res=400)
par(mfrow=c(3,3), mgp=c(1.5, .5, 0), tck=-.02, mar=c(3,3,1,1))
for(k in ind){
  ## yy <- seq(-10,10, len=100)
  yy <- seq(obs[k]-.5, obs[k]+.5, len=100)
  ## print(observation(k, yy[1]))
  ff <- function(y) {
    newobj$fn(observation(k, y))
  }
  plot(yy, sapply(yy, ff) , type='l', xlab='Residual',
       ylab='NLL', main=paste('Obs=',k))
  abline(v=obs[k])
  points(pred$mean[k], ff(pred$mean[k]), cex=1.5, col=2, pch=16)
}
dev.off()
