
## Copy of oneStepPredict modified to return key info for
## visualizing the residual method rather than the residuals
## themselves. -Cole
extractOSA <- function(obj, observation.name=NULL,
                       data.term.indicator=NULL,
                       method=c("oneStepGaussianOffMode", "fullGaussian",
                                "oneStepGeneric", "oneStepGaussian", "cdf"),
                       subset=NULL, conditional=NULL, discrete=NULL,
                       discreteSupport=NULL, range=c(-Inf, Inf), seed=123,
                       parallel=FALSE, trace=TRUE,
                       reverse=(method == "oneStepGaussianOffMode"), ...){
  if (missing(observation.name))  stop("'observation.name' must define a data component")
  if (!(observation.name %in% names(obj$env$data))) stop("'observation.name' must be in data component")
  method <- match.arg(method)
  if (is.null(data.term.indicator)) {
    if (method != "fullGaussian") {
      stop(paste0("method='", method, "' requires a 'data.term.indicator'"))
    }
  }
  if (!missing(discreteSupport) && !missing(range))
    stop("Cannot specify both 'discreteSupport' and 'range'")
  obs <- as.vector(obj$env$data[[observation.name]])
  if (is.null(discrete)) {
    ndup <- sum(duplicated(obs))
    if (ndup > 0) {
      warning("Observations do not look continuous. Number of duplicates = ", ndup)
      stop("Argument 'discrete' (TRUE/FALSE) must be specified.")
    }
    discrete <- FALSE
  }
  else {
    stopifnot(is.logical(discrete))
  }
  if (discrete) {
    if (!(method %in% c("oneStepGeneric", "cdf"))) {
      stop(paste0("method='", method, "' is not for discrete observations."))
    }
  }
  if (is.null(subset)) {
    subset <- 1:length(obs)
    subset <- setdiff(subset, conditional)
  }
  if (!is.null(conditional)) {
    if (length(intersect(subset, conditional)) > 0) {
      stop("'subset' and 'conditional' have non-empty intersection")
    }
  }
  unconditional <- setdiff(1:length(obs), union(subset, conditional))
  args <- as.list(obj$env)[intersect(names(formals(MakeADFun)),
                                     ls(obj$env))]
  if (length(obj$env$random))
    args$parameters <- obj$env$parList(par = obj$env$last.par.best)
  else args$parameters <- obj$env$parList(obj$env$last.par.best)
  names.random <- unique(names(obj$env$par[obj$env$random]))
  names.all <- names(args$parameters)
  fix <- setdiff(names.all, names.random)
  map <- lapply(args$parameters[fix], function(x) factor(x *
                                                         NA))
  ran.in.map <- names.random[names.random %in% names(args$map)]
  if (length(ran.in.map))
    map <- c(map, args$map[ran.in.map])
  args$map <- map
  args$random <- names.random
  args$regexp <- FALSE
  args$parameters[observation.name] <- args$data[observation.name]
  args$data[observation.name] <- NULL
  if (!is.null(data.term.indicator)) {
    one <- rep(1, length(obs))
    zero <- rep(0, length(obs))
    if (method == "cdf") {
      args$parameters[[data.term.indicator]] <- cbind(one, zero, zero)
    } else {
      args$parameters[[data.term.indicator]] <- cbind(one)
    }
  }
  if (length(unconditional) > 0) {
    if (is.null(data.term.indicator))
      stop("Failed to disable some data terms (because 'data.term.indicator' missing)")
    args$parameters[[data.term.indicator]][unconditional,
                                           1] <- 0
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
  }
  else {
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
  if (parallel) {
    nthreads.restore <- TMB::openmp()
    on.exit(TMB::openmp(nthreads.restore), add = TRUE)
    TMB::openmp(1)
    requireNamespace("parallel")
    lapply <- parallel::mclapply
  }
  tracefun <- function(k) if (trace) print(k)
  applyMethod <- function(oneStepMethod) {
    ord <- seq_along(subset)
    if (reverse)
      ord <- rev(ord)
    pred <- do.call("rbind", lapply(ord, oneStepMethod))
    pred <- as.data.frame(pred)[ord, ]
    pred$Fx <- 1/(1 + exp(pred$nlcdf.lower - pred$nlcdf.upper))
    pred$px <- 1/(exp(-pred$nlcdf.lower + pred$nll) +
                  exp(-pred$nlcdf.upper + pred$nll))
    if(discrete){
      if(!is.null(seed)){
        ## Restore RNG on exit:
        Random.seed <- .GlobalEnv$.Random.seed
        on.exit(.GlobalEnv$.Random.seed <- Random.seed)
        set.seed(seed)
      }
      intidx <- which(round(obs) == obs)
      U <- rep(0, nrow(pred))
      U[intidx] <- runif(length(intidx))
    } else {
      U <- 0
    }
    pred$residual <- qnorm(pred$Fx - U * pred$px)
    pred
  }
  if (method == "oneStepGaussian") {
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
      opt <- nlminb(obs[index], f, g)
      H <- optimHess(opt$par, f, g)
      c(observation = obs[index], mean = opt$par, sd = sqrt(1/H))
    }
    ord <- seq_along(subset)
    if (reverse) ord <- rev(ord)
    pred <- do.call("rbind", lapply(ord, oneStepGaussian))
    pred <- as.data.frame(pred)[ord, ]
    pred$residual <- (pred$observation - pred$mean)/pred$sd
    ## Cole added this bit
    LLcurves <- LLpoints <- list()
    for(k in 1:length(obs)){
      yy <- seq(obs[k]-2, obs[k]+2, len=500)
      ff <- function(y) newobj$fn(observation(k, y))
      LLcurves[[k]] <- data.frame(k=k, x=yy, y=sapply(yy, ff) )
      LLpoints[[k]] <- data.frame(k=k, x=pred$mean[k], y=ff(pred$mean[k]))
    }
    LLcurves <- do.call(rbind, LLcurves)
    LLpoints <- do.call(rbind, LLpoints)
    out <- list(method=method, LLcurves=LLcurves,
                LLpoints=LLpoints, pred=pred)
  }
  if (method == "oneStepGaussianOffMode") {
    p <- newobj$par
    newobj$fn(p)
    newobj$env$random.start <- expression({
      last.par[random]
    })
    oneStepGaussian <- function(k) {
      tracefun(k)
      index <- subset[k]
      f <- function(y) {
        newobj$fn(observation(k, y))
      }
      g <- function(y) {
        newobj$gr(observation(k, y))[obs.pointer[k]]
      }
      c(observation = obs[index], nll = f(obs[index]),
        grad = g(obs[index]))
    }
    ord <- seq_along(subset)
    if (reverse) ord <- rev(ord)
    pred <- do.call("rbind", lapply(ord, oneStepGaussian))
    pred <- as.data.frame(pred)[ord, ]
    W <- function(x) {
      rel.tol <- sqrt(.Machine$double.eps)
      logx <- log(x)
      fdivg <- function(y) (y - exp(logx - y))/(1 + y)
      y <- pmax(logx, 0)
      while (any(abs(logx - log(y) - y) > rel.tol, na.rm = TRUE)) {
        y <- y - fdivg(y)
      }
      y
    }
    getResid <- function(value, grad) {
      Rabs <- sqrt(W(exp(2 * (value - log(sqrt(2 * pi)) +
                              log(abs(grad))))))
      R <- sign(grad) * Rabs
      R
    }
    nll0 <- newobj$fn(observation(0))
    R <- getResid(diff(c(nll0, pred$nll)), pred$grad)
    M <- pred$observation - ifelse(pred$grad != 0, R * (R/pred$grad), 0)
    pred$mean <- M
    pred$residual <- R
  }
  if ((method == "oneStepGeneric") && missing(discreteSupport)) {
    p <- newobj$par
    newobj$fn(p)
    newobj$env$value.best <- -Inf
    nan2zero <- function(x) if (!is.finite(x)) 0 else x
    formals(tmbprofile)$ytol <- 10
    formals(tmbprofile)$ystep <- 0.5
    if (discrete) {
      formals(tmbprofile)$h <- 1
      integrate <- function(f, lower, upper, ...) {
        grid <- ceiling(lower):floor(upper)
        list(value = sum(f(grid)))
      }
    }
    oneStepGeneric <- function(k) {
      tracefun(k)
      ans <- try({
        index <- subset[k]
        f <- function(y) {
          newobj$fn(observation(k, y))
        }
        nll <- f(obs[index])
        newobj$env$last.par.best <- newobj$env$last.par
        slice <- tmbprofile(newobj, k, slice = TRUE, parm.range = range, ...)
        spline <- splinefun(slice[[1]], slice[[2]])
        spline.range <- range(slice[[1]])
        if (trace >= 2) {
          plotfun <- function(slice, spline) {
            plot(slice, type = "p", level = NULL)
            plot(spline, spline.range[1], spline.range[2],
                 add = TRUE)
            abline(v = obs[index], lty = "dashed")
          }
          if (trace >= 3) {
            slice$value <- exp(-(slice$value - nll))
            plotfun(slice, function(x) exp(-(spline(x) -
                                             nll)))
          }
          else plotfun(slice, spline)
        }
        F1 <- integrate(function(x)
          exp(-(spline(x) - nll)), spline.range[1], obs[index])$value
        F2 <- integrate(function(x)
          exp(-(spline(x) - nll)), obs[index] + discrete, spline.range[2])$value
        mean <- integrate(function(x)
          exp(-(spline(x) - nll)) * x, spline.range[1], spline.range[2])$value/(F1 + F2)
        nlcdf.lower = nll - log(F1)
        nlcdf.upper = nll - log(F2)
        c(nll = nll, nlcdf.lower = nlcdf.lower, nlcdf.upper = nlcdf.upper,
          mean = mean)
      })
      if (is(ans, "try-error"))
        ans <- NaN
      ans
    }
    pred <- applyMethod(oneStepGeneric)
  }
  if ((method == "oneStepGeneric") && !missing(discreteSupport)) {
    p <- newobj$par
    newobj$fn(p)
    obs <- as.integer(round(obs))
    if (is.null(discreteSupport)) {
      warning("Setting 'discreteSupport' to ", min(obs), ":", max(obs))
      discreteSupport <- min(obs):max(obs)
    }
    oneStepDiscrete <- function(k) {
      tracefun(k)
      ans <- try({
        index <- subset[k]
        f <- function(y) {
          newobj$fn(observation(k, y))
        }
        nll <- f(obs[index])
        F <- Vectorize(function(x) exp(-(f(x) - nll)))(discreteSupport)
        F1 <- sum(F[discreteSupport <= obs[index]])
        F2 <- sum(F[discreteSupport > obs[index]])
        nlcdf.lower = nll - log(F1)
        nlcdf.upper = nll - log(F2)
        c(nll = nll, nlcdf.lower = nlcdf.lower, nlcdf.upper = nlcdf.upper)
      })
      if (is(ans, "try-error")) ans <- NaN
      ans
    }
    pred <- applyMethod(oneStepDiscrete)
  }
  if (method == "fullGaussian") {
    browser()
    args2 <- args
    args2$random <- c(args2$random, observation.name)
    fix <- data.term.indicator
    args2$map[fix] <-
      lapply(args2$map[fix], function(x) factor(NA * unclass(x)))
    newobj2 <- do.call("MakeADFun", args2)
    newobj2$fn()
    mode <- newobj2$env$last.par
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
    h <- newobj2$env$spHess(mode, random = TRUE)
    i <- which(names(newobj2$env$par[newobj2$env$random]) == observation.name)
    Sigma <- solve(as.matrix(GMRFmarginal(h, i)))
    res <- obs[subset] - mode[i]
    L <- t(chol(Sigma))
    pred <- data.frame(residual = as.vector(solve(L, res)))
    out <- list(method=method, obs=obs[subset], mode=mode[i], Sigma=Sigma, pred=pred)
  }
  if (method == "cdf") {
    p <- newobj$par
    newobj$fn(p)
    cdf <- function(k) {
      tracefun(k)
      nll <- newobj$fn(observation(k))
      nlcdf.lower <- newobj$fn(observation(k, lower.cdf = TRUE))
      nlcdf.upper <- newobj$fn(observation(k, upper.cdf = TRUE))
      c(nll = nll, nlcdf.lower = nlcdf.lower, nlcdf.upper = nlcdf.upper)
    }
    pred <- applyMethod(cdf)
  }
  if(!exists('out')) stop("Function not setup to work with method ", method)
  return(out)
}




extractOSAcovar <- function(obj, observation.name){
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
  ### rotate then subset
  Sigma <- solve(as.matrix(newobj2$env$spHess(mode, random = TRUE)))
  i <- which(names(newobj2$env$par[newobj2$env$random]) == observation.name)
  return(Sigma[i,i])
}
