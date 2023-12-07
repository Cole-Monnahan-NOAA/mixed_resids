## Unit tests for functions from sim_data.R
source('../../R/sim_data.R')

## test cMatern using function defined in geoR::matern
context("cMatern unit test")
geoRmatern <- function (u, phi, kappa) 
{
  if (is.vector(u)) 
    names(u) <- NULL
  if (is.matrix(u)) 
    dimnames(u) <- list(NULL, NULL)
  uphi <- u/phi
  uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/ifelse(0, Inf, 
                                                    gamma(kappa))) * (uphi^kappa) * besselK(x = uphi, nu = kappa)), 
                 1)
  uphi[u > 600 * phi] <- 0
  return(uphi)
}
Range <- 20; Nu = 1
test_that("cMatern",{
  expect_equal(geoRmatern(0,Range,Nu), cMatern(0,Nu,1/Range))
  expect_equal(geoRmatern(5,Range,Nu), cMatern(5,Nu,1/Range))
  expect_equal(geoRmatern(1000,Range,Nu), cMatern(1000,Nu,1/Range))
})

## test sim_data.R
context("simdat linmod tests")
simulate_linmod <- function(seed, n){
  intercept <- 4
  slope <- -5
  set.seed(seed)
  x <- rnorm(n)
  set.seed(seed*2)
  eps <- rnorm(n)
  y0 <- eps + intercept+slope*x
  set.seed(seed*3)
  y1 <- y0*exp(rnorm(n))
  Data <- list(y0 = y0, y1 = y1, x=x)
  Par <- list(b0=intercept, b1=0, logsigma=0)
  return(list(Data=Data, Par=Par))
}

old.fn <- simulate_linmod(seed = 123, n = 50)
y0 <-  old.fn$Data$y0
y1 <-  old.fn$Data$y1
new.fn <- simdat(n=50, mod='linmod', cov.mod = "norm",
                 trueparms = list(theta=c(4,-5), sd.vec=c(1,1), fam = NULL, link = NULL), 
                 misp='overdispersion', seed=123)


test_that('linmod, misp=overdispersion',{
  expect_equal(old.fn$Data$x, new.fn$x[,2])
  expect_equal(y0, as.vector(new.fn$y0[[1]]))
  expect_equal(y1, as.vector(new.fn$y1[[1]]))
})
rm(y0,y1,old.fn,new.fn,simulate_linmod)

context("simdat randomwalk tests")
simulate_randomwalk <- function(seed, n, type, sd.vec, init_u){
  set.seed(seed)
  if(type == "LMM") mu <- 2
  if(type == "GLMM") mu <- .1
  sig_y <- sd.vec[1]
  sig_u <- sd.vec[2]
  huge <- 1e3
  ## simulate random track
  set.seed(seed)
  u <- rnorm(n-1,mean=0,sd=sig_u)
  Ypred <- rep(NA, n)
  Ypred[1] <- init_u
  for(t in 2:n) Ypred[t] <- Ypred[t-1]+u[t-1]+mu
  ## simulate random measurements
  set.seed(seed)
  if(type == "LMM"){
    Y <- rnorm(n, Ypred, sd=sig_y)
  } 
  if(type == "GLMM"){
    Y <- rgamma(n, 1/sig_y^2, scale = exp(Ypred)*sig_y^2)
  }
  out <- list(y=Y,u=Ypred)
}
old.fn <- simulate_randomwalk(123,50,"LMM", sd.vec = c(1,1), init_u = 5)
new.fn <- simdat(n=50, mod='randomwalk', 
                 trueparms = list(theta=2, sd.vec=c(1,1), init.u = 5, 
                 fam = "Gaussian", link = "identity"), 
                 misp=c('missre', 'normal-lognorm', 'mu0'), seed=123)
test_that('randomwalk, LMM',{
  expect_equal(old.fn$y, new.fn$y0)
  expect_equal(old.fn$y, new.fn$y1[[1]])
  expect_equal(old.fn$y, new.fn$y1[[2]])
  expect_equal(old.fn$y, new.fn$y1[[3]])
})
min.y <- rep(-999, 1000)
for(i in 1:1000){
  old.fn <- simulate_randomwalk(i,50,"LMM", sd.vec = c(1,1), init_u = 5)
  min.y[i] <- min(old.fn$y)
}
test_that('randomwalk, "min y',{
  expect_equal(TRUE, min(min.y) > 0 )
} )

old.fn <- simulate_randomwalk(123, 50, "GLMM", sd.vec = c(.5,.05), init_u = 1)
new.fn <- simdat(n=50, mod='randomwalk', 
                 trueparms = list(theta=.1, sd.vec=c(.5,.05), init.u = 1, fam = "Gamma", link = "log"), 
                 misp=c('missre', 'normal-lognorm', 'mu0'), seed=123)
test_that('randomwalk, misp=mu0',{
  expect_equal(old.fn$y, new.fn$y0)
  expect_equal(old.fn$y, new.fn$y1[[1]])
  expect_equal(old.fn$y, new.fn$y1[[2]])
  expect_equal(old.fn$y, new.fn$y1[[3]])
})
min.y <- max.y <- rep(-999,1000)
for(i in 1:1000){
  old.fn <- simulate_randomwalk(i,50,"GLMM", sd.vec = c(.5,.05), init_u = 1)
  min.y[i] <- min(old.fn$y)
  max.y[i] <- max(old.fn$y)
}
test_that('randomwalk, "min y',{
  expect_equal(TRUE, min(min.y) > 0 )
  expect_equal(TRUE, max(max.y) < 3000 )
} )


rm(old.fn,new.fn,min.y,max.y,simulate_randomwalk)


context("simdat simpleGLMM tests")
simulate_simpleGLMM <- function(seed, n=100, ng=5, 
  cov.mod=NULL, type){
  
  ## ngroups <- 5 #number of subjects
  ## nobs <- 100 #number of observations
  if(type == "LMM"){
    beta <- c(4,-8)
    sig.y <- 0.5
    sig.u <- 2
  }
  if(type == "GLMM"){
    beta <- log(.5)
    size <- .6
    sig.u <- .8 #sig_u
    theta <- size
  }

  #Set up cpvariate
  if(!is.null(cov.mod)){
    set.seed(seed)
    if(cov.mod == 'norm') X <- cbind(rep(1,n), rnorm(n))
    if(cov.mod == 'unif') X <- cbind(rep(1,n), runif(n, -.5,.5))
  } else {
    X <- matrix(1, n, 1)
  }
  mu <- X %*% beta

  #Simulate RE
  set.seed(seed)
  u <- rnorm(ng, 0, sig.u)
 
  y0 <- y1 <- matrix(0, n, ng)
  if(type == "LMM"){
    set.seed(seed*2)
    for(j in 1:ng){
      y0[,j] <- rnorm(n, mu + u[j], sig.y)
    }
    set.seed(seed*2)
    for(j in 1:ng){
      y1[,j] <- rnorm(n, mu + exp(u[j]), sig.y)
    }
  }
  

  if(type == "GLMM"){
    set.seed(seed*2)
    for(j in 1:ng){
      y0[,j] <- rnbinom(n, size, mu = exp(mu + u[j]))
    }
    set.seed(seed*2)
    for(j in 1:ng){
      y1[,j] <- rnbinom(n, size, mu = exp(mu + exp(u[j])))
    }
  }
  out <- list(y0 = as.vector(y0), y1 = as.vector(y1), u=u)
  return(out)
}
old.fn <- simulate_simpleGLMM(seed = 123, n = 100, ng = 5,
                              cov.mod = "unif",
                              type = "LMM")

new.fn <- simdat(n=100, ng=5, mod='simpleGLMM', cov.mod = "unif",
                 trueparms = list(beta=c(4,-8), sd.vec=c(.5,2), 
                                  fam = 'Gaussian', link = 'identity'), 
                 misp = c('missre', 'missunifcov', 'mispre'), 
                 seed = 123, type = "LMM")

test_that('simpleGLMM, type = LMM',{
  expect_equal(old.fn$y0, new.fn$y0)
  expect_equal(old.fn$y0, new.fn$y1[[1]])
  expect_equal(old.fn$y0, new.fn$y1[[2]])
  expect_equal(old.fn$y1, new.fn$y1[[3]])
})
old.fn <- simulate_simpleGLMM(seed = 123, n = 100, ng = 5,
                              type = "GLMM")
new.fn <- simdat(n=100, ng=5, mod='simpleGLMM',
                 trueparms = list(beta=log(.5), sd.vec=c(NA,.8),
                                  theta = .6, 
                                  fam = 'NB', link = 'log'), 
                 misp = c('missre', 'nb-pois', 'mispre'), 
                 seed = 123, type = "GLMM")

test_that('simpleGLMM, type = GLMM',{
  expect_equal(old.fn$y0, new.fn$y0)
  expect_equal(old.fn$y0, new.fn$y1[[1]])
  expect_equal(old.fn$y0, new.fn$y1[[2]])
  expect_equal(old.fn$y1, new.fn$y1[[3]])
})

min.y0 <- max.y0 <- rep(-999,1000)
min.y1 <- max.y1 <- rep(-999,1000)
for(i in 1:1000){
  old.fn <- simulate_simpleGLMM(seed = i, n = 100, ng = 5,
                              cov.mod = "unif",
                              type = "LMM")
  min.y0[i] <- min(old.fn$y0)
  max.y0[i] <- max(old.fn$y0)
  min.y1[i] <- min(old.fn$y1)
  max.y1[i] <- max(old.fn$y1)
}
test_that('simpleGLMM, "min/max y',{
  expect_equal(TRUE, max(max.y0 - min.y0) < 100 )
  expect_equal(TRUE, max(max.y1 - min.y1) < 3000 )
} )

max.y0 <- rep(-999,1000)
max.y1 <- rep(-999,1000)
for(i in 1:1000){
  old.fn <- simulate_simpleGLMM(seed = i, n = 100, ng = 5,
                              type = "GLMM")
  max.y0[i] <- max(old.fn$y0)
  max.y1[i] <- max(old.fn$y1)
}
test_that('simpleGLMM, "min/max y',{
  expect_equal(TRUE, max(max.y0) < 100 )
  #expect_equal(TRUE, max(max.y1) < 100 )
} )

rm(min.y,max.y,old.fn,new.fn,simulate_simpleGLMM)


context('simdat spatial tests')
setup_mesh <- function(n, seed, sp.parm){
  set.seed(seed)
  if(n < 500){
    nspobs <- 500
  } else {
    nspobs <- n
  }
  #simulate spatial
  Loc <- matrix(runif(nspobs*2,0,100),ncol=2)
  samp.idx <- sample(1:nspobs, n)
  loc <- Loc[samp.idx,]
  mesh <- R.utils::withTimeout( 
      fmesher::fm_mesh_2d(loc, max.edge = c(sp.parm/3,sp.parm), 
                          offset = c(sp.parm/10,sp.parm*2), min.angle = 26),
      timeout = 30, onTimeout = 'silent' )
  if(is.null(mesh)){
    system("Taskkill /IM fmesher.exe /F")
    # warning("mesh failed in rep=", ii)
    return(NULL)
  }
  out <- list(Loc = Loc, mesh = mesh, samp.idx = samp.idx)
  return(out)
}
class.vec <- rep(NA, 1000)
for(i in 1:1000){
  if(!is.null(setup_mesh(100,i,50)$mesh)  ){
    class.vec[i] <- 1
  }
}
expect_equal(1000, sum(class.vec, na.rm = TRUE))

simulate_spatial <- function(seed, n, misp=NULL, 
                             fam = NULL, type){
  sp.parm <- 50
  if(type == "LMM"){
    beta <- 20
    sp.var <- 1
    sig.y <- 1
  }

  if(type == "GLMM"){
    if(fam == "Poisson"){
      beta <- 0.5
      sp.var <- 0.25
    }
  }

  ## Simulate spatial random effects
  sp.dat <- setup_mesh(n,seed,sp.parm)
  dmat <- as.matrix(dist(sp.dat$Loc))
  set.seed(seed)
  Omega <- sim_omega(sp.parm,sp.var,dmat,method="R.matern")
  omega0 <- Omega[sp.dat$samp.idx]

  y0 <- rep(NA, n)
  if(type == "LMM"){
    set.seed(seed)
    y0 <- rnorm(n, beta + omega0, sig.y)
    set.seed(seed)
    y1 <- rnorm(n, beta + exp(omega0), sig.y)
  }
  if(type == "GLMM"){
    set.seed(seed)
    y0 <- rpois(n, exp(beta + omega0))
    y1 <- list()
    for(m in 1:length(misp)){
      if(misp[m] == "mispre"){
        set.seed(seed)
        y1[[m]] <- rpois(n, exp(beta + exp(omega0)))
      } else {
        y1[[m]] <- y0
      }
      if(misp[m] == "pois-zip"){
        set.seed(seed)
        y1[[m]] <- y1[[m]] * rbinom(n, 1, 0.7)
      }
    }
  }
  out <- list(y0 = y0, y1 = y1, omega = omega0)
  return(out)
}


old.fn <- simulate_spatial(123, 100, type = "LMM")
new.fn <- simdat(n=100, mod='spatial', type = "LMM", 
                      trueparms = list(beta=20, sd.vec=c(1,1), 
                                       sp.parm = 50, fam = 'Gaussian', 
                                       link = 'identity'),
                      misp = c("missre", "normal-gamma", "mispre"), seed=123)

test_that('spatial, LMM',{
  expect_equal(old.fn$y0, new.fn$y0)
  expect_equal(old.fn$y0, new.fn$y1[[1]])
  expect_equal(old.fn$y0, new.fn$y1[[2]])
  expect_equal(old.fn$y1, new.fn$y1[[3]])
})

# One misp fits a gamma model to normal data,
# make sure no values are below zero in simulations,
# also no extreme values
min.y0 <- max.y0 <- 
  min.y1 <- max.y1 <-rep(NA, 1000)
for(i in 1:1000){
  old.fn <- simulate_spatial(i, 100, type = "LMM")
  min.y0[i] <- min(old.fn$y0)
  max.y0[i] <- max(old.fn$y0)
  min.y1[i] <- min(old.fn$y1)
  max.y1[i] <- max(old.fn$y1)
}
test_that('spatial, LMM, bounds',{
  expect_gt(min(min.y0), 0)
  expect_gt(min(min.y1), 0)
  expect_lt(max(max.y0), 100)
  expect_lt(max(max.y1), 100)
})


old.fn <- simulate_spatial(123, 100, 
                           misp = c("missre", "pois-zip", "mispre"), 
                           fam = "Poisson", type = "GLMM")
new.fn <- simdat(n=100, mod='spatial', type = "GLMM", 
                      trueparms = list(beta=0.5, sd.vec=c(NA,sqrt(.25)), 
                                       sp.parm = 50, fam = 'Poisson', 
                                       link = 'log'),
                      misp = c("missre", "pois-zip", "mispre"), seed=123)


test_that('spatial, GLMM',{
  expect_equal(old.fn$omega, new.fn$random$omega0)
  expect_equal(old.fn$y0, new.fn$y0)
  expect_equal(old.fn$y1[[1]], new.fn$y1[[1]])
  expect_equal(old.fn$y1[[2]], new.fn$y1[[2]])
  expect_equal(old.fn$y1[[3]], new.fn$y1[[3]])
})

# test upper bound and %zero inflation
y0.max <- rep(NA, 1000)
y1.max <- y.zip <- list()
y1.max[[1]] <- y1.max[[2]] <- y1.max[[3]] <- rep(NA, 1000)
y.zip[[1]] <- y.zip[[2]] <- rep(NA, 0)
for(i in 1:1000){
  old.fn <- simulate_spatial(i, 100, 
                           misp = c("missre", "pois-zip", "mispre"), 
                           fam = "Poisson", type = "GLMM")
  y0.max[i] <- max(old.fn$y0)
  y1.max[[1]][i] <- max(old.fn$y1[[1]])
  y1.max[[2]][i] <- max(old.fn$y1[[2]])
  y1.max[[3]][i] <- max(old.fn$y1[[3]])
  y.zip[[1]][i] <- sum(old.fn$y0 == 0)
  y.zip[[2]][i] <- sum(old.fn$y1[[2]] == 0)
}

test_that('spatial, GLMM, bounds', {
  expect_lt(max(y0.max), 10000)
  expect_lt(max(y1.max[[3]]), 10000)
  expect_lt(max(y.zip[[1]]), 55)
  expect_gt(min(y.zip[[2]] - y.zip[[1]]), 0)
})
