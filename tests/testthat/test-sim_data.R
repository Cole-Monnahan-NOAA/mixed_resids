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
context("general simdat tests")
test_that('general simdat tests',{
  expect_error(new.fn <- simdat(n=50, mod='linmod', 
                                trueparms = list(theta=c(4,-5), sd.vec=1, sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                                misp='overdisp', seed=123))
  expect_error(new.fn <- simdat(n=50, mod='linmodel',
                                trueparms = list(theta=c(4,-5), sd.vec=1, sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                                misp='overdispersion', seed=123))
})

context("simdat linmod tests")
simulate_linmod <- function(seed, n){
  intercept <- 4
  slope <- -5
  sigma <- 1
  set.seed(seed)
  x <- rnorm(n)
  set.seed(seed*2)
  eps <- rnorm(n, 0, sigma)
  y <- eps + intercept+slope*x
  Data <- list(y=y, x=x)
  Par <- list(b0=intercept, b1=0, logsigma=0)
  return(list(Data=Data, Par=Par))
}

old.fn <- simulate_linmod(123,50)
y0 <-  old.fn$Data$y
new.fn <- simdat(n=50, mod='linmod', 
                 trueparms = list(theta=c(4,-5), sd.vec=1, sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                 misp='outliers', seed=123)

set.seed(123)
noutlier <- 5
y1 <- y0
ind <- sample(1:50, size=noutlier)
y1[ind] <- y0[ind]+rnorm(noutlier, 0, 4)

test_that('linmod, misp=outliers',{
  expect_equal(old.fn$Data$x, new.fn$x[,2])
  expect_equal(y0, as.vector(new.fn$y0))
  expect_equal(y1, as.vector(new.fn$y1))
})

new.fn <- simdat(n=50, mod='linmod', 
                 trueparms = list(theta=c(4,-5), sd.vec=1, sp.parm=0, sp.fam = NULL, sp.link = NULL),
                 misp='miss.cov', seed=123)
test_that('linmod, misp=miss.cov',{
  expect_equal(old.fn$Data$x, new.fn$x[,2])
  expect_equal(y0, as.vector(new.fn$y0))
  expect_equal(y0, as.vector(new.fn$y1))
})

set.seed(123*3)
y0 <- old.fn$Dat$y
y1 <- y0* exp(rnorm(50,0,log(4)))
new.fn <- simdat(n=50, mod='linmod', 
                 trueparms = list(theta=c(4,-5), sd.vec=c(1,log(4)), sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                 misp='overdispersion', seed=123)
test_that('linmod, misp=overdispersion',{
  expect_equal(old.fn$Data$x, new.fn$x[,2])
  expect_equal(y0, as.vector(new.fn$y0))
  expect_equal(y1, as.vector(new.fn$y1))
})

test_that('linmod, invalid misp',{
  expect_error(simdat(n=50, mod='linmod', 
                      trueparms = list(theta=c(4,-5), sd.vec=1, sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                      misp='mu0', seed=123))
  expect_error(simdat(n=50, mod='linmod', 
                      trueparms = list(theta=c(4,-5), sd.vec=1, sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                      misp='misp.omega', seed=123))
})
rm(y0,y1,old.fn,new.fn,ind,noutlier,simulate_linmod)

context("simdat randomwalk tests")
simulate_randomwalk <- function(seed,n){
  set.seed(seed)
  mu <- .75
  sigma <- 1
  tau <- 1
  huge <- 1e3
  ## simulate random track
  X <- rnorm(n,mean=0,sd=tau)
  Ypred <- rep(NA, n)
  Ypred[1] <- X[1]
  for(t in 2:n) Ypred[t] <- Ypred[t-1]+X[t]+mu
  ## simulate random measurements
  Y <- rnorm(n, Ypred, sd=sigma)
  out <- list(y=Y,x=X)
}
old.fn <- simulate_randomwalk(123,50)
new.fn <- simdat(n=50, mod='randomwalk', 
                 trueparms = list(theta=.75, sd.vec=c(1,1), sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                 misp='mu0', seed=123)
test_that('randomwalk, misp=mu0',{
  expect_equal(matrix(1,50,1), new.fn$x)
  expect_equal(old.fn$y, new.fn$y0)
  expect_equal(old.fn$y, new.fn$y1)
})

y0 <- old.fn$y
set.seed(123)
noutlier <- 5
y1 <- y0
ind <- sample(1:50, size=noutlier)
y1[ind] <- y0[ind]+rnorm(noutlier, 0, 4)
new.fn <- simdat(n=50, mod='randomwalk', 
                 trueparms = list(theta=.75, sd.vec=c(1,1), sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                 misp='outliers', seed=123)
test_that('randomwalk, misp=outliers',{
  expect_equal(y0, new.fn$y0)
  expect_equal(y1, new.fn$y1)
})
test_that('randomwalk, invalid misp',{
  expect_error(simdat(n=50, mod='randomwalk', 
                      trueparms = list(theta=.75, sd.vec=c(1,1), sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                      misp='miss.cov', seed=123))
  expect_error(simdat(n=50, mod='randomwalk', 
                      trueparms = list(theta=.75, sd.vec=c(1,1), sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                      misp='overdispersion', seed=123))
  expect_error(simdat(n=50, mod='randomwalk', 
                      trueparms = list(theta=.75, sd.vec=c(1,1), sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                      misp='misp.omega', seed=123))
})

rm(y0,y1,old.fn,new.fn,ind,noutlier,simulate_randomwalk)


context("simdat simpleGLMM tests")
simulate_simpleGLMM <- function(seed, n=10, ng=3, cov.mod='rnorm', misp){
  ## ngroups <- 3 #number of subjects
  ## nobs <- 10 #number of observations
  b0 <- 4
  sig2.y <- .5 #obs variance
  #groups are being simulated with high overlap - hard for model to differentiate
  #need to implement sum-to-zero constraint in u
  sig2.u <- 10 # between group variance
  set.seed(seed)

  u <- rnorm(ng, 0, sqrt(sig2.u))
  v <- rep(0,n)
  if(misp == 'overdispersion'){
    set.seed(seed*3)
    v <- rnorm(n)
  }
  y <- matrix(0, n, ng)
  set.seed(seed*2)
  for(j in 1:ng){
    y[,j] <- rnorm(n, b0 + u[j] + v, sqrt(sig2.y))
  }
  ## boxplot(y)
  Dat <- data.frame(y = as.vector(y), group = rep(1:ng, each = n))
  Data <- list(y = Dat[,1], group = Dat[,2]-1, u = u, sim_re = 0)
  Par <- list(b0=0, ln_sig_u=0, ln_sig_y=0, u=rep(0, ng))
  return(list(Data=Data, Par=Par))
}
old.fn <- simulate_simpleGLMM(123,10,5,misp = 'outliers')
new.fn <- simdat(n=10, ng=5, mod='simpleGLMM', 
                 trueparms = list(theta=4, sd.vec=sqrt(c(.5,10)), sp.parm=0, 
                                  fam = 'Gaussian', link = 'identity'), 
                 misp='outliers', seed=123)
y0 <- old.fn$Data$y
set.seed(123)
noutlier <- 5
y1 <- y0
ind <- sample(1:50, size=noutlier)
y1[ind] <- y0[ind]+rnorm(noutlier, 0, 4*sqrt(.5))

test_that('simpleGLMM, misp=outliers',{
  expect_equal(y0, new.fn$y0)
  expect_equal(y1, new.fn$y1)
})

old.fn <- simulate_simpleGLMM(123,10,5,misp = 'overdispersion')
y0 <- y1 <- old.fn$Data$y
new.fn <- simdat(n=10, ng=5, mod='simpleGLMM', 
                 trueparms = list(theta=4, sd.vec=c(sqrt(c(.5,10)),1), sp.parm = 0, 
                                  fam = 'Gaussian', link = 'identity'), 
                 misp='overdispersion', seed=123)
test_that('simpleGLMM, misp=overdispersion',{
  expect_equal(y0, as.vector(new.fn$y0))
  expect_equal(y1, as.vector(new.fn$y1))
})

set.seed(123)
X <- cbind(rep(1,10), rnorm(10))
mu <- X%*%c(4,-5)
y0 <- matrix(0, 10, 5)
set.seed(123*2)
for(j in 1:5){
  for(i in 1:10){
    y0[i,j] <- rnorm(1, mu[i] + old.fn$Data$u[j], sd=sqrt(0.5))
  }
}
y0 <- as.vector(y0)
y1 <- y0

new.fn <- simdat(n=10, ng=5, mod='simpleGLMM', 
                 trueparms = list(theta=c(4,-5), sd.vec=sqrt(c(.5,10)), sp.parm=0, fam = 'Gaussian', sp.link = 'identity'), 
                 misp='miss.cov', seed=123)
test_that('simpleGLMM, misp=miss.cov',{
  expect_equal(X[,2], as.vector(new.fn$x[,2]))
  expect_equal(y0, as.vector(new.fn$y0))
  expect_equal(y1, as.vector(new.fn$y0))
  expect_equal(y1, as.vector(new.fn$y1))
})

test_that('simpleGLMM, invalid misp',{
  expect_error(simdat(n=10, ng=5, mod='simpleGLMM', 
                      trueparms = list(theta=4, sd.vec=sqrt(c(.5,10)), sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                      misp='mu0', seed=123))
  expect_error(simdat(n=10, ng=5, mod='simpleGLMM', 
                      trueparms = list(theta=4, sd.vec=sqrt(c(.5,10)), sp.parm=0, sp.fam = NULL, sp.link = NULL), 
                                       misp='misp.omega', seed=123))
})
rm(y0,y1,old.fn,new.fn,ind,noutlier,simulate_simpleGLMM)


context('simdat spatial tests')
simulate_spatial <- function(seed,n){#},Fam,Link, misp){
  library(TMB)
  dyn.load(dynlib('../../src/spatial'))
  set.seed(seed)
  sp.var <- 0.5
  sd.obs <- .5
  Range <- 20
  ## Simulate spatial random effects
  Loc <- matrix(runif(n*2,0,n),ncol=2)
  dmat <- as.matrix(dist(Loc))
  mesh <- try(
    R.utils::withTimeout( INLA::inla.mesh.2d(Loc, max.edge = c(Range/3, Range), offset = c(2, Range*.75)),
                 timeout = 30, onTimeout = 'silent' ))
  if(is.character(mesh)){
    system("Taskkill /IM fmesher.exe /F")
   # warning("mesh failed in rep=", ii)
    return(NULL)
  }
  Omega <- sim_omega(Range,sp.var,dmat,method="TMB.spde",mesh=mesh)
  return(Omega[mesh$idx$loc])
}


old.fn <- simulate_spatial(123,50)#,'Poisson','log', misp='outliers')
X <- matrix(1,nrow=50)
mu <- X%*%2
set.seed(123)
y <- sim_y(Eta=mu, omega=old.fn,
            parm=0.5, fam='Poisson', link='log')
new.fn <- simdat(n=50, mod='spatial', 
                 trueparms = list(theta=2, sd.vec=c(.5,sqrt(0.5)), 
                                  sp.parm = 20, fam = 'Poisson', link = 'log'),
                 misp='outliers', seed=123)
y0 <- y
set.seed(123)
noutlier <- 5
y1 <- y0
ind <- sample(1:50, size=noutlier)
y1[ind] <- y0[ind]+rnorm(noutlier, 0, 4*.5)
test_that('spatial, misp=outliers',{
  expect_equal(X,new.fn$x)
  expect_equal(y0, new.fn$y0)
  expect_equal(y1, new.fn$y1)
})

set.seed(123)
X <- cbind(rep(1,50), rnorm(50))
mu <- X%*%c(1,2)
set.seed(123)
y <- sim_y(Eta=mu, omega=old.fn,
           parm=0.5, fam='Poisson', link='log')
new.fn <- simdat(n=50, mod='spatial', 
                 trueparms = list(theta=c(1,2), sd.vec=c(.5,sqrt(0.5)), 
                                  sp.parm = 20, fam = 'Poisson', link = 'log'),
                 misp='miss.cov', seed=123)
test_that('spatial, misp=miss.cov',{
  expect_equal(X,new.fn$x)
  expect_equal(y, new.fn$y0)
  expect_equal(y, new.fn$y1)
})

X <- matrix(1,nrow=50)
set.seed(123)
mu <- X%*%2 + rnorm(50,0,0.5*2)
set.seed(123)
y <- sim_y(Eta=mu, omega=old.fn,
           parm=0.5, fam='Poisson', link='log')

new.fn <- simdat(n=50, mod='spatial', 
                 trueparms = list(theta=2, sd.vec=c(.5,sqrt(0.5),0.5*2),
                                  sp.parm = 20, fam = 'Poisson', link = 'log'),
                 misp='overdispersion', seed=123)
test_that('spatial, misp=miss.cov',{
  expect_equal(X,new.fn$x)
  expect_equal(y, new.fn$y0)
  expect_equal(y, new.fn$y1)
})

X <- matrix(1,nrow=50)
mu <- X%*%2
set.seed(123)
y0 <- sim_y(Eta=mu, omega=old.fn,
           parm=0.5, fam='Poisson', link='log')
set.seed(123)
y1 <- sim_y(Eta=mu, omega=exp(old.fn),
            parm=0.5, fam='Poisson', link='log')
new.fn <- simdat(n=50, mod='spatial', 
                 trueparms = list(theta=2, sd.vec=c(.5,sqrt(0.5)), 
                                  sp.parm = 20, fam = 'Poisson', link = 'log'),
                 misp='misp.omega', seed=123)
test_that('spatial, misp=miss.cov',{
  expect_equal(X,new.fn$x)
  expect_equal(y0, new.fn$y0)
  expect_equal(y1, new.fn$y1)
})
test_that('spatial, invalid misp',{
  expect_error(simdat(n=50, mod='spatial', 
                      trueparms = list(theta=2, sd.vec=c(.5,sqrt(0.5)), sp.parm = 20, sp.fam = 'Poisson', sp.link = 'log'),
                      misp='mu0', seed=123))
})
