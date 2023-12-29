## Unit tests for functions from sim_data.R
source('../../R/sim_data.R')
source('../../R/model_fns.R')

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

linmod.par <- setup_trueparms(mod ='linmod', 
                              misp='overdispersion', 
                              fam = NULL, link = NULL)

old.fn <- simulate_linmod(seed = 123, n = 50)
y0 <-  old.fn$Data$y0
y1 <-  old.fn$Data$y1
new.fn <- simdat(n=50, mod='linmod', cov.mod = "norm",
                 trueparms = linmod.par, 
                 misp='overdispersion', seed=123)


test_that('linmod, misp=overdispersion',{
  expect_equal(old.fn$Data$x, new.fn$x[,2])
  expect_equal(y0, as.vector(new.fn$y0[[1]]))
  expect_equal(y1, as.vector(new.fn$y1[[1]]))
})


true.par <- list(
  h0 = list(beta = c(4,-5),
            ln_sig_y = 0 ),
  h1 = list(beta = c(4,-5),
            ln_sig_y = 0 )
)
est.par <- list(
  h0 = list(beta = c(0,0),
            ln_sig_y = 0 ),
  h1 = list(beta = c(0,0),
            ln_sig_y = 0 )
)

linmod.dat <- mkTMBdat(Dat = new.fn, Pars = linmod.pars, 
                         Mod="linmod", Misp ="overdispersion")
model.true.par <-  mkTMBpar(Pars = linmod.par, 
                               Dat = linmod.dat, 
                               Mod = "linmod", Misp = "overdispersion", 
                               doTrue = TRUE)
model.est.par <-  mkTMBpar(Pars = linmod.par, 
                               Dat = linmod.dat, 
                               Mod = "linmod", Misp = "overdispersion", 
                               doTrue = FALSE)
test_that("linmod, mkpar", {
  expect_equal(true.par, model.true.par)
  expect_equal(est.par, model.est.par)
})

model.map <- mkTMBmap(Pars = model.true.par, Mod = "linmod", 
                      Misp = "overdispersion", 
                      Type = "LMM")
test_that("linmod, map", {
  expect_equal(model.map, list(h0 = list(), h1 = list()))
})

model.random <- mkTMBrandom(Mod = "linmod", Misp = "overdispersion")
test_that("linmod, random", {
  expect_equal(model.random, list(h0 = NULL, h1 = NULL))
})
