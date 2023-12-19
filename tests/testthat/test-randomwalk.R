source('../../R/sim_data.R')
source('../../R/model_fns.R')

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

context("LMM randomwalk tests")
testsim <- simulate_randomwalk(123,50,"LMM", sd.vec = c(1,1), init_u = 5)
rw.par <- setup_trueparms(mod ='randomwalk', 
                          misp=c('missre', 'normal-lognorm', 'mu0'), 
                          fam = "Gaussian", link = "identity",
                          type = "LMM")
modelsim <- simdat(n=50, mod='randomwalk', 
                 trueparms = rw.par, 
                 misp=c('missre', 'normal-lognorm', 'mu0'), seed=123)

modeldat <- mkTMBdat(Dat = modelsim, 
                       Pars = rw.par,  
                       Mod = "randomwalk", 
                       Misp = c('missre', 'normal-lognorm', 'mu0'), 
                       Type = "LMM")
model.true.par <- mkTMBpar(Pars = rw.par, 
                           Dat = modelsim, 
                           Mod = "randomwalk", 
                           Misp = c('missre', 'normal-lognorm', 'mu0'), 
                           doTrue = TRUE)
model.est.par <- mkTMBpar(Pars = rw.par, 
                          Dat = modelsim, 
                          Mod = "randomwalk", 
                          Misp = c('missre', 'normal-lognorm', 'mu0'), 
                          doTrue = FALSE)
model.true.map <- mkTMBmap(Pars = model.true.par, 
                      Mod = "randomwalk", 
                      Misp = c('missre', 'normal-lognorm', 'mu0'), 
                      Type = "LMM")
model.est.map <- mkTMBmap(Pars = model.est.par, 
                      Mod = "randomwalk", 
                      Misp = c('missre', 'normal-lognorm', 'mu0'), 
                      Type = "LMM")
test_that('randomwalk, LMM, simdat',{
  expect_equal(testsim$y, modelsim$y0)
  expect_equal(testsim$y, modelsim$y1[[1]])
  expect_equal(testsim$y, modelsim$y1[[2]])
  expect_equal(testsim$y, modelsim$y1[[3]])
})
min.y <- rep(-999, 1000)
for(i in 1:1000){
  testsim <- simulate_randomwalk(i,50,"LMM", sd.vec = c(1,1), init_u = 5)
  min.y[i] <- min(testsim$y)
}
test_that('randomwalk, LMM, min y',{
  expect_equal(TRUE, min(min.y) > 0 )
} )

test_that('randomwalk, LMM, mkDat, y',{
    expect_equal(modelsim$y0, modeldat$h0$y)
    expect_equal(modelsim$y1[[1]], modeldat$h1[[1]]$y)
    expect_equal(modelsim$y1[[2]], modeldat$h1[[2]]$y)
    expect_equal(modelsim$y1[[3]], modeldat$h1[[3]]$y)
})

test_that('randomwalk, LMM, mkDat, mod and simre',{
    expect_equal(0, modeldat$h0$mod)
    expect_equal(0, modeldat$h0$sim_re)
    expect_equal(0, modeldat$h1[[1]]$mod)
    expect_equal(0, modeldat$h1[[1]]$sim_re)
    expect_equal(0, modeldat$h1[[2]]$mod)
    expect_equal(0, modeldat$h1[[2]]$sim_re)
    expect_equal(0, modeldat$h1[[3]]$mod)
    expect_equal(0, modeldat$h1[[3]]$sim_re)
})

test.true.par <- list(
  h0 = list(
    mu = 2, ln_sig_y = 0,
    ln_sig_u = 0,
    u = modelsim$random$u0
  ),
  h1 = list(
    list(
      mu = 2, ln_sig_y = 0,
      ln_sig_u = numeric(0),
      u = rep(0, length(modelsim$random$u0))
    ),
    list(
      mu = 2, ln_sig_y = 0,
      ln_sig_u = 0,
      u = modelsim$random$u1[[2]]
    ),
    list(
      mu = 0, ln_sig_y = 0,
      ln_sig_u = 0,
      u = modelsim$random$u1[[3]]
    )
  )
)
test.est.par <- test.true.par
test.est.par$h0$mu <- 
  test.est.par$h1[[1]]$mu <- 
    test.est.par$h1[[2]]$mu <- 0

test.est.par$h0$u <- 
  test.est.par$h1[[2]]$u <- 
    test.est.par$h1[[3]]$u <- 
      rep(1, length(modelsim$random$u0))

test.map <- list(
  h0 = list(),
  h1 = list(
    list(u = rep(factor(NA), length(modelsim$random$u0))),
    list(),
    list(mu = factor(NA))
  )
)

test_that("randomwalk, LMM, mkpar", {
  expect_equal(test.true.par, model.true.par)
  expect_equal(test.est.par, model.est.par)
})

test_that("randomwalk, LMM, mkmap", {
  expect_equal(test.map, model.true.map)
  expect_equal(test.map, model.est.map)
})

test.random <- (list(h0 = "u", h1 = list(NULL, "u", "u")))
model.random <- mkTMBrandom("randomwalk", c('missre', 'normal-lognorm', 'mu0'))
test_that("randomwalk, LMM, mkrandom", {
  expect_equal(test.random, model.random)
})


test.true.objpar <- list(
  h0 = test.true.par$h0[1:3],
  h1 = list(
    test.true.par$h1[[1]][1:2],
    test.true.par$h1[[2]][1:3],
    test.true.par$h1[[3]][2:3]
  )
) 
test.est.objpar <- list(
  h0 = test.est.par$h0[1:3],
  h1 = list(
    test.est.par$h1[[1]][1:2],
    test.est.par$h1[[2]][1:3],
    test.est.par$h1[[2]][2:3]
  )
)
dyn.load(TMB::dynlib("../../src/randomwalk"))
test_that("randomwalk, LMM, init.obj", {
for(h in 1:4){
  if(h == 1){
    model.true.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.true.par$h0, 
                     map = model.true.map$h0, random = model.random$h0, DLL = "randomwalk")
    model.est.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.est.par$h0, 
                     map = model.est.map$h0, random = model.random$h0, DLL = "randomwalk")
    expect_equal(unlist(test.true.objpar$h0), model.true.obj$par)
    expect_equal(unlist(test.est.objpar$h0), model.est.obj$par)
  }
  if(h>1){
    model.true.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.true.par$h1[[h-1]], 
                     map = model.true.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "randomwalk")
    model.est.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.est.par$h1[[h-1]], 
                     map = model.est.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "randomwalk")
    expect_equal(unlist(test.true.objpar$h1[[h-1]]), model.true.obj$par)
    expect_equal(unlist(test.est.objpar$h1[[h-1]]), model.est.obj$par)
  }
  rm(model.true.obj, model.est.obj)
}
})


context("GLMM randomwalk tests")
testsim <- simulate_randomwalk(123, 50, "GLMM", sd.vec = c(.5,.05), init_u = 1)
rw.par <- setup_trueparms(mod ='randomwalk', 
                          misp = c('missre', 'gamma-lognorm', 'mu0'), 
                          fam = "Gamma", link = "log",
                          type = "GLMM")
modelsim <- simdat(n=50, mod='randomwalk', 
                 trueparms = setup_trueparms(mod ='randomwalk', 
                                             misp=c('missre', 'gamma-lognorm', 'mu0'), 
                                             fam = "Gamma", link = "log",
                                             type = "GLMM"), 
                 misp=c('missre', 'normal-lognorm', 'mu0'), seed=123)
modeldat <- mkTMBdat(Dat = modelsim, 
                       Pars = setup_trueparms(mod ='randomwalk', 
                                              misp = c('missre', 'gamma-lognorm', 'mu0'), 
                                              fam = "Gamma", link = "log",
                                              type = "GLMM"),  
                       Mod = "randomwalk", 
                       Misp = c('missre', 'gamma-lognorm', 'mu0'), 
                       Type = "GLMM")
model.true.par <- mkTMBpar(Pars = rw.par, 
                           Dat = modelsim, 
                           Mod = "randomwalk", 
                           Misp = c('missre', 'gamma-lognorm', 'mu0'), 
                           doTrue = TRUE)
model.est.par <- mkTMBpar(Pars = rw.par, 
                          Dat = modelsim, 
                          Mod = "randomwalk", 
                          Misp = c('missre', 'gamma-lognorm', 'mu0'), 
                          doTrue = FALSE)
model.true.map <- mkTMBmap(Pars = model.true.par, 
                      Mod = "randomwalk", 
                      Misp = c('missre', 'gamma-lognorm', 'mu0'), 
                      Type = "GLMM")
model.est.map <- mkTMBmap(Pars = model.est.par, 
                      Mod = "randomwalk", 
                      Misp = c('missre', 'gamma-lognorm', 'mu0'), 
                      Type = "GLMM")
test_that('randomwalk, misp=mu0',{
  expect_equal(testsim$y, modelsim$y0)
  expect_equal(testsim$y, modelsim$y1[[1]])
  expect_equal(testsim$y, modelsim$y1[[2]])
  expect_equal(testsim$y, modelsim$y1[[3]])
})
min.y <- max.y <- rep(-999,1000)
for(i in 1:1000){
  testsim <- simulate_randomwalk(i,50,"GLMM", sd.vec = c(.5,.05), init_u = 1)
  min.y[i] <- min(testsim$y)
  max.y[i] <- max(testsim$y)
}
test_that('randomwalk, GLMM, min y',{
  expect_equal(TRUE, min(min.y) > 0 )
  expect_equal(TRUE, max(max.y) < 3000 )
} )

test_that('randomwalk, GLMM, mkDat, y',{
    expect_equal(modelsim$y0, modeldat$h0$y)
    expect_equal(modelsim$y1[[1]], modeldat$h1[[1]]$y)
    expect_equal(modelsim$y1[[2]], modeldat$h1[[2]]$y)
    expect_equal(modelsim$y1[[3]], modeldat$h1[[3]]$y)
})

test_that('randomwalk, GLMM, mkDat, mod and simre',{
    expect_equal(2, modeldat$h0$mod)
    expect_equal(0, modeldat$h0$sim_re)
    expect_equal(2, modeldat$h1[[1]]$mod)
    expect_equal(0, modeldat$h1[[1]]$sim_re)
    expect_equal(1, modeldat$h1[[2]]$mod)
    expect_equal(0, modeldat$h1[[2]]$sim_re)
    expect_equal(2, modeldat$h1[[3]]$mod)
    expect_equal(0, modeldat$h1[[3]]$sim_re)
})


test.true.par <- list(
  h0 = list(
    mu = 0.1, ln_sig_y = log(0.5),
    ln_sig_u = log(0.05),
    u = modelsim$random$u0
  ),
  h1 = list(
    list(
      mu = 0.1, ln_sig_y = log(0.5),
      ln_sig_u = numeric(0),
      u = rep(0, length(modelsim$random$u0))
    ),
    list(
      mu = 0.1, ln_sig_y = log(0.5),
      ln_sig_u = log(0.05),
      u = modelsim$random$u1[[2]]
    ),
    list(
      mu = 0, ln_sig_y = log(0.5),
      ln_sig_u = log(0.05),
      u = modelsim$random$u1[[3]]
    )
  )
)

test.est.par <- list(
  h0 = list(
    mu = 0, ln_sig_y = 0,
    ln_sig_u = 0,
    u = rep(1, length(modelsim$random$u0))
  ),
  h1 = list(
    list(
      mu = 0, ln_sig_y = 0,
      ln_sig_u = numeric(0),
      u = rep(0, length(modelsim$random$u0))
    ),
    list(
      mu = 0, ln_sig_y = 0,
      ln_sig_u = 0,
      u = rep(1, length(modelsim$random$u0))
    ),
    list(
      mu = 0, ln_sig_y = 0,
      ln_sig_u = 0,
      u = rep(1, length(modelsim$random$u0))
    )
  )
)

test.map <- list(
  h0 = list(),
  h1 = list(
    list(u = rep(factor(NA), length(modelsim$random$u0))),
    list(),
    list(mu = factor(NA))
  )
)

test_that("randomwalk, GLMM, mkpar", {
  expect_equal(test.true.par, model.true.par)
  expect_equal(test.est.par$h0, model.est.par$h0)
})

test_that("randomwalk, GLMM, mkmap", {
  expect_equal(test.map, model.true.map)
  expect_equal(test.map, model.est.map)
})

test.random <- (list(h0 = "u", h1 = list(NULL, "u", "u")))
model.random <- mkTMBrandom("randomwalk", c('missre', 'gamma-lognorm', 'mu0'))
test_that("randomwalk, GLMM, mkrandom", {
  expect_equal(test.random, model.random)
})


test.true.objpar <- list(
  h0 = test.true.par$h0[1:3],
  h1 = list(
    test.true.par$h1[[1]][1:2],
    test.true.par$h1[[2]][1:3],
    test.true.par$h1[[3]][2:3]
  )
) 
test.est.objpar <- list(
  h0 = test.est.par$h0[1:3],
  h1 = list(
    test.est.par$h1[[1]][1:2],
    test.est.par$h1[[2]][1:3],
    test.est.par$h1[[2]][2:3]
  )
)

test_that("randomwalk, GLMM, init.obj", {
for(h in 1:4){
  if(h == 1){
    model.true.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.true.par$h0, 
                     map = model.true.map$h0, random = model.random$h0, DLL = "randomwalk")
    model.est.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.est.par$h0, 
                     map = model.est.map$h0, random = model.random$h0, DLL = "randomwalk")
    expect_equal(unlist(test.true.objpar$h0), model.true.obj$par)
    expect_equal(unlist(test.est.objpar$h0), model.est.obj$par)
  }
  if(h>1){
    model.true.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.true.par$h1[[h-1]], 
                     map = model.true.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "randomwalk")
    model.est.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.est.par$h1[[h-1]], 
                     map = model.est.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "randomwalk")
    expect_equal(unlist(test.true.objpar$h1[[h-1]]), model.true.obj$par)
    expect_equal(unlist(test.est.objpar$h1[[h-1]]), model.est.obj$par)
  }
  rm(model.true.obj, model.est.obj)
}
})

dyn.unload(TMB::dynlib("../../src/randomwalk"))