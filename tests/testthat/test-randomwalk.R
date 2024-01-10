setwd("../..")
source('R/startup.R')

simulate_randomwalk <- function(seed, n, type, sd.vec, init_u){
  if(type == "LMM") mu <- 2
  if(type == "GLMM") mu <- .05
  sig_y <- sd.vec[1]
  sig_u <- sd.vec[2]
  huge <- 1e3
  ## simulate random track
  set.seed(seed)
  z <- cumsum( rnorm(n, 0, sig_u) ) 
  eta <- init_u + 1:n * mu + z
 
  
  set.seed(seed)
  if(type == "LMM"){
    Y <- rnorm(n, eta, sd=sig_y)
  } 
  if(type == "GLMM"){
    Y <- rgamma(n, 1/sig_y^2, scale = exp(eta)*sig_y^2)
  }

  out <- list(y = Y, u = z, eta = eta)
}

context("LMM randomwalk tests")
testsim <- simulate_randomwalk(123, 100, "LMM", sd.vec = c(1,1), init_u = 6)
rw.par <- setup_trueparms(mod ='randomwalk', 
                          misp=c('missre', 'normal-lognorm', 'mu0'), 
                          fam = "Gaussian", link = "identity",
                          type = "LMM")
modelsim <- simdat(n=100, mod='randomwalk', 
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
  testsimmin <- simulate_randomwalk(i,100,"LMM", sd.vec = c(1,1), init_u = 6)
  min.y[i] <- min(testsimmin$y)
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
    beta = 6,
    mu = 2, ln_sig_y = 0,
    ln_sig_u = 0,
    u = modelsim$random$u0
  ),
  h1 = list(
    list(
      beta = 6, 
      mu = 2, ln_sig_y = 0,
      ln_sig_u = numeric(0),
      u = rep(0, length(modelsim$random$u0))
    ),
    list(
      beta = 6, 
      mu = 2, ln_sig_y = 0,
      ln_sig_u = 0,
      u = modelsim$random$u1[[2]]
    ),
    list(
      beta = 6, 
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

test.est.par$h0$beta <- 
  test.est.par$h1[[1]]$beta <- 
  test.est.par$h1[[2]]$beta <-  
  test.est.par$h1[[3]]$beta <-0

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
  h0 = test.true.par$h0[1:4],
  h1 = list(
    test.true.par$h1[[1]][1:3],
    test.true.par$h1[[2]][1:4],
    test.true.par$h1[[3]][c(1,3:4)]
  )
) 
test.est.objpar <- list(
  h0 = test.est.par$h0[1:4],
  h1 = list(
    test.est.par$h1[[1]][1:3],
    test.est.par$h1[[2]][1:4],
    test.est.par$h1[[2]][c(1,3:4)]
  )
)
dyn.load(dynlib("src/randomwalk"))
test_that("randomwalk, LMM, init.obj", {

for(h in 1:4){
  if(h == 1){
    model.true.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.true.par$h0, 
                     map = model.true.map$h0, random = model.random$h0, DLL = "randomwalk")
    model.est.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.est.par$h0, 
                     map = model.est.map$h0, random = model.random$h0, DLL = "randomwalk")
    report.true <- model.true.obj$report()
    report.est <- model.est.obj$report()
    expect_equal(testsim$u, report.true$u)
    expect_equal(testsim$eta, report.true$eta)
    expect_equal(unlist(test.true.objpar$h0), model.true.obj$par)
    expect_equal(unlist(test.est.objpar$h0), model.est.obj$par)
    expect_equal(rep(1,100), report.est$u)
    expect_equal(rep(1,100), report.est$eta)
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
testsim <- simulate_randomwalk(123, 100, "GLMM", sd.vec = c(.5,.05), init_u = .1)
rw.par <- setup_trueparms(mod ='randomwalk', 
                          misp = c('missre', 'gamma-lognorm', 'mu0'), 
                          fam = "Gamma", link = "log",
                          type = "GLMM")
modelsim <- simdat(n=100, mod='randomwalk', 
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
  testsim <- simulate_randomwalk(i,100,"GLMM", sd.vec = c(.5,.05), init_u = .1)
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
    beta = 0.1,
    mu = 0.05, 
    ln_sig_y = log(0.5),
    ln_sig_u = log(0.05),
    u = modelsim$random$u0
  ),
  h1 = list(
    list(
      beta = 0.1,
      mu = 0.05, 
      ln_sig_y = log(0.5),
      ln_sig_u = numeric(0),
      u = rep(0, length(modelsim$random$u0))
    ),
    list(
      beta = 0.1,
      mu = 0.05, 
      ln_sig_y = log(0.5),
      ln_sig_u = log(0.05),
      u = modelsim$random$u1[[2]]
    ),
    list(
      beta = 0.1,
      mu = 0, 
      ln_sig_y = log(0.5),
      ln_sig_u = log(0.05),
      u = modelsim$random$u1[[3]]
    )
  )
)

test.est.par <- list(
  h0 = list(
    beta = 0,
    mu = 0, 
    ln_sig_y = 0,
    ln_sig_u = 0,
    u = rep(1, length(modelsim$random$u0))
  ),
  h1 = list(
    list(
      beta = 0,
      mu = 0, 
      ln_sig_y = 0,
      ln_sig_u = numeric(0),
      u = rep(0, length(modelsim$random$u0))
    ),
    list(
      beta = 0,
      mu = 0, 
      ln_sig_y = 0,
      ln_sig_u = 0,
      u = rep(1, length(modelsim$random$u0))
    ),
    list(
      beta = 0,
      mu = 0, 
      ln_sig_y = 0,
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
  expect_equal(test.est.par, model.est.par)
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
  h0 = test.true.par$h0[1:4],
  h1 = list(
    test.true.par$h1[[1]][1:3],
    test.true.par$h1[[2]][1:4],
    test.true.par$h1[[3]][c(1,3:4)]
  )
) 
test.est.objpar <- list(
  h0 = test.est.par$h0[1:4],
  h1 = list(
    test.est.par$h1[[1]][1:3],
    test.est.par$h1[[2]][1:4],
    test.est.par$h1[[2]][c(1,3:4)]
  )
)

test_that("randomwalk, GLMM, init.obj", {
  dharma.methods <- c('uncond', 'cond', 'uncond_nrot', 'cond_nrot' )
  osa.methods <- c('gen', 'cdf', 'mcmc','pears')
  osa.true.out <- osa.est.out <- osa.true.ks <- osa.est.ks <- list()
  test.mod <- run_iter(ii = 123, n = 100, ng = 0, mod = "randomwalk", cov.mod = NULL, 
                       misp =  c('missre', 'gamma-lognorm', 'mu0'), 
                       family = "Gamma", link = "log",
                       type = 'GLMM', do.true = TRUE, savefiles = FALSE)
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
    model.est.opt <- nlminb(model.est.obj$par, model.est.obj$fn, model.est.obj$gr)
    osa.true.out[[h]] <- list(
      gen.resid = oneStepPredict(model.true.obj, observation.name="y",
                    data.term.indicator='keep',
                    range = c(0,Inf),
                    method="oneStepGeneric", trace=FALSE,
                    discrete = FALSE)$residual,
      cdf.resid = oneStepPredict(model.true.obj, observation.name="y",
                    data.term.indicator='keep' ,
                    method="cdf", trace=FALSE,
                    discrete = FALSE)$residual
    )
    osa.true.ks[[h]] <- list(
      gen.pvalue = ks.test(osa.true.out[[h]]$gen.resid,'pnorm')$p.value,
      cdf.pvalue = ks.test(osa.true.out[[h]]$cdf.resid,'pnorm')$p.value
    )
    expect_equal(dplyr::filter(test.mod$pvals, 
                               test == "GOF.ks" & method == "gen" & 
                                 version == paste0("h", h-1))$pvalue,
                 osa.true.ks[[h]]$gen.pvalue)
    expect_equal(dplyr::filter(test.mod$pvals, 
                               test == "GOF.ks" & method == "cdf" & 
                                 version == paste0("h", h-1))$pvalue,
                 osa.true.ks[[h]]$cdf.pvalue)
    rm(model.true.obj, model.est.obj)
  }
})
