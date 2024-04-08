setwd("../..")
source('R/startup.R')

simulate_randomwalk <- function(seed, n, type, sd.vec, misp = NULL){
  if(type == "LMM") mu <- 2
  if(type == "GLMM") mu <- .02
  sig_y <- sd.vec[1]
  sig_u <- sd.vec[2]
  ## simulate random track
  set.seed(seed)
  u0 <- c(0, cumsum( rnorm(n-1, mu, sig_u) ))
  if(type == "LMM"){
    set.seed(seed*10)
    y0 <- rnorm(n, u0, sd=sig_y)
  }
  if(type == "GLMM"){
    set.seed(seed*10)
    y0 <- rgamma(n, 1/sig_y^2, scale = exp(u0)*sig_y^2)
  }
  
  u1 <- y1 <- list()
  for(m in 1:length(misp)){
    u1[[m]] <- u0
    y1[[m]] <- y0
    if(misp[m] == 'mispre'){
      if(type == "LMM"){
        set.seed(seed)
        u1[[m]] <- cumsum( exp(rnorm(n, 2, sig_u) ))
      }
      if(type == "GLMM"){
        set.seed(seed)
        u1[[m]] <- cumsum( rgamma(n, 0.5, 20))
      }
      
      if(type == "LMM"){
        set.seed(seed)
        y1[[m]] <- rnorm(n, u1[[m]], sd=sig_y)
      } 
      if(type == "GLMM"){
        set.seed(seed)
        y1[[m]] <- rgamma(n, 1/sig_y^2, scale = exp(u1[[m]])*sig_y^2)
      }
    }
    
    if(misp[m] == 'hsk'){
      set.seed(seed)
      y1[[m]] <- rnorm(n, u0, sd=sqrt((1:n)^2))
    }
    
    if(misp[m] == 'mu0'){
      set.seed(seed)
      u1[[m]] <- u0 <- c(0, cumsum( rnorm(n-1, 0, sig_u) ))
    }
  }
  
  out <- list(y0 = y0, y1 = y1, 
              u0 = u0, u1 = u1)
}

context("LMM randomwalk tests")
testsim <- simulate_randomwalk(123, 100, "LMM", sd.vec = c(1,1),
                               misp = c('missre', 'hsk', 'mu0', 'mispre'))
rw.par <- setup_trueparms(mod ='randomwalk', 
                          misp = c('missre', 'hsk', 'mu0', 'mispre'), 
                          fam = "Gaussian", link = "identity",
                          type = "LMM")
modelsim <- simdat(n=100, mod='randomwalk', type = "LMM",
                 trueparms = rw.par, 
                 misp = c('missre', 'hsk', 'mu0', 'mispre'), 
                 seed=123)

modeldat <- mkTMBdat(Dat = modelsim, 
                     Pars = rw.par,  
                     Mod = "randomwalk", 
                     Misp = c('missre', 'hsk', 'mu0', 'mispre'), 
                     Type = "LMM",
                     doTrue = 1)
model.true.par <- mkTMBpar(Pars = rw.par, 
                           Dat = modelsim, 
                           Mod = "randomwalk", 
                           Misp = c('missre', 'hsk', 'mu0', 'mispre'), 
                           doTrue = TRUE)
model.est.par <- mkTMBpar(Pars = rw.par, 
                          Dat = modelsim, 
                          Mod = "randomwalk", 
                          Misp = c('missre', 'hsk', 'mu0', 'mispre'), 
                          doTrue = FALSE)
model.true.map <- mkTMBmap(Pars = model.true.par, 
                      Mod = "randomwalk", 
                      Misp = c('missre', 'hsk', 'mu0', 'mispre'), 
                      Type = "LMM",
                      doTrue = TRUE)
model.est.map <- mkTMBmap(Pars = model.est.par, 
                      Mod = "randomwalk", 
                      Misp = c('missre', 'hsk', 'mu0', 'mispre'), 
                      Type = "LMM",
                      doTrue = FALSE)
test_that('randomwalk, LMM, simdat',{
  expect_equal(testsim$y0, modelsim$y0)
  expect_equal(testsim$y0, modelsim$y1[[1]])
  expect_equal(testsim$y1[[2]], modelsim$y1[[2]])
  expect_equal(testsim$y0, modelsim$y1[[3]])
  expect_equal(testsim$y1[[4]], modelsim$y1[[4]])
})
max.y <- matrix(-999, 1000, 4)
for(i in 1:1000){
  testsimmax <- simulate_randomwalk(i,100,"LMM", sd.vec = c(1,1),
                                    misp = c('missre', 'hsk', 'mu0', 'mispre'))
  max.y[i,1] <- max(testsimmax$y1[[1]])
  max.y[i,2] <- max(testsimmax$y1[[2]])
  max.y[i,3] <- max(testsimmax$y1[[3]])
  max.y[i,4] <- max(testsimmax$y1[[4]])
}
test_that('randomwalk, LMM, max y',{
  expect_equal(TRUE, max(max.y[,1]) < 500 )
  expect_equal(TRUE, max(max.y[,2]) < 600 )
  expect_equal(TRUE, max(max.y[,3]) < 500 )
  expect_equal(TRUE, max(max.y[,4]) < 2000 )
} )

test_that('randomwalk, LMM, mkDat, y',{
    expect_equal(modelsim$y0, modeldat$h0$y)
    expect_equal(modelsim$y1[[1]], modeldat$h1[[1]]$y)
    expect_equal(modelsim$y1[[2]], modeldat$h1[[2]]$y)
    expect_equal(modelsim$y1[[3]], modeldat$h1[[3]]$y)
    expect_equal(modelsim$y1[[4]], modeldat$h1[[4]]$y)
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
    expect_equal(0, modeldat$h1[[4]]$mod)
    expect_equal(0, modeldat$h1[[4]]$sim_re)
})

test.true.par <- list(
  h0 = list( #correct
    mu = 2, ln_sig_y = 0,
    ln_sig_u = 0,
    u = modelsim$random$u0
  ),
  h1 = list( #missre
    list(
      mu = 2, ln_sig_y = 0,
      ln_sig_u = numeric(0),
      u = rep(0, length(modelsim$random$u0))
    ),
    list( #hsk
      mu = 2, ln_sig_y = 0,
      ln_sig_u = 0,
      u = modelsim$random$u0
    ),
    list( #mu0
      mu = 0, ln_sig_y = 0,
      ln_sig_u = 0,
      u = modelsim$random$u1[[3]]
    ),
    list( #mispre
      mu = 2, ln_sig_y = 0,
      ln_sig_u = 0,
      u = modelsim$random$u0
    )
  )
)
test.est.par <- test.true.par
test.est.par$h0$mu <- 
  test.est.par$h1[[1]]$mu <- 
  test.est.par$h1[[2]]$mu <-  
  test.est.par$h1[[3]]$mu <-  
  test.est.par$h1[[4]]$mu <- 0

test.est.par$h0$u <- 
  test.est.par$h1[[2]]$u <- 
  test.est.par$h1[[3]]$u <- 
  test.est.par$h1[[4]]$u <- 
      rep(1, length(modelsim$random$u0))

test.map.true <- list(
  h0 = list(),
  h1 = list(
    list(u = rep(factor(NA), length(modelsim$random$u0))),
    list(),
    list(mu = factor(NA)),
    list()
  )
)
test.map.est <- list(
  h0 = list(),
  h1 = list(
    list(u = rep(factor(NA), length(modelsim$random$u0))),
    list(),
    list(mu = factor(NA)),
    list()
  )
)

test_that("randomwalk, LMM, mkpar", {
  expect_equal(test.true.par, model.true.par)
  expect_equal(test.est.par, model.est.par)
})

test_that("randomwalk, LMM, mkmap", {
  expect_equal(test.map.true, model.true.map)
  expect_equal(test.map.est, model.est.map)
})

test.random <- (list(h0 = "u", h1 = list(NULL, "u", "u", "u")))
model.random <- mkTMBrandom("randomwalk", c('missre', 'hsk', 'mu0', 'mispre'))
test_that("randomwalk, LMM, mkrandom", {
  expect_equal(test.random, model.random)
})


test.true.objpar <- list(
  h0 = test.true.par$h0[1:3],
  h1 = list(
    test.true.par$h1[[1]][1:2],
    test.true.par$h1[[2]][1:3],
    test.true.par$h1[[3]][2:3],
    test.true.par$h1[[4]][1:3]
  )
) 
test.est.objpar <- list(
  h0 = test.est.par$h0[1:3],
  h1 = list(
    test.est.par$h1[[1]][1:2],
    test.est.par$h1[[2]][1:3],
    test.est.par$h1[[3]][2:3],
    test.est.par$h1[[4]][1:3]
  )
)
dyn.load(dynlib("src/randomwalk"))
test_that("randomwalk, GLMM, init.obj", {

for(h in 1:4){
  if(h == 1){
    model.true.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.true.par$h0, 
                     map = model.true.map$h0, random = model.random$h0, DLL = "randomwalk")
    model.est.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.est.par$h0, 
                     map = model.est.map$h0, random = model.random$h0, DLL = "randomwalk")
    report.true <- model.true.obj$report()
    report.est <- model.est.obj$report()
    expect_equal(testsim$u0, report.true$u)
    expect_equal(unlist(test.true.objpar$h0), model.true.obj$par)
    expect_equal(unlist(test.est.objpar$h0), model.est.obj$par)
    expect_equal(rep(1,100), report.est$u)
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
testsim <- simulate_randomwalk(123, 100, "GLMM", sd.vec = c(.5,.2),
                               misp = c('missre', 'gamma-normal', 'mu0', 'mispre'))
rw.par <- setup_trueparms(mod ='randomwalk', 
                          misp = c('missre', 'gamma-normal', 'mu0', 'mispre'), 
                          fam = "Gamma", link = "log",
                          type = "GLMM")
modelsim <- simdat(n=100, mod='randomwalk', 
                 trueparms = setup_trueparms(mod ='randomwalk', 
                                             misp = c('missre', 'gamma-normal', 
                                                      'mu0', 'mispre'), 
                                             fam = "Gamma", link = "log",
                                             type = "GLMM"), 
                 misp = c('missre', 'gamma-normal', 'mu0', 'mispre'), 
                 type = "GLMM", seed=123)
modeldat <- mkTMBdat(Dat = modelsim, 
                       Pars = setup_trueparms(mod ='randomwalk', 
                                              misp = c('missre', 'gamma-normal', 
                                                       'mu0', 'mispre'), 
                                              fam = "Gamma", link = "log",
                                              type = "GLMM"),  
                       Mod = "randomwalk", 
                       Misp = c('missre', 'gamma-normal', 'mu0', 'mispre'), 
                       Type = "GLMM",
                     doTrue = 1)
model.true.par <- mkTMBpar(Pars = rw.par, 
                           Dat = modelsim, 
                           Mod = "randomwalk", 
                           Misp = c('missre', 'gamma-normal', 'mu0', 'mispre'), 
                           doTrue = TRUE)
model.est.par <- mkTMBpar(Pars = rw.par, 
                          Dat = modelsim, 
                          Mod = "randomwalk", 
                          Misp = c('missre', 'gamma-normal', 'mu0', 'mispre'), 
                          doTrue = FALSE)
model.true.map <- mkTMBmap(Pars = model.true.par, 
                      Mod = "randomwalk", 
                      Misp = c('missre', 'gamma-normal', 'mu0', 'mispre'), 
                      Type = "GLMM", doTrue = TRUE)
model.est.map <- mkTMBmap(Pars = model.est.par, 
                      Mod = "randomwalk", 
                      Misp = c('missre', 'gamma-normal', 'mu0', 'mispre'), 
                      Type = "GLMM", doTrue = FALSE)
test_that('randomwalk, misp=mu0',{
  expect_equal(testsim$y0, modelsim$y0)
  expect_equal(testsim$y0, modelsim$y1[[1]])
  expect_equal(testsim$y0, modelsim$y1[[2]])
  expect_equal(testsim$y0, modelsim$y1[[3]])
  expect_equal(testsim$y1[[4]], modelsim$y1[[4]])
})
max.y <- matrix(-999,1000,4)
for(i in 1:1000){
  testsim <- simulate_randomwalk(i,100,"GLMM", sd.vec = c(.5,.15),
                                 misp = c('missre', 'gamma-normal', 'mu0', 'mispre'))

  max.y[i,1] <- max(testsim$y1[[1]])
  max.y[i,2] <- max(testsim$y1[[2]])
  max.y[i,3] <- max(testsim$y1[[3]])
  max.y[i,4] <- max(testsim$y1[[4]])
}
test_that('randomwalk, GLMM, min y',{
  expect_equal(TRUE, max(max.y[,1]) < 3000 )
  expect_equal(TRUE, max(max.y[,2]) < 3000 )
  expect_equal(TRUE, max(max.y[,3]) < 3000 )
  expect_equal(TRUE, max(max.y[,4]) < 3000 )
} )

test_that('randomwalk, GLMM, mkDat, y',{
    expect_equal(modelsim$y0, modeldat$h0$y)
    expect_equal(modelsim$y1[[1]], modeldat$h1[[1]]$y)
    expect_equal(modelsim$y1[[2]], modeldat$h1[[2]]$y)
    expect_equal(modelsim$y1[[3]], modeldat$h1[[3]]$y)
    expect_equal(modelsim$y1[[4]], modeldat$h1[[4]]$y)
})

test_that('randomwalk, GLMM, mkDat, mod and simre',{
    expect_equal(2, modeldat$h0$mod)
    expect_equal(0, modeldat$h0$sim_re)
    expect_equal(2, modeldat$h1[[1]]$mod)
    expect_equal(0, modeldat$h1[[1]]$sim_re)
    expect_equal(0, modeldat$h1[[2]]$mod)
    expect_equal(0, modeldat$h1[[2]]$sim_re)
    expect_equal(2, modeldat$h1[[3]]$mod)
    expect_equal(0, modeldat$h1[[3]]$sim_re)
    expect_equal(2, modeldat$h1[[4]]$mod)
    expect_equal(0, modeldat$h1[[4]]$sim_re)
})


test.true.par <- list(
  h0 = list(
    mu = 0.01, 
    ln_sig_y = log(0.5),
    ln_sig_u = log(0.05),
    u = modelsim$random$u0
  ),
  h1 = list(
    list(
      mu = 0.01, 
      ln_sig_y = log(0.5),
      ln_sig_u = numeric(0),
      u = rep(0, length(modelsim$random$u0))
    ),
    list(
      mu = 0.01, 
      ln_sig_y = log(0.5),
      ln_sig_u = log(0.05),
      u = modelsim$random$u0
    ),
    list(
      mu = 0, 
      ln_sig_y = log(0.5),
      ln_sig_u = log(0.05),
      u = modelsim$random$u1[[3]]
    ),
    list(
      mu = 0.01, 
      ln_sig_y = log(0.5),
      ln_sig_u = log(0.05),
      u = modelsim$random$u0
    )
  )
)

test.est.par <- list(
  h0 = list(
    mu = 0, 
    ln_sig_y = 0,
    ln_sig_u = 0,
    u = rep(1, length(modelsim$random$u0))
  ),
  h1 = list(
    list(
      mu = 0, 
      ln_sig_y = 0,
      ln_sig_u = numeric(0),
      u = rep(0, length(modelsim$random$u0))
    ),
    list(
      mu = 0, 
      ln_sig_y = 0,
      ln_sig_u = 0,
      u = rep(1, length(modelsim$random$u0))
    ),
    list(
      mu = 0, 
      ln_sig_y = 0,
      ln_sig_u = 0,
      u = rep(1, length(modelsim$random$u0))
    ),
    list(
      mu = 0, 
      ln_sig_y = 0,
      ln_sig_u = 0,
      u = rep(1, length(modelsim$random$u0))
    )
  )
)

test.true.map <- list(
  h0 = list(),
  h1 = list(
    list(u = rep(factor(NA), length(modelsim$random$u0))),
    list(),
    list(mu = factor(NA)),
    list()
  )
)

test.est.map <- list(
  h0 = list(),
  h1 = list(
    list(u = rep(factor(NA), length(modelsim$random$u0))),
    list(),
    list(mu = factor(NA)),
    list()
  )
)

test_that("randomwalk, GLMM, mkpar", {
  expect_equal(test.true.par, model.true.par)
  expect_equal(test.est.par, model.est.par)
})

test_that("randomwalk, GLMM, mkmap", {
  expect_equal(test.true.map, model.true.map)
  expect_equal(test.est.map, model.est.map)
})

test.random <- (list(h0 = "u", h1 = list(NULL, "u", "u", "u")))
model.random <- mkTMBrandom("randomwalk", 
                            c('missre', 'gamma-normal', 'mu0', 'mispre'))
test_that("randomwalk, GLMM, mkrandom", {
  expect_equal(test.random, model.random)
})


test.true.objpar <- list(
  h0 = test.true.par$h0[1:3],
  h1 = list(
    test.true.par$h1[[1]][1:2],
    test.true.par$h1[[2]][1:3],
    test.true.par$h1[[3]][2:3],
    test.true.par$h1[[4]][1:3]
  )
) 
test.est.objpar <- list(
  h0 = test.est.par$h0[1:3],
  h1 = list(
    test.est.par$h1[[1]][1:2],
    test.est.par$h1[[2]][1:3],
    test.est.par$h1[[3]][2:3],
    test.est.par$h1[[4]][1:3]
  )
)

test_that("randomwalk, GLMM, init.obj", {
  dharma.methods <- c('uncond', 'cond', 'uncond_nrot', 'cond_nrot' )
  osa.methods <- c('gen', 'cdf', 'mcmc','pears')
  osa.true.out <- osa.est.out <- osa.true.ks <- osa.est.ks <- list()
  test.mod.true <- run_iter(ii = 123, n = 100, ng = 0, mod = "randomwalk", cov.mod = NULL, 
                       misp =  c('missre', 'gamma-normal', 'mu0', 'mispre'), 
                       family = "Gamma", link = "log",
                       type = 'GLMM', do.true = TRUE, savefiles = FALSE)
  test.mod.est <- run_iter(ii = 123, n = 100, ng = 0, mod = "randomwalk", cov.mod = NULL, 
                            misp =  c('missre', 'gamma-normal', 'mu0', 'mispre'), 
                            family = "Gamma", link = "log",
                            type = 'GLMM', do.true = FALSE, savefiles = FALSE)
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
    expect_equal(dplyr::filter(test.mod.true$pvals, 
                               test == "GOF.ks" & method == "gen" & 
                                 version == paste0("h", h-1))$pvalue,
                 osa.true.ks[[h]]$gen.pvalue)
    expect_equal(dplyr::filter(test.mod.true$pvals, 
                               test == "GOF.ks" & method == "cdf" & 
                                 version == paste0("h", h-1))$pvalue,
                 osa.true.ks[[h]]$cdf.pvalue)
    rm(model.true.obj, model.est.obj)
  }
})

context("test correct LMM model")
dyn.load(dynlib('src/randomwalk'))
n <- 500
pvals.fg <- pvals.uncond <- rep(0,n)
for(i in 1:n){
  set.seed(i)
  u0 <- c(0, cumsum( rnorm(99, 2, 1) ))
  set.seed(i*2)
  y0 <- rnorm(100, u0, sd=1)
  
  Dat <- list(y = y0,
              mod = 0,
              sim_re = 1)
  Par <- list(mu = 2,
              ln_sig_y = 0,
              ln_sig_u = 0,
              u = u0)
  obj <- MakeADFun(Dat, Par, random = "u", DLL = "randomwalk")
  res <- oneStepPredict(obj, observation.name = "y", method = "fullGaussian")
  pvals.fg[i] <- ks.test(res$residual, 'pnorm')$p.value
  sim.y <- replicate(1000,{obj$simulate()$y})
  res <- DHARMa::createDHARMa(sim.y, y0, rotation = "estimated")
  pvals.uncond[i] <- ks.test(res$scaledResiduals, 'punif')$p.value
}
test_that("type I error", {
expect_lt(sum(pvals.fg<0.05)/n, 0.05)
expect_lt(sum(pvals.uncond<0.05)/n, 0.07)
})
