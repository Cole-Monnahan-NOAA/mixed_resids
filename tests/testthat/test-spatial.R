source('../../R/sim_data.R')
source('../../R/model_fns.R')

context('spatial tests')
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
test_that("successful mesh", {
  expect_equal(1000, sum(class.vec, na.rm = TRUE))
})

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

context('LMM spatial tests')
testsim <- simulate_spatial(123, 100, type = "LMM")
spatial.par <- setup_trueparms(mod ='spatial', 
                          misp = c('missre', 'normal-gamma', 'mispre'), 
                          fam = "Gaussian", link = "identity",
                          type = "LMM")
modelsim <- simdat(n=100, mod='spatial', type = "LMM",
                 trueparms = setup_trueparms(mod ='spatial', 
                                             misp = c('missre', 'normal-gamma', 'mispre'), 
                                             fam = "Gaussian", link = "identity",
                                             type = "LMM"), 
                      misp = c("missre", "normal-gamma", "mispre"), seed=123)

modeldat <- mkTMBdat(Dat = modelsim, 
                       Pars = setup_trueparms(mod ='spatial', 
                                             misp = c("missre", "normal-gamma", "mispre"), 
                                             fam = "Gaussian", link = "identity",
                                             type = "LMM"),  
                       Mod = "spatial", 
                       Misp = c("missre", "normal-gamma", "mispre"), 
                       Type = "LMM")

model.true.par <- mkTMBpar(Pars = spatial.par, 
                           Dat = modelsim, 
                           Mod = "spatial", 
                           Misp = c("missre", "normal-gamma", "mispre"), 
                           Type = "LMM",
                           doTrue = TRUE)
model.est.par <- mkTMBpar(Pars = spatial.par, 
                          Dat = modelsim, 
                          Mod = "spatial", 
                          Misp = c("missre", "normal-gamma", "mispre"),
                          Type = "LMM", 
                          doTrue = FALSE)
model.true.map <- mkTMBmap(Pars = model.true.par, 
                      Mod = "spatial", 
                      Misp = c('missre', 'normal-gamma', 'mispre'), 
                      Type = "LMM")
model.est.map <- mkTMBmap(Pars = model.est.par, 
                      Mod = "spatial", 
                      Misp = c('missre', 'normal-gamma', 'mispre'), 
                      Type = "LMM")
test_that('spatial, LMM, simdata',{
  expect_equal(testsim$y0, modelsim$y0)
  expect_equal(testsim$y0, modelsim$y1[[1]])
  expect_equal(testsim$y0, modelsim$y1[[2]])
  expect_equal(testsim$y1, modelsim$y1[[3]])
})

# One misp fits a gamma model to normal data,
# make sure no values are below zero in simulations,
# also no extreme values
min.y0 <- max.y0 <- 
  min.y1 <- max.y1 <-rep(NA, 1000)
for(i in 1:1000){
  testsim <- simulate_spatial(i, 100, type = "LMM")
  min.y0[i] <- min(testsim$y0)
  max.y0[i] <- max(testsim$y0)
  min.y1[i] <- min(testsim$y1)
  max.y1[i] <- max(testsim$y1)
}
test_that('spatial, LMM, bounds',{
  expect_gt(min(min.y0), 0)
  expect_gt(min(min.y1), 0)
  expect_lt(max(max.y0), 100)
  expect_lt(max(max.y1), 100)
})
test_that('spatial, LMM, mkDat, y',{
    expect_equal(modelsim$y0, modeldat$h0$y)
    expect_equal(modelsim$y1[[1]], modeldat$h1[[1]]$y)
    expect_equal(modelsim$y1[[2]], modeldat$h1[[2]]$y)
    expect_equal(modelsim$y1[[3]], modeldat$h1[[3]]$y)
})

test_that('spatial, LMM, mkDat, fam, link and simre',{
    expect_equal(0, modeldat$h0$family)
    expect_equal(2, modeldat$h0$link)
    expect_equal(0, modeldat$h0$sim_re)
    expect_equal(0, modeldat$h1[[1]]$family)
    expect_equal(2, modeldat$h1[[1]]$link)
    expect_equal(0, modeldat$h1[[1]]$sim_re)
    expect_equal(00, modeldat$h1[[2]]$family)
    expect_equal(2, modeldat$h1[[2]]$link)
    expect_equal(0, modeldat$h1[[2]]$sim_re)
    expect_equal(0, modeldat$h1[[3]]$family)
    expect_equal(2, modeldat$h1[[3]]$link)
    expect_equal(0, modeldat$h1[[3]]$sim_re)
})

omega.true <- omega.h1 <- rep(0, modelsim$mesh$n)
omega.true[modelsim$mesh$idx$loc] <- as.vector(modelsim$random$omega0)
omega.h1[modelsim$mesh$idx$loc] <- as.vector(modelsim$random$omega1[[2]])
test.true.par <- list(
  h0 = list(
    beta = 20, 
    theta = log(1),
    ln_tau = log(1/(sqrt(4*pi)*sqrt(8)/50)),
    ln_kappa = log(sqrt(8)/50),
    omega = omega.true
  ),
  h1 = list(
    list(
      beta = 20, 
      theta = log(1),
      ln_tau = 0,
      ln_kappa = 0,
      omega = rep(0, length(omega.true))
    ),
    list(
      beta = 20, 
      theta = log(1),
      ln_tau = log(1/(sqrt(4*pi)*sqrt(8)/50)),
      ln_kappa = log(sqrt(8)/50),
      omega = omega.h1
    ),
    list(
      beta = 20, 
      theta = log(1),
      ln_tau = log(1/(sqrt(4*pi)*sqrt(8)/50)),
      ln_kappa = log(sqrt(8)/50),
      omega = omega.true
    )
  )
)
test.est.par <- list(
  h0 = list(
    beta = 0, 
    theta = 0,
    ln_tau = 0,
    ln_kappa = 0,
    omega = rep(0,length(omega.true))
  ),
  h1 = list(
    list(
      beta = 0, 
      theta = 0,
      ln_tau = 0,
      ln_kappa = 0,
      omega = rep(0,length(omega.true))
    ),
    list(
      beta = 0, 
      theta = 0,
      ln_tau = 0,
      ln_kappa = 0,
      omega = rep(0,length(omega.true))
    ),
    list(
      beta = 0, 
      theta = 0,
      ln_tau = 0,
      ln_kappa = 0,
      omega = rep(0,length(omega.true))
    )
  )
)

test.map <- list(
  h0 = list(),
  h1 = list(
    list(ln_tau = factor(NA),
         ln_kappa = factor(NA),
         omega = rep(factor(NA), length(omega.true))),
    list(),
    list()
  )
)

test_that("spatial, LMM, mkpar", {
  expect_equal(test.true.par, model.true.par)
  expect_equal(test.est.par, model.est.par)
})

test_that("spatial, LMM, mkmap", {
  expect_equal(test.map, model.true.map)
  expect_equal(test.map, model.est.map)
})

test.random <- (list(h0 = "omega", h1 = list(NULL, "omega", "omega")))
model.random <- mkTMBrandom("spatial", c('missre', 'norm-gamma', 'mispre'))
test_that("spatial, LMM, mkrandom", {
  expect_equal(test.random, model.random)
})

test.true.objpar <- list(
  h0 = test.true.par$h0[1:4],
  h1 = list(
    test.true.par$h1[[1]][1:2],
    test.true.par$h1[[2]][1:4],
    test.true.par$h1[[3]][1:4]
  )
) 
test.est.objpar <- list(
  h0 = test.est.par$h0[1:4],
  h1 = list(
    test.est.par$h1[[1]][1:2],
    test.est.par$h1[[2]][1:4],
    test.est.par$h1[[3]][1:4]
  )
)
dyn.load(TMB::dynlib("../../src/spatial"))
test_that("spatial, LMM, init.obj", {
for(h in 1:4){
  if(h == 1){
    model.true.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.true.par$h0, 
                     map = model.true.map$h0, random = model.random$h0, DLL = "spatial")
    model.est.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.est.par$h0, 
                     map = model.est.map$h0, random = model.random$h0, DLL = "spatial")
    test.true.objpar.h0 <- unlist(test.true.objpar$h0)
    test.est.objpar.h0 <- unlist(test.est.objpar$h0)
    expect_equal(test.true.objpar.h0, model.true.obj$par)
    expect_equal(test.est.objpar.h0, model.est.obj$par)
  }
  if(h>1){
    model.true.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.true.par$h1[[h-1]], 
                     map = model.true.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "spatial")
    model.est.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.est.par$h1[[h-1]], 
                     map = model.est.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "spatial")
    test.true.objpar.h1 <- unlist(test.true.objpar$h1[[h-1]])
    test.est.objpar.h1 <- unlist(test.est.objpar$h1[[h-1]])
    expect_equal(test.true.objpar.h1, model.true.obj$par)
    expect_equal(test.est.objpar.h1, model.est.obj$par)
  }
  rm(model.true.obj, model.est.obj)
}
})


context('GLMM spatial tests')
testsim <- simulate_spatial(123, 100, 
                           misp = c("missre", "pois-zip", "mispre"), 
                           fam = "Poisson", type = "GLMM")
spatial.par <- setup_trueparms(mod ='spatial', 
                               misp = c('missre', 'pois-zip', 'mispre'), 
                               fam = "Poisson", link = "log",
                               type = "GLMM")
modelsim <- simdat(n=100, mod='spatial', type = "GLMM", 
                 trueparms = setup_trueparms(mod ='spatial', 
                                             misp = c('missre', 'pois-zip', 'mispre'), 
                                             fam = "Poisson", link = "log",
                                             type = "GLMM"), 
                      misp = c("missre", "pois-zip", "mispre"), seed=123)

modeldat <- mkTMBdat(Dat = modelsim, 
                       Pars = setup_trueparms(mod ='spatial', 
                                             misp = c('missre', 'pois-zip', 'mispre'), 
                                             fam = "Poisson", link = "log",
                                             type = "GLMM"),  
                       Mod = "spatial", 
                       Misp = c('missre', 'pois-zip', 'mispre'), 
                       Type = "GLMM")

model.true.par <- mkTMBpar(Pars = spatial.par, 
                           Dat = modelsim, 
                           Mod = "spatial", 
                           Misp =  c("missre", "pois-zip", "mispre"), 
                           Type = "GLMM",
                           doTrue = TRUE)
model.est.par <- mkTMBpar(Pars = spatial.par, 
                          Dat = modelsim, 
                          Mod = "spatial", 
                          Misp = c("missre", "pois-zip", "mispre"),
                          Type = "GLMM", 
                          doTrue = FALSE)
model.true.map <- mkTMBmap(Pars = model.true.par, 
                      Mod = "spatial", 
                      Misp = c("missre", "pois-zip", "mispre"), 
                      Type = "GLMM")
model.est.map <- mkTMBmap(Pars = model.est.par, 
                      Mod = "spatial", 
                      Misp = c("missre", "pois-zip", "mispre"), 
                      Type = "GLMM")
test_that('spatial, GLMM',{
  expect_equal(testsim$omega, modelsim$random$omega0)
  expect_equal(testsim$y0, modelsim$y0)
  expect_equal(testsim$y1[[1]], modelsim$y1[[1]])
  expect_equal(testsim$y1[[2]], modelsim$y1[[2]])
  expect_equal(testsim$y1[[3]], modelsim$y1[[3]])
})

# test upper bound and %zero inflation
y0.max <- rep(NA, 1000)
y1.max <- y.zip <- list()
y1.max[[1]] <- y1.max[[2]] <- y1.max[[3]] <- rep(NA, 1000)
y.zip[[1]] <- y.zip[[2]] <- rep(NA, 0)
for(i in 1:1000){
  testsim <- simulate_spatial(i, 100, 
                           misp = c("missre", "pois-zip", "mispre"), 
                           fam = "Poisson", type = "GLMM")
  y0.max[i] <- max(testsim$y0)
  y1.max[[1]][i] <- max(testsim$y1[[1]])
  y1.max[[2]][i] <- max(testsim$y1[[2]])
  y1.max[[3]][i] <- max(testsim$y1[[3]])
  y.zip[[1]][i] <- sum(testsim$y0 == 0)
  y.zip[[2]][i] <- sum(testsim$y1[[2]] == 0)
}

test_that('spatial, GLMM, bounds', {
  expect_lt(max(y0.max), 10000)
  expect_lt(max(y1.max[[3]]), 10000)
  expect_lt(max(y.zip[[1]]), 55)
  expect_gt(min(y.zip[[2]] - y.zip[[1]]), 0)
})

test_that('spatial, GLMM, mkDat, y',{
    expect_equal(modelsim$y0, modeldat$h0$y)
    expect_equal(modelsim$y1[[1]], modeldat$h1[[1]]$y)
    expect_equal(modelsim$y1[[2]], modeldat$h1[[2]]$y)
    expect_equal(modelsim$y1[[3]], modeldat$h1[[3]]$y)
})

test_that('spatial, GLMM, mkDat, fam, link and simre',{
    expect_equal(200, modeldat$h0$family)
    expect_equal(0, modeldat$h0$link)
    expect_equal(0, modeldat$h0$sim_re)
    expect_equal(200, modeldat$h1[[1]]$family)
    expect_equal(0, modeldat$h1[[1]]$link)
    expect_equal(0, modeldat$h1[[1]]$sim_re)
    expect_equal(200, modeldat$h1[[2]]$family)
    expect_equal(0, modeldat$h1[[2]]$link)
    expect_equal(0, modeldat$h1[[2]]$sim_re)
    expect_equal(200, modeldat$h1[[3]]$family)
    expect_equal(0, modeldat$h1[[3]]$link)
    expect_equal(0, modeldat$h1[[3]]$sim_re)
})

omega.true <- omega.h1 <- rep(0, modelsim$mesh$n)
omega.true[modelsim$mesh$idx$loc] <- as.vector(modelsim$random$omega0)
omega.h1[modelsim$mesh$idx$loc] <- as.vector(modelsim$random$omega1[[2]])
test.true.par <- list(
  h0 = list(
    beta = 0.5, 
    theta = 0,
    ln_tau = log(1/(sqrt(4*pi)*sqrt(8)/50*sqrt(.25))),
    ln_kappa = log(sqrt(8)/50),
    omega = omega.true
  ),
  h1 = list(
    list(
      beta = 0.5, 
      theta = 0,
      ln_tau = 0,
      ln_kappa = 0,
      omega = rep(0, length(omega.true))
    ),
    list(
      beta = 0.5, 
      theta = 0,
      ln_tau = log(1/(sqrt(4*pi)*sqrt(8)/50*sqrt(.25))),
      ln_kappa = log(sqrt(8)/50),
      omega = omega.true
    ),
    list(
      beta = 0.5, 
      theta = 0,
      ln_tau = log(1/(sqrt(4*pi)*sqrt(8)/50*sqrt(.25))),
      ln_kappa = log(sqrt(8)/50),
      omega = omega.true
    )
  )
)
test.est.par <- list(
  h0 = list(
    beta = 0, 
    theta = 0,
    ln_tau = 0,
    ln_kappa = 0,
    omega = rep(0,length(omega.true))
  ),
  h1 = list(
    list(
      beta = 0, 
      theta = 0,
      ln_tau = 0,
      ln_kappa = 0,
      omega = rep(0,length(omega.true))
    ),
    list(
      beta = 0, 
      theta = 0,
      ln_tau = 0,
      ln_kappa = 0,
      omega = rep(0,length(omega.true))
    ),
    list(
      beta = 0, 
      theta = 0,
      ln_tau = 0,
      ln_kappa = 0,
      omega = rep(0,length(omega.true))
    )
  )
)

test.map <- list(
  h0 = list(theta = factor(NA)),
  h1 = list(
    list(theta = factor(NA),
         ln_tau = factor(NA),
         ln_kappa = factor(NA),
         omega = rep(factor(NA), length(omega.true))),
    list(theta = factor(NA)),
    list(theta = factor(NA))
  )
)

test_that("spatial, GLMM, mkpar", {
  expect_equal(test.true.par, model.true.par)
  expect_equal(test.est.par, model.est.par)
})

test_that("spatial, GLMM, mkmap", {
  expect_equal(test.map, model.true.map)
  expect_equal(test.map, model.est.map)
})

test.random <- (list(h0 = "omega", h1 = list(NULL, "omega", "omega")))
model.random <- mkTMBrandom("spatial", c('missre', 'norm-gamma', 'mispre'))
test_that("spatial, LMM, mkrandom", {
  expect_equal(test.random, model.random)
})

test.true.objpar <- list(
  h0 = test.true.par$h0[c(1,3,4)],
  h1 = list(
    test.true.par$h1[[1]][1],
    test.true.par$h1[[2]][c(1,3,4)],
    test.true.par$h1[[3]][c(1,3,4)]
  )
) 
test.est.objpar <- list(
  h0 = test.est.par$h0[c(1,3,4)],
  h1 = list(
    test.est.par$h1[[1]][1],
    test.est.par$h1[[2]][c(1,3,4)],
    test.est.par$h1[[3]][c(1,3,4)]
  )
)
test_that("spatial, GLMM, init.obj", {
for(h in 1:4){
  if(h == 1){
    model.true.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.true.par$h0, 
                     map = model.true.map$h0, random = model.random$h0, DLL = "spatial")
    model.est.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.est.par$h0, 
                     map = model.est.map$h0, random = model.random$h0, DLL = "spatial")
    test.true.objpar.h0 <- unlist(test.true.objpar$h0)
    test.est.objpar.h0 <- unlist(test.est.objpar$h0)
    expect_equal(test.true.objpar.h0, model.true.obj$par)
    expect_equal(test.est.objpar.h0, model.est.obj$par)
  }
  if(h>1){
    model.true.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.true.par$h1[[h-1]], 
                     map = model.true.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "spatial")
    model.est.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.est.par$h1[[h-1]], 
                     map = model.est.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "spatial")
    test.true.objpar.h1 <- unlist(test.true.objpar$h1[[h-1]])
    test.est.objpar.h1 <- unlist(test.est.objpar$h1[[h-1]])
    expect_equal(test.true.objpar.h1, model.true.obj$par)
    expect_equal(test.est.objpar.h1, model.est.obj$par)
  }
  rm(model.true.obj, model.est.obj)
}
})
