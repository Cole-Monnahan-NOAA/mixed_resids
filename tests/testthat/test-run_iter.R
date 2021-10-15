#Functional tests on setting up models correctly
## Unit tests for functions from resid_fns.R
source('../../R/sim_data.R')
source('../../R/model_fns.R')
source('../../R/resid_fns.R')
source('../../R/utils.R')
library(TMB)
library(INLA)
library(DHARMa)

## Test model in = model out
run_iter_test <- function(ii, n, ng=0, mod, cov.mod = 'norm', misp, fit.true = FALSE, savefiles=TRUE){
  dyn.load(dynlib(paste0('../../src/',mod)))
  true.parms <- setup_trueparms(mod,misp)
  
  ## simulate data with these parameters
  message(ii, ": Simulating data...")
  sim.dat <- simdat(n, ng, mod, cov.mod, true.parms, misp, ii)
  
  init.dat <- mkTMBdat(sim.dat, true.parms, mod, misp)
  init.par <- mkTMBpar(true.parms, sim.dat, mod, misp, fit.true)
  init.random <- mkTMBrandom(mod, misp, fit.true)
  init.map <- mkTMBmap(init.par, mod, misp, fit.true)
  mod.out <- list(h0 = NULL, h1 = NULL)
  for(h in 1:2){
    
    message(ii, ": Optimizing two competing models...")
    init.obj <- list(data = init.dat[[h]], parameters = init.par[[h]], map = init.map[[h]], random = init.random[[h]], DLL = mod)
    mod.out[[h]] <- fit_tmb(obj.args = init.obj, control = list(run.model = !fit.true, do.sdreport = TRUE))
  }
  
  return(mod.out)
}

context('overdispersion test')
mods <- c('linmod', 'simpleGLMM', 'spatial')
for(m in 1:length(mods)){
  test_that(paste0(mods[m], ' overdispersion'),{
    N <- 100; Ng <- 10 #Ng only used if mod == 'simpleGLMM'
    true.parms <- setup_trueparms(mods[m],'overdispersion')
    dat <- simdat(n=N, ng = 10, mod=mods[m], trueparms = true.parms, misp = 'overdispersion', seed=1)
    mod.true <- run_iter_test(ii=1, n=N, ng = Ng, mod = mods[m], misp = 'overdispersion', fit.true = TRUE)  
    
    h0.true.parm <- c(true.parms$theta, true.parms$sd.vec)
    h1.true.parm <- c(true.parms$theta, true.parms$sd.vec[1:2])
    
    if(!is.null(true.parms$fam)){
      if(true.parms$fam == 'Poisson'){
        h0.true.parm <- c(true.parms$theta, true.parms$sd.vec[-1])
        h1.true.parm <- c(true.parms$theta, true.parms$sd.vec[2])
      }
    } 
    
    h0.fit.parm <- c(mod.true$h0$obj$par[1], exp(mod.true$h0$obj$par[-1]))
    h1.fit.parm <- c(mod.true$h1$obj$par[1], exp(mod.true$h1$obj$par[-1]))
    
    if(mods[m] == 'linmod'){ 
      h0.true.parm <- h1.true.parm <- c(true.parms$theta, true.parms$sd.vec[1])
      h0.fit.parm <- c(mod.true$h0$obj$par[1:2], exp(mod.true$h0$obj$par[3]))
      h1.fit.parm <- c(mod.true$h1$obj$par[1:2], exp(mod.true$h1$obj$par[3]))
    }
    if(mods[m] == 'spatial'){
      Kappa <- sqrt(8)/true.parms$sp.parm
      Tau <- 1/(sqrt(4*pi)*Kappa*true.parms$sd.vec[2])
      # beta, theta, tau, kappa, sig_v
      h0.true.parm <- c(true.parms$theta, true.parms$sd.vec[1], Tau, Kappa, true.parms$sd.vec[3] )
      h1.true.parm <- c(true.parms$theta, true.parms$sd.vec[1], Tau, Kappa )
      if(true.parms$fam=="Poisson"){
        h0.true.parm[2] <- 0
        h1.true.parm[2] <- 0
      }
      # beta, theta, tau, kappa, sig_v
      h0.fit.parm <- c(mod.true$h0$obj$par[1:2], exp(mod.true$h0$obj$par[3:5]))
      h1.fit.parm <- c(mod.true$h1$obj$par[1:2], exp(mod.true$h1$obj$par[3:4]))
    }
    expect_equal(h0.true.parm, unname(h0.fit.parm))
    expect_equal(h1.true.parm, unname(h1.fit.parm))
    expect_equal(dat$y0, mod.true$h0$obj$env$data$y)
    expect_equal(dat$y1, mod.true$h1$obj$env$data$y)
    expect_equal(dat$x, mod.true$h0$obj$env$data$X)
    expect_equal(dat$x, mod.true$h1$obj$env$data$X)
    expect_equal(dat$random$u, mod.true$h0$obj$env$parList()$u)
    expect_equal(dat$random$omega,mod.true$h0$obj$env$parList()$omega)
    expect_equal(dat$random$u, mod.true$h1$obj$env$parList()$u)
    expect_equal(dat$random$omega, mod.true$h1$obj$env$parList()$omega)
    if(mods[m] == 'linmod'){
      expect_equal(NULL, mod.true$h0$obj$env$parList()$v)
      expect_equal(NULL, mod.true$h1$obj$env$parList()$v)
    } else {
      expect_equal(dat$random$v, mod.true$h0$obj$env$parList()$v)
      expect_error(expect_equal(dat$random$v, mod.true$h1$obj$env$parList()$v))
    }
  })
}

context('outliers test')
mods <- c('linmod', 'randomwalk', 'simpleGLMM', 'spatial')
for(m in 1:length(mods)){ 
  test_that(paste0(mods[m], ' outliers'),{
    N <- 100; Ng <- 10 #Ng only used if mod == 'simpleGLMM'
    true.parms <- setup_trueparms(mods[m],'outliers')
    dat <- simdat(n=N, ng = 10, mod=mods[m], trueparms = true.parms, misp = 'outliers', seed=1)
    mod.true <- run_iter_test(ii=1, n=N, ng = Ng, mod = mods[m], misp = 'outliers', fit.true = TRUE)  
    
    h0.true.parm <- c(true.parms$theta, true.parms$sd.vec)
    h1.true.parm <- c(true.parms$theta, true.parms$sd.vec[1:2])
    
    if(!is.null(true.parms$fam)){
      if(true.parms$fam == 'Poisson'){
        h0.true.parm <- c(true.parms$theta, true.parms$sd.vec[-1])
        h1.true.parm <- c(true.parms$theta, true.parms$sd.vec[2])
      }
    } 
    
    h0.fit.parm <- c(mod.true$h0$obj$par[1], exp(mod.true$h0$obj$par[-1]))
    h1.fit.parm <- c(mod.true$h1$obj$par[1], exp(mod.true$h1$obj$par[-1]))
    
    if(mods[m] == 'linmod'){ 
      h0.true.parm <- h1.true.parm <- c(true.parms$theta, true.parms$sd.vec)
      h0.fit.parm <- c(mod.true$h0$obj$par[1:2], exp(mod.true$h0$obj$par[3]))
      h1.fit.parm <- c(mod.true$h1$obj$par[1:2], exp(mod.true$h1$obj$par[3]))
    }
    if(mods[m] == 'spatial'){
      Kappa <- sqrt(8)/true.parms$sp.parm
      Tau <- 1/(sqrt(4*pi)*Kappa*true.parms$sd.vec[2])
      # beta, theta, tau, kappa
      h0.true.parm <- c(true.parms$theta, true.parms$sd.vec[1], Tau, Kappa )
      h1.true.parm <- c(true.parms$theta, true.parms$sd.vec[1], Tau, Kappa )
      if(true.parms$fam=="Poisson"){
        h0.true.parm[2] <- 0
        h1.true.parm[2] <- 0
      }
      # beta, theta, tau, kappa
      h0.fit.parm <- c(mod.true$h0$obj$par[1:2], exp(mod.true$h0$obj$par[3:4]))
      h1.fit.parm <- c(mod.true$h1$obj$par[1:2], exp(mod.true$h1$obj$par[3:4]))
    }
    expect_equal(h0.true.parm, unname(h0.fit.parm))
    expect_equal(h1.true.parm, unname(h1.fit.parm))
    expect_equal(dat$y0, mod.true$h0$obj$env$data$y)
    expect_equal(dat$y1, mod.true$h1$obj$env$data$y)
    if(mods[m] == 'randomwalk'){
      expect_equal(NULL, mod.true$h0$obj$env$data$X)
      expect_equal(NULL, mod.true$h1$obj$env$data$X)
    } else {
      expect_equal(dat$x, mod.true$h0$obj$env$data$X)
      expect_equal(dat$x, mod.true$h1$obj$env$data$X)
    }
    expect_equal(dat$random$u, mod.true$h0$obj$env$parList()$u)
    expect_equal(dat$random$omega,mod.true$h0$obj$env$parList()$omega)
    expect_equal(dat$random$u, mod.true$h1$obj$env$parList()$u)
    expect_equal(dat$random$omega, mod.true$h1$obj$env$parList()$omega)
    expect_equal(dat$random$v, mod.true$h0$obj$env$parList()$v)
    expect_equal(dat$random$v, mod.true$h1$obj$env$parList()$v)
    
  })
}

context('missing covariate test')
mods <- c('linmod', 'simpleGLMM', 'spatial')
for(m in 1:length(mods)){ 
  test_that(paste0(mods[m], ' miss.cov'),{
    N <- 100; Ng <- 10 #Ng only used if mod == 'simpleGLMM'
    true.parms <- setup_trueparms(mods[m],'miss.cov')
    dat <- simdat(n=N, ng = 10, mod=mods[m], trueparms = true.parms, misp = 'miss.cov', seed=1)
    mod.true <- run_iter_test(ii=1, n=N, ng = Ng, mod = mods[m], misp = 'miss.cov', fit.true = TRUE)  
    
    h0.true.parm <- c(true.parms$theta, true.parms$sd.vec)
    h1.true.parm <- c(true.parms$theta[1], true.parms$sd.vec)
    
    if(!is.null(true.parms$fam)){
      if(true.parms$fam == 'Poisson'){
        h0.true.parm <- c(true.parms$theta, true.parms$sd.vec[2])
        h1.true.parm <- c(true.parms$theta[1], true.parms$sd.vec[2])
      }
    } 
    
    if(mods[m] == 'spatial'){
      Kappa <- sqrt(8)/true.parms$sp.parm
      Tau <- 1/(sqrt(4*pi)*Kappa*true.parms$sd.vec[2])
      # beta, theta, tau, kappa
      h0.true.parm <- c(true.parms$theta, true.parms$sd.vec[1], Tau, Kappa )
      h1.true.parm <- c(true.parms$theta[1], true.parms$sd.vec[1], Tau, Kappa )
      if(true.parms$fam=="Poisson"){
        h0.true.parm[3] <- 0
        h1.true.parm[2] <- 0
      }
    }
    h0.fit.parm <- c(mod.true$h0$obj$par[1:2], exp(mod.true$h0$obj$par[-(1:2)]))
    h1.fit.parm <- c(mod.true$h1$obj$par[1], exp(mod.true$h1$obj$par[-1]))
    if(!is.null(true.parms$fam)){
      if(true.parms$fam=="Poisson"){
        h0.true.parm[3] <- 0
        h1.true.parm[2] <- 0
        h0.fit.parm[3] <- 0
        h1.fit.parm[2] <- 0
      }
    }

    expect_equal(h0.true.parm, unname(h0.fit.parm))
    expect_equal(h1.true.parm, unname(h1.fit.parm))
    expect_equal(dat$y0, mod.true$h0$obj$env$data$y)
    expect_equal(dat$y1, mod.true$h1$obj$env$data$y)
    expect_equal(dat$x, mod.true$h0$obj$env$data$X)
    expect_equal(as.matrix(dat$x[,1]), mod.true$h1$obj$env$data$X)
    expect_equal(dat$random$u, mod.true$h0$obj$env$parList()$u)
    expect_equal(dat$random$omega,mod.true$h0$obj$env$parList()$omega)
    expect_equal(dat$random$u, mod.true$h1$obj$env$parList()$u)
    expect_equal(dat$random$omega, mod.true$h1$obj$env$parList()$omega)
    expect_equal(dat$random$v, mod.true$h0$obj$env$parList()$v)
    expect_equal(dat$random$v, mod.true$h1$obj$env$parList()$v)
    
  })
}

context('mu0 test')
test_that('randomwalk mu0',{
  N <- 100
  true.parms <- setup_trueparms('randomwalk','mu0')
  dat <- simdat(n=N, ng = 10, mod='randomwalk', trueparms = true.parms, misp = 'mu0', seed=1)
  mod.true <- run_iter_test(ii=1, n=N, ng = Ng, mod ='randomwalk', misp = 'mu0', fit.true = TRUE)  
  
  expect_equal(c(true.parms$theta, true.parms$sd.vec), unname(c(mod.true$h0$obj$par[1], exp(mod.true$h0$obj$par[-1]))))
  expect_equal(true.parms$sd.vec, unname(exp(mod.true$h1$obj$par)))
  expect_equal(dat$y0, mod.true$h0$obj$env$data$y)
  expect_equal(dat$y1, mod.true$h1$obj$env$data$y)
  expect_equal(dat$random$u, mod.true$h0$obj$env$parList()$u)
  expect_equal(dat$random$u, mod.true$h1$obj$env$parList()$u)
  
})

context('misp.omega test')
test_that('spatial, misp.omega',{
  
  true.parms <- setup_trueparms('spatial','misp.omega')
  dat <- simdat(n=100, mod='spatial', trueparms = true.parms, misp = 'misp.omega', seed=1)
  mod.true <- run_iter_test(ii=1, n=100, ng=0, mod = 'spatial', misp = 'misp.omega', fit.true = TRUE) 
  Kappa <- sqrt(8)/true.parms$sp.parm
  Tau <- 1/(sqrt(4*pi)*Kappa*true.parms$sd.vec[2])
  # beta, theta, tau, kappa, sig_v
  h0.true.parm <- c(true.parms$theta, true.parms$sd.vec[1], Tau, Kappa, true.parms$sd.vec[3] )
  h1.true.parm <- c(true.parms$theta, true.parms$sd.vec[1], Tau, Kappa )
  if(true.parms$fam=="Poisson"){
    h0.true.parm[2] <- 0
    h1.true.parm[2] <- 0
  }
  # beta, theta, tau, kappa, sig_v
  h0.fit.parm <- c(mod.true$h0$obj$par[1:2], exp(mod.true$h0$obj$par[3:5]))
  h1.fit.parm <- c(mod.true$h1$obj$par[1:2], exp(mod.true$h1$obj$par[3:4]))
  
  expect_equal(h0.true.parm, unname(h0.fit.parm))
  expect_equal(h1.true.parm, unname(h1.fit.parm))
  expect_equal(dat$y0, mod.true$h0$obj$env$data$y)
  expect_equal(dat$y1, mod.true$h1$obj$env$data$y)
  expect_equal(dat$x, mod.true$h0$obj$env$data$X)
  expect_equal(dat$x, mod.true$h1$obj$env$data$X)
  expect_equal(dat$random$u, mod.true$h0$obj$env$parList()$u)
  expect_equal(dat$random$omega,mod.true$h0$obj$env$parList()$omega)
  expect_equal(dat$random$u, mod.true$h1$obj$env$parList()$u)
  expect_equal(dat$random$omega, mod.true$h1$obj$env$parList()$omega)
  expect_equal(dat$random$v, mod.true$h0$obj$env$parList()$v)
  expect_equal(dat$random$v, mod.true$h1$obj$env$parList()$v)
  
})

context('functional test on run_iter') #verify all combinations work: model x do.true x method x version
mod.misp <- list(linmod = c('overdispersion', 'outliers', 'miss.cov'),
                 randomwalk = c('mu0', 'outliers'),
                 simpleGLMM = c('overdispersion', 'outliers', 'miss.cov'),
                 spatial = c('overdispersion', 'outliers', 'miss.cov', 'misp.omega'))
osa.methods <- c('fg', 'osg', 'gen', 'cdf')
dharma.methods <- c('uncond', 'cond')
for(m in 1:4){
  if(m==3) N <- 20 else N <- 100
  if(names(mod.misp[m]) == 'spatial') osa.methods <- c('cdf')
  n.misp <- length(mod.misp[[m]])
  for(n in 1:n.misp){
    out.true <- run_iter(ii=1, n=N, ng=5, mod=names(mod.misp)[m], 
                                      misp=mod.misp[[m]][n], do.true = TRUE, savefiles=FALSE)
    out.false <- run_iter(ii=1, n=N, ng=5, mod=names(mod.misp)[m], 
                         misp=mod.misp[[m]][n], do.true = FALSE, savefiles=FALSE)
    test_that(paste0('model: ', names(mod.misp)[m], ', mis-specification: ', mod.misp[[m]][n]),{
      expect_equal(c('pvals', 'resids', 'mles'), names(out.true))
      expect_equal(c('pvals', 'resids', 'mles'), names(out.false))
      expect_true(all(osa.methods %in% out.true$pvals$method))
      expect_true(all(osa.methods %in% out.false$pvals$method))
      expect_true(all(dharma.methods %in% out.true$pvals$method))
      expect_true(all(dharma.methods %in% out.false$pvals$method))
      expect_equal(unique(out.true$pvals$model), names(mod.misp)[m])
      expect_equal(2*(length(osa.methods) + length(dharma.methods)), sum(out.true$pvals$test == 'GOF.ks')) #!cdf fails for linmod-overdispersion and miss.cov-h1
      expect_equal(2*(length(osa.methods) + length(dharma.methods)), sum(out.false$pvals$test == 'GOF.ks')) 
    })
    rm(out.true, out.false)
  }
}


# context('linmod tests')
# 
# test_that('linmod, misp=overdisp',{
#   linmod.opt <- run_iter_test(ii=1, n=100, mod = 'linmod', misp = 'overdispersion', fit.true = FALSE)  
#   expect_equal(rep(0,3), unname(linmod.opt$h0$obj$par))
#   expect_equal(rep(0,3), unname(linmod.opt$h1$obj$par))
#   expect_equal(c(lin.mod.true.parms$theta, log(lin.mod.true.parms$sd.vec[1])),  unname(linmod.opt$h0$opt$par), tolerance = 0.1)
#   expect_error(expect_equal(c(lin.mod.true.parms$theta, log(lin.mod.true.parms$sd.vec[1])),  unname(linmod.opt$h1$opt$par), tolerance = 0.1))
# 
#   })  
# 
# test_that('linmod, misp=miss.cov',{
#   
#   lin.mod.true.parms <- setup_trueparms('linmod','miss.cov')
#   lin.mod.simdat <- simdat(n=100, mod='linmod', trueparms = lin.mod.true.parms, misp = 'miss.cov', seed=1)
#   linmod.opt <- run_iter_test(ii=1, n=100, mod = 'linmod', misp = 'miss.cov', fit.true = FALSE)  
#   expect_equal(rep(0,3), unname(linmod.opt$h0$obj$par))
#   expect_equal(rep(0,2), unname(linmod.opt$h1$obj$par))
#   expect_equal(c(lin.mod.true.parms$theta, log(lin.mod.true.parms$sd.vec[1])),  unname(linmod.opt$h0$opt$par), tolerance = 0.1)
#   expect_error(expect_equal(c(lin.mod.true.parms$theta, log(lin.mod.true.parms$sd.vec[1])),  unname(linmod.opt$h1$opt$par), tolerance = 0.1))
#   
# }) 
# 
# test_that('linmod, misp=outliers',{
#   
#   lin.mod.true.parms <- setup_trueparms('linmod','outliers')
#   lin.mod.simdat <- simdat(n=100, mod='linmod', trueparms = lin.mod.true.parms, misp = 'outliers', seed=1)
#   linmod.opt <- run_iter_test(ii=1, n=100, mod = 'linmod', misp = 'outliers', fit.true = FALSE)  
#   expect_equal(rep(0,3), unname(linmod.opt$h0$obj$par))
#   expect_equal(rep(0,3), unname(linmod.opt$h1$obj$par))
#   expect_equal(c(lin.mod.true.parms$theta, log(lin.mod.true.parms$sd.vec[1])),  unname(linmod.opt$h0$opt$par), tolerance = 0.1)
#   expect_equal(c(lin.mod.true.parms$theta, log(lin.mod.true.parms$sd.vec[1])),  unname(linmod.opt$h1$opt$par), tolerance = 0.1)
#   expect_gte(unname(linmod.opt$h0$opt$par[1]), unname(linmod.opt$h1$opt$par[1]))
#   expect_gte(unname(linmod.opt$h0$opt$par[2]), unname(linmod.opt$h1$opt$par[2]))
#   expect_gte(unname(linmod.opt$h1$opt$par[3]), unname(linmod.opt$h0$opt$par[3]))
#   
# }) 
# 
# context('randomwalk tests')
# 
# test_that('randomwalk, misp=outliers',{
#   
#   rw.true.parms <- setup_trueparms('randomwalk','outliers')
#   rw.simdat <- simdat(n=100, mod='randomwalk', trueparms = rw.true.parms, misp = 'outliers', seed=1)
#   rw.opt <- run_iter_test(ii=1, n=100, mod = 'randomwalk', misp = 'outliers', fit.true = FALSE)  
#   expect_equal(rep(0,3), unname(rw.opt$h0$obj$par))
#   expect_equal(rep(0,3), unname(rw.opt$h1$obj$par))
#   expect_equal(c(rw.true.parms$theta, log(rw.true.parms$sd.vec[1])), unname(rw.opt$h0$opt$par[1:2]), tolerance = 0.15)
#   expect_equal(log(rw.true.parms$sd.vec[2]), unname(rw.opt$h0$opt$par[3]), tolerance = 0.25)
#   expect_error(expect_equal(c(rw.true.parms$theta, log(rw.true.parms$sd.vec[1])), unname(rw.opt$h1$opt$par[1:2]), tolerance = 0.15))
#   expect_equal(log(rw.true.parms$sd.vec[2]), unname(rw.opt$h0$opt$par[3]), tolerance = 0.25)
#   
# })
# 
# test_that('randomwalk, misp=mu0',{
#   
#   rw.true.parms <- setup_trueparms('randomwalk','mu0')
#   rw.simdat <- simdat(n=100, mod='randomwalk', trueparms = rw.true.parms, misp = 'mu0', seed=1)
#   rw.opt <- run_iter_test(ii=1, n=100, mod = 'randomwalk', misp = 'mu0', fit.true = FALSE)  
#   expect_equal(rep(0,3), unname(rw.opt$h0$obj$par))
#   expect_equal(rep(0,2), unname(rw.opt$h1$obj$par))
#   expect_equal(c(rw.true.parms$theta, log(rw.true.parms$sd.vec[1])), unname(rw.opt$h0$opt$par[1:2]), tolerance = 0.15)
#   expect_equal(log(rw.true.parms$sd.vec[2]), unname(rw.opt$h0$opt$par[3]), tolerance = 0.25)
#   expect_error(expect_equal(log(rw.true.parms$sd.vec[1]), unname(rw.opt$h1$opt$par[1]), tolerance = 0.15))
#   expect_equal(log(rw.true.parms$sd.vec[2]), unname(rw.opt$h0$opt$par[2]), tolerance = 0.25)
#   
# })
# 
# context('simpleGLMM tests')
# 
# test_that('simpleGLMM, misp=outliers',{
#   
#   sg.true.parms <- setup_trueparms('simpleGLMM','outliers')
#   sg.simdat <- simdat(n=30, ng=5, mod='simpleGLMM', trueparms = sg.true.parms, misp = 'outliers', seed=1)
#   sg.opt <- run_iter_test(ii=1, n=30, ng=5, mod = 'simpleGLMM', misp = 'outliers', fit.true = FALSE)  
#   if(sg.true.parms$fam == 'Poisson'){
#     expect_equal(rep(0,2), unname(sg.opt$h0$obj$par))
#     expect_equal(rep(0,2), unname(sg.opt$h1$obj$par))
#     #tolerance good for both correct and misspecified
#     expect_equal(c(sg.true.parms$theta, sg.true.parms$sd.vec[-1]),  unname(c(sg.opt$h0$opt$par[1], exp(sg.opt$h0$opt$par[-1]))), tolerance = 0.2)
#     expect_equal(c(sg.true.parms$theta, sg.true.parms$sd.vec[-1]), unname(c(sg.opt$h0$opt$par[1], exp(sg.opt$h0$opt$par[-1]))), tolerance = 0.2)
#   } else {
#     expect_equal(rep(0,3), unname(sg.opt$h0$obj$par))
#     expect_equal(rep(0,3), unname(sg.opt$h1$obj$par))
#     #tolerance good for both correct and misspecified
#     expect_equal(c(sg.true.parms$theta, sg.true.parms$sd.vec),  unname(c(sg.opt$h0$opt$par[1], exp(sg.opt$h0$opt$par[-1]))), tolerance = 0.15)
#     expect_equal(c(sg.true.parms$theta, sg.true.parms$sd.vec),  unname(c(sg.opt$h1$opt$par[1], exp(sg.opt$h1$opt$par[-1]))), tolerance = 0.15)
#   }
#   
# }) 
# 
# test_that('simpleGLMM, misp=overdisp',{
#   
#   sg.true.parms <- setup_trueparms('simpleGLMM','overdispersion')
#   sg.simdat <- simdat(n=30, ng=5, mod='simpleGLMM', trueparms = sg.true.parms, misp = 'overdispersion', seed=1)
#   sg.opt <- run_iter_test(ii=1, n=100, ng=10, mod = 'simpleGLMM', misp = 'overdispersion', fit.true = FALSE)  
#   if(sg.true.parms$fam == 'Poisson'){
#     expect_equal(rep(0,3), unname(sg.opt$h0$obj$par))
#     expect_equal(rep(0,2), unname(sg.opt$h1$obj$par))
#     expect_equal(sg.true.parms$theta, unname(sg.opt$h0$opt$par[1]), tol=.05)
#     expect_error(expect_equal(sg.true.parms$theta, unname(sg.opt$h1$opt$par[1]), tol=.05))
#     expect_equal(log(sg.true.parms$sd.vec[-1]),  unname(sg.opt$h0$opt$par[2:3]), tolerance = 0.2)
#     expect_equal(log(sg.true.parms$sd.vec[2]),  unname(sg.opt$h1$opt$par[2]), tolerance = 0.2)
#   } else {
#     expect_equal(rep(0,4), unname(sg.opt$h0$obj$par))
#     expect_equal(rep(0,3), unname(sg.opt$h1$obj$par))
#     # Normal-Normal case: identifiability issues between sig.y and sig.v?
#     # expect_equal(sg.true.parms$theta, unname(sg.opt$h0$opt$par[1]), tol=.05)
#     # expect_error(expect_equal(sg.true.parms$theta, unname(sg.opt$h1$opt$par[1]), tol=.05))
#     # expect_equal(sg.true.parms$sd.vec^2,  unname(exp(sg.opt$h0$opt$par[-1])^2), tolerance = 0.2)
#     # expect_error(expect_equal(sg.true.parms$sd.vec[1:2]^2,  unname(exp(sg.opt$h1$opt$par[-1])^2), tolerance = 0.2))
#     
#   }
# }) 
# 
# test_that('simpleGLMM, misp=miss.cov',{
#   
#   sg.true.parms <- setup_trueparms('simpleGLMM','miss.cov')
#   sg.simdat <- simdat(n=30, ng=5, mod='simpleGLMM', trueparms = sg.true.parms, misp = 'miss.cov', seed=1)
#   sg.opt <- run_iter_test(ii=1, n=30, ng=5, mod = 'simpleGLMM', misp = 'miss.cov', fit.true = FALSE)  
#   expect_equal(rep(0,4), unname(sg.opt$h0$obj$par))
#   expect_equal(rep(0,3), unname(sg.opt$h1$obj$par))
#   expect_equal(c(sg.true.parms$theta, log(sg.true.parms$sd.vec)),  unname(sg.opt$h0$opt$par), tolerance = 0.1)
#   expect_error(expect_equal(c(sg.true.parms$theta[1], log(sg.true.parms$sd.vec)),  unname(sg.opt$h1$opt$par), tolerance = 0.1))
#   
# })
# 
# context('spatial tests')
# 
# test_that('spatial, misp=outliers',{
#   
#   sp.true.parms <- setup_trueparms('spatial','outliers')
#   sp.simdat <- simdat(n=50, mod='spatial', trueparms = sp.true.parms, misp = 'outliers', seed=1)
#   sp.opt <- run_iter_test(ii=1, n=30, ng=5, mod = 'spatial', misp = 'outliers', fit.true = FALSE)  
#   expect_equal(rep(0,3), unname(sp.opt$h0$obj$par))
#   expect_equal(rep(0,3), unname(sp.opt$h1$obj$par))
#   expect_equal(c(sp.true.parms$theta, log(sp.true.parms$sd.vec)),  unname(sp.opt$h0$opt$par), tolerance = 0.2)
#   expect_equal(c(sp.true.parms$theta, log(sp.true.parms$sd.vec)),  unname(sp.opt$h1$opt$par), tolerance = 0.2)
#   
# }) 
# 
# test_that('spatial, misp=overdisp',{
#   
#   sp.true.parms <- setup_trueparms('spatial','overdispersion')
#   sp.simdat <- simdat(n=50, mod='spatial', trueparms = sp.true.parms, misp = 'overdispersion', seed=1)
#   sp.opt <- run_iter_test(ii=1, n=50, mod = 'spatial', misp = 'overdispersion', fit.true = FALSE)  
#   expect_equal(rep(0,4), unname(sp.opt$h0$obj$par))
#   expect_equal(rep(0,3), unname(sp.opt$h1$obj$par))
#   expect_equal(c(sp.true.parms$theta, log(sp.true.parms$sd.vec[1])),  unname(sp.opt$h0$opt$par), tolerance = 0.1)
#   expect_error(expect_equal(c(sp.true.parms$theta, log(sp.true.parms$sd.vec[1])),  unname(sp.opt$h1$opt$par), tolerance = 0.1))
#   
# }) 

