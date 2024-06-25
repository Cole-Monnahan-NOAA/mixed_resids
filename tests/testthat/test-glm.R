setwd("../..")
source('R/startup.R')

simulate_glm <- function(seed, n=100){
  
  ## ngroups <- 5 #number of subjects
  ## nobs <- 100 #number of observations
  beta = log(0.5)

  #Simulate RE
  set.seed(seed)
  y0 <- rpois(n, exp(beta))
 set.seed(seed)
  y1 <- rbinom(n,1,.7)*y0
  
  out <- list(y0 = y0, y1 = y1)
  return(out)
}

context("glm tests")
testsim <- simulate_glm(seed = 123, n = 100)
glm.par <- setup_trueparms(mod ='pois_glm', 
                          misp = c('pois-zip'), 
                          fam = "Poisson", link = "log")
modelsim <- simdat(n=100, mod='pois_glm', 
                 trueparms = setup_trueparms(mod ='pois_glm', 
                                             misp = 'pois-zip', 
                                             fam = "Poisson", link = "log"), 
                 misp = "pois-zip", 
                 seed = 123)

modeldat <- mkTMBdat(Dat = modelsim, 
                       Pars = setup_trueparms(mod ='pois_glm', 
                                             misp = 'pois-zip', 
                                             fam = "Poisson", link = "log"),  
                       Mod = "pois_glm", 
                       Misp = 'pois-zip')

model.true.par <- mkTMBpar(Pars = glm.par, 
                           Dat = modelsim, 
                           Mod = "pois_glm", 
                           Misp = 'pois-zip',
                           doTrue = TRUE)
model.est.par <- mkTMBpar(Pars = glm.par, 
                          Dat = modelsim, 
                          Mod = "pois_glm", 
                          Misp = 'pois-zip',
                          doTrue = FALSE)
model.true.map <- mkTMBmap(Pars = model.true.par, 
                           Mod = "pois_glm", 
                           Misp = 'pois-zip')
model.est.map <- mkTMBmap(Pars = model.est.par,
                          Mod = "pois_glm", 
                          Misp = 'pois-zip')
test_that('pois_glm, simdata',{
  expect_equal(testsim$y0, modelsim$y0)
  expect_equal(testsim$y1, modelsim$y1)
})

test_that('pois_glm, mkDat, y',{
    expect_equal(modelsim$y0, modeldat$h0$y)
    expect_equal(modelsim$y1, modeldat$h1$y)
})

test.true.par <- list(
  h0 = list(
    beta = log(.5) 
  ),
  h1 = list(
    beta = log(.5)
  ))

test.est.par <- list(
  h0 = list(
    beta = 0 
  ),
  h1 = list(
    beta = 0
  ))

test.map <- list(
  h0 = list(),
  h1 = list()
)

test_that("pois_glm, mkpar", {
  expect_equal(test.true.par, model.true.par)
  expect_equal(test.est.par, model.est.par)
})

test_that("pois_glm, mkmap", {
  expect_equal(test.map, model.true.map)
  expect_equal(test.map, model.est.map)
})

# 
# test.true.objpar <- list(
#   h0 = test.true.par$h0[c(1,2,4)],
#   h1 = list(
#     test.true.par$h1[[1]][1:2],
#     test.true.par$h1[[2]][c(1,2,4)],
#     test.true.par$h1[[3]][c(1,2,4)]
#   )
# ) 
# test.est.objpar <- list(
#   h0 = test.est.par$h0[c(1,2,4)],
#   h1 = list(
#     test.est.par$h1[[1]][1:2],
#     test.est.par$h1[[2]][c(1,2,4)],
#     test.est.par$h1[[3]][c(1,2,4)]
#   )
# )
# dyn.load(TMB::dynlib("src/simpleGLMM"))
# test_that("simpleGLMM, LMM, init.obj", {
# for(h in 1:4){
#   if(h == 1){
#     model.true.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.true.par$h0, 
#                      map = model.true.map$h0, random = model.random$h0, DLL = "simpleGLMM")
#     model.est.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.est.par$h0, 
#                      map = model.est.map$h0, random = model.random$h0, DLL = "simpleGLMM")
#     test.true.objpar.h0 <- unlist(test.true.objpar$h0)
#     test.est.objpar.h0 <- unlist(test.est.objpar$h0)
#     names(test.true.objpar.h0)[1:2] <- names(test.est.objpar.h0)[1:2] <- 'beta'
#     expect_equal(test.true.objpar.h0, model.true.obj$par)
#     expect_equal(test.est.objpar.h0, model.est.obj$par)
#   }
#   if(h>1){
#     model.true.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.true.par$h1[[h-1]], 
#                      map = model.true.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "simpleGLMM")
#     model.est.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.est.par$h1[[h-1]], 
#                      map = model.est.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "simpleGLMM")
#     test.true.objpar.h1 <- unlist(test.true.objpar$h1[[h-1]])
#     test.est.objpar.h1 <- unlist(test.est.objpar$h1[[h-1]])
#     if(h != 3) names(test.true.objpar.h1)[1:2] <- names(test.est.objpar.h1)[1:2] <- 'beta'
#     expect_equal(test.true.objpar.h1, model.true.obj$par)
#     expect_equal(test.est.objpar.h1, model.est.obj$par)
#   }
#   rm(model.true.obj, model.est.obj)
# }
# })
# 
# 
# context("GLMM simpleGLMM tests")
# testsim <- simulate_simpleGLMM(seed = 123, n = 100, ng = 5,
#                               type = "GLMM")
# simpleGLMM.par <- setup_trueparms(mod ='simpleGLMM', 
#                           misp = c('missre', 'nb-pois', 'mispre'), 
#                           fam = "NB", link = "log",
#                           type = "GLMM")
# modelsim <- simdat(n=100, ng=5, mod='simpleGLMM',
#                  trueparms = setup_trueparms(mod ='simpleGLMM', 
#                                              misp = c('missre', 'nb-pois', 'mispre'), 
#                                              fam = "NB", link = "log",
#                                              type = "GLMM"), 
#                  misp = c('missre', 'nb-pois', 'mispre'), 
#                  seed = 123, type = "GLMM")
# 
# modeldat <- mkTMBdat(Dat = modelsim, 
#                        Pars = setup_trueparms(mod ='simpleGLMM', 
#                                              misp=c('missre', 'nb-pois', 'mispre'), 
#                                              fam = "NB", link = "log",
#                                              type = "GLMM"),  
#                        Mod = "simpleGLMM", 
#                        Misp = c('missre', 'nb-pois', 'mispre'), 
#                        Type = "GLMM")
# 
# model.true.par <- mkTMBpar(Pars = simpleGLMM.par, 
#                            Dat = modelsim, 
#                            Mod = "simpleGLMM", 
#                            Misp = c('missre', 'nb-pois', 'mispre'), 
#                            Type = "GLMM",
#                            doTrue = TRUE)
# model.est.par <- mkTMBpar(Pars = simpleGLMM.par, 
#                           Dat = modelsim, 
#                           Mod = "simpleGLMM", 
#                           Misp = c('missre', 'nb-pois', 'mispre'), 
#                           Type = "GLMM",
#                           doTrue = FALSE)
# model.true.map <- mkTMBmap(Pars = model.true.par, 
#                       Mod = "simpleGLMM", 
#                       Misp = c('missre', 'nb-pois', 'mispre'), 
#                       Type = "GLMM")
# model.est.map <- mkTMBmap(Pars = model.est.par, 
#                       Mod = "simpleGLMM", 
#                       Misp = c('missre', 'nb-pois', 'mispre'), 
#                       Type = "GLMM")
# test_that('simpleGLMM, GLMM, simdata',{
#   expect_equal(testsim$y0, modelsim$y0)
#   expect_equal(testsim$y0, modelsim$y1[[1]])
#   expect_equal(testsim$y0, modelsim$y1[[2]])
#   expect_equal(testsim$y1, modelsim$y1[[3]])
# })
# 
# min.y0 <- max.y0 <- rep(-999,1000)
# min.y1 <- max.y1 <- rep(-999,1000)
# for(i in 1:1000){
#   testsim <- simulate_simpleGLMM(seed = i, n = 100, ng = 5,
#                               cov.mod = "unif",
#                               type = "LMM")
#   min.y0[i] <- min(testsim$y0)
#   max.y0[i] <- max(testsim$y0)
#   min.y1[i] <- min(testsim$y1)
#   max.y1[i] <- max(testsim$y1)
# }
# test_that('simpleGLMM, GLMM, min/max y',{
#   expect_equal(TRUE, max(max.y0 - min.y0) < 100 )
#   expect_equal(TRUE, max(max.y1 - min.y1) < 3000 )
# } )
# 
# max.y0 <- rep(-999,1000)
# max.y1 <- rep(-999,1000)
# for(i in 1:1000){
#   testsim <- simulate_simpleGLMM(seed = i, n = 100, ng = 5,
#                               type = "GLMM")
#   max.y0[i] <- max(testsim$y0)
#   max.y1[i] <- max(testsim$y1)
# }
# test_that('simpleGLMM, GLMM, min/max y',{
#   expect_equal(TRUE, max(max.y0) < 500 )
#   expect_equal(TRUE, max(max.y1) < 30000 )
# } )
# 
# test_that('simpleGLMM, GLMM, mkDat, y',{
#     expect_equal(modelsim$y0, modeldat$h0$y)
#     expect_equal(modelsim$y1[[1]], modeldat$h1[[1]]$y)
#     expect_equal(modelsim$y1[[2]], modeldat$h1[[2]]$y)
#     expect_equal(modelsim$y1[[3]], modeldat$h1[[3]]$y)
# })
# 
# test_that('simpleGLMM, GLMM, mkDat, fam, link and simre',{
#     expect_equal(600, modeldat$h0$family)
#     expect_equal(0, modeldat$h0$link)
#     expect_equal(0, modeldat$h0$sim_re)
#     expect_equal(600, modeldat$h1[[1]]$family)
#     expect_equal(0, modeldat$h1[[1]]$link)
#     expect_equal(0, modeldat$h1[[1]]$sim_re)
#     expect_equal(200, modeldat$h1[[2]]$family)
#     expect_equal(0, modeldat$h1[[2]]$link)
#     expect_equal(0, modeldat$h1[[2]]$sim_re)
#     expect_equal(600, modeldat$h1[[3]]$family)
#     expect_equal(0, modeldat$h1[[3]]$link)
#     expect_equal(0, modeldat$h1[[3]]$sim_re)
# })
# 
# test.true.par <- list(
#   h0 = list(
#     beta = log(2), 
#     ln_sig_y = 0,
#     theta = log(1),
#     ln_sig_u = log(1),
#     u = modelsim$random$u0
#   ),
#   h1 = list(
#     list(
#     beta = log(2), 
#     ln_sig_y = 0,
#     theta = log(1),
#     ln_sig_u = numeric(0),
#     u = rep(0, length(modelsim$random$u0))
#     ),
#     list(
#       beta = log(2), 
#       ln_sig_y = 0,
#       theta = 0,
#       ln_sig_u = log(1),
#       u = modelsim$random$u1[[2]]
#     ),
#     list(
#       beta = log(2), 
#       ln_sig_y = 0,
#       theta = log(1),
#       ln_sig_u = log(1),
#       u = modelsim$random$u1[[3]]
#     )
#   )
# )
# 
# test.est.par <- list(
#   h0 = list(
#     beta = 0, 
#     ln_sig_y = 0,
#     theta = 0,
#     ln_sig_u = 0,
#     u = rep(0, length(modelsim$random$u0))
#   ),
#   h1 = list(
#     list(
#       beta = 0, 
#       ln_sig_y = 0,
#       theta = 0,
#       ln_sig_u = numeric(0),
#       u = rep(0, length(modelsim$random$u0))
#     ),
#     list(
#       beta = 0, 
#       ln_sig_y = 0,
#       theta = 0,
#       ln_sig_u = 0,
#       u = rep(0, length(modelsim$random$u0))
#     ),
#     list(
#       beta = 0, 
#       ln_sig_y = 0,
#       theta = 0,
#       ln_sig_u = 0,
#       u = rep(0, length(modelsim$random$u0))
#     )
#   )
# )
# 
# test.map <- list(
#   h0 = list(ln_sig_y = factor(NA)),
#   h1 = list(
#     list(ln_sig_y = factor(NA),
#          u = rep(factor(NA), length(modelsim$random$u0))),
#     list(ln_sig_y = factor(NA), theta = factor(NA)),
#     list(ln_sig_y = factor(NA))
#   )
# )
# 
# test_that("simpleGLMM, GLMM, mkpar", {
#   expect_equal(test.true.par, model.true.par)
#   expect_equal(test.est.par, model.est.par)
# })
# 
# test_that("simpleGLMM, GLMM, mkmap", {
#   expect_equal(test.map, model.true.map)
#   expect_equal(test.map, model.est.map)
# })
# 
# test.random <- (list(h0 = "u", h1 = list(NULL, "u", "u")))
# model.random <- mkTMBrandom("simpleGLMM", c('missre', 'nb-pois', 'mispre'))
# test_that("simpleGLMM, GLMM, mkrandom", {
#   expect_equal(test.random, model.random)
# })
# 
# test.true.objpar <- list(
#   h0 = test.true.par$h0[c(1,3,4)],
#   h1 = list(
#     test.true.par$h1[[1]][c(1,3)],
#     test.true.par$h1[[2]][c(1,4)],
#     test.true.par$h1[[3]][c(1,3,4)]
#   )
# ) 
# test.est.objpar <- list(
#   h0 = test.est.par$h0[c(1,3,4)],
#   h1 = list(
#     test.est.par$h1[[1]][c(1,3)],
#     test.est.par$h1[[2]][c(1,4)],
#     test.est.par$h1[[3]][c(1,3,4)]
#   )
# )
# test_that("simpleGLMM, GLMM, init.obj", {
# for(h in 1:4){
#   if(h == 1){
#     model.true.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.true.par$h0, 
#                      map = model.true.map$h0, random = model.random$h0, DLL = "simpleGLMM")
#     model.est.obj <- TMB::MakeADFun(data = modeldat$h0, parameters = model.est.par$h0, 
#                      map = model.est.map$h0, random = model.random$h0, DLL = "simpleGLMM")
#     test.true.objpar.h0 <- unlist(test.true.objpar$h0)
#     test.est.objpar.h0 <- unlist(test.est.objpar$h0)
#     expect_equal(test.true.objpar.h0, model.true.obj$par)
#     expect_equal(test.est.objpar.h0, model.est.obj$par)
#   }
#   if(h>1){
#     model.true.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.true.par$h1[[h-1]], 
#                      map = model.true.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "simpleGLMM")
#     model.est.obj <- TMB::MakeADFun(data = modeldat$h1[[h-1]], parameters = model.est.par$h1[[h-1]], 
#                      map = model.est.map$h1[[h-1]], random = model.random$h1[[h-1]], DLL = "simpleGLMM")
#     test.true.objpar.h1 <- unlist(test.true.objpar$h1[[h-1]])
#     test.est.objpar.h1 <- unlist(test.est.objpar$h1[[h-1]])
#     expect_equal(test.true.objpar.h1, model.true.obj$par)
#     expect_equal(test.est.objpar.h1, model.est.obj$par)
#   }
#   rm(model.true.obj, model.est.obj)
# }
# })
# dyn.unload(TMB::dynlib("src/simpleGLMM"))
