## Placeholder script for analysis of residuals of spatial
## models. Started 11/2020 by Cole.

## Assumes working directory is root of the repo spatial_resids

### Startup
source("R/startup.R")
packageVersion('TMB')                   # 1.9.1
packageVersion('VAST')                  # 3.6.0
packageVersion('FishStatsUtils')        # 2.9.0
packageVersion('DHARMa')                # 0.3.3.0
### Some specific scripts to run. Each script runs a separate
### model, computes different residual types, saves these
### results into the results folder for post-hoc analysis, and
### makes plots comparing residual types

### Naming conventions:
## osa.cdf = one-step ahead (OSA) residual with CDF method
## osa.gen = OSA with generic method
## osa.osg = OSA with one-step Gaussian (failing)
## osa.fg  = OSA with full Guassian (failing)
## sim_cond = simulation residuals (DHARMa) conditioned on RE
## sim_uncon= simulation residual (DHARMa) unconditioned on RE
## sim_parcond = simulation residuals (DHARMa) conditioned on
##   random draws from the joint precision matrix (fixed + RE,
##   conditioned on data)

(cpus <- parallel::detectCores()-1)
reps <- 1:1000

do.true <- TRUE


### LM
### =======================================================================
### Simple linear model as sanity check. Some resid methods not
### applicable b/c no random effects
## possible mispecifications: overdispersion, outliers, misscov
## Set osa.methods and dharma.methods to NULL to turn off
osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'pears')
dharma.methods <- c('uncond_nrot', 'cond_nrot')
run_model(reps, mod='linmod', misp='overdispersion', cov.mod = 'norm',
          type = 'LMM', do.true = do.true)

### GLM
### =======================================================================
### Simple generalized linear model as a sanity check. Some resid methods not 
### applicable b/c no random effects
osa.methods <- c('gen', 'cdf', 'pears')
dharma.methods <- c('uncond_nrot', 'cond_nrot')
run_model(reps, mod='glm', misp='pois-zip', type = 'GLMM', do.true = do.true,
          family = "Poisson", link = "log")


### Random walk from the paper =================================================
## possible mispecifications: 'missre', 'normal-lognorm', 'gamma-lognorm', 'mu0'
#reps <- 1:441
a <- Sys.time()
for(dt in 1:2){
  if(dt == 1){do.true <- TRUE}
  if(dt == 2){do.true <- FALSE}
  dharma.methods <- c('uncond', 'cond', 'uncond_nrot', 'cond_nrot', 'process' )
  osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc','pears', 'process')
  run_model(reps, n = 100, mod='randomwalk', 
            misp =  c('missre', 'hsk', 'mu0', 'mispre'), 
            type = 'LMM', do.true = do.true)
  osa.methods <- c('gen', 'cdf', 'mcmc','pears')
  run_model(reps, n = 100, mod='randomwalk', 
            misp =  c('missre', 'gamma-normal', 'mu0', 'mispre'), 
            family = "Gamma", link = "log",
            type = 'GLMM', do.true = do.true)
}
b <- Sys.time()
c <- difftime(b, a, units = "mins")
save(c, file = "results/randomwalk_local_reps-441.RData")

### simpleGLMM with 5 groups===================================================
## possible mispecifications: 'missre', 'nb-pois', 'mispre', 'missunifcov', 'misscovnorm'
reps <- 1:441
a <- Sys.time()
for(dt in 1:2){
  if(dt == 1){do.true <- TRUE}
  if(dt == 2){do.true <- FALSE}
  dharma.methods <- c('uncond', 'cond', 'uncond_nrot', 'cond_nrot' )
  osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc','pears')
  run_model(reps, n = 20, ng = 5, mod='simpleGLMM', cov.mod = 'unif',
            misp = c('missre', 'missunifcov', 'mispre'), 
            type = 'LMM', do.true = do.true)
  osa.methods <- c('gen', 'cdf', 'mcmc','pears')
  run_model(reps, n = 20, ng = 5, mod='simpleGLMM',
            misp = c('missre', 'nb-pois', 'mispre'),
            type = 'GLMM',  family = "NB", link = "log",
            do.true = do.true)
}
b <- Sys.time()
save(c, file = "results/simpleGLMM_local_reps-441.RData")

### spatial ==================================================================                                                                                    ### Simple spatial SPDE model ==================================================
## possible mispecifications: 'missre', 'pois-zip', 'mispre', 'normal-gamma'
dharma.methods <- c('uncond', 'cond', 'uncond_nrot', 'cond_nrot' )
osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc', 'pears')
run_model(reps, n = 100, mod ='spatial',
          misp = c('missre', 'ln-error', 'mispre'), 
          type = 'LMM', do.true = do.true)
osa.methods <- c('gen', 'cdf', 'mcmc','pears')
run_model(reps, n = 100, mod ='spatial',
          misp = c('missre', 'pois-zip', 'mispre'),
          family = "Poisson", link = "log",
          type = 'GLMM', do.true = do.true)

## stop clusters
sfStop()


