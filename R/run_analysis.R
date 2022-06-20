## Placeholder script for analysis of residuals of spatial
## models. Started 11/2020 by Cole.

## Assumes working directory is root of the repo spatial_resids

### Startup
source("R/startup.R")
packageVersion('TMB')                   # 1.7.18
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

(cpus <- parallel::detectCores()-2)
reps <- 1:500

do.true <- TRUE
## Set osa.methods and dharma.methods to NULL to turn off
osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'pears')[-3]
dharma.methods <- c('uncond', 'cond')

## Simple linear model as sanity check. Some resid methods not
## applicable b/c no random effects
## possible mispecifications: overdispersion, outliers, misscov
run_model(reps, mod='linmod', misp='overdispersion', do.true = do.true)

## Random walk from the paper
## possible mispecifications: mu0, outliers
osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc', 're_mcmc', 'pears')
dharma.methods <- c('uncond', 'cond', 're_uncond', 
                    'cond_nrot', 'uncond_nrot', 're_uncond_nrot' )
run_model(reps, mod='randomwalk', misp='mu0', do.true = do.true)

## Andrea's simple GLMM with 5 groups
## possible mispecifications: overdispersion,  misscov 
##! outliers not set up correctly when distribution not normal (lognormal better misp?)
##! misp cannot be overdispersion when fam = Tweedie
osa.methods <- c('mcmc', 're_mcmc', 'pears')
dharma.methods <- c('uncond', 'uncond_nrot', 're_uncond_nrot', 'cond', 'cond_nrot' )
run_model(reps, ng = 5, mod='simpleGLMM', misp='deltagamma', do.true = do.true)

## Simple spatial SPDE model
## possible mispecifications: overdispersion, misscov, mispomega, dropRE - outliers not set up correctly when distribution not normal
#Turn off generic method - takes too long
osa.methods <- c('cdf', 'mcmc', 're_mcmc', 'pears', 'gen') #only 'cdf' and 'gen' suitable for discrete distributions
dharma.methods <- c('uncond', 'cond', 're_uncond', 
                    'uncond_nrot', 're_uncond_nrot' )
run_model(reps, n=200, mod='spatial', misp='overdispersion', do.true = do.true)
run_model(reps, n=200, mod='spatial', cov.mod = 'unif', misp='misscov', do.true = do.true)
run_model(reps, n=200, mod='spatial', misp='mispomega', do.true = do.true)
run_model(reps, n=200, mod='spatial', misp='dropRE', do.true = do.true)

# stop clusters
# sfStop()

