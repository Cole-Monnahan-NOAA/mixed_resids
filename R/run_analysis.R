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
test <- run_model(reps, mod='linmod', misp='overdispersion', do.true = do.true)



### Random walk from the paper =================================================
## possible mispecifications: mu0, outliers, normal
osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc','pears')
dharma.methods <- c('uncond', 'cond', 'uncond_nrot', 'cond_nrot' )
run_model(reps, mod='randomwalk', misp='mu0', do.true = do.true)

### simpleGLMM with 5 groups===================================================
## possible mispecifications: overdispersion,  missnormcov, missunifcov, deltagamma
##! outliers not set up correctly when distribution not normal (lognormal better misp?)
##! misp cannot be overdispersion when fam = Tweedie
osa.methods <- c('fg', 'osg', 'cdf', 'mcmc', 'pears')
dharma.methods <- c('uncond', 'cond', 'uncond_nrot', 'cond_nrot' )
run_model(reps, ng = 5, mod='simpleGLMM', misp='missnormcov', do.true = do.true)
run_model(reps, ng = 5, mod='simpleGLMM', misp='missunifcov', do.true = do.true)
osa.methods <- c('mcmc', 'pears')
dharma.methods <- c('uncond', 'cond', 'uncond_nrot', 'cond_nrot' )
run_model(reps, ng = 5, mod='simpleGLMM', misp='dropRE',
          family = "Tweedie", link = "log", do.true = do.true)

### Simple spatial SPDE model ==================================================
## possible mispecifications: overdispersion, misscov, mispomega,
## dropRE, aniso - outliers not set up correctly when
## distribution not normal
osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc', 'pears')
dharma.methods <- c('uncond', 'cond', 'uncond_nrot', 'cond_nrot' )
run_model(reps, n=100, mod='spatial', misp='mispomega', do.true = do.true)
#h1 mispomega models occasionally fail b/c the exp(omega) leads to convergence issues.
#Repeat failed model runs with new seed for h0 and h1:

pvals <- lapply(list.files('results', pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows
pvals.true <- dplyr::filter(pvals, model == "spatial" & misp == "mispomega" &
                         do.true == TRUE & version == "h1")
pvals.est <- dplyr::filter(pvals, model == "spatial" & misp == "mispomega" &
                             do.true == FALSE & version == "h1")
true.reps <- unique(pvals.true$replicate)
##bad reps
idx <- which(!(true.reps %in% unique(pvals.est$replicate) ))
#delete from spatial_mispomega_true...
metric.nms <- c("stats", "resids", "pvals", "mles")
for(i in seq_along(metric.nms)){
  unlink(paste0("results/spatial_true_mispomega_",
                metric.nms[i],
                "/",
                metric.nms[i],
                "_",
                idx,
                ".RDS"))
}
## Generate new seeds
reps <- (max(true.reps)+1):(max(true.reps)+length(idx)) 
rm(pvals, idx, i, metric.nms)
## Rerun new reps for do.true = TRUE and do.true = FALSE
run_model(reps, n=100, mod='spatial', misp='mispomega', do.true = TRUE)
run_model(reps, n=100, mod='spatial', misp='mispomega', do.true = FALSE)
run_model(reps, n=100, mod='spatial', misp='misscov', cov.mod = 'unif', do.true = do.true)

osa.methods <- c('gen', 'cdf', 'mcmc', 'pears') #only 'cdf' and 'gen' suitable for discrete distributions
run_model(reps, n=100, mod='spatial', misp='dropRE',
          family = "Poisson", link = "log", do.true = do.true)

## stop clusters
sfStop()


