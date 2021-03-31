## Placeholder script for analysis of residuals of spatial
## models. Started 11/2020 by Cole.

## Assumes working directory is root of the repo spatial_resids

### Startup
source("code/startup.R")
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

## ### Test a single iteration for debugging purposes
## TMB::compile('models/spatial.cpp')
## run.spatial.iter(1)
cpus <- parallel::detectCores()-2
Nreps <- cpus*50
source('code/run_spatial.R')


## Random walk from the paper
cpus <- parallel::detectCores()-2
Nreps <- cpus*100
source('code/run_randomwalk.R')

## Andrea's simple GLMM with 3 groups
cpus <- parallel::detectCores()-2
Nreps <- cpus*20
source('code/run_simpleGLMM.R')

## Simple linear model as sanity check. Some resid methods not
## applicable b/c no random effects
cpus <- parallel::detectCores()-2
Nreps <- cpus*200
source('code/run_linmod.R')


### Load results into workspace
source('code/load_results.R')
### Make quick plots of results
source("code/make_plots.R")



## Compare GOF across models
pvals <- bind_rows(randomwalk_pvals,
          linmod_pvals,
          spatial_pvals)
pvals <- pvals %>% filter(test=='GOF' )
pvals <- filter(pvals, (version=='m1' & model=='randomwalk') | version=='m0')
g <- ggplot(pvals, aes(pvalue)) + geom_histogram() +
  facet_grid(method~model, scales='free_y')
ggsave("plots/GOF_by_model.png", g, width=7, height=5)

## Quick test of nonparameteric tests on real normal samples
library(goftest)
Nsim <- 15000
Nreps <- c(10,100,1000)
mu <- 0; sigma <- 1
results <- list(); k <- 1
for(j in Nreps){
  for(i in 1:Nsim){
    x <- rnorm(j, mu, sigma)
    p.ks <- ks.test(x,'pnorm')$p.value
    p.ad <- ad.test(x, 'pnorm', estimated = TRUE)$p.value
    results[[k]] <- data.frame(Nreps=j, ks=p.ks, ad=p.ad)
    k <- k+1
  }
}
results.long <- bind_rows(results) %>%
  pivot_longer(-c(Nreps), names_to='test', values_to='pvalue')

g <- ggplot(results.long, aes(pvalue)) + geom_histogram() +
  facet_grid(Nreps~test, scales='free_y')
ggsave('plots/ks_vs_ad.png', g, width=7, height=5)
