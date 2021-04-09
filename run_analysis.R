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

(cpus <- parallel::detectCores()-2)
Nreps <- 1000

## Simple linear model as sanity check. Some resid methods not
## applicable b/c no random effects
source('code/run_linmod.R')
## Random walk from the paper
source('code/run_randomwalk.R')
## Andrea's simple GLMM with 3 groups
source('code/run_simpleGLMM.R')
## Simple spatial SPDE model
source('code/run_spatial.R')


### Load results into workspace
source('code/load_results.R')
### Make quick plots of results
source("code/make_plots.R")

## Loop over varying sample sizes to check that effect
sfInit( parallel=T, cpus=cpus )
mles <- pvals <- list(); k <- 1
sfExportAll()
for(nobs in c(10,20, 40, 80)){
  tmp <- sfLapply(1:Nreps, function(ii)
    run.linmod.iter(ii, nobs, FALSE))
  pvals[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, x$pvals)) %>% bind_rows
  mles[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, x$mles)) %>% bind_rows
  k <- k+1
}
mles <- bind_rows(mles)
pvals <- bind_rows(pvals) %>%
  filter(grepl('GOF', test)) %>%
  mutate(type=gsub('GOF.','',test))
g <- ggplot(pvals, aes(pvalue, fill=type)) + facet_grid(nobs~method) +
  geom_histogram(position='identity', alpha=.5)
ggsave("plots/linmod_pvals_by_dim.png", width=7, height=7)

## ## Quick test of nonparameteric tests on real normal samples
## library(goftest)
## Nsim <- 15000
## Nreps <- c(10,100,1000)
## mu <- 0; sigma <- 1
## results <- list(); k <- 1
## for(j in Nreps){
##   for(i in 1:Nsim){
##     x <- rnorm(j, mu, sigma)
##     x <- rt(j, 20)
##     p.ks <- ks.test(x,'pnorm')$p.value
##     p.ad <- ad.test(x, 'pnorm', estimated = TRUE)$p.value
##     results[[k]] <- data.frame(Nreps=j, ks=p.ks, ad=p.ad)
##     k <- k+1
##   }
## }
## results.long <- bind_rows(results) %>%
##   pivot_longer(-c(Nreps), names_to='test', values_to='pvalue')
## g <- ggplot(results.long, aes(pvalue)) + geom_histogram() +
##   facet_grid(Nreps~test, scales='free_y')
## ggsave('plots/ks_vs_ad.png', g, width=7, height=5)
