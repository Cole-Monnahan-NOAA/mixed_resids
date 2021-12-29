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

do.true <- FALSE
osa.methods <- NULL#c('fg', 'osg', 'gen', 'cdf', 'mcmc')[-3]
dharma.methods <- NULL#c('uncond', 'cond')

## Simple linear model as sanity check. Some resid methods not
## applicable b/c no random effects

## possible mispecifications: overdispersion, outliers, misscov

run_model(reps, mod='linmod', misp='overdispersion', do.true = do.true)
## Random walk from the paper
## possible mispecifications: mu0, outliers
run_model(reps, mod='randomwalk', misp='mu0', do.true = do.true)
## Andrea's simple GLMM with 5 groups
## possible mispecifications: overdispersion, outliers, misscov
run_model(reps, ng = 5, mod='simpleGLMM', misp='misscov', do.true = do.true)
## Simple spatial SPDE model
## possible mispecifications: overdispersion, outliers, misscov, mispomega
#Turn off generic method - takes too long
osa.methods <- c('cdf', 'mcmc') #only 'cdf' and 'gen' suitable for discrete distributions
run_model(reps, n=200, mod='spatial', misp='overdispersion', do.true = do.true)
run_model(reps, n=200, mod='spatial', cov.mod = 'unif', misp='misscov', do.true = do.true)
run_model(reps, n=200, mod='spatial', misp='mispomega', do.true = do.true)

# stop clusters
# sfStop()

pvals <- lapply(list.files('results', pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows

mles <- lapply(list.files('results', pattern='_mles.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows

resids <- lapply(list.files('results', pattern='_resids.RDS',
                          full.names=TRUE), readRDS) %>% bind_rows

stats <- lapply(list.files('results', pattern='_stats.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows

## Effect of do.true for true model
filter(pvals, version=='h0' &  test=='GOF.ks') %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(model~method, scales = "free_y") + geom_histogram(position='identity', alpha=.5)

filter(pvals, version=='h1' &  test=='GOF.ks') %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(model~method, scales = "free_y") + geom_histogram(position='identity', alpha=.5)

## Model version for KS test
pvals %>% filter(test== 'GOF.ks' & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=version)) +
  facet_grid(model~method, scales = "free_y") + geom_histogram(position='identity', alpha=.5)

## KS vs AD for true model - delete create new plots for each test
#pvals %>% filter(version=='h0'& test!= 'outlier' & do.true==TRUE) %>%
#  ggplot(aes(pvalue, fill=test)) +
#  facet_grid(model~method) + geom_histogram(position='identity', alpha=.5)
sp.pvals <- pvals[pvals$model=='spatial',]
sp.stats <- stats[stats$model=='spatial',]
nc.reps <- filter(sp.stats, !is.na(converge) & converge == 1)$replicate
# filter out models that didn't converge
sp.pvals <- sp.pvals[!(sp.pvals$replicate %in% nc.reps),]

sp.pvals %>% filter(do.true==TRUE & test == 'GOF.ks') %>%
  ggplot(aes(pvalue, fill=version)) +
  facet_grid(test~method) + geom_histogram(position='identity', alpha=.5)

sp.pvals %>% filter(do.true==FALSE & test == 'GOF.ks') %>%
  ggplot(aes(pvalue, fill=version)) +
  facet_grid(misp~method, scales = "free_y") + geom_histogram(position='identity', alpha=.5)

sp.pvals %>% filter(do.true==FALSE & test == 'SAC') %>%
  ggplot(aes(pvalue, fill=version)) +
  facet_grid(misp~method, scales = "free_y") + geom_histogram(position='identity', alpha=.5)

sp.pvals %>% filter(do.true==FALSE & test == 'outliers') %>%
  ggplot(aes(pvalue, fill=version)) +
  facet_grid(misp~method, scales = "free_y") + geom_histogram(position='identity', alpha=.5)


## Model version for KS test
pvals %>% filter(test== 'GOF.ks' & do.true==TRUE) %>%
  ggplot(aes(pvalue, fill=version)) +
  facet_grid(model~method) + geom_histogram(position='identity', alpha=.5)

#MLE plots
## mles %>% filter(model == 'spatial' & type == 'mle' &
##                   do.true == FALSE & h==0) %>%
##   ggplot(., aes(x=par, y=value)) + geom_violin() +
##   geom_hline(yintercept = -2) + geom_hline(yintercept = 1)

## Check for MLE consistency. True and estimated par names and lengths are
## different so need to hack this by assuming order is the same
## then calculating within a unique group
nc_id <- stats[stats$converge==1|stats$maxgrad>0.1,]$id
g <- filter(mles, h==0 & !(id %in% nc_id) & do.true == FALSE) %>%
  ggplot(., aes(x=par, y=rel_error)) + geom_violin() + 
  facet_wrap(~model, ncol=1, scales = 'free_y')
g

# g <- filter(mles, h==0 & !(id %in% nc_id)) %>%
#   mutate(value=if_else((grepl('ln', x=par)), exp(value), value)) %>%
#   filter(par!='sp.parm') %>%
#   group_by(type, h,  replicate, model, misp) %>%
#   mutate(parnum=1:n()) %>%
#   pivot_wider(id_cols=c(replicate, model, misp, parnum), names_from=type, values_from=value)%>%
#   mutate(abs_error=mle-true, rel_error=abs_error/true) %>% ungroup %>%
#   drop_na %>%
#   ggplot(aes(x=factor(parnum), rel_error)) +
#     facet_wrap(model~misp, scales='free_y') + geom_violin()
# g
# filter(g$data, rel_error > 1000)

saveRDS(mles, 'mles.RDS')
#! Not modified yet
# ## This script runs the same examples above but with varying
# ## sample sizes and tracks runtime and pvalues to see the effect
# ## of this dimension.
# (cpus <- parallel::detectCores()-2)
# Nreps <- 1000
# sfInit( parallel=TRUE, cpus=cpus )
# source('code/run_sample_sizes.R')

### Load results into workspace
source('code/load_results.R')
### Make quick plots of results
source("code/make_plots.R")




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
