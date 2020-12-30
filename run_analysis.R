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
## pearson = Pearson (naive) residuals
## osa = one-step ahead residual
## sim_cond = simulation residuals (DHARMa) conditioned on RE
## sim_uncon= simulation residual (DHARMa) unconditioned on RE

## Basic random walk time series model
Nreps <- 100
source('code/randomwalk.R')

## Add others. VAST, TMB spatial?

cpus <- parallel::detectCores()-1
Nreps <- 1000
source('code/spatial.R')


### Load results
source('code/load_results.R')


### Make plots of results
source("code/plot_results.R")

