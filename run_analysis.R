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

## ## Basic random walk time series model
## Nreps <- 1000
## source('code/randomwalk.R')

## Spatial model. For now this script does everything. We can
## restructure later.

## ### Test a single iteration for debugging purposes
## TMB::compile('models/spatial.cpp')
## run.spatial.iter(1)
cpus <- parallel::detectCores()-2
Nreps <- cpus*5
source('code/spatial.R')

## ### Load results
## source('code/load_results.R')
## ### Make plots of results
## source("code/plot_results.R")

