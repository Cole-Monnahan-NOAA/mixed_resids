## Placeholder script for analysis of residuals of spatial models



### Startup
source("code/startup.R")
packageVersion('TMB')
packageVersion('VAST')
packageVersion('FishStatsUtils')


### Some specific scripts to run. Each script runs a separate
### model, computes different residual types, and saves these
### results into the results folder for post-hoc analysis.

## Basic random walk time series model
source('code/randomwalk.R')

## Add others. VAST, TMB spatial?


### Load results
source('code/load_results.R')


### Make plots of results
source("code/plot_results.R")

