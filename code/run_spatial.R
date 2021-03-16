### Run simulation testing of spatial residuals for a simple
### spatial model.

## Clean up the old runs
unlink('results/spatial_resids', TRUE)
unlink('results/spatial_pvals', TRUE)
unlink('results/spatial_mles', TRUE)

## run.spatial.iter and other files are in startup.R
source("code/functions_spatial.R")
if(!exists('run.spatial.iter'))
  stop("Problem loading spatial functions")

message("Preparing workspace to run ", Nreps, " iterations in parallel...")
TMB::compile("models/spatial.cpp") # modified for simulation
sfInit( parallel=TRUE, cpus=cpus )
## sfExport('run.spatial.iter', 'sim.omega', 'cMatern', 'sim.data',
##          'rmvnorm_prec', 'add_aic')
sfExportAll()

message("Starting parallel runs...")
results <- sfLapply(1:Nreps, function(ii) run.spatial.iter(ii))

## Read results back in from file
fs <- list.files('results/spatial_pvals/', full.names=TRUE)
## Sometimes they fail to run for some unkonwn reason so try
## rerunning those ones once
if(length(fs)<Nreps){
  message("Rerunning some failed runs...")
  bad <- which.failed(Nreps)
  results <- sfLapply(bad, function(ii) run.spatial.iter(ii))
  fs <- list.files('results/spatial_pvals/', full.names=TRUE)
}
bad <- which.failed(Nreps)
if(length(bad)>0) warning(length(bad), " runs failed")

message("Spatial: processing and saving final results...")
## Read results back in from file
fs <- list.files('results/spatial_pvals/', full.names=TRUE)
pvals <- lapply(fs, readRDS) %>% bind_rows %>% filter(!is.na(pvalue))
saveRDS(pvals, file='results/spatial_pvals.RDS')
## Read in residuals
fs <- list.files('results/spatial_resids/', full.names=TRUE)
resids <- lapply(fs, readRDS) %>% bind_rows
saveRDS(resids, file='results/spatial_resids.RDS')
fs <- list.files('results/spatial_mles/', full.names=TRUE)
mles <- lapply(fs, readRDS) %>% bind_rows
saveRDS(mles, file='results/spatial_mles.RDS')

