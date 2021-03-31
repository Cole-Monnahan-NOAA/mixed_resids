## estimate and validate a random walk model with and without
## drift. based on code example provided by uffe h?gsbro thygesen
## and kasper kristensen, 2016, in the tmb package examples:
## randomwalkvalidation.r.

## modified starting 11/2020 by cole
TMB::compile("models/randomwalk.cpp")

## Clean up the old runs
unlink('results/randomwalk_resids', TRUE)
unlink('results/randomwalk_pvals', TRUE)
unlink('results/randomwalk_mles', TRUE)

message("Preparing workspace to run ", Nreps, " iterations in parallel...")
TMB::compile("models/randomwalk.cpp") # modified for simulation
sfInit( parallel=TRUE, cpus=cpus )
## sfExport('run.randomwalk.iter', 'sim.omega', 'cMatern', 'sim.data',
##          'rmvnorm_prec', 'add_aic')
sfExportAll()

message("Starting parallel runs...")
results <- sfLapply(1:Nreps, function(ii) run.randomwalk.iter(ii))

## ## Read results back in from file
## fs <- list.files('results/randomwalk_pvals/', full.names=TRUE)
## ## Sometimes they fail to run for some unkonwn reason so try
## ## rerunning those ones once
## if(length(fs)<Nreps){
##   message("Rerunning some failed runs...")
##   bad <- which.failed(Nreps)
##   results <- sfLapply(bad, function(ii) run.randomwalk.iter(ii))
##   fs <- list.files('results/randomwalk_pvals/', full.names=TRUE)
## }
## bad <- which.failed(Nreps)
## if(length(bad)>0) warning(length(bad), " runs failed")

message("Randomwalk: processing and saving final results...")
## Read results back in from file
fs <- list.files('results/randomwalk_pvals/', full.names=TRUE)
pvals <- lapply(fs, readRDS) %>% bind_rows %>% filter(!is.na(pvalue))
saveRDS(pvals, file='results/randomwalk_pvals.RDS')
## Read in residuals
fs <- list.files('results/randomwalk_resids/', full.names=TRUE)
resids <- lapply(fs, readRDS) %>% bind_rows
saveRDS(resids, file='results/randomwalk_resids.RDS')
fs <- list.files('results/randomwalk_mles/', full.names=TRUE)
mles <- lapply(fs, readRDS) %>% bind_rows
saveRDS(mles, file='results/randomwalk_mles.RDS')

