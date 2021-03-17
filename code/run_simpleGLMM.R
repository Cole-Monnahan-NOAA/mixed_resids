


## Clean up the old runs
unlink('results/simpleGLMM_resids', TRUE)
unlink('results/simpleGLMM_pvals', TRUE)
unlink('results/simpleGLMM_mles', TRUE)

message("Preparing workspace to run ", Nreps, " iterations in parallel...")
TMB::compile("models/simpleGLMM.cpp") # modified for simulation
sfInit( parallel=TRUE, cpus=cpus )
## sfExport('run.simpleGLMM.iter', 'sim.omega', 'cMatern', 'sim.data',
##          'rmvnorm_prec', 'add_aic')
sfExportAll()

message("Starting parallel runs...")
results <- sfLapply(1:Nreps, function(ii) run.simpleGLMM.iter(ii))

## ## Read results back in from file
## fs <- list.files('results/simpleGLMM_pvals/', full.names=TRUE)
## ## Sometimes they fail to run for some unkonwn reason so try
## ## rerunning those ones once
## if(length(fs)<Nreps){
##   message("Rerunning some failed runs...")
##   bad <- which.failed(Nreps)
##   results <- sfLapply(bad, function(ii) run.simpleGLMM.iter(ii))
##   fs <- list.files('results/simpleGLMM_pvals/', full.names=TRUE)
## }
## bad <- which.failed(Nreps)
## if(length(bad)>0) warning(length(bad), " runs failed")

message("SimpleGLMM: processing and saving final results...")
## Read results back in from file
fs <- list.files('results/simpleGLMM_pvals/', full.names=TRUE)
pvals <- lapply(fs, readRDS) %>% bind_rows %>% filter(!is.na(pvalue))
saveRDS(pvals, file='results/simpleGLMM_pvals.RDS')
## Read in residuals
fs <- list.files('results/simpleGLMM_resids/', full.names=TRUE)
resids <- lapply(fs, readRDS) %>% bind_rows
saveRDS(resids, file='results/simpleGLMM_resids.RDS')
fs <- list.files('results/simpleGLMM_mles/', full.names=TRUE)
mles <- lapply(fs, readRDS) %>% bind_rows
saveRDS(mles, file='results/simpleGLMM_mles.RDS')


