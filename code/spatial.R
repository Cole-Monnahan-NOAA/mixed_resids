### Run simulation testing of spatial residuals for a simple
### spatial model.

## run.spatial.iter and other files are in startup.R
if(!exists('run.spatial.iter'))
  stop("Source code/startup.R before running this")

message("Preparing workspace to run ", Nreps, " iterations in parallel...")
TMB::compile("models/spatial.cpp") # modified for simulation
sfInit( parallel=TRUE, cpus=cpus )
sfExport('run.spatial.iter', 'sim.omega', 'cMatern', 'sim.data')

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




message("Processing and saving final results...")
## Read results back in from file
pvals <- lapply(fs, readRDS) %>% bind_rows  %>%
  filter(!is.na(pvalue))
saveRDS(pvals, file='results/spatial_pvals.RDS')
## Read in residuals
fs <- list.files('results/spatial_resids/', full.names=TRUE)
resids <- lapply(fs, readRDS) %>% bind_rows
saveRDS(resids, file='results/spatial_resids.RDS')


message("Making spatial results plots...")
g <- ggplot(filter(pvals, test=='outlier') , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_outlier.png', g, width=5, height=5)
g <- ggplot(filter(pvals, test=='disp') , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_disp.png', g, width=5, height=5)
g <- ggplot(filter(pvals, test=='sac') , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_sac.png', g, width=5, height=5)
g <- ggplot(filter(results, test=='GOF'& (RE=='cond' | RE=='uncond')) , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_GOF_DHARMa.png', g, width=5, height=5)
g <- ggplot(filter(results, test=='GOF'& RE!='cond' & RE!='uncond') , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')

## resids.long <- pivot_longer(resids, c('osa', 'sim_cond', 'sim_uncond'),
##                             names_to='type',
##                             values_to='residual') %>% filter(replicate==1)
## ggplot(resids.long, aes(ytrue, residual, color=maxgrad>.01)) + geom_point() +
##   facet_grid(version~type) + scale_x_log10()
