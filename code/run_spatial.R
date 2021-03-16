### Run simulation testing of spatial residuals for a simple
### spatial model.

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


message("Processing and saving final results...")
## Read results back in from file
fs <- list.files('results/spatial_pvals/', full.names=TRUE)
pvals <- lapply(fs, readRDS) %>% bind_rows %>% filter(!is.na(pvalue))
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

g <- filter(pvals, test=='GOF' & grepl('osa', x=RE)) %>%
  ggplot(aes(pvalue)) + geom_histogram() +
         facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_GOF_osa.png', g, width=5, height=5)
g <- filter(pvals, test=='GOF' & !grepl('osa', x=RE)) %>%
  ggplot(aes(pvalue)) + geom_histogram() +
         facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_GOF_DHARMa.png', g, width=5, height=5)

g <- pivot_longer(resids, c('osa.cdf', 'osa.gen', 'sim_cond', 'sim_uncond', 'sim_parcond'),
                  names_to='type', values_to='residual') %>%
  filter(replicate<=5) %>%
  ggplot(aes(ytrue, residual, color=version)) +
  geom_point(alpha=.5) +
  facet_grid(replicate~type) + scale_x_log10()
ggsave('plots/spatial_residuals_examples.png', g, width=8, height=6)
