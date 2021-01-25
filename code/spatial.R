## estimate and validate a spatial model with and without
## drift. based on code example provided by uffe h?gsbro thygesen
## and kasper kristensen, 2016, in the tmb package examples:
## randomwalkvalidation.r.

## modified starting 11/2020 by cole
## modified to a spatial model by Andrea 12/2020

## run.spatial.iter is in startup.R

TMB::compile("models/spatial.cpp") # modified for simulation
sfInit( parallel=TRUE, cpus=cpus )
sfExport('run.spatial.iter', 'sim.omega', 'cMatern', 'sim.data')
results <- sfLapply(1:Nreps, function(ii) run.spatial.iter(ii))

## Read results back in from file
fs <- list.files('results/spatial_pvals', full.names=TRUE)
results <- lapply(fs, readRDS) %>% do.call(rbind, .) %>%
  filter(!is.na(pvalue))
## Did any fail to run? Try to rerun
bad <- which(!1:Nreps %in% results$replicate)
results <- sfLapply(bad, function(ii) run.spatial.iter(ii))

## Read results back in from file
fs <- list.files('results/spatial_pvals', full.names=TRUE)
results <- lapply(fs, readRDS) %>% do.call(rbind, .) %>%
  filter(!is.na(pvalue))
saveRDS(results, file='results/spatial_pvals.RDS')

g <- ggplot(filter(results, test=='outlier') , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_outlier.png', g, width=5, height=5)
g <- ggplot(filter(results, test=='disp') , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_disp.png', g, width=5, height=5)
g <- ggplot(filter(results, test=='sac') , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_sac.png', g, width=5, height=5)
g <- ggplot(filter(results, test=='GOF') , aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test, scales='free')
ggsave('plots/spatial_pvalues_GOF.png', g, width=5, height=5)


