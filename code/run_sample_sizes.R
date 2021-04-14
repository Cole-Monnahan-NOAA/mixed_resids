## Loop over varying sample sizes to check that effect on null
## distributions

source("code/startup.R")


## quick local functions to streamline code a bit
extract_runtime <- function(resids){
  select(resids, model, replicate, version, starts_with('runtime'))
}
process_results <- function(mles, runtimes, pvals, model){
  mles <- bind_rows(mles) %>% filter(version=='m0')
  runtimes <- bind_rows(runtimes) %>% filter(version=='m0') %>%
    pivot_longer(starts_with('runtime'), names_to='type',
                 values_to='runtime') %>%
    mutate(type=gsub('runtime.|runtime_', '', type)) %>%
    group_by(nobs, model, version, type) %>%
    summarize(med=median(runtime, na.rm=TRUE),
              lwr=quantile(runtime, .1, na.rm=TRUE),
              upr=quantile(runtime, .9, na.rm=TRUE),
              pct.na=sum(is.na(runtime)),
              n=length(runtime), .groups='drop')
  pvals <- bind_rows(pvals) %>%
    filter(version=='m0' & grepl('GOF', test)) %>%
    mutate(type=gsub('GOF.','',test))
  results <- list(mles=mles, pvals=pvals, runtimes=runtimes, model=model)
  saveRDS(results, file=paste0('results/',model,'_sample_sizes.RDS'))
  return(results)
}
plot_sample_sizes <- function(results){
  model <- results$model
  g <- ggplot(results$pvals, aes(pvalue, fill=type)) + facet_grid(nobs~method) +
    geom_histogram(position='identity', alpha=.5, bins=20)
  ggsave(paste0("plots/", model,"_pvals_by_dim.png"), g, width=7, height=7)
  g <- ggplot(results$runtimes,
              aes(nobs, med, ymin=lwr, ymax=upr,  color=type)) +
    geom_line()+
    geom_pointrange(fatten=2) + scale_y_log10()+labs(y='runtime (s)')
  ##  ggsave(paste0("plots/", model,"_runtimes_by_dim.png"), g, width=7, height=5)
  g <- ggplot(results$mles, aes(factor(nobs), mle-true)) +
    geom_violin() +
    geom_hline(yintercept=0, color='red') +
    facet_wrap('par') + labs(y='Absolute error')
  ggsave(paste0("plots/", model,"_mles_by_dim.png"),g, width=7, height=7)
}




### linmod
runtimes <- mles <- pvals <- list(); k <- 1
sfExportAll()
for(nobs in c(10, 50, 100)){
  sfExport('nobs')
  tmp <- sfLapply(1:Nreps, function(ii)
    run.linmod.iter(ii, nobs=nobs, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, x$pvals)) %>% bind_rows
  mles[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, x$mles)) %>% bind_rows
  runtimes[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, extract_runtime(x$resids))) %>% bind_rows
  k <- k+1
}
results.linmod <- process_results(mles, runtimes, pvals, model='linmod')
plot_sample_sizes(results.linmod)



runtimes <- mles <- pvals <- list(); k <- 1
sfExportAll()
for(nobs in c(10, 50, 100)){
  sfExport('nobs')
  tmp <- sfLapply(1:Nreps, function(ii)
    run.randomwalk.iter(ii, nobs=nobs, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, x$pvals)) %>% bind_rows
  mles[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, x$mles)) %>% bind_rows
  runtimes[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, extract_runtime(x$resids))) %>% bind_rows
  k <- k+1
}
results.randomwalk <- process_results(mles, runtimes, pvals, model='randomwalk')
plot_sample_sizes(results.randomwalk)


runtimes <- mles <- pvals <- list(); k <- 1
sfExportAll()
for(nobs in c(25, 50, 100)){
  sfExport('nobs')
  tmp <- sfLapply(1:Nreps, function(ii)
    run.spatial.iter(ii, nobs=nobs, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, x$pvals)) %>% bind_rows
  mles[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, x$mles)) %>% bind_rows
  runtimes[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, extract_runtime(x$resids))) %>% bind_rows
  k <- k+1
}
results.spatial <- process_results(mles, runtimes, pvals, model='spatial')
plot_sample_sizes(results.spatial)


## this one is a bit different since it has ngroups and nobs per
## group. Just increasing ngroups and leaving nobs the same (10), but
## from the residual standpoint I think it's ngroups*nobs that
## matters
runtimes <- mles <- pvals <- list(); k <- 1
sfExportAll()
for(ngroups in c(5, 10, 25)){
  nobs <- ngroups*10
  sfExport('ngroups', 'nobs')
  tmp <- sfLapply(1:Nreps, function(ii)
    run.simpleGLMM.iter(ii, ngroups=ngroups, nobs=10,  savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, x$pvals)) %>% bind_rows
  mles[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, x$mles)) %>% bind_rows
  runtimes[[k]] <- lapply(tmp, function(x) cbind(nobs=nobs, extract_runtime(x$resids))) %>% bind_rows
  k <- k+1
}
results.simpleGLMM <- process_results(mles, runtimes, pvals, model='simpleGLMM')
plot_sample_sizes(results.simpleGLMM)



## Combine together to make runtime plots
results.simpleGLMM <- readRDS('results/simpleGLMM_sample_sizes.RDS')
results.linmod <- readRDS('results/linmod_sample_sizes.RDS')
results.randomwalk <- readRDS('results/randomwalk_sample_sizes.RDS')
results.spatial <- readRDS('results/spatial_sample_sizes.RDS')
runtimes.all <- rbind(results.simpleGLMM$runtimes,
                      results.linmod$runtimes,
                      results.randomwalk$runtimes,
                      results.spatial$runtimes)
g <- ggplot(runtimes.all,
            aes(nobs, med, ymin=lwr, ymax=upr,  color=type)) +
  geom_line()+
  geom_pointrange(fatten=2) + scale_y_log10()+ scale_x_log10()+
  facet_wrap('model', scales='free')+
  labs(y='runtime (s)')
ggsave('plots/runtimes.png', g, width=7, height=5)

