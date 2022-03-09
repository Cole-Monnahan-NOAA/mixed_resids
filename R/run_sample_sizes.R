## Loop over varying sample sizes to check that effect on null
## distributions

source("R/startup.R")


## quick local functions to streamline code a bit
extract_runtime <- function(stats){
  select(stats, model, replicate, version, starts_with('runtime'))
}
process_results <- function(mles, runtimes, pvals, model){
  if(!is.null(mles))
    mles <- bind_rows(mles) %>% filter(h==0) #filter(version=='h0')
  runtimes <- bind_rows(runtimes) %>% filter(version=='h0') %>%
    pivot_longer(starts_with('runtime'), names_to='type',
                 values_to='runtime') %>%
    mutate(type=gsub('runtime.|runtime_', '', type)) %>%
    group_by(nobs, model, version, type) %>%
    summarize(med=median(runtime, na.rm=TRUE),
              lwr=quantile(runtime, .1, na.rm=TRUE),
              upr=quantile(runtime, .9, na.rm=TRUE),
              pct.na=sum(is.na(runtime)),
              n=length(runtime), .groups='drop')
  if(!is.null(pvals))
    pvals <- bind_rows(pvals) %>%
      filter(version=='h0' & grepl('GOF', test)) %>%
      mutate(type=gsub('GOF.','',test))
  results <- list(mles=mles, pvals=pvals, runtimes=runtimes, model=model)
  saveRDS(results,
          file=paste0('results/',model,'_sample_sizes.RDS'))
  return(results$runtimes)
}
plot_sample_sizes <- function(results){
  model <- results$model
  g <- ggplot(filter(results$pvals, !is.na(pvalue)), aes(pvalue, fill=type)) + facet_grid(nobs~method) +
    geom_histogram(position='identity', alpha=.5, bins=20)
  ggsave(paste0("plots/", model,"_pvals_by_dim.png"), g, width=8, height=5)
  g <- ggplot(results$runtimes,
              aes(nobs, med, ymin=lwr, ymax=upr,  color=type)) +
    geom_line()+
    geom_pointrange(fatten=2) + scale_y_log10()+labs(y='runtime (s)')
  ##  ggsave(paste0("plots/", model,"_runtimes_by_dim.png"), g, width=7, height=5)
  g <- ggplot(results$mles, aes(factor(nobs), mle-true)) +
    geom_violin() +
    geom_hline(yintercept=0, color='red') +
    facet_wrap('par') + labs(y='Absolute error')
  ggsave(paste0("plots/", model,"_mles_by_dim.png"),g, width=8, height=5)
}
get.value <- function(x, val, nobs){
  if(is.null(x)) return(NULL)
  if(val!='runtimes')
    data.frame(nobs=nobs, x[[val]])
  else
    data.frame(nobs=nobs, extract_runtime(x[['stats']]))
}



cpus <- parallel::detectCores()-1
Nreps <- 10
do.true <- FALSE
osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc', 'pears')[-c(2,3,5)]
dharma.methods <- c('uncond', 'cond')

## bunch of machinery here to look at more than runtimes which is
## turned off for now

### linmod
runtimes <- mles <- pvals <- list(); k <- 1
(nobsvec <- 2^c(4:8))
for(nobs in nobsvec){
  sfInit( parallel=cpus>1, cpus=cpus )
  sfExportAll()
  tmp <- sfLapply(1:Nreps, function(ii)
    run_iter(ii, n=nobs, ng, mod='linmod', cov.mod='norm', misp='overdispersion', do.true=do.true, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nobs))
  mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nobs))
  runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'runtimes', nobs))
  k <- k+1
}
results.linmod <- process_results(mles, runtimes, pvals, model='linmod')
## plot_sample_sizes(results.linmod)

runtimes <- mles <- pvals <- list(); k <- 1
(nobsvec <- 2^c(4:8))
for(nobs in nobsvec){
  sfExportAll()
  tmp <- sfLapply(1:Nreps, function(ii)
    run_iter(ii, n=nobs, ng, mod='randomwalk', cov.mod='norm',
             misp='mu0', do.true=do.true, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nobs))
  mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nobs))
  runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nobs))
  k <- k+1
}
results.randomwalk <- process_results(mles, runtimes, pvals, model='randomwalk')
plot_sample_sizes(results.randomwalk)


runtimes <- mles <- pvals <- list(); k <- 1
(nobsvec <- 2^c(6:7))
for(nobs in nobsvec){
  sfExportAll()
  tmp <- sfLapply(1:Nreps, function(ii)
    run_iter(ii, n=nobs, ng, mod='spatial', cov.mod='norm',
             misp='mispomega', do.true=do.true, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nobs))
  mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nobs))
  runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nobs))
  k <- k+1
}
results.spatial <- process_results(mles, runtimes, pvals, model='spatial')
#plot_sample_sizes(results.spatial)


## this one is a bit different since it has ngroups and nobs per
## group. Just increasing ngroups and leaving nobs the same (10), but
## from the residual standpoint I think it's ngroups*nobs that
## matters
runtimes <- mles <- pvals <- list(); k <- 1
(ngroupsvec <- 2^c(4:6))
for(ngroups in ngroupsvec){
  nobs <- ngroups*10
  sfExportAll()
  tmp <- sfLapply(1:Nreps, function(ii)
    run_iter(ii, n=nobs, ng=10, mod='simpleGLMM', cov.mod='norm',
             misp='misscov', do.true=do.true, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nobs))
  mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nobs))
  runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nobs))
  k <- k+1
}
results.simpleGLMM <- process_results(mles, runtimes, pvals, model='simpleGLMM')
##plot_sample_sizes(results.simpleGLMM)



## Combine together to make runtime plots
results.simpleGLMM <- readRDS('results/simpleGLMM_sample_sizes.RDS')
results.linmod <- readRDS('results/linmod_sample_sizes.RDS')
results.randomwalk <- readRDS('results/randomwalk_sample_sizes.RDS')
results.spatial <- readRDS('results/spatial_sample_sizes.RDS')

runtimes.all <- rbind(results.simpleGLMM$runtimes,
                      results.linmod$runtimes,
                      results.randomwalk$runtimes,
                      results.spatial$runtimes)

runtimes.all <- rbind(results.linmod, results.randomwalk,
                      results.spatial, results.simpleGLMM) %>% filter(!is.na(med))
g <- ggplot(runtimes.all,
            aes(nobs, med, ymin=lwr, ymax=upr,  color=type)) +
  geom_line()+
  geom_pointrange(fatten=2) + scale_y_log10()+ scale_x_log10()+
  facet_wrap('model', scales='free')+
  labs(y='runtime (s)')
g


ggsave('plots/runtimes.png', g, width=7, height=5)

