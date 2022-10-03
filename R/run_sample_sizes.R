## Loop over varying sample sizes to check that effect on null
## distributions

source("R/startup.R")


## quick local functions to streamline code a bit
extract_runtime <- function(stats){
  select(stats, model, misp, replicate, version, starts_with('runtime'))
}
process_results <- function(mles, runtimes, pvals, model, misp, vary = NULL){
  if(!is.null(mles))
    mles <- bind_rows(mles) %>% filter(h==0) #filter(version=='h0')
  runtimes <- bind_rows(runtimes) %>% filter(version=='h0') %>%
    pivot_longer(starts_with('runtime'), names_to='type',
                 values_to='runtime') %>%
    mutate(type=gsub('runtime.|runtime_', '', type)) %>%
    group_by(nobs, model, version, type) %>%
    summarize(med=median(runtime, na.rm=TRUE),
              lwr=quantile(runtime, .25, na.rm=TRUE),
              upr=quantile(runtime, .75, na.rm=TRUE),
              pct.na=sum(is.na(runtime)),
              n=length(runtime), .groups='drop')
  if(!is.null(pvals))
    pvals <- bind_rows(pvals) %>%
      filter(#version=='h0' & 
               grepl('GOF', test)) %>%
      mutate(type=gsub('GOF.','',test))
  results <- list(mles=mles, pvals=pvals, runtimes=runtimes, model=model, misp = misp)
  
  if(model == 'simpleGLMM'){
    filename <- paste0('results/',model,'_', misp,'_', vary, '_sample_sizes.RDS')
  } else {
    filename <- paste0('results/',model,'_', misp, '_sample_sizes.RDS')
  }
  saveRDS(results, file = filename)
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

Nreps <- 50
do.true <- FALSE
# osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc', 'pears')[-c(2,3,6)]
# dharma.methods <- c('uncond', 'cond')

## bunch of machinery here to look at more than runtimes which is
## turned off for now

### linmod -- probably doesn't make sense to include these?
osa.methods <- c('fg', 'osg', 'gen', 'cdf')
dharma.methods <- c('uncond_nrot', 'cond_nrot')
runtimes <- mles <- pvals <- list(); k <- 1
(nobsvec <- 2^c(4:12))
for(nobs in nobsvec){
  sfInit( parallel=cpus>1, cpus=cpus )
  sfExportAll()
  tmp <- sfLapply(1:Nreps, function(ii)
    run_iter(ii, n=nobs, mod='linmod', cov.mod='norm', 
             misp='overdispersion', family = "Gaussian", 
             link = "identity", do.true=do.true, savefiles=FALSE))
  
  pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nobs))
  mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nobs))
  runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nobs))
  k <- k+1
  results.linmod <- process_results(mles, runtimes, pvals, 
                                    model='linmod', misp = 'overdispersion')
  sfStop()
}
## plot_sample_sizes(results.linmod)

### randomwalk
osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc')
dharma.methods <- c('uncond', 'cond', 
                    'uncond_nrot', 'cond_nrot' )
runtimes <- mles <- pvals <- list(); k <- 1
(nobsvec <- 2^c(4:12))
for(nobs in nobsvec){
  if(nobs>600) osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc')[-3] #don't run gen for large models
  sfInit( parallel=cpus>1, cpus=cpus )
  sfExportAll()
  tmp <- sfLapply(1:Nreps, function(ii)
    run_iter(ii, n=nobs, mod='randomwalk', cov.mod='norm',
             misp='mu0', family = "Gaussian", link = "identity", 
             do.true=do.true, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nobs))
  mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nobs))
  runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nobs))
  k <- k+1
  results.randomwalk <- process_results(mles, runtimes, pvals, 
                                        model='randomwalk', misp = "mu0")
  sfStop()
}
## plot_sample_sizes(results.randomwalk)

###simpleGLMM
osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc')
dharma.methods <- c('uncond', 'cond', 
                    'uncond_nrot', 'cond_nrot' )
## this one is a bit different since it has ngroups and nobs per
## group. Just increasing ngroups and leaving nobs the same (10), but
## from the residual standpoint I think it's ngroups*nobs that
## matters
runtimes <- mles <- pvals <- list(); k <- 1
(ngroupsvec <- 2^c(4:11))
for(ngroups in ngroupsvec){
  sfInit( parallel=cpus>1, cpus=cpus )
  sfExportAll()
  tmp <- sfLapply(1:Nreps, function(ii)
    run_iter(ii, n=10, ng=ngroups, mod='simpleGLMM', cov.mod='norm',
             misp='missunifcov', family = "Gaussian", link = "identity",
             do.true=do.true, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nobs))
  mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nobs))
  runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nobs))
  k <- k+1
  results.simpleGLMM <- process_results(mles, runtimes, pvals, 
                                        model='simpleGLMM', misp='missunifcov',
                                        vary = "grps")
  sfStop()
}

(nobsvec <- 2^c(4:11))
for(nobs in nobsvec){
  sfInit( parallel=cpus>1, cpus=cpus )
  sfExportAll()
  tmp <- sfLapply(1:Nreps, function(ii)
    run_iter(ii, n=nobs, ng=5, mod='simpleGLMM', cov.mod='norm',
             misp='missunifcov', family = "Gaussian", link = "identity",
             do.true=do.true, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nobs))
  mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nobs))
  runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nobs))
  k <- k+1
  results.simpleGLMM <- process_results(mles, runtimes, pvals, 
                                        model='simpleGLMM', misp='missunifcov',
                                        vary = "obs")
  sfStop()
}
##plot_sample_sizes(results.simpleGLMM)

### Spatial
runtimes <- mles <- pvals <- list(); k <- 1
(nobsvec <- 2^c(6:10))
for(nobs in nobsvec){
  sfInit( parallel=cpus>1, cpus=cpus )
  sfExportAll()
  tmp <- sfLapply(1:Nreps, function(ii)
    run_iter(ii, n=nobs, ng, mod='spatial', cov.mod='norm',
             misp='mispomega', family = "Gaussian", link = "identity",
             do.true=do.true, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nobs))
  mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nobs))
  runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nobs))
  k <- k+1
  results.spatial <- process_results(mles, runtimes, pvals, 
                                     model='spatial',  misp='mispomega')
}
#plot_sample_sizes(results.spatial)





## Combine together to make runtime plots
results.simpleGLMM <- readRDS('results/simpleGLMM_sample_sizes.RDS')
results.linmod <- readRDS('results/linmod_overdispersion_sample_sizes.RDS')
results.randomwalk <- readRDS('results/randomwalk_sample_sizes.RDS')
results.spatial <- readRDS('results/spatial_sample_sizes.RDS')
runtimes.all <- rbind(results.simpleGLMM$runtimes,
                     ## results.linmod$runtimes,
                      results.randomwalk$runtimes,
                      results.spatial$runtimes)
## runtimes.all <- rbind(results.linmod, results.randomwalk,
##                       results.spatial, results.simpleGLMM)
runtimes.all <- runtimes.all %>% filter(!is.na(med))

g <- ggplot(runtimes.all,
            aes(nobs, med, ymin=lwr, ymax=upr,  color=type)) +
  geom_line()+
  geom_pointrange(fatten=2) + scale_y_log10()+ scale_x_log10()+
  facet_wrap('model', scales='free', ncol=1)+
  labs(y='runtime (s)')
g


ggsave('plots/runtimes.png', g, width=5, height=7)

#Type I error
results.linmod$pvals %>% 
  filter(version == "h0") %>% 
  group_by(nobs, method) %>% 
  summarise(t1_err = sum(pvalue<0.05)/5) %>%
  ggplot(aes(x = nobs, y = t1_err, color = method)) + geom_point()
pow <- results.linmod$pvals %>% 
  filter(version == "h1") %>% 
  group_by(nobs, method) %>% 
  summarise(power = sum(pvalue<=0.05)/5)
pow %>%
  ggplot(aes(x = nobs, y = power, color = method)) + geom_point()
results.linmod$pvals %>% ggplot(aes(x=nobs, y = pvalue)) + geom_point() + facet_grid(version~method)
