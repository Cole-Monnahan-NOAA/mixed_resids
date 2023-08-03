## Loop over varying sample sizes to check that effect on null
## distributions

source("R/startup.R")

cpus <- parallel::detectCores()-1

Nreps <- 100
do.true <- FALSE
osa.methods <- c('fg', 'osg', 'gen', 'cdf', 'mcmc')
dharma.methods <- c('uncond', 'cond',
                    'uncond_nrot', 'cond_nrot' )

### randomwalk
runtimes <- mles <- pvals <- list(); k <- 1
(nobsvec <- 2^c(5:11))
for(nobs in nobsvec){
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
## this one is a bit different since it has ngroups and nobs per
## group. Just increasing ngroups and leaving nobs the same (10), but
## from the residual standpoint I think it's ngroups*nobs that
## matters
runtimes <- mles <- pvals <- list(); k <- 1
# fixing group size or fixing obs produce the same results, just need to run one
# (ngroupsvec <- 2^c(5:11))
# nobs <- 8
# for(nxng in ngroupsvec){
#   sfInit( parallel=cpus>1, cpus=cpus )
#   sfExportAll()
#   tmp <- sfLapply(1:Nreps, function(ii)
#     run_iter(ii, n=nobs, ng=nxng/nobs, mod='simpleGLMM', cov.mod='norm',
#              misp='missunifcov', family = "Gaussian", link = "identity",
#              do.true=do.true, savefiles=FALSE))
#   pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nxng))
#   mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nxng))
#   runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nxng))
#   k <- k+1
#   results.simpleGLMM <- process_results(mles, runtimes, pvals,
#                                         model='simpleGLMM', misp='missunifcov',
#                                         vary = "grps")
#   sfStop()
# }

(nobsvec <- 2^c(5:11))
ngroups <- 4
for(nxng in nobsvec){
  sfInit( parallel=cpus>1, cpus=cpus )
  sfExportAll()
  tmp <- sfLapply(1:Nreps, function(ii)
    run_iter(ii, n=nxng/ngroups, ng=ngroups, mod='simpleGLMM', cov.mod='norm',
             misp='missunifcov', family = "Gaussian", link = "identity",
             do.true=do.true, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nxng))
  mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nxng))
  runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nxng))
  k <- k+1
  results.simpleGLMM <- process_results(mles, runtimes, pvals,
                                        model='simpleGLMM', misp='missunifcov',
                                        vary = "obs")
  sfStop()
}
##plot_sample_sizes(results.simpleGLMM)

### Spatial
runtimes <- mles <- pvals <- list(); k <- 1
(nobsvec <- 2^c(5:11))
for(nobs in nobsvec){
  sfInit( parallel=cpus>1, cpus=cpus )
  sfExportAll()
  tmp <- sfLapply(1:Nreps, function(ii)
    run_iter(ii, n=nobs, ng = 0, mod='spatial', cov.mod='norm',
             misp='mispomega', family = "Gaussian", link = "identity",
             do.true=do.true, savefiles=FALSE))
  pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nobs))
  mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nobs))
  runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nobs))
  k <- k+1
  results.spatial <- process_results(mles, runtimes, pvals,
                                     model='spatial',  misp='mispomega')
  sfStop()
}
plot_sample_sizes(results.spatial)

## bunch of machinery here to look at more than runtimes which is
## turned off for now

### linmod -- probably doesn't make sense to include these?
# osa.methods <- c('fg', 'osg', 'gen', 'cdf')
# dharma.methods <- c('uncond_nrot', 'cond_nrot')
# runtimes <- mles <- pvals <- list(); k <- 1
# (nobsvec <- 2^c(5:11))
# for(nobs in nobsvec){
#   sfInit( parallel=cpus>1, cpus=cpus )
#   sfExportAll()
#   tmp <- sfLapply(1:Nreps, function(ii)
#     run_iter(ii, n=nobs, mod='linmod', cov.mod='norm',
#              misp='overdispersion', family = "Gaussian",
#              link = "identity", do.true=do.true, savefiles=FALSE))
#
#   pvals[[k]] <- lapply(tmp, function(x) get.value(x, 'pvals', nobs))
#   mles[[k]] <- lapply(tmp, function(x)  get.value(x, 'mles', nobs))
#   runtimes[[k]] <- lapply(tmp, function(x)  get.value(x, 'stats', nobs))
#   k <- k+1
#   results.linmod <- process_results(mles, runtimes, pvals,
#                                     model='linmod', misp = 'overdispersion')
#   sfStop()
# }
## plot_sample_sizes(results.linmod)


#
# Output moved to Fig-Table Rmarkdown
#
# ## Combine together to make runtime plots
# results.simpleGLMM.obs <- readRDS('results/simpleGLMM_missunifcov_obs_sample_sizes.RDS')
# results.simpleGLMM.grps <- readRDS('results/simpleGLMM_missunifcov_grps_sample_sizes.RDS')
#
# results.linmod <- readRDS('results/linmod_overdispersion_sample_sizes.RDS')
# results.randomwalk <- readRDS('results/randomwalk_mu0_sample_sizes.RDS')
# results.spatial <- readRDS('results/spatial_mispomega_sample_sizes.RDS')
# runtimes.all <- rbind(## results.linmod$runtimes,
#                       results.randomwalk$runtimes,
#                       results.simpleGLMM.obs$runtimes,
#                      results.simpleGLMM.grps$runtimes,
#                       results.spatial$runtimes)
# ## runtimes.all <- rbind(results.linmod, results.randomwalk,
# ##                       results.spatial, results.simpleGLMM)
# runtimes.all <- runtimes.all %>% filter(!is.na(med))
#
# g <- runtimes.all %>%
#  # dplyr::filter(type == "cdf" | type == "cond" | type == "gen" |
#   #                type == "osg") %>%
#   ggplot(.,aes(nobs, med, ymin=lwr, ymax=upr,  color=type)) +
#   geom_line()+
#   geom_pointrange(fatten=2) + scale_y_log10()+ scale_x_log10()+
#   facet_wrap('model', scales='free', ncol=1)+
#   labs(y='runtime (s)') +
#   theme_classic() +
#   scale_colour_viridis_d()
# g
#
#
# ggsave('plots/runtimes.png', g, width=5, height=7)
#
# # Type I error
# results.linmod$pvals %>%
#   filter(version == "h0") %>%
#   group_by(nobs, method) %>%
#   summarise(t1_err = sum(pvalue<0.05)/sum(pvalue>=0)) %>%
#   ggplot(aes(x = nobs, y = t1_err, color = method)) +
#   geom_point() +
#   facet_wrap(~method) +
#   theme_classic()
#
# # Power
# pow <- results.linmod$pvals %>%
#   filter(version == "h1") %>%
#   group_by(nobs, method) %>%
#   summarise(power = sum(pvalue<=0.05)/sum(pvalue>=0))
# pow %>%
#   ggplot(aes(x = nobs, y = power, color = method)) +
#   geom_line()  +
#   facet_wrap(~method) +
#   theme_classic()
#
# ## randomwalk
# # Type I error
# results.randomwalk$pvals %>%
#   filter(version == "h0") %>%
#   group_by(nobs, method) %>%
#   summarise(t1_err = sum(pvalue<0.05)/sum(pvalue>=0)) %>%
#   ggplot(aes(x = nobs, y = t1_err, color = method)) +
#   geom_line() +
#   facet_wrap(~method) +
#   theme_classic()
#
# # Power
# pow <- results.randomwalk$pvals %>%
#   filter(version == "h1") %>%
#   group_by(nobs, method) %>%
#   summarise(power = sum(pvalue<=0.05)/sum(pvalue>=0))
# pow %>%
#   ggplot(aes(x = nobs, y = power, color = method)) +
#   geom_line()  +
#   facet_wrap(~method) +
#   theme_classic()
#
# ## simpleGLMM
# # Type I error
# results.simpleGLMM.obs$pvals %>%
#   filter(version == "h0") %>%
#   group_by(nobs, method) %>%
#   summarise(t1_err = sum(pvalue<0.05)/sum(pvalue>=0)) %>%
#   ggplot(aes(x = nobs, y = t1_err, color = method)) +
#   geom_line() +
#   facet_wrap(~method) +
#   theme_classic()
# results.simpleGLMM.grps$pvals %>%
#   filter(version == "h0") %>%
#   group_by(nobs, method) %>%
#   summarise(t1_err = sum(pvalue<0.05)/sum(pvalue>=0)) %>%
#   ggplot(aes(x = nobs, y = t1_err, color = method)) +
#   geom_line() +
#   facet_wrap(~method) +
#   theme_classic()
#
# # Power
# pow <- results.simpleGLMM.obs$pvals %>%
#   filter(version == "h1") %>%
#   group_by(nobs, method, misp) %>%
#   summarise(power = sum(pvalue<=0.05)/sum(pvalue >=0))
# pow %>%
#   ggplot(aes(x = nobs, y = power, color = method)) +
#   geom_line()  +
#   facet_wrap(~method) +
#   theme_classic()
# pow <- results.simpleGLMM.grps$pvals %>%
#   filter(version == "h1") %>%
#   group_by(nobs, method, misp) %>%
#   summarise(power = sum(pvalue<=0.05)/sum(pvalue >=0))
# pow %>%
#   ggplot(aes(x = nobs, y = power, color = method)) +
#   geom_line()  +
#   facet_wrap(~method) +
#   theme_classic()
