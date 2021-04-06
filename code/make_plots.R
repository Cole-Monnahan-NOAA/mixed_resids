bins <- 20

if(exists("spatial_pvals")){
message("Spatial: making results plots...")
## g <- ggplot(filter(spatial_pvals, test=='outlier') , aes(pvalue, )) + geom_histogram(bins=bins) +
##   facet_grid(method~test+version, scales='free')
## ggsave('plots/spatial_pvalues_outlier.png', g, width=5, height=5)
## g <- ggplot(filter(spatial_pvals, test=='disp') , aes(pvalue, )) + geom_histogram(bins=bins) +
##   facet_grid(method~test+version, scales='free')
## ggsave('plots/spatial_pvalues_disp.png', g, width=5, height=5)
g <- ggplot(filter(spatial_pvals, test=='sac') , aes(pvalue, )) + geom_histogram(bins=bins) +
  facet_grid(method~test+version, scales='free')
ggsave('plots/spatial_pvalues_sac.png', g, width=5, height=5)
g <- filter(spatial_pvals, grepl('GOF', test) & grepl('osa', x=method)) %>%
  mutate(test=gsub('GOF.', '', x=test), method=gsub('osa.', '', x=method))  %>%
  ggplot(aes(pvalue, fill=test)) + geom_histogram(bins=bins) +
  facet_grid(method~version, scales='free')
ggsave('plots/spatial_pvalues_GOF_osa.png', g, width=5, height=5)
g <- filter(spatial_pvals, grepl('GOF', test) & !grepl('osa', x=method)) %>%
  mutate(test=gsub('GOF.', '', x=test), method=gsub('osa.', '', x=method))  %>%
  ggplot(aes(pvalue, fill=test)) + geom_histogram(bins=bins) +
  facet_grid(method~version, scales='free')
ggsave('plots/spatial_pvalues_GOF_DHARMa.png', g, width=5, height=5)
## g <- pivot_longer(resids, c('osa.cdf', 'osa.gen', 'sim_cond', 'sim_uncond', 'sim_parcond'),
##                   names_to='type', values_to='residual') %>%
##   filter(replicate<=5) %>%
##   ggplot(aes(ytrue, residual, color=version)) +
##   geom_point(alpha=.5) +
##   facet_grid(replicate~type) + scale_x_log10()
## ggsave('plots/spatial_residuals_examples.png', g, width=8, height=6)
spatial_mles <- spatial_mles %>% mutate(abserror=true-mle,
                                        relerror=abserror/true)
g <- ggplot(spatial_mles, aes(x=version, abserror)) + geom_violin() +
  facet_wrap('par', scales='free') + geom_hline(yintercept=0,
  color='red') + labs(y='absolute error (MLE-truth)')
ggsave('plots/spatial_mles.png', g, width=ggwidth, height=ggheight)
}

if(exists('randomwalk_pvals')){
message("Randomwalk: making results plots...")
## g <- ggplot(filter(randomwalk_pvals, test=='outlier') , aes(pvalue, )) + geom_histogram(bins=bins) +
##   facet_grid(method~test+version, scales='free')
## ggsave('plots/randomwalk_pvalues_outlier.png', g, width=5, height=5)
## g <- ggplot(filter(randomwalk_pvals, test=='disp') , aes(pvalue, )) + geom_histogram(bins=bins) +
##   facet_grid(method~test+version, scales='free')
## ggsave('plots/randomwalk_pvalues_disp.png', g, width=5, height=5)

g <- filter(randomwalk_pvals, grepl('GOF', test) & grepl('osa', x=method)) %>%
  mutate(test=gsub('GOF.', '', x=test), method=gsub('osa.', '', x=method))  %>%
  ggplot(aes(pvalue, fill=test)) + geom_histogram(bins=bins) +
  facet_grid(method~version, scales='free')
ggsave('plots/randomwalk_pvalues_GOF_osa.png', g, width=5, height=5)
g <- filter(randomwalk_pvals, grepl('GOF', test) & !grepl('osa', x=method)) %>%
  mutate(test=gsub('GOF.', '', x=test), method=gsub('osa.', '', x=method))  %>%
  ggplot(aes(pvalue, fill=test)) + geom_histogram(bins=bins) +
  facet_grid(method~version, scales='free')
ggsave('plots/randomwalk_pvalues_GOF_DHARMa.png', g, width=5, height=5)
## g <- pivot_longer(randomwalk_resids, c('osa.cdf', 'osa.gen', 'sim_cond', 'sim_uncond', 'sim_parcond'),
##                   names_to='type', values_to='residual') %>%
##   filter(replicate<=5) %>%
##   ggplot(aes(ytrue, residual, color=version)) +
##   geom_point(alpha=.5) +
##   facet_grid(replicate~type) + scale_x_log10()
## ggsave('plots/randomwalk_residuals_examples.png', g, width=8, height=6)
randomwalk_mles <- randomwalk_mles %>% mutate(abserror=true-mle,
                                        relerror=abserror/true)
g <- ggplot(randomwalk_mles, aes(x=version, abserror)) + geom_violin() +
  facet_wrap('par', scales='free') + geom_hline(yintercept=0,
  color='red') + labs(y='absolute error (MLE-truth)')
ggsave('plots/randomwalk_mles.png', g, width=ggwidth, height=ggheight)
}

if(exists('simpleGLMM_pvals')){
message("SimpleGLMM: making results plots...")
## g <- ggplot(filter(simpleGLMM_pvals, test=='outlier') , aes(pvalue, )) + geom_histogram(bins=bins) +
##   facet_grid(method~test+version, scales='free')
## ggsave('plots/simpleGLMM_pvalues_outlier.png', g, width=5, height=5)
## g <- ggplot(filter(simpleGLMM_pvals, test=='disp') , aes(pvalue, )) + geom_histogram(bins=bins) +
##   facet_grid(method~test+version, scales='free')
## ggsave('plots/simpleGLMM_pvalues_disp.png', g, width=5, height=5)
g <- filter(simpleGLMM_pvals, grepl('GOF', test) & grepl('osa', x=method)) %>%
  mutate(test=gsub('GOF.', '', x=test), method=gsub('osa.', '', x=method))  %>%
  ggplot(aes(pvalue, fill=test)) + geom_histogram(bins=bins) +
  facet_grid(method~version, scales='free')
ggsave('plots/simpleGLMM_pvalues_GOF_osa.png', g, width=5, height=5)
g <- filter(simpleGLMM_pvals, grepl('GOF', test) & !grepl('osa', x=method)) %>%
  mutate(test=gsub('GOF.', '', x=test), method=gsub('osa.', '', x=method))  %>%
  ggplot(aes(pvalue, fill=test)) + geom_histogram(bins=bins) +
  facet_grid(method~version, scales='free')
ggsave('plots/simpleGLMM_pvalues_GOF_DHARMa.png', g, width=5, height=5)
## g <- pivot_longer(resids, c('osa.cdf', 'osa.gen', 'sim_cond', 'sim_uncond', 'sim_parcond'),
##                   names_to='type', values_to='residual') %>%
##   filter(replicate<=5) %>%
##   ggplot(aes(ytrue, residual, color=version)) +
##   geom_point(alpha=.5) +
##   facet_grid(replicate~type) + scale_x_log10()
## ggsave('plots/simpleGLMM_residuals_examples.png', g, width=8, height=6)
simpleGLMM_mles <- simpleGLMM_mles %>% mutate(abserror=true-mle,
                                        relerror=abserror/true)
g <- ggplot(simpleGLMM_mles, aes(x=version, abserror)) + geom_violin() +
  facet_wrap('par', scales='free') + geom_hline(yintercept=0,
  color='red') + labs(y='absolute error (MLE-truth)')
ggsave('plots/simpleGLMM_mles.png', g, width=ggwidth, height=ggheight)
}

if(exists('linmod_pvals')){
message("linmod: making results plots...")
## g <- ggplot(filter(linmod_pvals, test=='outlier') , aes(pvalue, )) + geom_histogram(bins=bins) +
##   facet_grid(method~test+version, scales='free')
## ggsave('plots/linmod_pvalues_outlier.png', g, width=5, height=5)
## g <- ggplot(filter(linmod_pvals, test=='disp') , aes(pvalue, )) + geom_histogram(bins=bins) +
##   facet_grid(method~test+version, scales='free')
## ggsave('plots/linmod_pvalues_disp.png', g, width=5, height=5)
g <- filter(linmod_pvals, grepl('GOF', test) & grepl('osa', x=method)) %>%
  mutate(test=gsub('GOF.', '', x=test), method=gsub('osa.', '', x=method))  %>%
  ggplot(aes(pvalue, fill=test)) + geom_histogram(bins=bins) +
  facet_grid(method~version, scales='free')
ggsave('plots/linmod_pvalues_GOF_osa.png', g, width=5, height=5)
g <- filter(linmod_pvals, grepl('GOF', test) & !grepl('osa', x=method)) %>%
  mutate(test=gsub('GOF.', '', x=test), method=gsub('osa.', '', x=method))  %>%
  ggplot(aes(pvalue, fill=test)) + geom_histogram(bins=bins) +
  facet_grid(method~version, scales='free')
ggsave('plots/linmod_pvalues_GOF_DHARMa.png', g, width=5, height=5)
g <- pivot_longer(linmod_resids, c('pearsons', 'sim_cond', 'sim_parcond'),
                  names_to='type', values_to='residual') %>%
  filter(replicate<=5) %>%
  ggplot(aes(x, residual, color=version)) +
  geom_point(alpha=.5) +
  facet_grid(replicate~type)
ggsave('plots/linmod_residuals_examples.png', g, width=8, height=6)
linmod_mles <- linmod_mles %>% mutate(abserror=true-mle,
                                        relerror=abserror/true)
g <- filter(linmod_mles, !is.na(mle)) %>% ggplot(aes(x=version, abserror)) + geom_violin() +
  facet_wrap('par', scales='free') + geom_hline(yintercept=0,
  color='red') + labs(y='absolute error (MLE-truth)')
ggsave('plots/linmod_mles.png', g, width=ggwidth, height=ggheight)
}
