


message("Spatial: making results plots...")
g <- ggplot(filter(spatial_pvals, test=='outlier') , aes(pvalue, )) + geom_histogram() +
  facet_grid(method~test+version, scales='free')
ggsave('plots/spatial_pvalues_outlier.png', g, width=5, height=5)
g <- ggplot(filter(spatial_pvals, test=='disp') , aes(pvalue, )) + geom_histogram() +
  facet_grid(method~test+version, scales='free')
ggsave('plots/spatial_pvalues_disp.png', g, width=5, height=5)
g <- ggplot(filter(spatial_pvals, test=='sac') , aes(pvalue, )) + geom_histogram() +
  facet_grid(method~test+version, scales='free')
ggsave('plots/spatial_pvalues_sac.png', g, width=5, height=5)
g <- filter(spatial_pvals, test=='GOF' & grepl('osa', x=method)) %>%
  ggplot(aes(pvalue)) + geom_histogram() +
  facet_grid(method~test+version, scales='free')
ggsave('plots/spatial_pvalues_GOF_osa.png', g, width=5, height=5)
g <- filter(spatial_pvals, test=='GOF' & !grepl('osa', x=method)) %>%
  ggplot(aes(pvalue)) + geom_histogram() +
  facet_grid(method~test+version, scales='free')
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



message("Randomwalk: making results plots...")
g <- ggplot(filter(randomwalk_pvals, test=='outlier') , aes(pvalue, )) + geom_histogram() +
  facet_grid(method~test+version, scales='free')
ggsave('plots/randomwalk_pvalues_outlier.png', g, width=5, height=5)
g <- ggplot(filter(randomwalk_pvals, test=='disp') , aes(pvalue, )) + geom_histogram() +
  facet_grid(method~test+version, scales='free')
ggsave('plots/randomwalk_pvalues_disp.png', g, width=5, height=5)
g <- filter(randomwalk_pvals, test=='GOF' & grepl('osa', x=method)) %>%
  ggplot(aes(pvalue)) + geom_histogram() +
  facet_grid(method~test+version, scales='free')
ggsave('plots/randomwalk_pvalues_GOF_osa.png', g, width=5, height=5)
g <- filter(randomwalk_pvals, test=='GOF' & !grepl('osa', x=method)) %>%
  ggplot(aes(pvalue)) + geom_histogram() +
  facet_grid(method~test+version, scales='free')
ggsave('plots/randomwalk_pvalues_GOF_DHARMa.png', g, width=5, height=5)
## g <- pivot_longer(resids, c('osa.cdf', 'osa.gen', 'sim_cond', 'sim_uncond', 'sim_parcond'),
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
