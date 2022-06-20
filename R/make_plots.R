### Startup
source("R/startup.R")

pvals <- lapply(list.files('results', pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows

mles <- lapply(list.files('results', pattern='_mles.RDS',
                          full.names=TRUE), readRDS) %>% bind_rows

resids <- lapply(list.files('results', pattern='_resids.RDS',
                            full.names=TRUE), readRDS) %>% bind_rows

stats <- lapply(list.files('results', pattern='_stats.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows

## check on if runs worked out
check_runs <- group_by(pvals, model, misp, version, do.true) %>% summarize(count=length(unique(replicate)))
check_runs[1:8,]
check_runs[9:12,]
check_runs[13:16,]

## Check for MLE consistency.
nc_id <- stats[stats$converge==1|stats$maxgrad>0.1,]$id
#only use spatial misp = dropRE
mles <- mles[!(mles$misp == "mispomega"),]
mles <- mles[!(mles$model == "spatial" & mles$misp == "overdispersion"),]
mles <- mles[!(mles$model == "spatial" & mles$misp == "misscov"),]
g <- filter(mles, h==0 & !(id %in% nc_id) & do.true == FALSE) %>%
  ggplot(., aes(x=par, y=bias)) + geom_violin() +
  facet_wrap(~model, nrow=1, scales = 'free_x')
g

#filter before plotting
pvals <- pvals[!(pvals$model=='linmod'&pvals$method=='uncond'),]
#only use spatial misp == dropRE for main plot and drop method == 'gen'
pvals <- pvals[pvals$method != "gen",]
pvals <- pvals[pvals$method != "pears.adj",]
pvals <- pvals[!(pvals$misp == "mispomega"),]
pvals <- pvals[!(pvals$model == "spatial" & pvals$misp == "overdispersion"),]
pvals <- pvals[!(pvals$model == "spatial" & pvals$misp == "misscov"),]

goodreps <- filter(pvals, model=='spatial' & version=='h1') %>%
  summarize(replicate=unique(replicate)) %>% pull(replicate)
which(! (1:500 %in% goodreps))
library(viridis)
pvals$method <- factor(pvals$method, 
                       levels = c('fg', 'osg','cdf', 'mcmc',  're_mcmc', 
                                  're_obs_mcmc','cond','uncond', 
                                  're_uncond', 
                                  'uncond_nrot', 're_uncond_nrot',
                                  'pears'))
## Effect of do.true for true model
png(filename='plots/h0_TF.png', width=2600, height=1500,res=300)
filter(pvals, version=='h0' &  test=='GOF.ks') %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(model~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))# + 
#scale_x_continuous(labels = scales::number_format(accuracy = 0.1))
dev.off()

png(filename='plots/h1_TF.png', width=2600, height=1500,res=300)
filter(pvals, version=='h1' &  test=='GOF.ks'& do.true == FALSE) %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(model~method, scales = "free_y") + geom_histogram(position='identity', alpha=.5)+
  scale_fill_viridis_d(labels = c('estimated'), name = 'GOF p-value')+
  scale_color_viridis_d(labels = c('estimated'), name = 'GOF p-value') +
  theme(axis.text=element_text(size=6))
dev.off()

## Model version for KS test
png(filename='plots/h0h1_F.png', width=2600, height=1500,res=300)
pvals %>% filter(test== 'GOF.ks' & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=version)) +
  facet_grid(model~method, scales = "free_y") + geom_histogram(position='identity', alpha=.5) +
  theme(axis.text=element_text(size=6))
dev.off()

## KS vs AD for true model - delete create new plots for each test
#pvals %>% filter(version=='h0'& test!= 'outlier' & do.true==TRUE) %>%
#  ggplot(aes(pvalue, fill=test)) +
#  facet_grid(model~method) + geom_histogram(position='identity', alpha=.5)


pvals <- lapply(list.files('results', pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows
sp.pvals <- pvals[pvals$model=='spatial',]
sp.stats <- stats[stats$model=='spatial',]
nc.reps <- filter(sp.stats, !is.na(converge) & converge == 1)$replicate
# filter out models that didn't converge
sp.pvals <- sp.pvals[!(sp.pvals$replicate %in% nc.reps),]
sp.pvals <- sp.pvals[sp.pvals$method != 'pears.adj',]

# sp.pvals %>% filter(do.true==TRUE) %>%
#   ggplot(aes(pvalue, fill=version)) +
#   facet_grid(test~method) + geom_histogram(position='identity', alpha=.5)
# 
# sp.pvals %>% filter(do.true==FALSE) %>%
#   ggplot(aes(pvalue, fill=version)) +
#   facet_grid(test~method, scales = "free_y") + geom_histogram(position='identity', alpha=.5)

sp.pvals %>% filter(do.true==FALSE & test == 'GOF.ks') %>%
  ggplot(aes(pvalue, fill=version)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d()

# sp.pvals %>% filter(do.true==FALSE & test == 'outliers') %>%
#   ggplot(aes(pvalue, fill=version)) +
#   facet_grid(misp~method, scales = "free_y") + 
#   geom_histogram(position='identity', alpha=.5) +
#   scale_fill_viridis_d()

sp.pvals %>% filter(do.true==FALSE & test == 'SAC') %>%
  ggplot(aes(pvalue, fill=version)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d()

sp.pvals %>% filter(do.true==FALSE & test == 'disp') %>%
  ggplot(aes(pvalue, fill=version)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d()



#MLE plots
## mles %>% filter(model == 'spatial' & type == 'mle' &
##                   do.true == FALSE & h==0) %>%
##   ggplot(., aes(x=par, y=value)) + geom_violin() +
##   geom_hline(yintercept = -2) + geom_hline(yintercept = 1)

## Check for MLE consistency.
nc_id <- stats[stats$converge==1|stats$maxgrad>0.1,]$id

g <- filter(mles, h==0 & !(id %in% nc_id) & do.true == FALSE) %>%
  ggplot(., aes(x=par, y=bias)) + geom_violin() +
  facet_wrap(~model, nrow=1, scales = 'free_x')
g

saveRDS(mles, 'mles.RDS')
#! Not modified yet
# ## This script runs the same examples above but with varying
# ## sample sizes and tracks runtime and pvalues to see the effect
# ## of this dimension.
# (cpus <- parallel::detectCores()-2)
# Nreps <- 1000
# sfInit( parallel=TRUE, cpus=cpus )
# source('code/run_sample_sizes.R')

### Load results into workspace
source('code/load_results.R')
### Make quick plots of results
source("code/make_plots.R")




## ## Quick test of nonparameteric tests on real normal samples
## library(goftest)
## Nsim <- 15000
## Nreps <- c(10,100,1000)
## mu <- 0; sigma <- 1
## results <- list(); k <- 1
## for(j in Nreps){
##   for(i in 1:Nsim){
##     x <- rnorm(j, mu, sigma)
##     x <- rt(j, 20)
##     p.ks <- ks.test(x,'pnorm')$p.value
##     p.ad <- ad.test(x, 'pnorm', estimated = TRUE)$p.value
##     results[[k]] <- data.frame(Nreps=j, ks=p.ks, ad=p.ad)
##     k <- k+1
##   }
## }
## results.long <- bind_rows(results) %>%
##   pivot_longer(-c(Nreps), names_to='test', values_to='pvalue')
## g <- ggplot(results.long, aes(pvalue)) + geom_histogram() +
##   facet_grid(Nreps~test, scales='free_y')
## ggsave('plots/ks_vs_ad.png', g, width=7, height=5)

#Single instance plots =====================================
#randomwalk
pvals <- lapply(list.files('results', pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows
pvals <- pvals[pvals$model == 'randomwalk' & pvals$misp == 'mu0',]
pvals <- pvals[pvals$method == 'cond' | pvals$method == 'uncond',]
filter(pvals, version=='h0' &  test=='GOF.ks') %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_wrap(~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))


#simpleGLMM
pvals <- lapply(list.files('results', pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows
pvals <- pvals[pvals$model == 'simpleGLMM' & pvals$misp == 'deltagamma',]
pvals <- pvals[pvals$method == 'cond' | pvals$method == 'uncond',]
filter(pvals, version=='h0' &  test=='GOF.ks') %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_wrap(~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))
