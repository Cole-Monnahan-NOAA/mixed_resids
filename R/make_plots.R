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
pvals %>% dplyr::filter(., model == "linmod" & test == "GOF.ks") %>% 
  group_by(., method, misp, version, do.true) %>% 
  summarize(count=length(unique(replicate))) %>%
  print(n=28)

unique(pvals[pvals$model=="randomwalk",]$misp)
pvals %>% dplyr::filter(., model == "randomwalk" & test == "GOF.ks") %>% 
  group_by(., method, misp, version, do.true) %>% 
  summarize(count=length(unique(replicate))) %>%
  print(n=48) # method = fg & version = h1: 317 out of 500 sucessfully ran

unique(pvals[pvals$model=="simpleGLMM",]$misp)
pvals %>% dplyr::filter(., model == "simpleGLMM" & test == "GOF.ks" & misp == "dropRE") %>% 
  group_by(., method, misp, version, do.true) %>% 
  summarize(count=length(unique(replicate))) %>%
  print(n=14)
pvals %>% dplyr::filter(., model == "simpleGLMM" & test == "GOF.ks" & misp == "missnormcov") %>% 
  group_by(., method, misp, version, do.true) %>% 
  summarize(count=length(unique(replicate))) %>%
  print(n=32)
pvals %>% dplyr::filter(., model == "simpleGLMM" & test == "GOF.ks" & misp == "missunifcov") %>% 
  group_by(., method, misp, version, do.true) %>% 
  summarize(count=length(unique(replicate))) %>%
  print(n=32)

unique(pvals[pvals$model=="spatial",]$misp)
pvals %>% dplyr::filter(., model == "spatial" & test == "GOF.ks" & misp == "dropRE") %>% 
  group_by(., method, misp, version, do.true) %>% 
  summarize(count=length(unique(replicate))) %>%
  print(n=40)
pvals %>% dplyr::filter(., model == "spatial" & test == "GOF.ks" & misp == "mispomega") %>% 
  group_by(., method, misp, version, do.true) %>% 
  summarize(count=length(unique(replicate))) %>%
  print(n=28)

## Check for MLE consistency.
nc_id <- stats[!is.na(stats$converge) & (stats$converge==1|stats$maxgrad>0.1),]$id
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
#linear regression
pvals <- lapply(list.files('results', pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows
pvals <- pvals[pvals$do.true == TRUE,]
pvals.dharma <- pvals[pvals$model == 'linmod' & 
                        pvals$type == "sim",]

pvals.osa <- pvals[pvals$model == 'linmod' & 
                     pvals$type == "osa",]
filter(pvals.dharma, test=='GOF.ks' & version == "h0") %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))
filter(pvals.dharma, test=='GOF.ks' & version == "h1" & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))

#linear model - OSA residuals
filter(pvals.osa, test=='GOF.ks' & version == "h0") %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))
filter(pvals.osa, test=='GOF.ks' & version == "h1" & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))

#randomwalk - DHARMa residuals
pvals <- lapply(list.files('results', pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows
pvals.dharma <- pvals[pvals$model == 'randomwalk' & 
                        pvals$type == "sim",]

pvals.osa <- pvals[pvals$model == 'randomwalk' & 
                     pvals$type == "osa",]

filter(pvals.dharma, test=='GOF.ks' & version == "h0") %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))
filter(pvals.dharma, test=='GOF.ks' & version == "h1" & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))

#randomwalk - OSA residuals
filter(pvals.osa, test=='GOF.ks' & version == "h0") %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))
filter(pvals.osa, test=='GOF.ks' & version == "h1" & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))

#simpleGLMM - DHARMa residuals
pvals <- lapply(list.files('results', pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows
pvals <- pvals[pvals$model == 'simpleGLMM',]
pvals.dharma <- pvals[pvals$model == 'simpleGLMM' & 
                        pvals$type == "sim",]
pvals.osa <- pvals[pvals$model == 'simpleGLMM' & 
                     pvals$type == "osa",]

filter(pvals.dharma, test=='GOF.ks' & version == "h0") %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))
filter(pvals.dharma, test=='GOF.ks' & version == "h1" & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))

#simpleGLMM - OSA residuals
filter(pvals.osa, test=='GOF.ks' & version == "h0") %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))
filter(pvals.osa, test=='GOF.ks' & version == "h1" & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))

#calculate power
power.df <- filter(pvals, test=='GOF.ks' & version == "h1" & 
         do.true==FALSE & misp == "missunifcov" )  %>% 
  pivot_wider(., id_cols = !type, names_from = "method", 
              values_from = "pvalue") %>% 
  dplyr::select(unique(pvals$method))
apply(power.df, 2, function(x) sum(x < 0.05, na.rm = TRUE))/
apply(power.df, 2, function(x) sum(x > 0, na.rm = TRUE))

power.df <- filter(pvals, test=='GOF.ks' & version == "h1" & 
         do.true==FALSE & misp == "missnormcov" )  %>% 
  pivot_wider(., id_cols = !type, names_from = "method", 
              values_from = "pvalue") %>% 
  dplyr::select(unique(pvals$method))
apply(power.df, 2, function(x) sum(x < 0.05, na.rm = TRUE))/
apply(power.df, 2, function(x) sum(x > 0, na.rm = TRUE))

#calculate Type I error
tIerr.df <- filter(pvals, test=='GOF.ks' & version == "h0" & 
                     do.true==FALSE & misp == "missunifcov" )  %>% 
  pivot_wider(., id_cols = !type, names_from = "method", 
              values_from = "pvalue") %>% 
  dplyr::select(unique(pvals$method))
apply(tIerr.df, 2, function(x) sum(x < 0.05, na.rm = TRUE))/
apply(tIerr.df, 2, function(x) sum(x > 0, na.rm = TRUE))

tIerr.df <- filter(pvals, test=='GOF.ks' & version == "h0" & 
                     do.true==FALSE & misp == "missnormcov" )  %>% 
  pivot_wider(., id_cols = !type, names_from = "method", 
              values_from = "pvalue") %>% 
  dplyr::select(unique(pvals$method))
apply(tIerr.df, 2, function(x) sum(x < 0.05, na.rm = TRUE))/
apply(tIerr.df, 2, function(x) sum(x > 0, na.rm = TRUE))

#spatial - DHARMa residuals
pvals <- lapply(list.files('results', pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% bind_rows
# pvals <- pvals[pvals$method != "cond",]
pvals.dharma <- pvals[pvals$model == 'spatial'  & pvals$type == "sim",]
pvals.osa <- pvals[pvals$model == 'spatial' & pvals$type == "osa",]

filter(pvals.dharma, test=='GOF.ks' & version == "h0") %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))
filter(pvals.dharma, test=='GOF.ks' & version == "h1" & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))

#spatial - OSA residuals
filter(pvals.osa, test=='GOF.ks' & version == "h0") %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))
filter(pvals.osa, test=='GOF.ks' & version == "h1" & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))

#SAC
filter(pvals.dharma, test=='SAC' & version == "h0") %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))
filter(pvals.dharma, test=='SAC' & version == "h1" & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))

#spatial - OSA residuals
filter(pvals.osa, test=='SAC' & version == "h0") %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))
filter(pvals.osa, test=='SAC' & version == "h1" & do.true==FALSE) %>%
  ggplot(aes(pvalue, fill=do.true, color=do.true)) +
  facet_grid(misp~method, scales = "free_y") + 
  geom_histogram(position='identity', alpha=.5) +
  scale_fill_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value') + 
  scale_color_viridis_d(labels = c('estimated','theoretical'), name = 'GOF p-value')  +
  theme(axis.text=element_text(size=6))


#simpleGLMM and spatial power and Type I error
name.map <- c(
  simpleGLMM_dropRE = "GLMM Tweedie",
  simpleGLMM_missnormcov = "Mixed Model Norm X",
  simpleGLMM_missunifcov = "Mixed Model Unif X",
  spatial_dropRE = "GLMM Spatial Poisson",
  spatial_mispomega = "Mixed Model Spatial"
)
method.level <- c(
  uncond_nrot = "Unconditional",
  uncond = "Unconditional",
  pears = "Conditional",
  osg = "Conditional",
  mcmc = "Unconditional",
  gen = "Conditional",
  fg = "Unconditional",
  cond_nrot = "Conditional",
  cond = "Conditional",
  cdf = "Conditional"
)
method.name <- c(
  uncond_nrot = "Not Rotated ecdf",
  uncond = "Rotated ecdf",
  pears = "Pearson",
  osg = "One-step Gaussian",
  mcmc = "MCMC",
  gen = "One-step Generic",
  fg = "Full Gaussian",
  cond_nrot = "Not Rotated ecdf",
  cond = "Rotated ecdf",
  cdf = "One-step cdf"
)
#Type I error
pvals.err <- lapply(list.files('results', pattern='_pvals.RDS',
                                 full.names=TRUE), readRDS) %>% 
  bind_rows %>%
  dplyr::filter((model == "spatial" | model == "simpleGLMM") & 
                  misp != "deltagamma" & test=='GOF.ks' & version == "h0" & 
                  do.true==FALSE & method != "re_uncond_nrot" &
                  method != "re_uncond" &
                  method != "re_mcmc" &
                  method != "re_fg") 
pvals.err$mod <- paste0(pvals.err$model, "_", pvals.err$misp)
pvals.err$mod.name <- unname(name.map[pvals.err$mod])
pvals.err$method.name <- unname(method.name[pvals.err$method])
pvals.power$method.level <- unname(method.level[pvals.err$method])

pvals.err %>% group_by(mod.name, method.name, method.level) %>%
  summarize(typeIerror = sum(pvalue <= 0.05)/sum(pvalue >= 0)) %>%
  ggplot(., aes(x = typeIerror, y = method.name)) + geom_point() + 
  facet_grid(method.level~mod.name, scales = "free_y") + 
  geom_vline(xintercept = 0.05)

pvals.err <- lapply(list.files('results', pattern='_pvals.RDS',
                           full.names=TRUE), readRDS) %>% 
  bind_rows %>%
  dplyr::filter((model == "spatial" | model == "simpleGLMM") & 
                  misp != "deltagamma" & test=='GOF.ks' & version == "h0" & 
                  do.true==FALSE)
#Power
pvals.power <- lapply(list.files('results', pattern='_pvals.RDS',
                                 full.names=TRUE), readRDS) %>% 
  bind_rows %>%
  dplyr::filter((model == "spatial" | model == "simpleGLMM") & 
                  misp != "deltagamma" & test=='GOF.ks' & version == "h1" & 
                  do.true==FALSE & method != "re_uncond_nrot" &
                  method != "re_uncond" &
                  method != "re_mcmc" &
                  method != "re_fg") 
pvals.power$mod <- paste0(pvals.power$model, "_", pvals.power$misp)
pvals.power$mod.name <- unname(name.map[pvals.power$mod])
pvals.power$method.name <- unname(method.name[pvals.power$method])
pvals.power$method.level <- unname(method.level[pvals.power$method])

pvals.power %>% group_by(mod.name, method.name, method.level) %>%
  summarize(power = sum(pvalue <= 0.05)/sum(pvalue >= 0)) %>%
  ggplot(., aes(x = power, y = method.name)) + geom_point() + 
  facet_grid(method.level~mod.name, scales = "free_y") + 
  geom_vline(xintercept = 0.95)

pvals.power %>% group_by(mod.name, method.name, method.level) %>%
  summarize(power = sum(pvalue <= 0.05)/sum(pvalue >= 0)) %>%
  ggplot(., aes(x = mod.name, y = power, color = method.name, group = method.name)) + geom_line() + 
  # geom_point(stat='summary', fun=sum) +
  # stat_summary(fun=sum, geom="line") +
  facet_wrap(~method.level, scales = "free_y") + 
  geom_vline(xintercept = 0.95)

