


## Clean up the old runs
unlink('results/simpleGLMM_resids', TRUE)
unlink('results/simpleGLMM_pvals', TRUE)
unlink('results/simpleGLMM_mles', TRUE)

message("Preparing workspace to run ", Nreps, " iterations in parallel...")
TMB::compile("models/simpleGLMM.cpp") # modified for simulation
sfInit( parallel=TRUE, cpus=cpus )
## sfExport('run.simpleGLMM.iter', 'sim.omega', 'cMatern', 'sim.data',
##          'rmvnorm_prec', 'add_aic')
sfExportAll()

message("Starting parallel runs...")
results <- sfLapply(1:Nreps, function(ii) run.simpleGLMM.iter(ii))

## ## Read results back in from file
## fs <- list.files('results/simpleGLMM_pvals/', full.names=TRUE)
## ## Sometimes they fail to run for some unkonwn reason so try
## ## rerunning those ones once
## if(length(fs)<Nreps){
##   message("Rerunning some failed runs...")
##   bad <- which.failed(Nreps)
##   results <- sfLapply(bad, function(ii) run.simpleGLMM.iter(ii))
##   fs <- list.files('results/simpleGLMM_pvals/', full.names=TRUE)
## }
## bad <- which.failed(Nreps)
## if(length(bad)>0) warning(length(bad), " runs failed")

message("SimpleGLMM: processing and saving final results...")
## Read results back in from file
fs <- list.files('results/simpleGLMM_pvals/', full.names=TRUE)
pvals <- lapply(fs, readRDS) %>% bind_rows %>% filter(!is.na(pvalue))
saveRDS(pvals, file='results/simpleGLMM_pvals.RDS')
## Read in residuals
fs <- list.files('results/simpleGLMM_resids/', full.names=TRUE)
resids <- lapply(fs, readRDS) %>% bind_rows
saveRDS(resids, file='results/simpleGLMM_resids.RDS')
fs <- list.files('results/simpleGLMM_mles/', full.names=TRUE)
mles <- lapply(fs, readRDS) %>% bind_rows
saveRDS(mles, file='results/simpleGLMM_mles.RDS')



stop("Temp code beyond here")
temp <- resids %>%   filter(version=='m0' & replicate==1) %>%
  select(c(y, ypred, version, osa.fg))

ggplot(temp, aes(y, osa.fg, color=version)) + geom_point()
ggplot(temp, aes(sample=osa.fg, color=version)) + geom_qq() + facet_wrap('version')

ad.test(temp$osa.fg, "pnorm", mean=mean(temp$osa.fg),
        sd=sd(temp$osa.fg), estimated=TRUE)
ad.test(temp$osa.fg, "pnorm")


pvals %>% filter(version=='m0' & method=='osa.fg' & grepl('GOF', x=test)) %>%
  ggplot(aes(pvalue, fill=test)) +
  geom_histogram(position='identity', alpha=.5)


## Quick exploration of cholesky rotation
set.seed(2355)
N <- 5
Sigma <- rWishart(1, df=N, Sigma=diag(N))[,,1]
cov2cor(Sigma)
(ses <- sqrt(diag(Sigma)))

Y <- rmvnorm(500, sigma=Sigma)
pairs(Y)

L <- t(chol(Sigma))
L %*% t(L)-Sigma # check


z <- rnorm(N)
y <- L %*% z
pairs(rbind(Y,t(y)), col=c(rep(1,500), 2), pch=16)

Z <- solve(L) %*% t(Y) %>% t
pairs(Z)
