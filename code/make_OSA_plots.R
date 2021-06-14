
## Pull out pieces of OSA calcs for FullGaussian and recreate
## minimally to explore. Using the function extractOSA which is a
## copy of oneStepPredict modified to return the information we
## need to visualize.
rm(list=ls())
source('code/functions_simpleGLMM.R')
source('code/extractOSA.R')
## Build object and optimize
TMB::compile('models/simpleGLMM.cpp')
dyn.load(dynlib('models/simpleGLMM'))
out <- simulate.simpleGLMM(seed=123, ngroups=15, nobs=6)
obj <- MakeADFun(out$Data, out$Par, random = 'u', DLL = 'simpleGLMM')
opt <- nlminb(obj$par, obj$fn, obj$gr)

### ------------------------------------------------------------
### First try visualizing the Full Guassian method
FG <- extractOSA(obj, observation.name='y', method='fullGaussian')
FG$Sigma[1:3, 1:3] %>% cov2cor
## Hack way to generate random draws and plot the real data vs
## them.
library(mvtnorm)
nsim <- 1000
x <- rmvnorm(n=nsim, mean=FG$mode, sigma=FG$Sigma)
x <- rbind(x, FG$obs, FG$mode)
## The first group
## ind <- c(1,2,3, 11,12,13, 21,22,23)
ind <- which(out$Data$group==1)
ind <- 1:6
labs <- paste0('y[',ind,']\n', round(FG$pred[ind,1],2))
## png('plots/glmm_fullGaussian_mvn.png', width=7, height=7,
##     units='in', res=500)
pairs(x[,ind], upper.panel=NULL, labels=labs,
      col=c(rep(rgb(0,0,0,.1),nsim), rgb(1,0,0), rgb(0,1,0)),
      pch=16, gap=0)
## dev.off()
## Now look at Cholesky rotated version
xseq <- seq(-5,5, len=1000); yseq <- dnorm(xseq)
png('plots/glmm_fullGaussian_znorm.png', width=9, height=5,
    units='in', res=500)
par(mfrow=c(1,2), mgp=c(1.5,.5,0), tck=-.02, mar=c(3,3,.5,.5))
plot(xseq, yseq, type='l', xlab='Residual', ylab='Density')
mycol <- rep(c(1,2,3), each=10)
text(FG$pred[,1], y=dnorm(FG$pred[,1]), labels=out$Data$group,#paste0('y[', 1:30,']'),
     cex=.7, col=mycol)
boxplot(FG$pred[,1]~out$Data$group, xlab='Group', ylab='Residual')
dev.off()


qqnorm(FG$pred$residual); abline(0,1)
resids

### ------------------------------------------------------------
## oneStepGaussian (OSG)
OSG <- extractOSA(obj, 'y', 'keep', method='oneStepGaussian')
png('plots/glmm_oneStepGauss_marginals.png', width=7, height=5,
    units='in', res=400)
par(mfrow=c(2,1), mgp=c(1.5, .5, 0), tck=-.02, mar=c(3,3,1,1))
ind <-  c(1,11,21)#c(1,2,3, 11,12,13, 21,22,23)
plot(0,0, type='n', xlab='Y',
     ylab='Probability', xlim=c(2,8), ylim=c(0,1e-12))
mycols <- 1:3
for(k in ind){
  g <- out$Data$group[k]+1
  pp <- OSG$LLpoints[OSG$LLpoints$k==k,]
  ll <- OSG$LLcurve[OSG$LLcurve$k==k,]
  lines(ll$x, exp(-ll$y), col=mycols[g])
  ## abline(v=OSG$pred$obs[k], col=mycols[g])
  rug(x=OSG$pred$obs[k], col=mycols[g])
  points(pp$x, exp(-pp$y), cex=1.5, col=mycols[g], pch=16)
}
ind <-  1:23# c(1,2,3, 11,12,13, 21,22,23)
plot(0,0, type='n', xlab='Y',
     ylab='Scaled Probability', xlim=c(2,8), ylim=c(0,1))
mycols <- 1:3
for(k in ind){
  g <- out$Data$group[k]+1
  pp <- OSG$LLpoints[OSG$LLpoints$k==k,]
  ll <- OSG$LLcurve[OSG$LLcurve$k==k,]
  lines(ll$x, exp(-ll$y+pp$y), col=mycols[g])
  ## abline(v=OSG$pred$obs[k], col=mycols[g])
  rug(x=OSG$pred$obs[k], col=mycols[g])
  points(pp$x, exp(-pp$y+pp$y), cex=1.5, col=mycols[g], pch=16)
}
dev.off()


### ------------------------------------------------------------
### Compare the different residuals
resids.all <- data.frame(id=1:length(Dat$y), obs=Dat$y, group=Dat$group,
                         FG=FG$pred$residual,
                         OSG=OSG$pred$residual)
resids.long <- resids.all %>%
  pivot_longer(-c(obs, group, id), names_to='method',
               values_to='residual')
g <- ggplot(resids.long, aes(id, residual, color=method)) +
  geom_line()
ggsave('plots/glmm_compare_residuals.png', g, width=7, height=5)
