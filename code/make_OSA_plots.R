
## Pull out pieces of OSA calcs for FullGaussian and recreate
## minimally to explore. Using the function extractOSA which is a
## copy of oneStepPredict modified to return the information we
## need to visualize.
rm(list=ls())
source('code/run_simpleGLMM.R')
source('code/extractOSA.R')
## Rebuild obj just in case
obj <- MakeADFun(Data, Par, random = 'u', DLL = 'simpleGLMM')
opt <- nlminb(obj$par, obj$fn, obj$gr)

### ------------------------------------------------------------
### First try visualizing the Full Guassian method
FG <- extractOSA(obj, observation.name='y', method='fullGaussian')
## Hack way to generate random draws and plot the real data vs
## them.
library(mvtnorm)
nsim <- 1000
x <- rmvnorm(n=nsim, mean=FG$mode, sigma=FG$Sigma)
x <- rbind(x, FG$obs, FG$mode)
ind <- c(1,2,3, 11,12,13, 21,22,23)
labs <- paste0('y[',ind,']\n', round(FG$pred[ind,1],2))
png('plots/glmm_fullGaussian_mvn.png', width=7, height=7,
    units='in', res=500)
pairs(x[,ind], upper.panel=NULL, labels=labs,
      col=c(rep(rgb(0,0,0,.1),nsim), rgb(1,0,0), rgb(0,1,0)),
      pch=16, gap=0)
dev.off()
## Now look at Cholesky rotated version
xseq <- seq(-5,5, len=1000); yseq <- dnorm(xseq)
png('plots/glmm_fullGaussian_znorm.png', width=9, height=5,
    units='in', res=500)
par(mfrow=c(1,2), mgp=c(1.5,.5,0), tck=-.02, mar=c(3,3,.5,.5))
plot(xseq, yseq, type='l', xlab='Residual', ylab='Density')
mycol <- rep(c(1,2,3), each=10)
text(FG$pred[,1], y=dnorm(FG$pred[,1]), labels=paste0('y[', 1:30,']'),
     cex=.7, col=mycol)
boxplot(FG$pred[,1]~Dat[,2], xlab='Group', ylab='Residual')
dev.off()



### ------------------------------------------------------------
## oneStepGaussian (OSG)
OSG <- extractOSA(obj, 'y', 'keep', method='oneStepGaussian')
png('plots/glmm_oneStepGauss_marginals.png', width=7, height=7,
    units='in', res=400)
par(mfrow=c(3,3), mgp=c(1.5, .5, 0), tck=-.02, mar=c(3,3,1,1))
ind <- c(1,2,3, 11,12,13, 21,22,23)
for(k in ind){
  ll <- OSG$LLcurve[OSG$LLcurve$k==k,]
  pp <- OSG$LLpoints[OSG$LLpoints$k==k,]
  plot(ll$x, ll$y , type='l', xlab='Residual',
       ylab='NLL', main=paste('Obs=',k))
  abline(v=OSG$pred$obs[k])
  points(pp$x, pp$y, cex=1.5, col=2, pch=16)
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
