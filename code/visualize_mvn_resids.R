library(mvtnorm)
library(plotly)
library(goftest)
set.seed(1)
C <- diag(2)
C[1,2] <- C[2,1] <- 0.7
obs <- as.numeric(rmvnorm(1, sigma=C)) 

#Plot MVN
x.seq <- y.seq <- seq(-4,4,.1)
dmvn <- mpdf <- mcdf <- matrix(0,length(x.seq),length(y.seq))
for(i in 1:length(x.seq)){
  for(j in 1:length(y.seq)){
    dmvn[i,j] <- dmvnorm(c(x.seq[i],y.seq[j]), sigma = C)
    mpdf[i,j] <- pmvnorm(lower = -Inf, upper = c(x.seq[i], y.seq[j]), mean = c(0,0), corr = C)
    mcdf[i,j] <- qmvnorm(mpdf[i,j], sigma = C)$quantile
  }
}
plot_ly(x = ~x.seq, y = ~y.seq, z = ~dmvn, type = 'surface') 
plot_ly(x = ~x.seq, y = ~y.seq, z = ~mpdf, type = 'surface') 
plot_ly(x = ~x.seq, y = ~y.seq, z = ~mcdf, type = 'surface') 

#theoretical cdf:
qmvnorm(pmvnorm(lower = -Inf, upper = obs, corr = C), sigma = C)$quantile
#marginal cdfs:
qnorm(pnorm(obs)) #same as original obs

#marginal cdfs after rotation:
## L * L^T = C
L <- t(chol(C))
## OSA residuals
r1 <- solve(L, obs)
obs;r1

#Generate n-dimensional mvn data
set.seed(1)
C <- exp(-as.matrix(dist(seq(0,50,by=.1))))
obs <- as.numeric(rmvnorm(1, sigma=C)) 
hist(obs, breaks=50)
n.obs <- length(obs)

## L * L^T = C
L <- t(chol(C))

## OSA residuals
r1 <- solve(L, obs)
qqnorm(r1, main="OSA")
abline(0,1)
ad.test(r1, 'pnorm')$p.value

n.sim <- 1000
y.sim <- t(rmvnorm(n.sim,sigma=C))
sim.rot <- solve(L, y.sim)

#y.sim is correlated within each simulation
pairs(t(y.sim[1:2,]), labels = c('obs1', 'obs2'))
#rotated simulation is uncorrelated
pairs(t(sim.rot[1:2,]), labels = c('obs1', 'obs2'))

#but between simulations are iid for both
pairs(y.sim[,1:2], labels = c('sim1', 'sim2'))
pairs(sim.rot[,1:2], labels = c('sim1', 'sim2'))

#Correlation within simulations, correlated obs
res.cc <- DHARMa::createDHARMa(y.sim,obs)$scaledResiduals
gap::qqunif(res.cc,logscale = FALSE);abline(0,1)
ks.test(res.cc, 'punif')$p.value
plot(1:n.obs, res.cc)


#Correlation within simulations, rotated obs
res.ci <- DHARMa::createDHARMa(y.sim,r1)$scaledResiduals 
gap::qqunif(res.ci, logscale = FALSE);abline(0,1)
ks.test(res.ci, 'punif')$p.value
plot(1:n.obs, res.ci)

#iid within simulations, correlated obs
res.ic <- DHARMa::createDHARMa(sim.rot,obs)$scaledResiduals 
gap::qqunif(res.ic, logscale = FALSE);abline(0,1)
ks.test(res.ic, 'punif')$p.value
plot(1:n.obs, res.ic)

#iid within simulations, rotated obs
res.ii <- DHARMa::createDHARMa(sim.rot,r1)$scaledResiduals 
gap::qqunif(res.ii, logscale = FALSE);abline(0,1)
ks.test(res.ii, 'punif')$p.value
plot(1:n.obs, res.ii)
