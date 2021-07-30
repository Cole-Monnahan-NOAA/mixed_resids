library(mvtnorm)
library(plotly)
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