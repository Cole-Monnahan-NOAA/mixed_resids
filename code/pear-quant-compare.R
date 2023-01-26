## IID =======================================================================
set.seed(123)
# Approx. normal gamma
y1 <- rgamma(1000, shape = 100, scale = .1)# mean=10; var=1
hist(y1)
# Skewed gamma
y2 <- rgamma(1000, shape = .5, scale = 3)# mean=1.5; var = 4.5
hist(y2)

output.iid <- function(y){
  hist(y)
  mu <- mean(y)
  sig2. <- var(y)
  pear.res <- (y-mu)/sqrt(sig2.)
  Fx <- pgamma(y, shape = mu^2/sig2., scale = sig2./mu)
  quant.res <- qnorm(Fx)
  pear.gof <- ks.test(pear.res, "pnorm")$p.value
  quant.gof <- ks.test(quant.res, "pnorm")$p.value
  (qqnorm(pear.res,
          main = "Pearson Normal Q-Q Plot",
          xlab = paste("Ks Test", round(pear.gof, 3))));abline(0,1, col = "red", lwd=1.5)
  (qqnorm(quant.res,
          main = "Quantile Normal Q-Q Plot",
          xlab = paste("Ks Test", round(quant.gof, 3))));abline(0,1, col = "red", lwd=1.5)
  
  
  return(invisible(y))
  
}
# iid
par(mfrow = c(2,3))
output.iid(y1)
output.iid(y2)


## Random grouping term ==========================================================
set.seed(123)
eta.mat <- as.matrix(exp(rnorm(5, 2.3, .5)))
# Approx. normal gamma
sig2 <- 1
y3 <- apply(eta.mat, 1, function(x) rgamma(200, x^2/sig2, scale = sig2/x ))
hist(y3)
# Skewed gamma
eta.mat <- as.matrix(exp(rnorm(5, 0.4, .5)))
sig2 <- 4.5
y4 <- apply(eta.mat, 1, function(x) rgamma(200, x^2/sig2, scale = sig2/x ))
hist(y4)

output.groupRE <- function(y){
  boxplot(y, horizontal=TRUE)
  #re-scaling by group mean and sd
  Fx <- apply(y, 2, function(x) pgamma(x, mean(x)^2/var(x), scale = var(x)/mean(x)))
  pear.res <- apply(y, 2, function(x) (x - mean(x))/sd(x))
  quant.res <- qnorm(Fx)
  pear.gof <- ks.test(pear.res, "pnorm")$p.value
  quant.gof <- ks.test(quant.res, "pnorm")$p.value
  (qqnorm(pear.res,
          main = "Pearson Normal Q-Q Plot",
          xlab = paste("Ks Test", round(pear.gof, 3))));abline(0,1, col = "red", lwd=1.5)
  (qqnorm(quant.res,
          main = "Quantile Normal Q-Q Plot",
          xlab = paste("Ks Test", round(quant.gof, 3))));abline(0,1, col = "red", lwd=1.5)
  
  
  return(invisible(y))
  
}
par(mfrow = c(2,3))
#approx normal gamma with group RE term
output.iid(as.vector(y3))
#apply rotation - fixes problem
y3rot <- output.groupRE(y3)

par(mfrow = c(2,3))
#skewed gamma with group RE term
output.iid(as.vector(y4))
output.groupRE(y4)

dev.off()

## Banded Correlation Structure ====================================================
library(mvtnorm)
library(moments)

#Copula example
set.seed(1)
C <- exp(-as.matrix(dist(seq(0,1,by=.1))))
L <- t(chol(C))
n.sim <- 500
eta <- t(rmvnorm(n.sim,sigma=C))

#Use Copula method
y5 <- qgamma(pnorm(eta), shape = 100, scale = .1)
y6 <- qgamma(pnorm(eta), shape = .5, scale = 3)
#Standard Normal distribution: skewness = 0, kurtosis = 3
moments::skewness(as.vector(y5)); moments::kurtosis(as.vector(y5))
moments::skewness(as.vector(y6)); moments::kurtosis(as.vector(y6))


#pgamma removes skewness and lowers kurtosis
r5 <- qnorm(pgamma(y5, shape = 100, scale = .1))
r6 <- qnorm(pgamma(y6, shape = .5, scale = 3))
moments::skewness(as.vector(r5)); moments::kurtosis(as.vector(r5))
moments::skewness(as.vector(r6)); moments::kurtosis(as.vector(r6))

qqnorm(r5);abline(0,1); ks.test(r5, 'pnorm')
qqnorm(r6);abline(0,1); ks.test(r6, 'pnorm')

#rotation removes correlation
r5.rot <- solve(L, r5)
r6.rot <- solve(L, r6)
qqnorm(r5.rot); abline(0,1); ks.test(r5.rot, "pnorm")
qqnorm(r6.rot); abline(0,1); ks.test(r6.rot, "pnorm")

# Correlation in mean process
# Normal case
set.seed(1)
C <- exp(-as.matrix(dist(seq(0,10,length.out = 500))))
sig <- 1.2
S <- diag(rep(sig, ncol(C))) %*% C %*%  diag(rep(sig, ncol(C)))
eta <- as.numeric(rmvnorm(1, rep(0, ncol(C)), sigma=S)) 
L <- t(chol(S))

y7 <- eta
mode <- y7 - mean(eta)
r7.nrot <- mode/sd(eta)# this is essentially the Pearson residual (y7-mean)/sd
r7.rot <- t(solve(L, mode))
gof.nrot <- ks.test(r7.nrot, "pnorm")$p.value
gof.rot <- ks.test(r7.rot, "pnorm")$p.value

par(mfrow = c(2,3))
hist(y7)
(qqnorm(r7.nrot,
        main = "Unrotated Normal Q-Q Plot",
        xlab = paste("Ks Test", round(gof.nrot, 3))));abline(0,1, col = "red", lwd=1.5)
(qqnorm(r7.rot,
        main = "Rotated Normal Q-Q Plot",
        xlab = paste("Ks Test", round(gof.rot, 3))));abline(0,1, col = "red", lwd=1.5)

# Gamma case
sig2 <- 1
y8 <- rgamma( length(eta), exp(eta)^2/sig2, scale = sig2/exp(eta) )
hist(y8)
moments::skewness(y8); moments::kurtosis(y8)

#Standard Normal distribution: skewness = 0, kurtosis = 3
#pgamma first rotation second results in high kurtosis
mu.x <- mean(y8)
var.x <- var(y8)
r8.unif <- pgamma(y8, mu.x^2/var.x, scale = var.x/mu.x)
r8.norm <- qnorm(r8.unif)
r8.gamma.rot <- solve(L, r8.norm)
hist(r8.gamma.rot)
moments::skewness(r8.gamma.rot); moments::kurtosis(r8.gamma.rot)

#rotating first gives same results as r8.gamma.rot
r8.rot <- solve(L, y8)
hist(r8.rot)


#Generate a DHARMa rotation scenario using the simulated RE, not obs to calculate Sigma for rotation
#y8.sim.rot <- matrix(0, 1000, 500)
y8.sim.mat <- eta.sim.mat <- matrix(0, 500, 1000)
for(i in 1:1000){
  eta <- as.numeric(rmvnorm(1, rep(0, ncol(C)), sigma=S))
  eta.sim.mat[,i] <- eta
  #use simulated fitted predicted response in transformed space to calculate covariance
  y8.sim <- rgamma(length(eta), exp(eta)^2/sig2, scale = sig2/exp(eta) )
  y8.sim.mat[,i] <- y8.sim
}

ks.test(createDHARMa(y8.sim.mat, as.vector(y8), rotation = NULL)$scaledResiduals, "punif")
ks.test(createDHARMa(y8.sim.mat, as.vector(y8), rotation = "estimated")$scaledResiduals, "punif")
ks.test(createDHARMa(y8.sim.mat, as.vector(y8), rotation = S)$scaledResiduals, "punif")
#rotate with just the correlation matrix? - need to test under simulation
ks.test(createDHARMa(y8.sim.mat, as.vector(y8), rotation = C)$scaledResiduals, "punif")
#rotate with just the correlation matrix of the simulation?
ks.test(createDHARMa(y8.sim.mat, as.vector(y8), 
                     rotation = as.matrix(Matrix::nearPD(cor(t(y8.sim.mat)))$mat))$scaledResiduals, "punif")

covar <- Matrix::nearPD(cov(t(eta.sim.mat)))$mat
L <- t(chol(covar))
y8.sim.rot <- apply(y8.sim.mat, 2, function(x) solve(L, x))
y8.rot <- solve(L, y8)
r8.sim <- ecdf(y8.sim.rot)(y8.rot)
hist(r8.sim)
gap::qqunif(r8.sim, logscale = FALSE)
hist(qnorm(r8.sim))
ks.test(r8.sim, "punif")
moments::skewness(qnorm(r8.sim)); moments::kurtosis(qnorm(r8.sim))

