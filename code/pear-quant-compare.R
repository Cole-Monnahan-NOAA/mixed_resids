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
y7 <- eta
r7.nrot <- y7 #mean=0; sd=1, this is essentially the Pearson residual (y7-0)/1
r7.rot <- solve(L, y7)
gof.nrot <- ks.test(r7.nrot, "pnorm")$p.value
gof.rot <- ks.test(r7.rot, "pnorm")$p.value

par(mfrow = c(2,3))
boxplot(t(y7), horizontal=TRUE)
(qqnorm(r7.nrot,
        main = "Unrotated Normal Q-Q Plot",
        xlab = paste("Ks Test", round(gof.nrot, 3))));abline(0,1, col = "red", lwd=1.5)
(qqnorm(r7.rot,
        main = "Rotated Normal Q-Q Plot",
        xlab = paste("Ks Test", round(gof.rot, 3))));abline(0,1, col = "red", lwd=1.5)

# Gamma case
sig2 <- 1
y8 <- t(apply(exp(eta), 1:2, function(x) rgamma(1, x^2/sig2, scale = sig2/x )))
hist(y8)
moments::skewness(as.vector(y8)); moments::kurtosis(as.vector(y8))

#pgamma removes skewness and lowers kurtosis but returns zeros - cannot convert with qnorm, returns -Inf values
r8.unif <- apply(y8, 2, function(x) pgamma(x, mean(x)^2/var(x), scale = var(x)/mean(x)))
hist(r8.unif)
moments::skewness(as.vector(r8.unif)); moments::kurtosis(as.vector(r8.unif))

#rotating first removes correlation but retains kurtosis - cannot use pgamma because negative values
S <- cov(y8)
round(S,2)
L <- t(chol(S))
mode <- apply(y8, 2, function(x) x - mean(x))
r8.rot <- t(solve(L, t(mode)))
hist(r8.rot)
mean(r8.rot) #mean of 0
var(as.vector(r8.rot)) # var of 1
moments::skewness(as.vector(r8.rot))
moments::kurtosis(as.vector(r8.rot))
