## IID =======================================================================
set.seed(123)
# Approx. normal gamma
y1 <- rgamma(1000, shape = 100, scale = .1)# mean=10; var=1
hist(y1)
# Skewed gamma
y2 <- rgamma(1000, shape = .1, scale = sqrt(10))# mean =approx. .316; var = 1
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
eta.mat <- as.matrix(exp(2 + rnorm(5, 0, .5)))
# Approx. normal gamma
alpha <- 100
y3 <- apply(eta.mat, 1, function(x) rgamma(200, shape = alpha, scale = x/alpha ))
hist(y3)
# Skewed gamma
#eta.mat <- as.matrix(exp(2 + rnorm(5, 0.4, .5)))
alpha <- 0.1
y4 <- apply(eta.mat, 1, function(x) rgamma(200, shape = alpha, scale = x/alpha ))
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
#residual by group specific mean and sd
y3rot <- output.groupRE(y3)

par(mfrow = c(2,3))
#skewed gamma with group RE term
output.iid(as.vector(y4))
output.groupRE(y4)

dev.off()

## Banded Correlation Structure ====================================================
library(mvtnorm)

# Normal case
set.seed(1)
C <- exp(-as.matrix(dist(seq(0,10,length.out = 500))))
sig <- 3
S <- diag(rep(sig, ncol(C))) %*% C %*%  diag(rep(sig, ncol(C)))
eta <- as.numeric(rmvnorm(1, rep(0, ncol(C)), sigma=S)) 
L <- t(chol(S))

y7 <- eta + 2
mode <- y7 - mean(y7)
r7.pears <- mode/sd(y7)
r7.rot <- t(solve(L, mode))
gof.pears <- ks.test(r7.pears, "pnorm")$p.value
gof.rot <- ks.test(r7.rot, "pnorm")$p.value

par(mfrow = c(1,3))
hist(y7)
(qqnorm(r7.pears,
        main = "Pearson Normal Q-Q Plot",
        xlab = paste("Ks Test", round(gof.pears, 3))));abline(0,1, col = "red", lwd=1.5)
(qqnorm(r7.rot,
        main = "Rotated Normal Q-Q Plot",
        xlab = paste("Ks Test", round(gof.rot, 3))));abline(0,1, col = "red", lwd=1.5)

set.seed(123)
# Gamma case
alpha = .1
y8 <- rgamma( length(eta), shape = alpha, scale = exp(eta+2)/alpha )
mu.x <- mean(y8)
var.x <- var(y8)
mode <- y8 - mu.x
r8.pears <- mode/sqrt(var.x)
r8 <- qnorm(pgamma(y8, mu.x^2/var.x, scale = var.x/mu.x))
r8.rot <- solve( t(chol(C)), r8)
hist(r8);hist(r8.rot)
ks.test(r8, "pnorm");ks.test(r8.rot, "pnorm")

gof.pears <- ks.test(r8.pears, "pnorm")$p.value
gof.quant <- ks.test(r8, "pnorm")$p.value
gof.rot <- ks.test(r8.rot, "pnorm")$p.value

par(mfrow = c(2,4))
hist(r8.pears);hist(r8);hist(r8.rot);hist(y8)
(qqnorm(r8.pears,
        main = "Pearson Skewed Q-Q Plot",
        xlab = paste("Ks Test", round(gof.pears, 3))));abline(0,1, col = "red", lwd=1.5)
(qqnorm(r8,
        main = "Rotated Skewed Q-Q Plot",
        xlab = paste("Ks Test", round(gof.quant, 3))));abline(0,1, col = "red", lwd=1.5)
(qqnorm(r8.rot,
        main = "Rotated Skewed Q-Q Plot",
        xlab = paste("Ks Test", round(gof.rot, 3))));abline(0,1, col = "red", lwd=1.5)



