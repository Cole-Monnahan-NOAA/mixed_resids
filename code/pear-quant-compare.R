set.seed(123)
y1 <- rgamma(1000, shape = 100, scale = .1)
hist(y1)
y2 <- rgamma(1000, shape = .5, scale = 3)
hist(y2)
eta.mat <- as.matrix(exp(rnorm(5, 0, .5)))
sig2 <- .001
y3 <- apply(eta.mat, 1, function(x) rgamma(200, x^2/sig2, scale = sig2/x ))
hist(y3)
sig2 <- 1
y4 <- apply(eta.mat, 1, function(x) rgamma(200, x^2/sig2, scale = sig2/x ))
hist(y4)

output <- function(y){

  if(is.matrix(y)){
    #sig2. <- var(as.vector(y))
    Sig <- cov(y)
    L <- t(chol(Sig))
    mode <- apply(y, 2, function(x) x-mean(x))
    y <- as.vector(apply(mode, 1, function(x) solve(L,x)))
    quant.res <- qnorm(pnorm(y))
    }
  mu <- mean(y)
  sig2. <- var(y)
  pear.res <- (y-mu)/sqrt(sig2.)
  if(!exists("quant.res")){
    Fx <- pgamma(y, shape = mu^2/sig2., scale = sig2./mu)
    quant.res <- qnorm(Fx)
  }
  pear.gof <- ks.test(pear.res, "pnorm")$p.value
  quant.gof <- ks.test(quant.res, "pnorm")$p.value
  print(hist(y))
  par(mfrow = c(1,2))
  print(qqnorm(pear.res, 
               main = "Pearson Normal Q-Q Plot", 
               xlab = paste("Ks Test", round(pear.gof, 3))));abline(0,1, col = "red", lwd=1.5)
  print(qqnorm(quant.res, 
               main = "Quantile Normal Q-Q Plot", 
               xlab = paste("Ks Test", round(quant.gof, 3))));abline(0,1, col = "red", lwd=1.5)
}
output(y1)
output(y2)
#approx normal gamma with group RE term
output(as.vector(y3))
#apply rotation - fixes problem
output(y3)
#skewed gamma with group RE term
output(as.vector(y4))
#apply rotation - residuals still correlated
output(y4)
