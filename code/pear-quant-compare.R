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
    boxplot(y, horizontal=TRUE)
    #sig2. <- var(as.vector(y))
    Sig <- cov(y)
    L <- t(chol(Sig))
    mode <- apply(y, 2, function(x) x-mean(x))
    yrot <- t(apply(mode, 1, function(x) solve(L,x)))
    y <- as.vector(yrot)
    quant.res <- qnorm(pnorm(y))
  } else {
    hist(y)
    yrot <- y
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
  (qqnorm(pear.res,
               main = "Pearson Normal Q-Q Plot",
               xlab = paste("Ks Test", round(pear.gof, 3))));abline(0,1, col = "red", lwd=1.5)
  (qqnorm(quant.res,
               main = "Quantile Normal Q-Q Plot",
               xlab = paste("Ks Test", round(quant.gof, 3))));abline(0,1, col = "red", lwd=1.5)
  return(invisible(yrot))
}


par(mfrow = c(2,3))
output(y1)
output(y2)

par(mfrow = c(2,3))
#approx normal gamma with group RE term
output(as.vector(y3))
#apply rotation - fixes problem
y3rot <- output(y3)

par(mfrow = c(2,3))
#skewed gamma with group RE term
output(as.vector(y4))
#apply rotation - residuals still correlated
y4rot <- output(y4)

## What is the rotation doing?
pairs(rbind(y3,y3rot), col=rep(c(1,2), each=nrow(y3)))
pairs(rbind(y4,y4rot), col=rep(c(1,2), each=nrow(y4)))
