library(tidyverse)
library(mvtnorm)
library(ellipse)



## Cole's experiments to explore why marginals aren't working
set.seed(1) # 1 is good
n <- 3
N <- 5000
S <- toeplitz((n:1)/n)
C <- rWishart::rWishart(1,n,Sigma=S)[,,1]
L <- t(chol(C))
mean <- rep(0,n)
obs <- as.numeric(rmvnorm(1, mean=mean, sigma=C))
draws <- rmvnorm(n=N, mean=mean,sigma=C)
## Correctly rotated space
draws.rotated <- t(solve(L, t(draws)) )
obs.rotated <- solve(L, obs)
cols <- c(rgb(0,0,0,.5), rgb(1,0,0,.5))
png('plots/demo_rotation.png', width=5, height=3.5, units='in', res=1000)
par(mfcol=c(n,n), mar=0*c(.65,.65,.65,.65), oma=c(.5,.5,.5,0),
    mgp=c(.5,.1,0), tck=-.01, cex.axis=.6 )
for(i in 1:n){
  for(j in 1:n){
    xlim <- range(c(draws[, c(i,j)], draws.rotated[,c(i,j)]))
    ylim <- range(c(draws[, c(i,j)], draws.rotated[,c(i,j)]))
    if(i==j){
      x <- seq(xlim[1], xlim[2], len=500)
      y1 <- dnorm(x, 0, sd=sqrt(C[i,i]))
      y2 <- dnorm(x,0,sd=1)
      plot(x,y1, ylim=c(0, max(c(y1,y2))*1.4), type='l',
           col=cols[1], axes=FALSE)
      lines(x,y2, col=cols[2])
      ## hist(draws[,j], xlab='', main='', ylab='', col=cols[1], border=cols[1],
      ##      ylim=c(0,.5), xlim=xlim, freq=FALSE)
      ## abline(v=obs[j], col=cols[1])
      points(c(obs[j], obs[j])[2],
             c(0,dnorm(obs[j],0,sqrt(C[j,j])))[2], col=cols[1], pch=15)
      lines(c(obs[j], obs[j]), c(0,dnorm(obs[j],0,sqrt(C[j,j]))), col=cols[1])
      mtext(line=-1, col=cols[1], paste("Marginal percentile=", round(mean(obs[j] > draws[,j]),3)), cex=.7)
      ## hist(draws.rotated[,j], xlab='', main='', ylab='',
      ##      add=TRUE, col=cols[2], freq=FALSE, border=cols[2])
      ## abline(v=obs.rotated[j], col=cols[2])
      points(c(obs.rotated[j], obs.rotated[j])[2],
             c(0,dnorm(obs.rotated[j],0,sqrt(C[j,j])))[2], col=cols[2], pch=15)
      lines(c(obs.rotated[j], obs.rotated[j]), c(0,dnorm(obs.rotated[j],0,sqrt(C[j,j]))), col=cols[2])
      mtext(line=-2, col=cols[2], paste("Marginal percentile=", round(mean(obs.rotated[j] > draws.rotated[,j]),3)), cex=.7)
      box()
    }
    if(i<j){
      ## plot(draws[,i], draws[,j], ann=FALSE, pch=16, cex=.25,
      ##      col=cols[1])
      plot(obs[i], obs[j],col=4, cex=1, pch=16, xlim=xlim,
           ylim=ylim, axes=FALSE)
      ## points(draws.rotated[,i], draws.rotated[,j], ann=FALSE,
      ##        pch=16, cex=.25,  col=cols[2])
      ## points(obs[i], obs[j],col=3, cex=2, pch=16)
      arrows(obs[i], obs[j], obs.rotated[i], obs.rotated[j],
             length=.05, lwd=1.5, col=4)
      lines(ellipse(C[c(i,j),c(i,j)], centre=mean[c(i,j)]), col=cols[1])
      lines(ellipse(diag(2), centre=mean[c(i,j)]),  col=cols[2])
      box()
    }
    if(i>j) {plot(1,1, type='n', axes=FALSE, ann=FALSE)}
  }
}
dev.off()




## ## Cole's experiments to explore why marginals aren't working
## set.seed(122)
## n <- 4
## N <- 1000
## S <- toeplitz((n:1)/n)
## ## C <- rWishart::rWishart(1,n,Sigma=(1:n)*diag(n))[,,1]
## C <- rWishart::rWishart(1,n,Sigma=S)[,,1]
## L <- t(chol(C))
## Lbad <- diag(sqrt(diag(C)))
## obs <- as.numeric(rmvnorm(1,sigma=C))
## draws <- rmvnorm(n=N, sigma=C)
## ## Correctly rotated space
## draws.rotated <- t(solve(L, t(draws)) )
## obs.rotated <- solve(L, obs)
## ## incorrectly rotated
## draws.bad <- t(solve(Lbad, t(draws)))
## obs.bad <- solve(Lbad, obs)
## tmp <- 10
## for(step in 1:3){
##   png(paste0('plots/demo_pairs',step,'.png'), width=5, height=3.5, units='in', res=1000)
##   par(mfcol=c(n,n), mar=c(.65,.65,.65,.65), oma=c(.5,.5,.5,0),
##       mgp=c(.5,.1,0), tck=-.01, cex.axis=.6 )
##   for(i in 1:n){
##     for(j in 1:n){
##       if(i==j){
##         if(step==1){
##           hist(draws[,j], xlab='', main='', ylab='')
##           abline(v=obs[j], col='red')
##           mtext(paste("Marginal percentile=", round(mean(obs[j] > draws[,j]),3)), cex=.5)
##         }
##         if(step==2){
##           hist(draws.bad[,j], xlab='', main='', ylab='')
##           abline(v=obs.bad[j], col='red')
##           mtext(paste("Marginal percentile=", round(mean(obs.bad[j] > draws.bad[,j]),3)), cex=.5)
##         }
##         if(step==3){
##           hist(draws.rotated[,j], xlab='', main='', ylab='')
##           abline(v=obs.rotated[j], col='red')
##           mtext(paste("Marginal percentile=", round(mean(obs.rotated[j] > draws.rotated[,j]),3)), cex=.5)
##         }
##         box()
##       }
##       if(i<j){
##         if(step==1){
##           plot(draws[,i], draws[,j], ann=FALSE, pch=16, cex=.25, col=rgb(0,0,0,.5))
##         }
##         if(step==2){
##           plot(draws.bad[,i], draws.bad[,j], ann=FALSE, pch=16, cex=.25, col=rgb(0,0,0,.5))
##           arrows(draws[1:tmp, i], draws[1:tmp,j], draws.bad[1:tmp,i],
##                  draws.bad[1:tmp,j], length=.05, col=4, lwd=1.5)
##           arrows(obs[i], obs[j], obs.bad[i], obs.bad[j], length=.05,
##                  col=2, lwd=1.5)
##         }
##         if(step==3){
##           plot(draws.rotated[,i], draws.rotated[,j], ann=FALSE, pch=16, cex=.25, col=rgb(0,0,0,.5))
##           arrows(draws[1:tmp, i], draws[1:tmp,j], draws.rotated[1:tmp,i],
##                  draws.rotated[1:tmp,j], length=.05, col=4,
##                  lwd=1.5)
##           arrows(obs[i], obs[j], obs.rotated[i], obs.rotated[j],
##                  length=.05, lwd=1.5, col=2)
##         }
##         points(draws[1:tmp,i], draws[1:tmp,j], pch=21, col=4, bg='white')
##         points(obs[i], obs[j], pch=21, col=2, bg='white')
##       }
##       if(i>j) {plot(1,1, type='n', axes=FALSE, ann=FALSE)}
##     }
##   }
## dev.off()
## }
## ## Massage for ggplot
## xx <- rbind(data.frame(type='OSA', draws.rotated),
##             data.frame(type='dharma', draws.bad)) %>%
##   pivot_longer(-type)
## yy <- rbind(data.frame(type='OSA', t(obs.rotated)),
##             data.frame(type='dharma', t(obs.bad))) %>%
##   pivot_longer(-type)
## g <- ggplot(xx, aes(value, fill=type)) +
##   geom_histogram(bins=50, position='identity', alpha=.5) +
##   facet_wrap('name') +
##   geom_vline(data=yy, aes(xintercept=value, color=type))
## g
## ggsave('plots/demo_marginals.png', g, width=7, height=5)
## ## Try to visualize the rotations by conneting points
## nshow <- 7
## png('plots/demo_arrows.png', width=7, height=5, units='in', res=500)
## par(mfcol=c(2,2), mgp=c(1,.5,0), mar=c(2,2,.5,.5))
## for(i in 1:ncol(draws)){
##   xlim <- range(c(draws[1:nshow,i], draws.rotated[1:nshow,],
##                   draws.bad[1:nshow,i]))
##    xlim <- range(c(draws[1:nshow,], draws.rotated[1:nshow,],
##                    draws.bad[1:nshow,]))
##   xlim <- c(-2.5,2.5)
##   plot(0, type='n', xlim=xlim,
##        ylim=c(-1.05,1.05), ann=FALSE)
##   ## arrows(obs[1], obs[2], obs.rotated[1], obs.rotated[2])
##   ## arrows(draws[,1], draws[,2], draws.rotated[,2],
##   ##        draws.rotated[,2])
##   arrows(draws[1:nshow,i], 0, draws.rotated[1:nshow,i], -1, length=.05)
##   arrows(draws[1:nshow,i], 0, draws.bad[1:nshow,i], 1, length=.05, col=2)
##   points(draws[1:nshow,i], y=rep(0, nshow), pch=16)
## }
## dev.off()
