library(tidyverse)
library(mvtnorm)

## Cole's experiments to explore why marginals aren't working
set.seed(122)
n <- 4
S <- toeplitz((n:1)/n)
C <- rWishart::rWishart(1,n,Sigma=(1:n)*diag(n))[,,1]
C <- rWishart::rWishart(1,n,Sigma=S)[,,1]
L <- t(chol(C))
Lbad <- diag(sqrt(diag(C)))
obs <- as.numeric(rmvnorm(1,sigma=C))
draws <- rmvnorm(n=1000, sigma=C)
## Correctly rotated space
draws.rotated <- t(solve(L, t(draws)) )
obs.rotated <- solve(L, obs)
## incorrectly rotated
draws.bad <- t(solve(Lbad, t(draws)))
obs.bad <- solve(Lbad, obs)
png('plots/demo_pairs1.png', width=7, height=5, units='in', res=500)
pairs(draws)
dev.off()
png('plots/demo_pairs2.png', width=7, height=5, units='in', res=500)
pairs(rbind(draws, draws.bad),
      col=rep(1:2, each=1000))
dev.off()
png('plots/demo_pairs3.png', width=7, height=5, units='in', res=500)
pairs(rbind(draws.bad, draws.rotated),
      col=rep(2:3, each=1000))
dev.off()


## Massage for ggplot
xx <- rbind(data.frame(type='OSA', draws.rotated),
            data.frame(type='dharma', draws.bad)) %>%
  pivot_longer(-type)
yy <- rbind(data.frame(type='OSA', t(obs.rotated)),
            data.frame(type='dharma', t(obs.bad))) %>%
  pivot_longer(-type)
g <- ggplot(xx, aes(value, fill=type)) +
  geom_histogram(bins=50, position='identity', alpha=.5) +
  facet_wrap('name') +
  geom_vline(data=yy, aes(xintercept=value, color=type))
g
ggsave('plots/demo_marginals.png', g, width=7, height=5)

## Try to visualize the rotations by conneting points
nshow <- 7
png('plots/demo_arrows.png', width=7, height=5, units='in', res=500)
par(mfcol=c(2,2), mgp=c(1,.5,0), mar=c(2,2,.5,.5))
for(i in 1:ncol(draws)){
  xlim <- range(c(draws[1:nshow,i], draws.rotated[1:nshow,],
                  draws.bad[1:nshow,i]))
   xlim <- range(c(draws[1:nshow,], draws.rotated[1:nshow,],
                   draws.bad[1:nshow,]))
  xlim <- c(-2.5,2.5)
  plot(0, type='n', xlim=xlim,
       ylim=c(-1.05,1.05), ann=FALSE)
  ## arrows(obs[1], obs[2], obs.rotated[1], obs.rotated[2])
  ## arrows(draws[,1], draws[,2], draws.rotated[,2],
  ##        draws.rotated[,2])
  arrows(draws[1:nshow,i], 0, draws.rotated[1:nshow,i], -1, length=.05)
  arrows(draws[1:nshow,i], 0, draws.bad[1:nshow,i], 1, length=.05, col=2)
  points(draws[1:nshow,i], y=rep(0, nshow), pch=16)
}
dev.off()
