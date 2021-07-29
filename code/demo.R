library(tidyverse)
library(mvtnorm)

## Cole's experiments to explore why marginals aren't working
set.seed(12)
n <- 9
C <- rWishart::rWishart(1,n,diag(n))[,,1]
L <- t(chol(C))
Lbad <- diag(sqrt(diag(C)))
obs <- as.numeric(rmvnorm(1,sigma=C))
draws <- rmvnorm(n=10000, sigma=C)

## Correctly rotated space
draws.rotated <- t(solve(L, t(draws)) )
obs.rotated <- solve(L, obs)
## incorrectly rotated
draws.bad <- t(solve(Lbad, t(draws)))
obs.bad <- solve(Lbad, obs)

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
ggsave('plots/demo_marginals.png', g, width=7, height=5)

## Try to visualize the rotations by conneting points
nshow <- 15
png('plots/demo_arrows.png', width=7, height=5, units='in', res=500)
par(mfrow=c(3,3), mgp=c(1,.5,0), mar=c(2,2,.5,.5))
for(i in 1:ncol(draws)){
  plot(0, type='n', xlim=range(draws[1:nshow,i]),
       ylim=c(-1.05,1.05), ann=FALSE)
  ## arrows(obs[1], obs[2], obs.rotated[1], obs.rotated[2])
  ## arrows(draws[,1], draws[,2], draws.rotated[,2],
  ##        draws.rotated[,2])
  arrows(draws[1:nshow,i], 0, draws.rotated[1:nshow,i], 1, length=.05)
  arrows(draws[1:nshow,i], 0, draws.bad[1:nshow,i], -1, length=.05, col=2)
}
dev.off()
