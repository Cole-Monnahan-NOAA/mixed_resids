## estimate and validate a random walk model with and without
## drift. based on code example provided by uffe høgsbro thygesen
## and kasper kristensen, 2016, in the tmb package examples:
## randomwalkvalidation.r.

## modified starting 11/2020 by cole

## randomwalk is conditional, randomwalk2 is unconditional
compile("models/randomwalk.cpp") # modified for simulation
dyn.load(dynlib("models/randomwalk"))

## simulate data with these parameters
set.seed(2342)
mu <- 0.75
sigma <- 1
s <- 1
huge <- 1e3
## simulate random track
nt <- 100
X <- c(0,cumsum(rnorm(nt-1,mean=mu,sd=sigma)))
## simulate random measurements
Y <- X + rnorm(nt,sd=s)
data <- list(y=Y,huge=huge)
parameters <- list(x=X, mu=0, logsigma=log(sigma), logs=log(s))

## estimate states and parameters under h0: mu=0
obj0 <- MakeADFun(data, parameters, random=c("x"),
                  dll="randomwalk", map=list(mu=factor(NA)))
opt0 <- do.call("optim",obj0)
sdr0 <- sdreport(obj0)
estX0 <- summary(sdr0,"random")
## estimate states and parameters under h1: mu != 0
obj1 <- MakeADFun(data, parameters, random=c("x"), dll="randomwalk")
opt1 <- do.call("optim",obj1)
sdr1 <- sdreport(obj1)
estX1 <- summary(sdr1,"random")

### Naive normalized residuals
resid0 <- Y-estX0[,1]
resid1 <- Y-estX1[,1]
Norm.resid0 <- resid0 / estX0[,2]
Norm.resid1 <- resid1 / estX1[,2]

### OSA residuals
## Generate one step predictions with the models fitted under H0 and H1
predict0  <- oneStepPredict(obj0, observation.name="y", method="fullGaussian")$residual
predict1  <- oneStepPredict(obj1, observation.name="y", method="fullGaussian")$residual

### DHARMa resids
library(DHARMa)
sim0 <- replicate(1000, {obj0$simulate()$y})
sim0.d <- createDHARMa(sim0, Y, integerResponse=FALSE)
sim1 <- replicate(1000, {obj1$simulate()$y})
sim1.d <- createDHARMa(sim1, Y, integerResponse=FALSE)
## I don't know the right way to do this yet
sim0.resid <- qnorm(sim0.d$scaledResiduals)
sim1.resid <- qnorm(sim1.d$scaledResiduals)
plot(sim0.d, quantreg = FALSE)
plot(sim1.d, quantreg = FALSE)
## matplot(sim0[,1:5])
## matplot(sim1[,1:5])

resids <- rbind(data.frame(model='m0', Pearson=Norm.resid0, OSA=predict0, simulation=sim0.resid),
      data.frame(model='m1', Pearson=Norm.resid1, OSA=predict1,
                 simulation=sim1.resid)) %>% cbind(x=1:100)
resids.long <- resids %>% pivot_longer(-c(x, model))

ggplot(resids.long, aes(x, value, color=name)) +
  geom_jitter(width=.2, alpha=.6) + facet_wrap('model', ncol=1)
## ggplot(resids, aes(x=name, value, color=name)) + geom_violin() + facet_wrap('model')

library(GGally)
ggpairs(resids, columns=2:4, mapping=aes(color=model))
