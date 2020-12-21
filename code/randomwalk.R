## estimate and validate a random walk model with and without
## drift. based on code example provided by uffe høgsbro thygesen
## and kasper kristensen, 2016, in the tmb package examples:
## randomwalkvalidation.r.

## modified starting 11/2020 by cole

## randomwalk is conditional, randomwalk2 is unconditional
compile("models/randomwalk.cpp") # modified for simulation
dyn.load(dynlib("models/randomwalk"))

sim_pvalues_list <- list()

for(ii in 1:Nreps){
## simulate data with these parameters
message("Simulating data...")
set.seed(ii)
mu <- 0.75
sigma <- 1
s <- 1
huge <- 1e3
## simulate random track
nt <- 50
X <- c(0,cumsum(rnorm(nt-1,mean=mu,sd=sigma)))
## simulate random measurements
Y <- X + rnorm(nt,sd=s)
data <- list(y=Y,huge=huge)
parameters <- list(x=X, mu=0, logsigma=log(sigma), logs=log(s))

## estimate states and parameters under h0: mu=0
message("Optimizing two competing models...")
obj0 <- MakeADFun(data, parameters, random=c("x"),
                  dll="randomwalk", map=list(mu=factor(NA)))
trash <- obj0$env$beSilent()
opt0 <- do.call("optim",obj0)
sdr0 <- sdreport(obj0)
estX0 <- summary(sdr0,"random")
## estimate states and parameters under h1: mu != 0
obj1 <- MakeADFun(data, parameters, random=c("x"), dll="randomwalk")
trash <- obj1$env$beSilent()
opt1 <- do.call("optim",obj1)
sdr1 <- sdreport(obj1)
estX1 <- summary(sdr1,"random")

message("Calculating residuals..")
### Naive normalized residuals
pearson0 <- (Y-estX0[,1])/ estX0[,2]
pearson1 <- (Y-estX1[,1]) / estX1[,2]
### OSA residuals
## Generate one step predictions with the models fitted under H0
## and H1
osa0 <- tryCatch(oneStepPredict(obj0, observation.name="y",
                           method="fullGaussian")$residual,
                 error=function(e) 'error')
osa1 <- tryCatch(oneStepPredict(obj1, observation.name="y",
                                method="fullGaussian")$residual,
                 error=function(e) 'error')
  if(is.character(osa0) | is.character(osa1)){
    warning("OSA failed in rep=", ii)
    next
  }

  ### Validation using one sample from the posterior

C0 <- solve(obj0$env$spHess(random=TRUE))  ## Covariance matrix of random effects
Xr0 <- mvrnorm(1,estX0[,1],C0)             ## Generate one sample of random effects
W0 <- diff(Xr0)/exp(opt0$par["logsigma"])  ## Compute corresponding driving process noise

C1 <- solve(obj1$env$spHess(random=TRUE))  ## .. repeat for H1
Xr1 <- mvrnorm(1,estX1[,1],C1)
W1 <- (diff(Xr1)-opt1$par["mu"])/exp(opt1$par["logsigma"])


### DHARMa resids, both conditional and unconditional
tmp <- replicate(1000, {obj0$simulate()$y})
dharma0_cond <- createDHARMa(tmp, Y, integerResponse=FALSE)
tmp <- replicate(1000, {obj0$simulate()$y2})
dharma0_uncond <- createDHARMa(tmp, Y, integerResponse=FALSE)
tmp <- replicate(1000, {obj1$simulate()$y})
dharma1_cond <- createDHARMa(tmp, Y, integerResponse=FALSE)
tmp <- replicate(1000, {obj1$simulate()$y2})
dharma1_uncond <- createDHARMa(tmp, Y, integerResponse=FALSE)

  ## warning("don't know the right way to calculate DHARMa resids")
sim0_cond <- qnorm(dharma0_cond$scaledResiduals)
sim0_uncond <- qnorm(dharma0_uncond$scaledResiduals)
sim1_cond <- qnorm(dharma0_cond$scaledResiduals)
sim1_uncond <- qnorm(dharma0_uncond$scaledResiduals)

### Combine together in tidy format for analysis and plotting
d0 <- data.frame(x=1:nt, version='m0', pearson=pearson0,
                 osa=osa0, sim_cond=sim0_cond,
                 sim_uncond=sim0_uncond)
d1 <- data.frame(x=1:nt, version='m1', pearson=pearson1,
                 osa=osa1, sim_cond=sim1_cond,
                 sim_uncond=sim1_uncond)
resids <- rbind(d0, d1)
resids.long <- resids %>% pivot_longer(-c(x, version))

### Extract p-values calculated by DHARMa
## Note: Type binomial for continuous, if integer be careful. Not
## sure if we want two-sided for dispersion? Using defaults for
## now.
disp0_uncond <- testDispersion(dharma0_uncond, plot=FALSE)
outlier0_uncond <- testOutliers(dharma0_uncond, type='binomial', plot=FALSE)
disp1_uncond <- testDispersion(dharma1_uncond, plot=FALSE)
outlier1_uncond <- testOutliers(dharma1_uncond, type='binomial', plot=FALSE)
disp0_cond <- testDispersion(dharma0_cond, plot=FALSE)
outlier0_cond <- testOutliers(dharma0_cond, type='binomial', plot=FALSE)
disp1_cond <- testDispersion(dharma1_cond, plot=FALSE)
outlier1_cond <- testOutliers(dharma1_cond, type='binomial', plot=FALSE)
pvals <- rbind(
data.frame(version='m0', RE='cond', test='outlier', pvalue=outlier0_cond$p.value),
data.frame(version='m0', RE='uncond', test='outlier', pvalue=outlier0_uncond$p.value),
data.frame(version='m0', RE='cond', test='disp', pvalue=disp0_cond$p.value),
data.frame(version='m0', RE='uncond', test='disp', pvalue=disp0_uncond$p.value),
data.frame(version='m1', RE='cond', test='outlier', pvalue=outlier1_cond$p.value),
data.frame(version='m1', RE='uncond', test='outlier', pvalue=outlier1_uncond$p.value),
data.frame(version='m1', RE='cond', test='disp', pvalue=disp1_cond$p.value),
data.frame(version='m1', RE='uncond', test='disp', pvalue=disp1_uncond$p.value))
pvals$replicate <- ii
sim_pvalues_list[[ii]] <- pvals


### Exploratory plots for first replicate
  if(ii==1){
    g <- ggplot(resids.long, aes(x, value, color=name)) +
      geom_jitter(width=.2, alpha=.6) + facet_wrap('version', ncol=1) +
      labs(x='Time Step', y='Residual Value', title='Random Walk')
    ggsave('plots/randomwalk_resids_by_time.png', g, width=7, height=5)
    ## ggplot(resids, aes(x=name, value, color=name)) + geom_violin() + facet_wrap('version')

    g <- ggpairs(resids, columns=3:6, mapping=aes(color=version), title='Random Walk')
    ggsave('plots/randomwalk_resids_pairs.png', g, width=7, height=5)

    tmp <- data.frame(x=1:nt, Xtrue=X, X0=estX0[,1],
                      X1=estX1[,1])
    tmp2 <- data.frame(x=1:nt, X0_diff=estX0[,1]-X,
                      X1_diff=estX1[,1]-X)
    g <- pivot_longer(tmp2, -x) %>%
      ggplot(aes(x, value, color=name)) + geom_line() +
                      labs(y='Relative diff to true X')
    ggsave('plots/randomwalk_fits.png', g, width=7, height=5)

  ## What do the simulated data look like?
  ff <- function(x, v, re) data.frame(t=1:nt, version=v, RE=re, x$simulatedResponse[,1:20])
  g <- rbind(ff(dharma0_cond, 'm0', 'cond'),
             ff(dharma0_uncond, 'm0', 'uncond'),
             ff(dharma1_cond, 'm1', 'cond'),
             ff(dharma1_uncond, 'm1', 'uncond')) %>%
    pivot_longer(cols=c(-t, -version, -RE), names_prefix="X",
                 names_to='replicate', values_to='y') %>%
    mutate(replicate=as.numeric(replicate)) %>%
    ggplot(aes(t, y, group=replicate)) + geom_line(alpha=.5) +
    facet_grid(version~RE)
  ggsave('plots/randomwalk_simdata.png', g, width=7, height=5)
  }

}

sim_pvalues <- do.call(rbind, sim_pvalues_list)
g <- ggplot(sim_pvalues, aes(pvalue, )) + geom_histogram() +
  facet_grid(version+RE~test)
ggsave('plots/randomwalk_pvalues.png', g, width=7, height=5)
