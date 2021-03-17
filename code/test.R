library(DHARMa)
library(mvtnorm)
library(magrittr)
library(tidyr)
library(ggplot2)
library(gridExtra)
n.sim <- 500
n.samp <- 100
Range <- c(20,60)
Kappa <- sqrt(8)/Range
CV <- c(0.5,1,2)
sig2 <- 1
B0 <- 1
n <- 100
results <- list()



for(r in 1:length(Range)){
  for(c in 1:length(CV)){
    nm <- paste0('R',r,'CV',c)
    p.val.omega <- matrix(NA,n.sim,3)
    p.val.y <- matrix(NA,n.sim,10)
    for(ii in 1:n.sim){
      set.seed(ii) 
      Loc <- matrix(runif(n*2,0,100), ncol=2)
      dmat <- as.matrix(dist(Loc, upper = TRUE, diag = TRUE))
      mu <- rep(0,n)
      
      Sigma <- diag(n)*sig2
      cnt <- 1
      for(j1 in 2:n){
        for(j2 in 1:(j1-1)){
          Sigma[j1,j2] <- Sigma[j2,j1] <- sig2 * 
            (Kappa[r] * dmat[j1,j2]) * besselK(Kappa[r] * dmat[j1,j2], 1) 
        }
      }
      set.seed(ii)
      u <- replicate(1000,{rnorm(n, 0, sqrt(sig2))  })
      set.seed(ii)
      omega1 <- t(mvtnorm::rmvnorm(1000,rep(0,n),Sigma))
      res1 <- createDHARMa(u[,2:1000], u[,1])$scaledResiduals
      res2 <- createDHARMa(omega1[,2:1000],omega1[,1])$scaledResiduals
      omega2 <-  solve(t(chol(Sigma)), omega1)
      ## correction factor
      res3 <- createDHARMa(omega2[,2:1000], omega2[,1])$scaledResiduals
      p.val.omega[ii,1] <- suppressWarnings(ks.test(res1, 'punif')$p.value)
      p.val.omega[ii,2] <- suppressWarnings(ks.test(res2, 'punif')$p.value) #matern
      p.val.omega[ii,3] <- suppressWarnings(ks.test(res3, 'punif')$p.value)
      
      #Simulate obs
      y0 <- rgamma(n, shape = 1/CV[c]^2, scale = CV[c]^2*exp(B0 + u[,1]))
      y1 <- rgamma(n, shape = 1/CV[c]^2, scale = CV[c]^2*exp(B0 + u[,1] + omega1[,1]))
      
      #Simulate testing data
      ysim0.c <- replicate(1000, {rgamma(n, shape = 1/CV[c]^2, scale = CV[c]^2*exp(B0 + u[,1]))})
      ysim1.c <- replicate(1000, {rgamma(n, shape = 1/CV[c]^2, scale = CV[c]^2*exp(B0 + u[,1] + omega1[,1]))})
      ysim0.u <- ysim1.u <- ysim1.u.adj <- matrix(0,n,1000)
      for(i in 1:1000){
        ysim0.u[,i] <- rgamma(n, shape = 1/CV[c]^2, scale = CV[c]^2*exp(B0 + u[,i]))  
        ysim1.u[,i] <- rgamma(n, shape = 1/CV[c]^2, scale = CV[c]^2*exp(B0 + u[,i] + omega1[,i])) 
      }
      res0.c.y0 <- createDHARMa(ysim0.c,y0)
      res0.u.y0 <- createDHARMa(ysim0.u,y0)
      res0.c.y1 <- createDHARMa(ysim0.c,y1)
      res0.u.y1 <- createDHARMa(ysim0.u,y1)
      res1.c.y0 <- createDHARMa(ysim1.c,y0)
      res1.u.y0 <- createDHARMa(ysim1.u,y0)
      res1.u.adj.y0 <- createDHARMa(exp(solve(t(chol(Sigma)), log(ysim1.u))), exp(solve(t(chol(Sigma)), log(y0))))
      res1.c.y1 <- createDHARMa(ysim1.c,y1)
      res1.u.y1 <- createDHARMa(ysim1.u,y1)
      res1.u.adj.y1 <- createDHARMa(exp(solve(t(chol(Sigma)), log(ysim1.u))), exp(solve(t(chol(Sigma)), log(y1))))
      
      p.val.y[ii,1] <- suppressWarnings(ks.test(res0.c.y0$scaledResiduals, 'punif')$p.value)
      p.val.y[ii,2] <- suppressWarnings(ks.test(res0.u.y0$scaledResiduals, 'punif')$p.value)
      p.val.y[ii,3] <- suppressWarnings(ks.test(res1.c.y0$scaledResiduals, 'punif')$p.value)
      p.val.y[ii,4] <- suppressWarnings(ks.test(res1.u.y0$scaledResiduals, 'punif')$p.value)
      p.val.y[ii,5] <- suppressWarnings(ks.test(res1.u.adj.y0$scaledResiduals, 'punif')$p.value)
      p.val.y[ii,6] <- suppressWarnings(ks.test(res0.c.y1$scaledResiduals, 'punif')$p.value)
      p.val.y[ii,7] <- suppressWarnings(ks.test(res0.u.y1$scaledResiduals, 'punif')$p.value)
      p.val.y[ii,8] <- suppressWarnings(ks.test(res1.c.y1$scaledResiduals, 'punif')$p.value)
      p.val.y[ii,9] <- suppressWarnings(ks.test(res1.u.y1$scaledResiduals, 'punif')$p.value)
      p.val.y[ii,10] <- suppressWarnings(ks.test(res1.u.adj.y1$scaledResiduals, 'punif')$p.value)
    }
    results[[nm]] <- list(p.val.omega, p.val.y)
  }
}

c.names <-  c("sim0 Conditional y0", "sim0 Unconditional y0",
              "sim1 Conditional y0", "sim1 Unconditional y0", "sim1adj Unconditional y0",
              "sim0 Conditional y1", "sim0 Unconditional y1",
              "sim1 Conditional y1", "sim1 Unconditional y1", "sim1adj Unconditional y1")

plot.fun <- function(i, obs){
  p.val <- as.data.frame(results[[i]][[2]])
  colnames(p.val) <- c.names
  p.val.df <- pivot_longer(p.val, starts_with("sim"), names_to = c("H", "Sim", "Obs"), 
                           names_sep = " ", values_to = "pval")
  p <- p.val.df %>% dplyr::filter(Obs == obs) %>% 
    ggplot(., aes(x=pval)) + facet_grid(H ~ Sim, scales = 'free_y') + geom_histogram(binwidth = .01)
  return(p)
}

p0 <- lapply(1:6, function(x) plot.fun(x,"y0"))
p1 <- lapply(1:6, function(x) plot.fun(x,"y1"))

col.nm <- c("CV = 0.5", "CV = 1", "CV = 2")
row.nm <- c("Range = 20", "Range = 60")

png(file = "plots/sp+re_y0.png", width = 1200, height = 800)
grid.arrange(arrangeGrob(p0[[1]], left = row.nm[1], top = col.nm[1]),
             arrangeGrob(p0[[2]], top = col.nm[2]),
             arrangeGrob(p0[[3]], top = col.nm[3]),
             arrangeGrob(p0[[4]], left = row.nm[2]),
             arrangeGrob(p0[[5]]),
             arrangeGrob(p0[[6]]), ncol = 3
)
dev.off()

png(file = "plots/sp+re_y1.png", width = 1200, height = 800)
grid.arrange(arrangeGrob(p1[[1]], left = row.nm[1], top = col.nm[1]),
             arrangeGrob(p1[[2]], top = col.nm[2]),
             arrangeGrob(p1[[3]], top = col.nm[3]),
             arrangeGrob(p1[[4]], left = row.nm[2]),
             arrangeGrob(p1[[5]]),
             arrangeGrob(p1[[6]]), ncol = 3
)
dev.off()

p.val.omega <- as.data.frame(results[[1]][[1]])
colnames(p.val.omega) <- c('u', 'omega', 'omega.adj')
p.val.omega <- pivot_longer(p.val.omega, cols = 1:3, names_to = 'sim', values_to = 'pval')
png(file = "plots/re.png", width = 800, height = 300)
p.val.omega %>% 
  ggplot(., aes(x = pval)) + facet_wrap(~sim, nrow=1) + geom_histogram()
dev.off()