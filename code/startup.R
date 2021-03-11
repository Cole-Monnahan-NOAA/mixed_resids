library(MASS) # rmvnorm
library(TMB)
library(INLA)
library(DHARMa)
## library(VAST)
library(ggplot2)
library(dplyr)
library(tidyr)
## Has ggpairs for comparing residual types
library(GGally)
library(snowfall)
library(R.utils)
library(Matrix)

## Some global settings
ggwidth <- 7
ggheight <- 5
theme_set(theme_bw())

message("Loading functions...")
## Function to simulate parameters from the joint precisions
## matrix (fixed + random effects). Modified from
## FishStatsUtils::simulate_data
rmvnorm_prec <- function(mu, prec ) {
  ##set.seed( random_seed )
  z = matrix(rnorm(length(mu)), ncol=1)
  L = Matrix::Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.vector(z)
  return(mu + z)
}

## Quick fn to check for failed runs by looking at results output
## that doesn't exist
which.failed <- function(Nreps){
  success <- gsub('results/spatial_pvals/pvals_|.RDS', "", x=fs) %>%
    as.numeric()
  fail <- which(! 1:Nreps %in% success)
  fail
}


add_aic <- function(opt,n){
  opt$AIC <- TMBhelper::TMBAIC(opt, n=Inf)
  opt$AICc <- TMBhelper::TMBAIC(opt, n=n)
  opt$BIC <- TMBhelper::TMBAIC(opt, p=log(n))
  opt
}
