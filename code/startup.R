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
source("functions.R")
