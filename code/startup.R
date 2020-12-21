library(MASS) # rmvnorm
library(TMB)
library(DHARMa)
library(VAST)
library(ggplot2)
library(dplyr)
library(tidyr)
## Has ggpairs for comparing residual types
library(GGally)

## Some global settings
ggwidth <- 7
ggheight <- 5
theme_set(theme_bw())
