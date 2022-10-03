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
library(goftest)


## Some global settings
ggwidth <- 7
ggheight <- 5
theme_set(theme_bw())
source('R/sim_data.R')
source('R/utils.R')
source('R/model_fns.R')
source('R/resid_fns.R')


message("Compiling models if needed...")
TMB::compile('src/linmod.cpp')
TMB::compile('src/randomwalk.cpp')
TMB::compile('src/simpleGLMM.cpp')
TMB::compile('src/spatial.cpp', framework = "TMBad")
