library(MASS) # rmvnorm
library(TMB)
library(fmesher)
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
library(tweedie)
library(car)
library(ape)
library(phylolm)

## Some global settings
ggwidth <- 7
ggheight <- 5
theme_set(theme_bw())
source('R/sim_data.R')
source('R/utils.R')
source('R/model_fns.R')
source('R/resid_fns.R')


message("Compiling models if needed...")
TMB::compile('src/linmod.cpp', framework = "TMBad")
TMB::compile('src/randomwalk.cpp', framework = "TMBad")
TMB::compile('src/simpleGLMM.cpp', framework = "TMBad")
TMB::compile('src/spatial.cpp', framework = "TMBad")
TMB::compile('src/pois_glm.cpp', framework = "TMBad")
TMB::compile('src/phylo.cpp', framework = "TMBad")

