

if(file.exists('results/spatial_resids.RDS')){
  message("Loading spatial results")
  spatial_mles <- readRDS('results/spatial_mles.RDS')
  spatial_resids <- readRDS('results/spatial_resids.RDS')
  spatial_pvals <- readRDS('results/spatial_pvals.RDS')
}
if(file.exists('results/randomwalk_resids.RDS')){
  message("Loading randomwalk results")
  randomwalk_mles <- readRDS('results/randomwalk_mles.RDS')
  randomwalk_resids <- readRDS('results/randomwalk_resids.RDS')
  randomwalk_pvals <- readRDS('results/randomwalk_pvals.RDS')
}

if(file.exists('results/simpleGLMM_resids.RDS')){
  message("Loading simpleGLMM results")
  simpleGLMM_mles <- readRDS('results/simpleGLMM_mles.RDS')
  simpleGLMM_resids <- readRDS('results/simpleGLMM_resids.RDS')
  simpleGLMM_pvals <- readRDS('results/simpleGLMM_pvals.RDS')
}

if(file.exists('results/linmod_resids.RDS')){
  message("Loading linmod results")
  linmod_mles <- readRDS('results/linmod_mles.RDS')
  linmod_resids <- readRDS('results/linmod_resids.RDS')
  linmod_pvals <- readRDS('results/linmod_pvals.RDS')
}

