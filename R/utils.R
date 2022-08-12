#' Setup C++ TMB code
#'
#' @param comp boolean that determines whether or not to compile the C++ file
#'
#' @return DLL is compiled and loaded
#' @export
setupTMB <- function(dll.name, comp=FALSE){
  if(!(paste0(dll.name, '.dll') %in% list.files('src')) |
     comp==TRUE){
    try(dyn.unload(dynlib(paste0('src/', dll.name))))
    TMB::compile(paste0('src/', dll.name,'.cpp'))
  }
  suppressWarnings(dyn.load(dynlib(paste0('src/', dll.name))))
}

#' Helper function to run TMB code and return output

#' @param obj.args List of arguments for TMB MakeADFun() function
#' @param opt.args List of arguments for nlminb() function
#' @param control List controlling model runs and standard error reporting
#'
#' @return Fitted objective function, nlminb output, reported values from model, sdreport if true
#'
#' @noRd
fit_tmb <- function(obj.args, opt.args = list(control = list(iter = 800, eval = 800), 
                                              scale = 1,
                                              lower = -Inf, upper = Inf ),
                    control = list(run.model = TRUE, do.sdreport = TRUE)){
  obj <- do.call(MakeADFun, obj.args)
  if(control$run.model){
    opt <- with(obj, do.call(nlminb,  c(list(par, fn, gr), opt.args) ))
    report <- obj$report(obj$env$last.par.best)
    aic <- add_aic(opt, n=length(obj$env$data$y))
    if(control$do.sdreport){
      sdr <- sdreport(obj)
      fit.results <- list(obj = obj, opt = opt, report = report, sdr = sdr, aic = aic)
    } else {
      fit.results <- list(obj = obj, opt = opt, report = report, aic = aic)
    }
  } else {
    fit.results <- list(obj = obj, report = obj$report())
  }
  return(fit.results)
}

## Quick fn to check for failed runs by looking at results output
## that doesn't exist
which.failed <- function(reps){
  success <- gsub('results/spatial_pvals/pvals_|.RDS', "", x=fs) %>%
    as.numeric()
  fail <- which(! reps %in% success)
  fail
}


add_aic <- function(opt,n){
  opt$AIC <- TMBhelper::TMBAIC(opt, n=Inf)
  opt$AICc <- TMBhelper::TMBAIC(opt, n=n)
  opt$BIC <- TMBhelper::TMBAIC(opt, p=log(n))
  opt
}


run_model <- function(reps, n=100, ng=0, mod, cov.mod = 'norm', 
                      misp, family = "Gaussian", link = "identity",
                      do.true = FALSE, savefiles=TRUE){

  ## Clean up the old runs
  unlink(paste0('results/',mod,'_resids', TRUE))
  unlink(paste0('results/',mod,'_pvals', TRUE))
  unlink(paste0('results/', mod,'_mles', TRUE))

  message(mod,": Preparing workspace to run ", length(reps), " iterations in parallel...")
  ## TMB::compile(paste0("src/",mod,".cpp")) # modified for simulation
  sfInit( parallel=cpus>1, cpus=cpus )
  ## sfExport('run.linmod.iter', 'sim.omega', 'cMatern', 'sim.data',
  ##          'rmvnorm_prec', 'add_aic')
  sfExportAll()

  message("Starting parallel runs...")
  results <- sfLapply(reps, function(ii) run_iter(ii, n, ng, mod, cov.mod, 
                                                  misp, family, link, 
                                                  do.true, savefiles))
  #results <- sapply(reps, function(ii) run_iter(ii, n, ng, mod, cov.mod, misp, do.true, savefiles))

  ## ## Read results back in from file
  ## fs <- list.files('results/linmod_pvals/', full.names=TRUE)
  ## ## Sometimes they fail to run for some unkonwn reason so try
  ## ## rerunning those ones once
  ## if(length(fs)<length(reps)){
  ##   message("Rerunning some failed runs...")
  ##   bad <- which.failed(reps)
  ##   results <- sfLapply(bad, function(ii) run.linmod.iter(ii))
  ##   fs <- list.files('results/linmod_pvals/', full.names=TRUE)
  ## }
  ## bad <- which.failed(reps)
  ## if(length(bad)>0) warning(length(bad), " runs failed")

  sfStop()
  message("Processing and saving final results...")
  if(do.true) mod <- paste0(mod, "_true")
  ## Read results back in from file
  fs <- list.files(paste0('results/',mod, '_', misp, '_pvals/'), full.names=TRUE)
  pvals <- lapply(fs, readRDS) %>% bind_rows %>% filter(!is.na(pvalue))
  saveRDS(pvals, file=paste0('results/',mod, '_', misp, '_pvals.RDS'))
  ## Read in residuals
  fs <- list.files(paste0('results/',mod, '_', misp, '_resids/'), full.names=TRUE)
  resids <- lapply(fs, readRDS) %>% bind_rows
  saveRDS(resids, file=paste0('results/',mod, '_', misp,'_resids.RDS'))
  fs <- list.files(paste0('results/',mod, '_', misp,'_mles/'), full.names=TRUE)
  mles <- lapply(fs, readRDS) %>% bind_rows
  saveRDS(mles, file=paste0('results/',mod, '_', misp,'_mles.RDS'))
  fs <- list.files(paste0('results/',mod, '_', misp,'_stats/'), full.names=TRUE)
  stats <- lapply(fs, readRDS) %>% bind_rows
  saveRDS(stats, file=paste0('results/',mod, '_', misp,'_stats.RDS'))
  return(results)
}


#' Extract marginal precision matrix for subset
#'
#' @param Q spHess object from TMB
#' @param i index of subset
#' @param ... additional arguments
#'
#' @return marginal Precision matrix
#' @export
#'
#' @examples
GMRFmarginal <- function(Q, i, ...) {
  ind <- 1:nrow(Q)
  i1 <- (ind)[i]
  i0 <- setdiff(ind, i1)
  if (length(i0) == 0)
    return(Q)
  Q0 <- as(Q[i0, i0, drop = FALSE], "symmetricMatrix")
  L0 <- Cholesky(Q0, ...)
  ans <- Q[i1, i1, drop = FALSE] - Q[i1, i0, drop = FALSE] %*%
    solve(Q0, Q[i0, i1, drop = FALSE])
  ans
}

