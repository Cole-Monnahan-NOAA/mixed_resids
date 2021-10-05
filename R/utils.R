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
                                              hessian = NULL, scale = 1,
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


run_model <- function(reps, n=100, ng=0, mod, cov.mod = 'norm', misp, do.true = FALSE, savefiles=TRUE){

  ## Clean up the old runs
  unlink(paste0('results/',mod,'_resids', TRUE))
  unlink(paste0('results/',mod,'_pvals', TRUE))
  unlink(paste0('results/', mod,'_mles', TRUE))

  message("Preparing workspace to run ", length(reps), " iterations in parallel...")
  ## TMB::compile(paste0("src/",mod,".cpp")) # modified for simulation
  sfInit( parallel=TRUE, cpus=cpus )
  ## sfExport('run.linmod.iter', 'sim.omega', 'cMatern', 'sim.data',
  ##          'rmvnorm_prec', 'add_aic')
  sfExportAll()

  message("Starting parallel runs...")
  results <- sfLapply(reps, function(ii) run_iter(ii, n, ng, mod, cov.mod, misp, do.true, savefiles))
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

  message("Processing and saving final results...")
  ## Read results back in from file
  fs <- list.files(paste0('results/',mod,'_pvals/'), full.names=TRUE)
  pvals <- lapply(fs, readRDS) %>% bind_rows %>% filter(!is.na(pvalue))
  saveRDS(pvals, file=paste0('results/',mod,'_pvals.RDS'))
  ## Read in residuals
  fs <- list.files(paste0('results/',mod,'_resids/'), full.names=TRUE)
  resids <- lapply(fs, readRDS) %>% bind_rows
  saveRDS(resids, file=paste0('results/',mod,'_resids.RDS'))
  fs <- list.files(paste0('results/',mod,'_mles/'), full.names=TRUE)
  ##! mles need to be standardize, so far they are saved as list
  #mles <- lapply(fs, readRDS) %>% bind_rows
  #saveRDS(mles, file=paste0('results/',mod,'_mles.RDS'))

}
