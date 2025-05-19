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


run_model <- function(reps, n=100, ng=0, mod, cov.mod = NULL,
                      misp, type, family = "Gaussian", link = "identity",
                      do.true = FALSE, savefiles=TRUE){
  if(do.true){
    mod.name <- paste0(mod, "_true")
  } else {
    mod.name <- mod
  }
  res.name <- paste0('results/', mod.name, '_', type)
  ## Clean up the old runs
  # unlink(paste0('results/', res.name))
  # unlink(paste0('results/', res.name))
  # unlink(paste0('results/', res.name))

  message(mod,": Preparing workspace to run ", length(reps), " iterations in parallel...")
  ## TMB::compile(paste0("src/",mod,".cpp")) # modified for simulation
  sfInit( parallel=cpus>1, cpus=cpus )
  ## sfExport('run.linmod.iter', 'sim.omega', 'cMatern', 'sim.data',
  ##          'rmvnorm_prec', 'add_aic')
  sfExportAll()
  sfLibrary(TMB)
  sfLibrary(DHARMa)
  sfLibrary(fmesher)
  sfLibrary(dplyr)
  sfLibrary(tidyr)
  sfLibrary(R.utils)
  sfLibrary(goftest)
  sfLibrary(tweedie)


  message("Starting parallel runs...")
  results <- sfLapply(reps, function(ii) run_iter(ii, n, ng, mod, cov.mod,
                                                  misp, type, family, link,
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
  ## Read results back in from file
  fs <- list.files(paste0(res.name, '_pvals/'), full.names=TRUE)
  pvals <- lapply(fs, readRDS) %>% bind_rows %>% filter(!is.na(pvalue))
  saveRDS(pvals, file=paste0(res.name, '_pvals.RDS'))
  ## Read in residuals
  fs <- list.files(paste0(res.name, '_resids/'), full.names=TRUE)
  resids <- lapply(fs, readRDS) %>% bind_rows
  saveRDS(resids, file=paste0(res.name,'_resids.RDS'))
  fs <- list.files(paste0(res.name,'_mles/'), full.names=TRUE)
  mles <- lapply(fs, readRDS) %>% bind_rows
  saveRDS(mles, file=paste0(res.name,'_mles.RDS'))
  fs <- list.files(paste0(res.name,'_stats/'), full.names=TRUE)
  stats <- lapply(fs, readRDS) %>% bind_rows
  saveRDS(stats, file=paste0(res.name,'_stats.RDS'))
  return(invisible(results))
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

### functions for running the speed testing code
## quick local functions to streamline code a bit
extract_runtime <- function(stats){
  select(stats, model, misp, replicate, version, starts_with('runtime'))
}
process_results <- function(mles, runtimes, pvals, type, 
                            model, misp, vary = NULL){
  if(!is.null(mles)){
    mles <- bind_rows(mles) %>% filter(h==0)} #filter(version=='h0')
  runtimes <- bind_rows(runtimes) %>% filter(version=='h0') %>%
    pivot_longer(starts_with('runtime'), names_to='method',
                 values_to='runtime') %>%
    mutate(type=gsub('runtime.|runtime_', '', method)) %>%
    group_by(nobs, model, version, method) %>%
    summarize(med=median(runtime, na.rm=TRUE),
              lwr=quantile(runtime, .25, na.rm=TRUE),
              upr=quantile(runtime, .75, na.rm=TRUE),
              pct.na=sum(is.na(runtime)),
              n=length(runtime), .groups='drop')
  if(!is.null(pvals)){
    if(length(pvals)!=0){
      pvals <- bind_rows(pvals) %>%
        filter(#version=='h0' &
          grepl('GOF', test)) %>%
        mutate(type=gsub('GOF.','',test))
    }
  }
  results <- list(mles=mles, pvals=pvals, runtimes=runtimes, 
                  model=model, misp = misp)

  if(model == 'simpleGLMM'){
    filename <- paste0('results/',model,'_', misp,'_', type, '_',
                       vary, '_sample_sizes.RDS')
  } else {
    filename <- paste0('results/',model,'_', misp, '_', type,
                       '_sample_sizes.RDS')
  }
  saveRDS(results, file = filename)
  return(results$runtimes)
}
plot_sample_sizes <- function(results){
  model <- results$model
  g <- ggplot(filter(results$pvals, !is.na(pvalue)), aes(pvalue, fill=type)) + facet_grid(nobs~method) +
    geom_histogram(position='identity', alpha=.5, bins=20)
  ggsave(paste0("plots/", model,"_pvals_by_dim.png"), g, width=8, height=5)
  g <- ggplot(results$runtimes,
              aes(nobs, med, ymin=lwr, ymax=upr,  color=type)) +
    geom_line()+
    geom_pointrange(fatten=2) + scale_y_log10()+labs(y='runtime (s)')
  ##  ggsave(paste0("plots/", model,"_runtimes_by_dim.png"), g, width=7, height=5)
  g <- ggplot(results$mles, aes(factor(nobs), mle-true)) +
    geom_violin() +
    geom_hline(yintercept=0, color='red') +
    facet_wrap('par') + labs(y='Absolute error')
  ggsave(paste0("plots/", model,"_mles_by_dim.png"),g, width=8, height=5)
}
get.value <- function(x, val, nobs){
  if(is.null(x)) return(NULL)
  new.df <- x[[val]]
  if(nrow(new.df) > 0){
    if(val!='runtimes') {
      return(data.frame(nobs=nobs, new.df))
    } else {
      return(data.frame(nobs=nobs, extract_runtime(x[['stats']])))
    }
  }
}
# get.value <- function(x, val, nobs){
#   if(is.null(x)) return(NULL)
#   if(val!='runtimes')
#     data.frame(nobs=nobs, x[[val]])
#   else
#     data.frame(nobs=nobs, extract_runtime(x[['stats']]))
# }

get.bad.reps <- function(df, Mod, Misp, doTrue){
  reps <- dplyr::filter(df, model == Mod & misp == Misp &
                  do.true == doTrue) %>%
    dplyr::pull(replicate) %>% unique() 
  idx <- which(!(1:1000 %in% reps ))
  return(idx)
}


#' 'squeeze' transform 
#'
#' @param u input
#'
#' @return [0,1] -> (0,1) to machine tolerance 
#'
squeeze <- function(u){
  eps <- .Machine$double.eps
  u = (1.0 - eps) * (u - .5) + .5
  return(u)
}
