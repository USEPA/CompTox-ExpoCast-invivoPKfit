get_error_sd <- function(pk_fit,
                            newdata,
                            DTXSID_in,
                            Species_in,
                            Analysis_Type_in,
                         Studies.Analyzed_in,
                            model_in,
                         modelfun = "analytic",
                         fit_conc_dose = TRUE,
                         optimx_args = list(
                           "method" = "bobyqa",
                           "itnmax" = 1e6,
                           "control" = list("kkt" = FALSE)
                         )){

  newdata <- copy(newdata)

  #get appropriate subset of pk_fit
  pk_sub <- unique(pk_fit[DTXSID %in% DTXSID_in &
                            Species %in% Species_in &
                            Analysis_Type %in% Analysis_Type_in &
                            Studies.Analyzed %in% Studies.Analyzed_in, ])

  #this should be a one-row subset
  if(nrow(pk_sub)>1) {
    stop("Selected subset of pk_fit has more than one unique row!")
  }

  #get parameter names for this model
  param_names <- get_model_paramnames(model = model_in)
  #get appropriate columns of pk_sub -- named [param].[model]
  param_names_dot <- paste(param_names, model_in, sep = ".")
  #extract appropriate columns of pk_sub
  params <- pk_sub[, .SD, .SDcols = intersect(param_names_dot,
                                              names(pk_sub))]
  #rename to strip the ".[model]" part
  setnames(params,
           param_names_dot,
           param_names,
           skip_absent = TRUE)
  #set Rblood2plasma to 1, if it is there and NA
  if("Rblood2plasma" %in% names(params)){
    if(is.na(params$Rblood2plasma)){
      params$Rblood2plasma <- 1
    }
  }
  #convert to list
  params <- as.list(params)
  #remove any NA params
  params <- params[!(sapply(params,
                            is.na))]



  #To get starting value for params: First get model-predicted concentrations
  newdata_pred <- get_predictions(pk_fit = pk_sub,
                                  newdata = newdata,
                                  DTXSID_in = DTXSID_in,
                                  Species_in = Species_in,
                                  Analysis_Type_in = Analysis_Type_in,
                                  Studies.Analyzed_in = Studies.Analyzed_in,
                                  model_in = model_in)
  #get residuals
  newdata_pred[Dose > 0 & Detect %in% "Detect",
               resid := Value - pred_conc]
  #get dose-normalized residuals
  newdata_pred[Dose > 0 & Detect %in% "Detect",
               resid_Dose := Value_Dose - (pred_conc/Dose)]

  #get SD of residuals
  if(fit_conc_dose %in% TRUE){
  sd_start <- sd(newdata_pred$resid_Dose, na.rm = TRUE)
  }else{
    sd_start <- sd(newdata_pred$resid, na.rm = TRUE)
  }


  #optimize for sigma only

  #Hold model params constant & optimize only error SD
  sigma_fit <- tryCatch(do.call(optimx::optimx,
                       args = c(
                         #
                         list(par = c("sigma" = sd_start),
                              fn = log_likelihood,
                              lower = c("sigma" = 1e-8),
                              upper = c("sigma" = Inf)),
                         #method and control
                         optimx_args,
                         list(
                           const_params = params,
                           DF = newdata,
                           modelfun = modelfun,
                           model = model_in,
                           fit_conc_dose = fit_conc_dose,
                           force_finite = (optimx_args$method %in% "L-BFGS-B"),
                           negative = TRUE
                         )
                       )
  ),
  error = function(err) {
    return(err$message)
  })

  if(inherits(sigma_fit, "optimx")){
  sigma_est <- stats::coef(sigma_fit)

  out_DF <- data.frame("error_SD" = sigma_est)
  #Keep information about optimization
  #number of function evals
  out_DF$fevals <- as.integer(sigma_fit$fevals)
  #convergence code
  out_DF$convcode <- as.integer(sigma_fit$convcode)
  #number of iterations
  out_DF$niter <- as.integer(sigma_fit$niter)
  }else{
    sigma_est <- NA_real_
    out_DF <- data.frame("error_SD" = sigma_est)
    #Keep information about optimization
    #number of function evals
    out_DF$fevals <- NA_integer_
    #convergence code
    out_DF$convcode <-  NA_integer_
    #number of iterations
    out_DF$niter <-  NA_integer_
  }

  #optimziation method used
  out_DF$method <- optimx_args$method
  #optimx control parameters
  #record control params
  out_DF[paste0("control_",
                names(optimx_args$control))] <- optimx_args$control

  return(as.data.table(out_DF))
}
