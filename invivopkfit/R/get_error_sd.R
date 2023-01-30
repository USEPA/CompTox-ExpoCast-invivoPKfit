#'Estimate pooled error standard deviation with fixed parameters
#'
#'Given a table of model parameters (e.g. fitted using [fit_all()]), estimate
#'the residual error standard deviation assuming residuals are iid following a
#'zero-mean normal distribution
#'
#'@param pk_fit A table of fitted PK parameters, as produced by
#'  `postprocess_data()`. Must contain variables `DTXSID`, `Species`,
#'  `Analysis_Type`, and `Studies.Analyzed`. Must also contain variables
#'  corresponding to the parameters of the model specified in argument
#'  `model_in`, as required for the route and media of data in `newdata` (as
#'  given by [get_model_paramnames()]). These variables must be named as
#'  `[param].[model]`. For example, for the 1-compartment model, if `newdata`
#'  contains only IV-dosing data measured in plasma, then `pk_fit` must contain
#'  variables named `kelim.1compartment` and `Vdist.1compartment`.
#'@param newdata A `data.table` which must contain variables named `Time`,
#'  `Dose`, `iv`, and `Media`.
#'@param DTXSID_in The DSSTox Substance ID for which to evaluate
#'@param Species_in The species for which to evaluate
#'@param Analysis_Type_in The analysis type that produced the fit being
#'  evaluated: one of 'Joint', 'Separate', or 'Pooled'
#'@param Studies.Analyzed_in The comma-separated string of studies included in
#'  the fit being analyzed.
#'@param model_in The model to evaluate: one of 'flat', '1compartment', or
#'  '2compartment'
#'@param optimx_args A named list of additional arguments to [optimx::optimx()],
#'  other than `par`, `fn`, `lower`, and `upper`. Default is:
#'
#'    ```
#'     list(
#'           "method" = "bobyqa",
#'           "itnmax" = 1e6,
#'          "control" = list("kkt" = FALSE))
#'    ```
#'  Briefly:  `"method"` allows you to select an optimization algorithm.
#'  `"itnmax"` controls the maximum number of iterations allowed in an attempt
#'  to optimize. `"control"` is itself a named list of control parameters
#'  relevant to the selected method. See documentation for [optimx::optimx()]
#'  for more details and more options. Note lower and upper bounds (box
#'  constraints) will be supplied; if you want them to be respected, please
#'  choose a method that allows box constraints (e.g. "bobyqa" or "L-BFGS-B").
#'
#'@return A one-row `data.table` with a variable `sigma` reporting the estimated
#'  residual error standard deviation, and other variables: `fevals` reporting
#'  the number of fucntion evaluations performed as returned by
#'  [optimx::optimx()]; `convcode` reporting the integer convergence code
#'  returned by [optimx::optimx()]; `niter`reporting the number of iterations
#'  performed by [optimx::optimx()]; `method` reporting the optimization method
#'  used (as specified in input argument `optimx_args`); and variables reporting
#'  each control parameter specified in `optimx_args`, prefixed with `control_`.
#'  (For example, using the default value of `optimx_args`, there will be a
#'  variable named `control_kkt` whose value will be FALSE.)
#' @author Caroline Ring
get_error_sd <- function(pk_fit,
                            newdata,
                            DTXSID_in,
                            Species_in,
                            Analysis_Type_in,
                         Studies.Analyzed_in,
                            model_in,
                         modelfun = "analytic",
                         fit_conc_dose = TRUE,
                         fit_log_conc = FALSE,
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

  if(fit_log_conc %in% FALSE){
  #get residuals
  newdata_pred[Dose > 0 & Detect %in% "Detect",
               resid := Value - pred_conc]
  #get dose-normalized residuals
  newdata_pred[Dose > 0 & Detect %in% "Detect",
               resid_Dose := Value_Dose - (pred_conc/Dose)]
  }else{
    #get residuals
    newdata_pred[Dose > 0 & Detect %in% "Detect",
                 resid := log(Value) - log(pred_conc)]
    #get dose-normalized residuals
    newdata_pred[Dose > 0 & Detect %in% "Detect",
                 resid_Dose := log(Value_Dose) - log(pred_conc/Dose)]
  }

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
                           fit_log_conc = fit_log_conc,
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

  out_DF <- data.frame("sigma" = sigma_est)
  #Keep information about optimization
  #number of function evals
  out_DF$fevals <- as.integer(sigma_fit$fevals)
  #convergence code
  out_DF$convcode <- as.integer(sigma_fit$convcode)
  #number of iterations
  out_DF$niter <- as.integer(sigma_fit$niter)
  }else{
    sigma_est <- NA_real_
    out_DF <- data.frame("sigma" = sigma_est)
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
