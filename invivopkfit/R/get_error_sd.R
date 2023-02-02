#'Estimate pooled error standard deviation with fixed parameters
#'
#'Given a table of model parameters (e.g. fitted using [fit_all()]), estimate
#'the residual error standard deviation assuming residuals are iid following a
#'zero-mean normal distribution
#'
#'@param pk_fit_row One row from a table of fitted PK parameters, as produced by
#'  [merge_fits()]. Must contain variables `DTXSID`, `Species`,
#'  `Analysis_Type`, and `Studies.Analyzed`. Must also contain variables
#'  corresponding to the parameters of the model specified in argument
#'  `model`, as required for the route and media of data in `newdata` (as
#'  given by [get_model_paramnames()]). These variables must be named as
#'  `[param].[model]`. For example, for the 1-compartment model, if `newdata`
#'  contains only IV-dosing data measured in plasma, then `pk_fit` must contain
#'  variables named `kelim.1compartment` and `Vdist.1compartment`.
#'@param newdata A `data.table` which must contain variables named `Time`,
#'  `Dose`, `iv`, and `Media`.
#'@param model The model to evaluate: one of 'flat', '1compartment', or
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
get_error_sd <- function(pk_fit_row,
                            newdata,
                            model,
                         modelfun = "analytic",
                         log_trans = TRUE,
                         dose_norm = FALSE,
                         optimx_args = list(
                           "method" = "bobyqa",
                           "itnmax" = 1e6,
                           "control" = list("kkt" = FALSE)
                         )){
#Make a copy of newdata so it behaves as though passed by value, not by reference
  newdata <- copy(newdata)

  #To get starting value: Get SD of residuals

  #First get predictions
  newdata[, pred := predict_conc(pk_fit = pk_fit_row,
                                  newdata = newdata,
                                  model = model)]
#Then calc residuals under conditions log_trans and dose_norm
  if(log_trans %in% FALSE){
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
  if(dose_norm %in% TRUE){
  sd_start <- sd(newdata_pred$resid_Dose, na.rm = TRUE)
  }else{
    sd_start <- sd(newdata_pred$resid, na.rm = TRUE)
  }

  #Handle the case when PO data was dropped because there were not enough data in
  #absorption phase: in this case, "kgutabs" will NOT be in the fitted params even
  #though PO data exist.
  params <- get_fitted_params(pk_fit = pk_fit_row,
                              model = model)
  if(any(newdata$Route %in% "po") &
     !any(names(params) %in% "kgutabs") &
     !(model %in% "flat")){
    if(suppress.messages %in% FALSE){
      message("newdata contains oral data, but fit was IV-only, likely because of insufficient absorption-phase data. Optimizing pooled error SD for IV data only.")
    }
    #keep only IV data
    newdata <- newdata[Route %in% "iv"]
  }

  #Hold model params constant & optimize only a pooled sigma
  sigma_fit <- tryCatch(do.call(optimx::optimx,
                       args = c(
                         #
                         list(par = c("sigma" = sd_start),
                              fn = log_likelihood,
                              lower = c("sigma" = 1e-8),
                              upper = c("sigma" = 1e8)),
                         #method and control
                         optimx_args,
                         #other args to log_likelihood()
                         list(
                           const_params = params, #fitted model params
                           DF = newdata,
                           modelfun = modelfun,
                           model = model,
                           fit_conc_dose = dose_norm,
                           fit_log_conc = log_trans,
                           force_finite = (optimx_args$method %in% "L-BFGS-B"),
                           negative = TRUE
                         )
                       )
  ),
  error = function(err) {
    return(err$message)
  })

  #if fit was successful, get fitted sigma value
  if(inherits(sigma_fit, "optimx")){
  sigma_est <- stats::coef(sigma_fit)

  out_DF <- data.frame("sigma" = sigma_est)
  #Keep information about optimization
  #number of function evals
  out_DF$fevals <- sigma_fit$fevals
  #convergence code
  out_DF$convcode <- sigma_fit$convcode
  #number of iterations
  out_DF$niter <- sigma_fit$niter
  }else{
    out_DF <- data.frame("sigma" = NA_real_)
    #Keep information about optimization
    #number of function evals
    out_DF$fevals <- NA_real_
    #convergence code
    out_DF$convcode <-  NA_real_
    #number of iterations
    out_DF$niter <-  NA_real_
  }

  #optimziation method used
  out_DF$method <- optimx_args$method
  #optimx control parameters
  #record control params
  out_DF[paste0("control_",
                names(optimx_args$control))] <- optimx_args$control

  return(as.data.table(out_DF))
}
