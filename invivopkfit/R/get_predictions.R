#' Get model predicted concentration and AUC
#'
#' Given a table of fitted PK parameters and a table of input data (time, dose,
#' route, and medium), get predicted concentrations and AUCs.
#'
#' @param pk_fit A table of fitted PK parameters, as produced by
#'   `postprocess_data()`. Must contain variables `DTXSID`, `Species`,
#'   `Analysis_Type`, and `Studies.Analyzed`. Must also contain variables
#'   corresponding to the parameters of the model specified in argument
#'   `model_in`, as required for the route and media of data in `newdata` (as
#'   given by [get_model_paramnames()]). These variables must be named as
#'   `[param].[model]`. For example, for the 1-compartment model, if `newdata`
#'   contains only IV-dosing data measured in plasma, then `pk_fit` must contain
#'   variables named `kelim.1compartment` and `Vdist.1compartment`.
#' @param newdata A `data.table` which must contain variables named `Time`,
#'   `Dose`, `iv`, and `Media`.
#' @param DTXSID_in The DSSTox Substance ID for which to evaluate
#' @param Species_in The species for which to evaluate
#' @param Analysis_Type_in The analysis type that produced the fit being
#'   evaluated: one of 'Joint', 'Separate', or 'Pooled'
#' @param Studies.Analyzed_in The comma-separated string of studies included in
#'   the fit being analyzed.
#' @param model_in The model to evaluate: one of 'flat', '1compartment', or
#'   '2compartment'
#' @return A `data.table`, the same number of rows as `newdata`, with variables
#'   `pred_conc` and `pred_auc` corresponding to the model-predicted
#'   concentration and AUC.
#' @author Caroline Ring
#' @export
get_predictions <- function(pk_fit,
                            newdata,
                            DTXSID_in,
                            Species_in,
                            Analysis_Type_in,
                            Studies.Analyzed_in,
                            model_in){

  newdata <- copy(newdata)
if(!(grepl(x = model_in,
           pattern = "None"))){
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
  #get names of appropriate columns of pk_sub -- named [param].[model]
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
  #get model function
  cpfun <- get_model_function(model = model_in)
  #evaluate model function
  #if any error (e.g. no params for this model & dataset) return NAs
  preds <- tryCatch(do.call(cpfun,
                   args = list("params" = params,
                               "time" = newdata$Time,
                               "dose" = newdata$Dose,
                               "iv.dose" = newdata$iv,
                               "medium" = newdata$Media)),
                   error = function(err) rep(NA_real_, nrow(newdata)))

  # AUCS #
  #get model AUC function
  aucfun <- get_model_auc_function(model = model_in)
  #evaluate model AUC function at each time point
  aucs <- tryCatch(do.call(aucfun,
                  args = list("params" = params,
                              "time" = newdata$Time,
                              "dose" = newdata$Dose,
                              "iv.dose" = newdata$iv,
                              "medium" = newdata$Media)),

                  error = function(err) rep(NA_real_, nrow(newdata)))
  newdata[, pred_conc := preds]
  newdata[, pred_auc := aucs]
}else{
  newdata[, pred_conc := NA_real_]
  newdata[, pred_auc := NA_real_]
}

  return(newdata)
}
