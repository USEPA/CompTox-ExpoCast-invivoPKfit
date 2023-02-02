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
#' @param model The model to evaluate: one of 'flat', '1compartment', or
#'   '2compartment'
#' @return A `data.table`, the same number of rows as `newdata`, with variables
#'   `pred_conc` and `pred_auc` corresponding to the model-predicted
#'   concentration and AUC.
#' @author Caroline Ring
#' @export
predict_conc <- function(pk_fit_row,
                            newdata,
                            model){
  #get fitted parameters
params <- get_fitted_params(pk_fit = pk_fit_row,
                            model_in = model)

#Handle the case when PO data was dropped because there were not enough data in
#absorption phase: in this case, "kgutabs" will NOT be in the fitted params even
#though PO data exist.
newdata[, get_pred := TRUE]
if(any(newdata$Route %in% "po") &
   !any(names(params) %in% "kgutabs") &
   !(model %in% "flat")){
  if(suppress.messages %in% FALSE){
    message("newdata contains oral data, but fit was IV-only, likely because of insufficient absorption-phase data. Returning predictions for IV data only, with NAs inserted for oral data.")
  }
  newdata[Route %in% "po", get_pred := FALSE]
}

  #get model function
  cpfun <- get_model_function(model = model)

  #get predictions only for rows of newdata where get_pred == TRUE
  #this will insert NAs where get_pred == FALSE
 newdata[get_pred %in% TRUE,
                   preds := tryCatch(do.call(cpfun,
                                    args = list("params" = params,
                                                "time" = Time,
                                                "dose" = Dose,
                                                "iv.dose" = iv,
                                                "medium" = Media)),
                            error = function(err){
                              rep(NA_real_, .N)
                            } )
                   ]

  preds <- newdata$pred

  return(preds)
}
