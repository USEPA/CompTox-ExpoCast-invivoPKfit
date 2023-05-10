#' Get model function
#'
#' Get the function that predicts area under the concentration-time curve for a
#' given model name
#'
#' @param model The name of a model. One of 'flat', '1compartment', or
#'   '2compartment'.
#' @return The analytic AUC function corresponding to `model`: `auc_flat()` for
#'   'flat', `auc_1comp()` for '1compartment', `auc_2comp()` for '2compartment'.
#'   Returned as a function object, not as a character string.
#' @author Caroline Ring

get_model_auc_function <- function(model){
  if(model %in% "flat"){
    modelfun <- auc_flat
  }else if(model %in% "1compartment"){
    modelfun <- auc_1comp
  }else if(model %in% "2compartment"){
    modelfun <- auc_2comp
  }else{
    stop(paste("invivopkfit::get_model_auc_function():",
               "model is not one of 'flat', '1compartment' or '2compartment'.",
               "Only these three models are currently implemented."))
  }

  return(modelfun)
}
