#' Get model function
#'
#' Get the function that predicts concentration for a given model name
#'
#' @param model The name of a model. One of 'flat', '1compartment', or
#'   '2compartment'.
#' @return The function corresponding to `model`: `cp_flat()` for 'flat',
#'   `cp_1comp()` for '1compartment', `cp_2comp()` for '2compartment'. Returned
#'   as a function object, not as a character string.
#' @author Caroline Ring

get_model_function <- function(model){
  if(model %in% "flat"){
    modelfun <- cp_flat
  }else if(model %in% "1compartment"){
    modelfun <- cp_1comp
  }else if(model %in% "2compartment"){
    modelfun <- cp_2comp
  }else{
    stop(paste("invivopkfit::get_model_function():",
               "model is not one of 'flat', '1compartment' or '2compartment'.",
               "Only these three models are currently implemented."))
  }

  return(modelfun)
}
