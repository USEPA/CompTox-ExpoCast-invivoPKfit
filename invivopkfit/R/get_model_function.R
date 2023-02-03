#' Get model function
#'
#' Get the function that predicts concentration for a given model name
#'
#' @param model The name of a model. One of 'flat', '1compartment', or
#'   '2compartment'.
#' @param model.type "analytic" or "full". Default "analytic" to use the
#'   analytic model solution. "full" to solve the ODE numerically.
#' @return Character: The name of the function corresponding to `model`:
#'   `cp_flat()` for 'flat', `cp_1comp()` for '1compartment', `cp_2comp()` for
#'   '2compartment'.
#' @author Caroline Ring

get_model_function <- function(model, model.type){
  if(model %in% "flat"){
    modelfun <- "cp_flat"
  }else if(model %in% "1compartment"){
    if(model.type %in% "analytic"){
    modelfun <- "cp_1comp"
    }else{
      modelfun <- "httk::solve_1comp"
    }
  }else if(model %in% "2compartment"){
    if(model.type %in% "analytic"){
    modelfun <- "cp_2comp"
    }else{
      modelfun <- stop(paste("invivopkfit::get_model_function():",
                             "model is '2compartment' but model.type is not 'analytic'.",
                             "Currently, only the analytic solution is implemented for the 2-compartment model."))
    }
  }else{
    stop(paste("invivopkfit::get_model_function():",
               "model is not one of 'flat', '1compartment' or '2compartment'.",
               "Only these three models are currently implemented."))
  }

  return(modelfun)
}
