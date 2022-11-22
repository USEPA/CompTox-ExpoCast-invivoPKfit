#' Get model parameter names
#'
#' Get parameter names required by a specified model
#'
#'
#' @param model The name of the model. Currently only "flat", "1compartment",
#'   and "2compartment" are implemented.
#' @return A character vector of parameter names for the specified model.

get_model_paramnames <- function(model){
  # Get parameter names for the specified model
  if(model %in% "flat"){
    param_names <- "A"
  }else if(model %in% "1compartment"){
    #all possible 1-compartment params
    param_names <- c("kelim", "Vdist",
                     "Fgutabs", "kgutabs", "Fgutabs_Vdist")
  }else if(model %in% "2compartment"){
    #all possible 2-compartment params
    param_names <- c("kelim", "V1", "k21", "k12",
                     "kgutabs", "Fgutabs", "Fgutabs_V1")
  }else{
    stop(paste("get_model_paramnames(): Model",
                   model,
                   "is not recognized. Model should be one of",
                   paste("flat",
                   "1compartment",
                   "2compartment",
                   sep = "; ")))
  }

  return(param_names)
}
