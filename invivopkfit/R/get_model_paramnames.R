#' Get model parameter names
#'
#' Get parameter names required by a specified model
#'
#' @param model The name of the model

get_model_paramnames <- function(model){
  # Get parameter names for the specified model
  if(model %in% "flat"){
    param_names <- "A"
  }else if(model %in% "1compartment"){
    param_names <- c("kelim", "Vdist", "Fgutabs", "kgutabs")
  }else if(model %in% "2compartment"){
    param_names <- c("kelim", "V1", "Ralphatokelim", "Fbetaofalpha")
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
