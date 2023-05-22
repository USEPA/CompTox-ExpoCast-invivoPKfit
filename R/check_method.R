#' Check methods
#'
#' Check methods for validity
#'
#' Helper function to ensure that a list of methods specified by the user matches the methods available in the fitted [pk()] object.
#'
#' @param obj A [pk()] object
#' @param method A user-supplied `character` vector of method names
#' @return `TRUE` if all `method %in% obj$settings_optimx$method`; otherwise stops with an error
#' @author Caroline Ring
check_method <- function(obj, method){
  #check that all methods are valid
  if(!(all(method %in% obj$settings_optimx$method))){
    stop(paste("All values in `method` must be found in `obj$settings_optimx$method.",
               paste0("`method` = ", paste(method, sep = ", ")),
               paste0("`obj$settings_optimx$method` = ", paste(obj$settings_optimx$method, collapse = ", ")),
               sep = "\n"))

  }else{
    return(TRUE)
  }
}

#' Check models
#'
#' Check models for validity
#'
#' Helper function to ensure that a list of models specified by the user matches the models available in the fitted [pk()] object.
#'
#' @param obj A [pk()] object
#' @param model A user-supplied `character` vector of model names
#' @return `TRUE` if all `model %in% names(obj$stat_model)`; otherwise stops with an error
#' @author Caroline Ring
check_model <- function(obj, model){
  #check that all models are valid
  if(!(all(model %in% names(obj$stat_model)))){
    stop(paste("All values in `model` must be found in `names(obj$stat_model).",
               paste0("`model` = ", paste(model, sep = ", ")),
               paste0("`names(obj$stat_model)` = ",
                      paste(names(obj$stat_model), collapse = ", ")
               ),
               sep = "\n"))

  }else{
    return(TRUE)
  }
}
