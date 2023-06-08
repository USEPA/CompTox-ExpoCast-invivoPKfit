#'Get winning model
#'
#'Get winning model for a fitted `pk` object
#'
#'Get the winning model (i.e. the model with the lowest value of the criterion
#'specified in `criterion`) for a fitted `pk` object, for a specified method,
#'and optionally for a specified new dataset.
#'
#'@param obj A [pk()] object
#'@param newdata Optional: A `data.frame` containing new data to plot. Must
#'  contain at least variables `Chemical`, `Species`, `Route`, `Media`, `Dose`,
#'  `Time`, `Time.Units`, `Conc`, `Detect`, `Conc_SD`.  Default `NULL`, to use
#'  the data in `obj$data`.
#'@param method Character: One or more of the [optimx::optimx()] methods used in
#'  fitting. The winning model will be determined for each of these methods.
#'  Default `NULL` to get the winning model for each method in
#'  `obj$settings_optimx$method`.
#'@param criterion The name of a criterion function to use for model comparison.
#'  Default "AIC". Must be the name of a function that (as for `AIC`) accepts
#'  arguments `obj`, `newdata`, `method` and `model` (may accept other
#'  arguments, specified in `...`) and returns output as for `AIC`: a named list
#'  of numeric vectors (named for each of the model names in `model`), where
#'  each vector has elements named for each of the method names in `method`,
#'  containing the criterion value calculated for that model fitted using that
#'  method.
#'@param ... Optional: Other arguments to `criterion` function.
#'@return A named character array with one element for each item in `method`,
#'  giving the name of the winning model. The return value has attribute
#'  `criterion` giving the name of the criterion function used to compare
#'  models. (If the return value is named `x`, you can access this attribute
#'  using `attr(x, 'criterion')`.
#'@export
#' @author Caroline Ring
get_winning_model.pk <- function(obj,
                                 newdata = NULL,
                                 method = NULL,
                                 criterion = "AIC",
                                 ...){
  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(method)) method <- obj$settings_optimx$method
  if(is.null(newdata)) newdata <- obj$data

  #check that all methods are valid
  if(!(all(method %in% obj$settings_optimx$method))){
    stop(paste("All values in `method` must be found in `obj$settings_optimx$method.",
               paste0("`method` = ", paste(method, sep = ", ")),
               paste0("`obj$settings_optimx$method` = ", paste(obj$settings_optimx$method)),
               sep = "\n"))
  }

  model_compare <- compare_models(obj = obj,
                                  newdata = newdata,
                                  model = NULL,
                                  method = method,
                                  criterion = criterion,
                                  ...)

 #return the winning model for each method
  mc_sort <- model_compare[order(model_compare$method,
                      model_compare[[criterion]],
                      decreasing = FALSE,
                      na.last =TRUE),]

  winmodels <- tapply(mc_sort$model,
         mc_sort$method,
         FUN = function(x) x[1])


  attr(winmodels, "criterion") <- criterion

  return(winmodels)

}
