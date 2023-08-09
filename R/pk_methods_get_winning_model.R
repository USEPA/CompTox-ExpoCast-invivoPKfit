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
#'  arguments, specified in `...`) and returns output as for `AIC`: a data.frame
#'  with a column with the same name as `criterion` that has calculated values
#'  for model comparison. The "winning" value will be the smallest value.
#'@param ... Optional: Other arguments to `criterion` function.
#'@return A data.frame with one row for each `data_group`, `model` and `method` and
#'  The return value has attribute `criterion` giving the name of the criterion function used to compare
#'  models.
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
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  if (is.null(method)) method <- obj$settings_optimx$method
  # AIC and logLik will take newdata = NULL


  #check that all methods are valid
  if (!(all(method %in% obj$settings_optimx$method))) {
    stop(paste("All values in `method` must be found in `obj$settings_optimx$method.",
               paste0("`method` = ", paste(method, sep = ", ")),
               paste0("`obj$settings_optimx$method` = ", paste(obj$settings_optimx$method)),
               sep = "\n"))
  }

  # N.B. this updated method makes compare_models() superfluous

  model_compare <- do.call(criterion,
                           args = list(obj = obj,
                                       newdata = newdata,
                                       method = method))

 #return the winning model for each method
 winmodels <- model_compare %>%
   dplyr::group_by(!!!obj$data_group, method) %>%
   dplyr::filter(if_any(contains(criterion), ~ . == min(.)),
          method == method) %>%
    dplyr::ungroup() %>%
    dplyr::select(!!!obj$data_group,
                  method, model)

  return(winmodels)

}
