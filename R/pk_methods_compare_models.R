#' Model comparison for [pk()] objects
#'
#' Perform model comparison for a fitted [pk()] object.
#'
#' Models are compared according to the goodness-of-fit criterion named in
#' "criterion", and the name of the winning model is returned.
#'
#' @param obj A [pk()] model object. Must be fitted, or the function will exit
#'   with an error.
#' @param newdata Optional: A `data.frame` containing new data for which to
#'   compute the TK stats. Must contain at least variables `Chemical`,
#'   `Species`, `Route`, `Media`, `Dose`, and any other variables named in
#'   `tk_grouping`. Default `NULL`, to use the data in `obj$data`.
#' @param model Character: One or more of the models fitted. Default `NULL` to
#'   return TK stats for all models.
#' @param method Character: One or more of the [optimx::optimx()] methods used.
#'   Default `NULL` to return TK stats for all methods.
#' @param criterion The name of a criterion function to use for model
#'   comparison. Default "AIC". Must be the name of a function that (as for
#'   `AIC`) accepts arguments `obj`, `newdata`, `method` and `model` (may accept
#'   other arguments, specified in `...`) and returns output as for `AIC`: a
#'   named list of numeric vectors (named for each of the model names in
#'   `model`), where each vector has elements named for each of the method names
#'   in `method`, containing the criterion value calculated for that model
#'   fitted using that method.
#' @param ... Optional: Other arguments to `criterion` function.
#' @return A `data.frame` with variables
#' - `model`: The name of each model
#' - `method`: The name of each method
#' - A variable named for `criterion` (e.g. if `criterion = "AIC"` then the result will have a
#'   variable named `AIC`): The criterion value for each model/method
#' @export
#' @author Caroline Ring
compare_models.pk <- function(obj,
                             newdata = NULL,
                             model = NULL,
                             method = NULL,
                             criterion = "AIC",
                             ...){
  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$settings_optimx$method
  if(is.null(newdata)) newdata <- obj$data

  #check that all methods are valid
  if(!(all(method %in% obj$settings_optimx$method))){
    stop(paste("All values in `method` must be found in `obj$settings_optimx$method.",
               paste0("`method` = ", paste(method, sep = ", ")),
               paste0("`obj$settings_optimx$method` = ", paste(obj$settings_optimx$method)),
               sep = "\n"))
  }

  crit_list <- do.call(criterion,
          args = c(list(obj = obj),
                      list(newdata = newdata,
                      model = model,
                      method = method),
                   list(...)))

  critDF_list <- sapply(model, #loop over models
         function(this_model){
           #convert criterion values to a data.frame with variables model, method, and (criterion)
          foo <- data.frame(model = this_model,
                     method = names(crit_list[[this_model]]),
                     criterion = crit_list[[this_model]]
                     )
          names(foo) <- c("model", "method", criterion)
          foo

         },
         simplify = FALSE)

  do.call(rbind, critDF_list)
}
