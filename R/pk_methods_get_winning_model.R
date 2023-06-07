get_winning_model.pk <- function(obj,
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

  model_compare <- compare_models(obj = obj,
                                  newdata = newdata,
                                  model = model,
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
