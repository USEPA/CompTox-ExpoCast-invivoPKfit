#' Get predictions
#'
#' Extract predictions from a fitted `pk` object.
#'
#' @param obj A [pk] object.
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions. If NULL (the default), then predictions will be made for the
#'   data in `obj$data`. `newdata` is required to contain at least the following
#'   variables: `Time`, `Dose`, `Route`, and `Media`. `Time` will be transformed
#'   according to the transformation in `obj$scales$time` before predictions are
#'   made.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions. If NULL (the default), predictions will be returned for
#'   all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions. If NULL (the default), predictions will be
#'   returned for all of the models in `obj$optimx_settings$method`.
#' @param type Either `"conc"` (the default) or `"auc"`. `type = "conc"`
#'   predicts concentrations; `type = "auc"` predicts area under the
#'   concentration-time curve (AUC).
#' @param exclude Logical: `TRUE` to return `NA_real_` for any observations in
#'   the data marked for exclusion (if there is a variable `exclude` in the
#'   data, an observation is marked for exclusion when `exclude %in% TRUE`).
#'   `FALSE` to return the prediction for each observation, regardless of
#'   exclusion. Default `TRUE`.
#' @return A named list of numeric matrixes. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, *i.e.*, each PK model that was
#'   fitted to the data. Each list element is a matrix with the same number of
#'   rows as the data in `obj$data` (corresponding to the rows in `obj$data`),
#'   and as many columns as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The column names are the method names. Each column
#'   contains the predictions of the model fitted by the corresponding method.
#'   These predictions are concentrations in the same units as
#'   `obj$data$Conc.Units`; any concentration transformations (in
#'   `obj$scale$conc`) are *not* applied.
#' @export
#' @author Caroline Ring
predict.pk <- function(obj,
                       newdata = NULL,
                       model = NULL,
                       method = NULL,
                       type = "conc",
                       exclude = TRUE,
                       ...
){

  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method

  coefs <- coef(obj = obj,
                model = model,
                method = method)

  if(is.null(newdata)) newdata <- obj$data

  #check that newdata has the required variables
  req_vars <- c("Dose",
                "Route",
                "Media")
  time_vars <- c("Time",
                 "Time_trans")
  if(!(all(req_vars %in% names(newdata))) |
     !(sum(time_vars %in% names(newdata)) > 0)){
    stop(paste("predict.pk(): newdata is missing one or more required variables:",
               "Dose, Route, Media, and either Time or Time_trans"))
  }

  if(!is.numeric(newdata$Dose) |
     !is.numeric(newdata$Time) |
     !is.character(newdata$Route)|
     !is.character(newdata$Media)){
    stop(paste("predict.pk(): One or more variables in newdata is not of the required type:",
               "Dose (should be numeric); Time (should be numeric);",
               "Route (should be character); Media (should be character)"))
  }

  if(!("Time_trans" %in% names(newdata))){
    #transform time if needed
    #first, default to identity transformation if none is specified
    if(is.null(obj$scales$time$new_units)){
      obj$scales$time$new_units <- "identity"
    }

    from_units <- unique(newdata$Time.Units)
    to_units <- ifelse(obj$scales$time$new_units %in% "identity",
                       from_units,
                       obj$scales$time$new_units)


    if(obj$scales$time$new_units %in% "auto"){
      to_units <- auto_units(y = newdata$Time,
                             from = from_units)
    }
    if(!obj$settings_preprocess$suppress.messages){
      message(paste("Converting time from",
                    from_units,
                    "to",
                    to_units))
    }
    newdata$Time_trans <- tryCatch(convert_time(x = newdata$Time,
                                                from = from_units,
                                                to = to_units,
                                                inverse = FALSE),
                                   error = function(err){
                                     warning(paste("invivopkfit::predict.pk():",
                                                   "Error in transforming time in `newdata` using convert_time():",
                                                   err$message))
                                     return(NA_real_)
                                   })
  }


  #loop over models
  sapply(model,
         function(this_model){
           this_coef_mat <- coefs[[this_model]]
           apply(this_coef_mat,
                 1,
                 function(this_coef_row){
                   #get coefficients
                   this_coef <- as.list(this_coef_row)
                   #get model function to be evaluated
                   this_model_fun <- ifelse(type %in% "conc",
                                            obj$stat_model[[this_model]]$conc_fun,
                                            ifelse(type %in% "auc",
                                                   obj$stat_model[[this_model]]$auc_fun,
                                                   NULL)
                   )

                   #evaluate model function
                   preds <- do.call(this_model_fun,
                                    args = list(params = this_coef,
                                                dose = newdata$Dose,
                                                time = newdata$Time_trans,
                                                route = newdata$Route,
                                                medium = newdata$Media))
                   if(exclude %in% TRUE){
                   #set NA for excluded data, if any
                   preds[obj$data$exclude %in% TRUE] <- NA_real_
                   }

                   preds

                 })
         },
         simplify = FALSE,
         USE.NAMES = TRUE)
}
