#' Get predictions
#'
#' Extract predictions from a fitted `pk` object.
#'
#' @param obj A [pk] object.
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions. If NULL (the default), then predictions will be made for the
#'   data in `obj$data`. `newdata` is required to contain at least the following
#'   variables: `Time`, `Time.Units`, `Dose`, `Route`, and `Media`. If variable
#'   `Time_trans` is not present, then `Time` will be transformed according to
#'   the transformation in `obj$scales$time` before making predictions;
#'   otherwise, `Time_trans` will be used to make predictions.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions. If NULL (the default), predictions will be returned for
#'   all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions. If NULL (the default), predictions will be
#'   returned for all of the models in `obj$settings_optimx$method`.
#' @param type Either `"conc"` (the default) or `"auc"`. `type = "conc"`
#'   predicts concentrations; `type = "auc"` predicts area under the
#'   concentration-time curve (AUC).
#' @param exclude Logical: `TRUE` to return `NA_real_` for any observations in
#'   the data marked for exclusion (if there is a variable `exclude` in the
#'   data, an observation is marked for exclusion when `exclude %in% TRUE`).
#'   `FALSE` to return the prediction for each observation, regardless of
#'   exclusion. Default `TRUE`.
#' @param use_scale_conc Possible values: `TRUE`, `FALSE`, or a named list with
#'   elements `dose_norm` and `log10_trans` which themselves should be either
#'   `TRUE` or `FALSE`. If `use_scale_conc = TRUE`, then the concentration
#'   scaling/transformations in `obj` will be applied to both predicted and
#'   observed concentrations before the log-likelihood is computed. If
#'   `use_scale_conc = FALSE` (the default for this function), then no
#'   concentration scaling or transformation will be applied before the
#'   log-likelihood is computed. If `use_scale_conc = list(dose_norm = ...,
#'   log10_trans = ...)`, then the specified dose normalization and/or
#'   log10-transformation will be applied.
#' @return A named list of numeric matrixes. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, *i.e.*, each PK model that
#'   was fitted to the data. Each list element is a matrix with the same number
#'   of rows as the data in `obj$data` (corresponding to the rows in
#'   `obj$data`), and as many columns as there were [optimx::optimx()] methods
#'   (specified in [settings_optimx()]). The column names are the method names.
#'   Each column contains the predictions of the model fitted by the
#'   corresponding method. If `use_scale_conc %in% FALSE`, these predictions are
#'   un-transformed concentrations in the same units as `obj$data$Conc.Units`.
#'   If `use_scale_conc %in% TRUE`, the predictions are transformed
#'   concentrations in the same units as `obj$data$Conc_trans.Units`.
#' @export
#' @author Caroline Ring
#' @family methods for fitted pk objects
predict.pk <- function(obj,
                       newdata = NULL,
                       model = NULL,
                       method = NULL,
                       type = "conc",
                       exclude = TRUE,
                       use_scale_conc = FALSE,
                       ...
){

  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$settings_optimx$method

  method_ok <- check_method(obj = obj, method = method)
  model_ok <- check_model(obj = obj, model = model)

  coefs <- coef(obj = obj,
                model = model,
                method = method)

  if(is.null(newdata)) newdata <- obj$data

  newdata_ok <- check_newdata(newdata = newdata,
                              olddata = obj$data,
                              req_vars = c("Time",
                                  "Time.Units",
                                           "Dose",
                                           "Route",
                                           "Media"),
                              exclude = exclude)

  #scale time if needed
  if(!("Time_trans" %in% names(newdata))){
  newdata$Time_trans <- convert_time(x = newdata$Time,
                                     from = newdata$Time.Units,
                                     to = obj$scales$time$new_units)
  }

  #apply transformations if so specified
  conc_scale <- conc_scale_use(obj = obj,
                               use_scale_conc = use_scale_conc)


  #loop over models
  sapply(model,
         function(this_model){
           this_coef_mat <- coefs[[this_model]][method, ]
           if(!is.matrix(this_coef_mat)){
             #in case there is only 1 method and the matrix is therefore 1-row and gets converetd into a vector,
             #convert it back
             this_coef_mat <- matrix(this_coef_mat,
                                     nrow = 1,
                                     ncol = length(this_coef_mat))
             colnames(this_coef_mat) <- colnames(coefs[[this_model]])
             rownames(this_coef_mat) <- method
           }
           apply(this_coef_mat,
                 1,
                 function(this_coef_row){
                   #get coefficients
                   this_coef <- this_coef_row
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
                     if("exclude" %in% names(newdata)){
                   #set NA for excluded data, if any
                   preds[newdata$exclude %in% TRUE] <- NA_real_
                     }
                   }

                  if(conc_scale$dose_norm %in% TRUE){
                    preds <- preds/newdata$Dose
                  }

                   if(conc_scale$log10_trans %in% TRUE){
                     preds <- log10(preds)
                   }

                   preds

                 })
         },
         simplify = FALSE,
         USE.NAMES = TRUE)
}
