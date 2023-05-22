#' Get residuals
#'
#' Extract residuals from a fitted `pk` object.
#'
#' Residuals are `obs - pred` in general, where `obs` is the observed
#' concentration value and `pred` is the predicted concentration value.
#'
#' For non-detect observations, residual is zero if `pred` is also below the
#' LOQ. Otherwise, the residual is the difference between the LOQ and `pred`.
#'
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions and compute residuals. If NULL (the default), then residuals
#'   will be computed for the data in `obj$data`. `newdata` is required to
#'   contain at least the following variables: `Time`, `Time.Units`, `Dose`,
#'   `Route`, `Media`, `Conc`, `Detect`. If variable `Time_trans` is not
#'   present, then `Time` will be transformed according to the transformation in
#'   `obj$scales$time` before making predictions; otherwise, `Time_trans` will
#'   be used to make predictions.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions and calculate residuals. If NULL (the default), residuals
#'   will be returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions and calculate residuals. If NULL (the
#'   default), residuals will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @param exclude Logical: `TRUE` to return `NA_real_` for any observations in
#'   the data marked for exclusion (if there is a variable `exclude` in the
#'   data, an observation is marked for exclusion when `exclude %in% TRUE`).
#'   `FALSE` to return the residual for each observation, regardless of
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
#' @return A named list of numeric matrices. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a matrix with the same number of
#'   rows as the data in `obj$data` (corresponding to the rows in `obj$data`),
#'   and as many columns as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The column names are the method names.  Each column
#'   contains the residuals (observed - predicted) of the model fitted by the
#'   corresponding method.  If `use_scale_conc %in% FALSE`, these residuals are
#'   in the same units as `obj$data$Conc.Units`. If `use_scale_conc %in% TRUE`,
#'   the residuals are in the same units as `obj$data$Conc_trans.Units`. If
#'   `use_scale_conc` was a named list, then the residuals are in units of
#'   `obj$data$Conc.Units` transformed as specified in `use_scale_conc`.
#' @export
#' @author Caroline Ring
#' @family methods for fitted pk objects
residuals.pk <- function(obj,
                         newdata = NULL,
                         model = NULL,
                         method = NULL,
                         exclude = TRUE,
                         use_scale_conc = FALSE,
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

  newdata_ok <- check_newdata(newdata = newdata,
                              olddata = obj$data,
                              req_vars = c("Time",
                                                 "Time.Units",
                                  "Dose",
                                  "Route",
                                  "Media",
                                  "Conc",
                                  "Detect"),
                              exclude = exclude)

  #Get predictions
  preds <- predict(obj,
                   newdata = newdata,
                   model = model,
                   method = method,
                   type = "conc",
                   exclude = exclude,
                   use_scale_conc = use_scale_conc)

  obs <- newdata$Conc

  #apply transformations if so specified
  conc_scale <- conc_scale_use(obj = obj,
                               use_scale_conc = use_scale_conc)
  if(conc_scale$dose_norm %in% TRUE){
    obs <- obs/newdata$Dose
  }

  if(conc_scale$log10_trans %in% TRUE){
    obs <- log10(obs)
  }


  resids <- sapply(preds,
         function(this_pred){
           apply(this_pred,
                 2,
                 function(x){
                   restmp <- ifelse(newdata$Detect %in% FALSE &
                            x <= obs,
                          0,
                          obs - x)
                   if(exclude %in% TRUE){
                     if("exclude" %in% names(newdata)){
                     restmp <- ifelse(newdata$exclude %in% TRUE,
                                      NA_real_,
                                      restmp)
                     }
                   }
                   restmp
                 }
           )
         },
         simplify = FALSE,
         USE.NAMES = TRUE)



  resids
}
