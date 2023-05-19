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
#'   contain at least the following variables: `Time`, `Dose`, `Route`, `Media`,
#'   `Conc`, `Detect`, `N_Subjects`, `Conc_SD`. `Time` will be transformed
#'   according to the transformation in `obj$scales$time`. Residuals will be
#'   calculated on un-transformed concentration data (i.e., ignoring any
#'   scalings or transformations in `obj$scales$conc`).
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
#' @param scale_conc Logical: `TRUE` to apply the concentration
#'   transformation(s) defined in `obj$scales$conc` to both predictions and
#'   observations before calculating residuals. `FALSE` to calculate residuals
#'   for un-transformed predictions and observations. Default `FALSE`.
#' @return A named list of numeric matrices. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a matrix with the same number of
#'   rows as the data in `obj$data` (corresponding to the rows in `obj$data`),
#'   and as many columns as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The column names are the method names.  Each column
#'   contains the residuals (observed - predicted) of the model fitted by the
#'   corresponding method.  If `scale_conc %in% FALSE`, these residuals are in
#'   the same units as `obj$data$Conc.Units`. If `scale_conc %in% TRUE`, the
#'   residuals are in the same units as `obj$data$Conc_trans.Units`.
#' @export
#' @author Caroline Ring
residuals.pk <- function(obj,
                         newdata = NULL,
                         model = NULL,
                         method = NULL,
                         exclude = TRUE,
                         scale_conc = FALSE,
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

  preds <- predict(obj,
                   newdata = newdata,
                   model = model,
                   method = method,
                   type = "conc",
                   exclude = exclude,
                   scale_conc = scale_conc)

  obs <- newdata$Conc

  if(scale_conc %in% TRUE){
    obs <- rlang::eval_tidy(obj$scales$conc$expr, data = cbind(newdata,
                                                                 data.frame(".conc" = obs)))
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
                     restmp[newdata$exclude %in% TRUE] <- NA_real_
                   }
                   restmp
                 }
           )
         },
         simplify = FALSE,
         USE.NAMES = TRUE)



  resids
}
