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
#' @param type Either `"conc"` (the default) or `"auc"`. `type = "conc"`
#'   predicts concentrations; `type = "auc"` predicts area under the
#'   concentration-time curve (AUC). Currently, only `type = "conc"` is
#'   implemented.
#' @return A named list of numeric matrices. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a matrix with the same number of
#'   rows as the data in `obj$data` (corresponding to the rows in `obj$data`),
#'   and as many columns as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The column names are the method names.  Each column
#'   contains the residuals (observed - predicted) of the model fitted by the
#'   corresponding method. These residuals are concentrations in the same units
#'   as `obj$data$Conc.Units`; any concentration transformations (in
#'   `obj$scale$conc`) are *not* applied.
#' @export
#' @author Caroline Ring
residuals.pk <- function(obj,
                         newdata = NULL,
                         model = NULL,
                         method = NULL){

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data

  preds <- predict(obj,
                   newdata = newdata,
                   model = model,
                   method = method,
                   type = "conc")

  obs <- newdata$Conc

  sapply(preds,
         function(this_pred){
           apply(this_pred,
                 2,
                 function(x){
                   ifelse(newdata$Detect %in% FALSE &
                            x <= obs,
                          0,
                          obs - x)
                 }
           )
         },
         simplify = FALSE,
         USE.NAMES = TRUE)
}
