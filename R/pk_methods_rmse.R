#' Root mean squared error
#'
#' Extract root mean squared error of a fitted `pk` object
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions and compute RMSsE. If NULL (the default), then RMSEs will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Dose`, `Route`, and `Media`. `Time`
#'   will be transformed according to the transformation in `obj$scales$time`
#'   before RMSEs are calculated.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions and calculate RMSEs. If NULL (the default), RMSEs will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions and calculate RMSEs. If NULL (the default),
#'   RMSEs will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @param type Either `"conc"` (the default) or `"auc"`. `type = "conc"`
#'   predicts concentrations; `type = "auc"` predicts area under the
#'   concentration-time curve (AUC). Currently, only `type = "conc"` is
#'   implemented.
#' @return A named list of numeric vectors. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a numeric vector with as many
#'   elements as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The vector names are the method names.  Each vector
#'   element contains the root mean squared error (`sqrt(mean(observed -
#'   predicted)^2)`) of the model fitted by the corresponding method, using the
#'   data in `newdata`. These RMSEs are concentrations in the same units as
#'   `obj$data$Conc.Units`; any concentration transformations (in
#'   `obj$scale$conc`) are *not* applied.
#' @export
#' @author Caroline Ring
rmse.pk <- function(obj,
                    newdata = NULL,
                    model = NULL,
                    method = NULL,
                    type = "conc"){
  if(!(type %in% "conc")) stop(paste("Error in residuals.pk():",
                                     "only type = 'conc' is currently implemented;",
                                     "residuals for type = 'auc' are not available yet"))
  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data

  resids <- residuals(obj,
                      newdata = newdata,
                      model = model,
                      method = method,
                      type = type)

  sapply(resids,
         function(this_resid){
           apply(this_resid,
                 2,
                 function(x) sqrt(mean(x*x)))
         },
         simplify = FALSE,
         USE.NAMES = TRUE)
}
