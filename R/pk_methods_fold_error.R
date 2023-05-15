#' Fold errors
#'
#' Calculate fold errors for a fitted `pk` object.
#'
#' Here, fold error is defined as `observed/predicted`.
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to compute
#'   fold errors. If NULL (the default), then fold errors will be computed for
#'   the data in `obj$data`. `newdata` is required to contain at least the
#'   following variables: `Time`, `Dose`, `Route`, and `Media`. `Time` will be
#'   transformed according to the transformation in `obj$scales$time` before
#'   fold errors are calculated.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions and calculate fold errors. If NULL (the default), fold errors will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions and calculate RMSEs. If NULL (the default),
#'   fold errors will be returned for all of the models in
#'   `obj$settings_optimx$method`.
#' @return A named list of numeric matrices. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a matrix with the same number of
#'   rows as the data in `obj$data` (corresponding to the rows in `obj$data`),
#'   and as many columns as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The column names are the method names.  Each column
#'   contains the fold errors (observed/predicted) of the model fitted by the
#'   corresponding method. These residuals are concentrations in the same units
#'   as `obj$data$Conc.Units`; any concentration transformations (in
#'   `obj$scale$conc`) are *not* applied.
#' @export
#' @author Caroline Ring
fold_errors.pk <- function(obj,
                          newdata = NULL,
                          model = NULL,
                          method = NULL){
  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$settings_optimx$method
  if(is.null(newdata)) newdata <- obj$data

  #get predicted concentrations
  preds <- predict(obj,
                   newdata = newdata,
                   model = model,
                   method = method,
                   type = "conc")

  sapply(preds,
         function(this_pred){
           apply(this_pred,
                 2,
                 function(x){
                   x <- ifelse(newdata$Detect %in% FALSE &
                                 (x <= newdata$Conc) %in% TRUE,
                               newdata$Conc,
                               x)
                   newdata$Conc/x
                 }
           )
         },
         simplify = FALSE,
         USE.NAMES = TRUE
  )
}
