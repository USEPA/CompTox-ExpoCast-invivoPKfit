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
#' @return  A data.frame with one row for each `data_group`, `model` and `method`.
#'   A column contains the fold errors (observed/predicted) of the model fitted by the
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
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  if (is.null(model)) model <- names(obj$stat_model)
  if (is.null(method)) method <- obj$settings_optimx$method

  if (!is.null(newdata)) {
    check_newdata(newdata = newdata,
                  olddata = obj$data,
                  req_vars = "Conc_trans")
  }

  #get predicted concentrations
  preds <- predict(obj,
                   newdata = newdata,
                   model = model,
                   method = method,
                   use_scale_conc = FALSE,
                   type = "conc") %>%
    mutate(Fold_Error = Conc_est/Conc)

  return(preds)

}
