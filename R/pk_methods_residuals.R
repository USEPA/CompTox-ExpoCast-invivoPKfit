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
#'   `Route`, `Media`, `Conc`, `Detect`.
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
#'   observed concentrations before the residuals are computed (i.e., the
#'   residuals will be computed on the same scale as the model was originally
#'   fitted). If `use_scale_conc = FALSE` (the default for this function), then
#'   no concentration scaling or transformation will be applied before the
#'   residuals are computed (i.e., the residuals will be computed on natural
#'   scale concentration data). If `use_scale_conc = list(dose_norm = ...,
#'   log10_trans = ...)`, then the specified dose normalization and/or
#'   log10-transformation will be applied.
#' @param ... Additional arguments not currently used.
#' @return A data.frame with the final column being calculated residuals.
#'   There is one row per each [optimx::optimx()] methods (specified in
#'   [settings_optimx()]), and `data_group`.  The final column
#'   contains the residuals (observed - predicted) of the model fitted by the
#'   corresponding method.  If `use_scale_conc %in% FALSE`, these residuals are
#'   in the same units as `obj$data$Conc.Units`. If `use_scale_conc %in% TRUE`,
#'   the residuals are in the same units as `obj$data$Conc_trans.Units`. If
#'   `use_scale_conc` was a named list, then the residuals are in units of
#'   `obj$data$Conc.Units` transformed as specified in `use_scale_conc`.
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
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

  method_ok <- check_method(obj = obj, method = method)
  model_ok <- check_model(obj = obj, model = model)

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

  preds <- predict(obj,
                   newdata = newdata,
                   model = model,
                   method = method,
                   type = "conc",
                   exclude = exclude,
                   use_scale_conc = use_scale_conc)


  #remove any excluded observations & corresponding predictions, if so specified
  if (exclude %in% TRUE) {
    if ("exclude" %in% names(newdata)) {
      newdata <- newdata %>% dplyr::filter(exclude %in% FALSE)
    }
  }

  req_vars <- c(names(preds),
                "Conc",
                "Detect",
                "exclude")

  new_preds <- suppressMessages(dplyr::left_join(preds, newdata) %>%
    dplyr::select(dplyr::all_of(req_vars)) %>%
    dplyr::ungroup())


  # Conc_trans columns will contain transformed values,
  conc_scale <- conc_scale_use(obj = obj,
                               use_scale_conc = use_scale_conc)

  message("residuals.pk(): Residuals calculated using the following transformations: \n",
          "Dose-normalization ", conc_scale$dose_norm, "\n",
          "log-transformation ", conc_scale$log10_trans)

  #apply dose-normalization if specified
  # conditional mutate ifelse
  resids <- new_preds %>%
    dplyr::mutate(
      Conc_set = ifelse(rep(conc_scale$dose_norm, NROW(Dose)),
                        ifelse(rep(conc_scale$log10_trans, NROW(Dose)),
                               log10(Conc/Dose),
                               Conc / Dose),
                        ifelse(rep(conc_scale$log10_trans, NROW(Dose)),
                               log10(Conc),
                               Conc)
                        ),
      Residuals = ifelse(Detect %in% FALSE & Conc_est <= Conc_set,
                         0,
                         Conc_est - Conc_set),
      .after = Conc_est) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

  return(resids)
}
