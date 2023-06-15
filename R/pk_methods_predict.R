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
                       ...){

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

  coefs <- my_coef(obj = obj) %>%
    dplyr::filter(model %in% model,
                  method %in% method)

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

  req_vars <- c("Time",
                "Time.Units",
                "Dose",
                "Route",
                "Media",
                "Value",
                "Value.Units")

  trans_vars <- c("Time_trans",
                  "Time_trans.Units",
                  "Conc_trans",
                  "Conc_trans.Units")

  newdata <- newdata %>%
    dplyr::select(!!!obj$data_group,
                  all_of(req_vars),
                  any_of(trans_vars)) %>%
    group_by(!!!obj$data_group,
             Route, Media) %>%
    nest(.key = "observations")

  newdata <- left_join(coefs, newdata,
                       relationship = "many-to-many")
  # From here, trying to call map within mutate to add the prediction column by
  # calling the appropriate model function and feeding it the parameters in coefs_vector

  newdata <- newdata %>%
    dplyr::group_by(Route, Media,
                    .add = TRUE) %>%
    mutate(model_fun = case_when(
      type == "conc" ~obj$stat_model[[model]]$conc_fun,
      type == "auc"  ~obj$stat_model[[model]]$auc_fun,
      .default = obj$stat_model[[model]]$conc_fun)) %>%
    reframe(predictions =
              map(observations,
                  .f = \(x) {
                    x %>%
                      group_by(Time_trans) %>%
                      mutate(Estimate = sapply(coefs_vector,
                                               FUN = model_fun,
                                               time = Time_trans,
                                               dose = ifelse(obj$scales$conc$dose_norm,
                                                             1, Dose),
                                               route = Route,
                                               medium = Media,
                                               simplify = TRUE,
                                               USE.NAMES = TRUE))
                  })) %>%
    unnest(predictions)

  if (type == "conc") newdata <- rename(newdata, Conc_est = "Estimate")
  if (type == "auc") newdata <- rename(newdata, AUC_est = "Estimate")


  return(newdata)
}
