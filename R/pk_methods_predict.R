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
#' @param include_NAs Logical: `FALSE` by default. Determines whether to include
#'  aborted fits which have NAs as coefficients.
#' @return A data.frame with one row for each `data_group`, `model` and `method`.
#'   A column that contains the predicted concentration or AUC at that timepoint
#'   given the TK parameters for that `model` and `method` specified in [coefs()].
#'   If `use_scale_conc %in% FALSE`, these predictions are
#'   un-transformed concentrations in the same units as `obj$data$Conc.Units`.
#'   If `use_scale_conc %in% TRUE`, the predictions are transformed
#'   concentrations in the same units as `obj$data$Conc_trans.Units`.
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @family methods for fitted pk objects
predict.pk <- function(obj,
                       newdata = NULL,
                       model = NULL,
                       method = NULL,
                       type = "conc",
                       exclude = TRUE,
                       use_scale_conc = FALSE,
                       suppress_messages = TRUE,
                       include_NAs = FALSE,
                       ...) {
  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  if (is.null(model))
    model <- names(obj$stat_model)
  if (is.null(method))
    method <- obj$settings_optimx$method

  method_ok <- check_method(obj = obj, method = method)
  model_ok <- check_model(obj = obj, model = model)

  coefs <- coef(
    obj = obj,
    model = model,
    method = method,
    drop_sigma = TRUE,
    include_NAs = include_NAs
  )

  other_vars <- NULL

  if (is.null(newdata)) {
    newdata <- obj$data

    other_vars <- ggplot2::vars(
      Value,
      Value.Units,
      Time_trans.Units,
      Conc,
      Conc.Units,
      Conc_trans,
      Conc_trans.Units,
      Detect,
      exclude
    )
  }

  newdata_ok <- check_newdata(
    newdata = newdata,
    olddata = obj$data,
    req_vars = c("Time",
                 "Time.Units",
                 "Dose",
                 "Route",
                 "Media"),
    exclude = exclude
  )

  #scale time if needed
  if (!("Time_trans" %in% names(newdata))) {
    newdata$Time_trans <- convert_time(
      x = newdata$Time,
      from = newdata$Time.Units,
      to = obj$scales$time$new_units
    )
  }

  #apply transformations if so specified
  conc_scale <- conc_scale_use(obj = obj,
                               use_scale_conc = use_scale_conc)

  req_vars <- ggplot2::vars(Time,
                            Time.Units,
                            Time_trans,
                            Dose,
                            Route,
                            Media)

  newdata <- newdata %>%
    dplyr::select(!!!union(obj$data_group, req_vars),
                  !!!other_vars) %>%
    dplyr::group_by(!!!union(obj$data_group,
                             ggplot2::vars(Route, Media))) %>%
    tidyr::nest(.key = "observations") %>%
    dplyr::ungroup()

  newdata <- tidyr::expand_grid(expand_grid(model, method),
                                newdata)

  newdata <- dplyr::left_join(coefs, newdata)


  newdata <- newdata %>%
    rowwise() %>%
    filter(!is.null(observations)) %>%
    ungroup()

  # After join it is joined by model, method, Chemical, Species
  newdata <- newdata %>%
    dplyr::group_by(model, method,
                    !!!obj$data_group,
                    Route, Media) %>%
    dplyr::mutate(
      model_fun = dplyr::case_when(
        type == "conc" ~ obj$stat_model[[model]]$conc_fun,
        type == "auc"  ~ obj$stat_model[[model]]$auc_fun,
        .default = obj$stat_model[[model]]$conc_fun
      )
    )

  newdata <- newdata %>%
    dplyr::reframe(predictions =
                     purrr::map(observations,
                                .f = \(x) {
                                  x %>%
                                    dplyr::rowwise() %>%
                                    dplyr::mutate(
                                      Estimate = tryCatch(
                                        sapply(
                                          coefs_vector,
                                          FUN = model_fun,
                                          time = Time_trans,
                                          dose = ifelse(conc_scale$dose_norm,
                                                        1, Dose),
                                          route = Route,
                                          medium = Media,
                                          simplify = TRUE,
                                          USE.NAMES = TRUE
                                        ),
                                        error = function(err) {
                                          if (!suppress_messages) {
                                            message(paste("Unable to run",
                                                          model_fun, "for",
                                                          Chemical, Species,
                                                          "data grouping.",
                                                          "Likely an aborted fit,",
                                                          "it is missing estimated parameters."))
                                          }
                                          # Return Value
                                          NA
                                        })
                                    )
                                })) %>%
    tidyr::unnest(predictions)

  if (type == "conc")
    newdata <- dplyr::rename(newdata, Conc_est = "Estimate")
  if (type == "auc")
    newdata <- dplyr::rename(newdata, AUC_est = "Estimate")

  if (conc_scale$dose_norm) {
    message("Note that the estimated values are Dose-normalized")
  } else {
    message("Note these values scale with Dose! (Not Dose-normalized)")
  }
  return(newdata)
}
