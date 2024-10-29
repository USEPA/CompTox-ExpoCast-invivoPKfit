#' Get predictions
#'
#' Extract predictions from a fitted `pk` object.
#'
#' @param obj A [pk] object.
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions. If NULL (the default), then predictions will be made for the
#'   data in `obj$data`. `newdata` is required to contain at least the following
#'   variables: `Time`, `Time.Units`, `Dose`, `Route`, and `Media`.
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
#' @param suppress_messages Whether to suppress warnings and other
#'  diagnostic messages during plotting
#' @param include_NAs Logical: `FALSE` by default. Determines whether to include
#'  aborted fits which have NAs as coefficients.
#' @param ... Additional arguments.
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
  if (is.null(model)) model <- names(obj$stat_model)
  if (is.null(method)) method <- obj$settings_optimx$method

  # From here it needs to output a named numeric vector coefs_vector
  # for the model functions
  # Responsibility for some checks done in coefs
  coefs <- coef(
    obj = obj,
    model = model,
    method = method,
    drop_sigma = TRUE,
    include_NAs = include_NAs,
    suppress_messages = suppress_messages
  ) %>%
    dplyr::select(-c(Time.Units, Time_trans.Units))

  # This setup allows for a more stable call to the model functions later on
  fun_models <- data.frame(
    model_name = unname(vapply(obj$stat_model, \(x) {x$name}, character(1))),
    model_fun = if (type == "auc") {
      unname(sapply(obj$stat_model, \(x) {x$auc_fun}))
    } else {
      unname(sapply(obj$stat_model, \(x) {x$conc_fun}))
    }
  )

  data_group_vars <- sapply(obj$data_group,
                            rlang::as_label)

  req_vars <- union(obj$data_group,
                    ggplot2::vars(
                      Conc.Units,
                      Time,
                      Time.Units,
                      Dose,
                      Route,
                      Media))

  if (is.null(newdata)) {
    newdata <- obj$data
  }

  # Check if there are other variables

  newdata_ok <- check_newdata(
    newdata = newdata,
    olddata = obj$data,
    req_vars = sapply(req_vars,
                      rlang::as_label),
    exclude = exclude
  )


  #apply transformations if so specified
  conc_scale <- conc_scale_use(obj = obj,
                               use_scale_conc = use_scale_conc)

  # Make observations into nested list-column
  newdata <- newdata %>%
    dplyr::group_by(!!!obj$data_group) %>%
    tidyr::nest(.key = "observations") %>%
    dplyr::ungroup()

  # Set each observation per data group with a model and method used
  newdata <- tidyr::expand_grid(tidyr::expand_grid(model, method),
                                newdata)

  # Add the coefs
  newdata <- dplyr::left_join(coefs, newdata,
                              by = c(data_group_vars,
                                     "model", "method"))

  # Remove any NULL observations
  newdata <- newdata %>%
    dplyr::rowwise() %>%
    dplyr::filter(!is.null(observations)) %>%
    dplyr::ungroup()


  # After join it is joined by model, method, Chemical, Species
  # Set a new column for the model function
  newdata <- newdata %>%
    dplyr::left_join(fun_models,
                     join_by(model == model_name)) %>%
    dplyr::distinct()

  # Get predictions
  # Note that the model functions only need Time, Dose, Route, and Medium
  newdata <- newdata %>% # Rowwise drops
    dplyr::rowwise(model, method, !!!obj$data_group) %>% # Needs to include columns outside the nest
    dplyr::summarise(predictions = list(
      observations %>%
        dplyr::mutate(
          Dose_tmp = dplyr::if_else(rep(conc_scale$dose_norm,
                                        NROW(Dose)),
                                    1.0,
                                    Dose),
          Estimate = tryCatch(
            do.call(model_fun,
                    list(coefs_vector,
                         time = Time,
                         dose = Dose_tmp,
                         route = Route,
                         medium = Media
                         )),
            error = function(err) {
              if (!suppress_messages) {
                message(paste("predict.pk(): Unable to run",
                              model_fun, "for",
                              Chemical, Species,
                              "data grouping.",
                              "Likely an aborted fit,",
                              "it is missing estimated parameters."))
              }
              # Return Value
              NA_real_
            }), #end tryCatch
          .after = Conc.Units)  %>%  #end dplyr::mutate
        dplyr::select(-Dose_tmp)
    )) %>%
    tidyr::unnest(predictions)



  #If log10 transformation was specified, then apply it now
  if (type %in% "conc") {
    newdata <- dplyr::rename(newdata, Conc_est = "Estimate")
  #apply log10-trans to predicted conc, if so specified
  if (conc_scale$log10_trans %in% TRUE) {
    newdata <- newdata %>%
      dplyr::mutate(Conc_est = log10(Conc_est))
  }
  } else if (type %in% "auc") {
    newdata <- dplyr::rename(newdata, AUC_est = "Estimate")
  #note that it doesn't make sense to log10-trans AUC
    if (suppress_messages %in% FALSE) {
      message("predict.pk(): Log10 transformation was specified, but was not used because `type == 'AUC'`.")
    }
  }

  if (suppress_messages %in% FALSE) {
    if(conc_scale$dose_norm) {
    message("predict.pk(): Note that the predicted values are for dose 1.0 (dose-normalized)")
  } else {
    message("predict.pk(): Note that the predicted values are not dose-normalized")
  }
  }

  message("predict.pk(): These predictions have been made using un-scaled Time and 1/hour rate constants from coefs()")
  return(newdata)
}
