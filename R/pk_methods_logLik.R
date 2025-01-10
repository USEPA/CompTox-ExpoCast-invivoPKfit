#' Log-likelihood
#'
#' Extract log-likelihood(s) from a fitted `pk` object
#'
#' For details on how the log-likelihood is calculated, see [log_likelihood()].
#'
#' # New levels in `newdata`
#'
#' The log-likelihood requires an error variance for each observation. Depending
#' on the error model used for the model fit, each observation may have a
#' different corresponding error variance, based on its unique combinations of
#' levels of factor variables as specified in [stat_error_model()] when setting
#' up the [pk()] object. (The error model specified in the `pk` object can be
#' viewed using [get_error_group()]).
#'
#' If you are supplying new data in `newdata`, then the unique combinations of
#' the factor levels for the new observations will be used to find the matching
#' error hyperparameter. If the new data contains new combinations of factor
#' levels not found in the data used to fit the model, then the following
#' procedure will be used to calculate the log-likelihood for the observations
#' with new levels:
#'
#' If there are \eqn{i = 1, 2, ... G} groups with unique combinations of the
#' factor levels in the original data, then there are corresponding error
#' standard deviations \eqn{\sigma_i = \sigma_1, \sigma_2, ..., \sigma_G}.
#'
#' Each observation with a new combination of factor levels will have \eqn{G}
#' differerent log-likelihoods computed, as though it were part of each of the
#' \eqn{G} existing groups. Then, the average of these \eqn{G} log-likelihoods
#' will be taken and assigned to the observation. In effect, each observation
#' with a new level is treated as though it is equally likely to belong to any
#' of the existing groups.
#'
#' # Scaling and transformation of concentration variables in `newdata`
#'
#' This function differs from some of the other methods for a fitted [pk()]
#' object that accept `newdata`, in that there is no `use_scale_conc` argument
#' for [logLik.pk()]. You cannot specify a different scaling/transformation for
#' concentration variables -- you have to use the same scaling/transformation
#' that was used to fit the model. This is because the log-likelihood depends on
#' the fitted values of the error variance hyperparameters, and those are valid
#' only for the transformation of the concentration data that was used to fit
#' the model.
#'
#'
#' @param object A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to compute
#'   log-likelihood. If NULL (the default), then log-likelihoods will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Time.Units`, `Dose`, `Route`,`Media`, `Conc`,
#'   `Detect`, `Conc_SD`, `N_Subjects`. `Time` will be transformed according to
#'   the transformation in `obj$scales$time` before making predictions. `Conc`
#'   will be transformed according to the transformation in `obj$scales$conc`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   calculate log-likelihood. If NULL (the default), log-likelihoods will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to calculate log-likelihoods. If NULL (the default),
#'   log-likelihoods will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @param force_finite Logical: Whether to force return of a finite value (e.g.
#'  as required by method `L-BFGS-B` in [optimx::optimx()]). Default FALSE. If
#'  TRUE, then if the log-likelihood works out to be non-finite, then it will be
#'  replaced with `.Machine$double.xmax`.
#' @param negative Logical: Whether to return the *negative* log-likelihood
#'  (i.e., the log-likelihood multiplied by negative 1). Default `FALSE`.
#' @param exclude Logical: `TRUE` to compute the log-likelihood excluding any
#'   observations in the data marked for exclusion (if there is a variable
#'   `exclude` in the data, an observation is marked for exclusion when `exclude
#'   %in% TRUE`). `FALSE` to include all observations when calculating the
#'   log-likelihood, regardless of exclusion status. Default `TRUE`.
#' @param drop_obs Logical: `TRUE` to drop the observations column after calculating log-likelihood.
#' @param ... Additional arguments. Not in use currently.
#' @return A data.frame with coefficients and log-likelihood values calculated using `newdata`.
#'   There is one row for each model in `obj`'s [stat_model()] element and
#'   each [optimx::optimx()] method (specified in [settings_optimx()]).
#' @export
#' @importFrom stats logLik
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @family fit evaluation metrics
#' @family log likelihood functions
#' @family methods for fitted pk objects
logLik.pk <- function(object,
                      newdata = NULL,
                      model = NULL,
                      method = NULL,
                      negative = FALSE,
                      force_finite = FALSE,
                      exclude = TRUE,
                      drop_obs = TRUE, ...) {

  suppress.messages <- object$settings_preprocess$suppress.messages

  # ensure that the model has been fitted
  check <- check_required_status(obj = object,
                                 required_status = status_fit)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }
  other_vars <- NULL
  if (is.null(model)) model <- names(object$stat_model)
  if (is.null(method)) method <- object$settings_optimx$method
  if (is.null(newdata)) {
    newdata <- object$data

    other_vars <- ggplot2::vars(
      Conc,
      Conc.Units,
      Time_trans.Units,
      Conc_trans,
      Conc_trans.Units,
      data_sigma_group,
      exclude
    )
  }

  method_ok <- check_method(obj = object, method = method)
  model_ok <- check_model(obj = object, model = model)

  # check variables in newdata
# including error grouping variables
  err_grp_vars <- sapply(eval(get_error_group(object)),
                         rlang::as_label)
  data_grp_vars <- sapply(eval(get_data_group(object)),
                         rlang::as_label)

  newdata_ok <- check_newdata(newdata = newdata,
                              olddata = object$data,
                              req_vars = union(c("Time",
                                           "Dose",
                                           "Route",
                                           "Media",
                                           "Conc",
                                           "Detect",
                                           "Conc_SD",
                                           "N_Subjects"),
                                err_grp_vars),
                              exclude = exclude)

  # time_scale_check
  if (!all(newdata$Time_trans.Units %in% "hours")) {
    message("logLik.pk(): Scaling these transformed time units back into hours ",
            "for log-likelihood calculation, to match time units of coefficients")
    # scale time if needed
    if (!("Time_trans" %in% names(newdata))) {
      newdata$Time_trans <- convert_time(x = newdata$Time,
                                         from = newdata$Time.Units,
                                         to = "hours")
      newdata$Time_trans.Units <- rep("hours", nrow(newdata))
    }
    if (!suppress.messages & (object$status < 5)) {
      print(newdata %>%
              dplyr::select(!!!object$data_group, Time.Units, Time_trans.Units) %>%
              dplyr::filter(Time.Units != Time_trans.Units) %>%
              dplyr::distinct())
    }
  }


  # get transformations to apply
  conc_scale <- conc_scale_use(obj = object,
                               use_scale_conc = TRUE)

  # remove any excluded observations & corresponding predictions, if so specified
  if (exclude %in% TRUE && "exclude" %in% names(newdata)) {
      newdata <- subset(newdata, exclude %in% FALSE)
  }

  # get coefs data.frame for each model and method
  # must include sigma value
  coefs <- coef(
    obj = object,
    model = model,
    method = method,
    drop_sigma = FALSE,
    suppress_message = TRUE)

  req_vars <- ggplot2::vars(Time,
                            Time.Units,
                            Time_trans,
                            Time_trans.Units,
                            Dose,
                            Route,
                            Media,
                            Conc,
                            Conc_SD,
                            N_Subjects,
                            Detect,
                            pLOQ)



  newdata <- suppressMessages(newdata %>%
    dplyr::select(!!!union(object$data_group, req_vars),
                  !!!other_vars))




    newdata <- newdata %>%
      dplyr::group_by(!!!object$data_group) %>%
      tidyr::nest(.key = "observations") %>%
      dplyr::ungroup()


  newdata <- tidyr::expand_grid(tidyr::expand_grid(model, method),
                                newdata)
  # This setup allows for a more stable call to the model functions later on
  fun_models <- data.frame(
    model_name = unname(sapply(object$stat_model, \(x) {x$name})),
    model_fun = unname(sapply(object$stat_model, \(x) {x$conc_fun}))
  )

  newdata <- dplyr::left_join(coefs, newdata,
                              by = c("model", "method",
                                      data_grp_vars))


  newdata <- newdata %>%
    dplyr::rowwise() %>%
    dplyr::filter(!is.null(observations)) %>%
    dplyr::left_join(fun_models, join_by(model == model_name)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

  newdata <- newdata %>%
    dplyr::rowwise() %>%
    dplyr::mutate(log_likelihood = log_likelihood(
      par = coefs_vector,
      data = observations,
      data_sigma_group = observations$data_sigma_group,
      modelfun = model_fun,
      dose_norm = conc_scale$dose_norm,
      log10_trans = conc_scale$log10_trans,
      negative = negative,
      force_finite = force_finite,
      suppress.messages = suppress.messages))

 if (drop_obs == TRUE) {
   newdata <- newdata %>%
     dplyr::select(-observations)
 }

  return(newdata)
}
