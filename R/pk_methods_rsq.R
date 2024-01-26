#' Calculate R-squared for observed vs. predicted values
#'
#' Calculate the square of the Pearson correlation coefficient (r) between
#' observed and model-predicted values
#'
#' Calculate the square of the Pearson correlation coefficient (r) between
#' observed and model-predicted values, when observed data may be left-censored
#' (non-detect) or may be reported in summary form (as sample mean, sample
#' standard deviation, and sample number of subjects). Additionally, handle the
#' situation when observed data and predictions need to be log-transformed
#' before RMSE is calculated.
#'
#' \eqn{r^2} is calculated according to the following formula, to properly
#' handle multi-subject observations reported in summary format:
#'
#' \deqn{ r^2 = \left(
#' \frac{
#' \sum_{i=1}^G \mu_i n_i \bar{y}_i -
#' (\bar{\mu} + \bar{y}) \sum_{i=1}^G n_i \mu_i +
#' (\bar{mu} \bar{y}) \sum_{i=1}^G n_i
#' }
#' { \sqrt{
#' \sum_{i=1}^G (n_i - 1) s_i^2 +
#' \sum_{i=1}^G n_i \bar{y}_i^2 -
#' 2 \bar{y} \sum_{i=1}^G n_i \bar{y}_i +
#' N + \bar{y}^2
#'   }
#' \sqrt{
#' \sum_{i=1}^G n_i \mu_i^2 -
#' 2 \bar{y} \sum_{i=1}^G n_i \mu_i +
#' N + \bar{y}^2
#' }
#' } \right)^2
#' }
#'
#' In this formula, there are \eqn{G} groups (reported observations). (For CvTdb
#' data, a "group" is a specific combination of chemical, species, route,
#' medium, dose, and timepoint.) \eqn{n_i} is the number of subjects for group
#' \eqn{i}. \eqn{\bar{y}_i} is the sample mean for group \eqn{i}. \eqn{s_i} is
#' the sample standard deviation for group \eqn{i}.\eqn{\mu_i} is the
#' model-predicted value for group \eqn{i}.

#' \eqn{\bar{y}} is the grand mean of observations:
#'
#' \deqn{ \bar{y} = \frac{ \sum_{i=1}^G n_i \bar{y}_i }{\sum_{i=1}^G n_i} }
#'
#' \eqn{\bar{\mu}} is the grand mean of predictions:
#'
#' \deqn{ \bar{\mu} = \frac{ \sum_{i=1}^G n_i \mu_i }{\sum_{i=1}^G n_i} }
#'
#' \eqn{N} is the grand total of subjects:
#'
#' \deqn{N = \sum_{i=1}^G n_i}
#'
#' For the non-summary case (\eqn{N} single-subject observations, with all
#' \eqn{n_i = 1}, \eqn{s_i = 0}, and \eqn{\bar{y}_i = y_i}), this formula
#' reduces to the familiar formula
#'
#' \deqn{ r^2 = \left( \frac{\sum_{i=1}^N (y_i - \bar{y}) (\mu_i - \bar{\mu})}
#' {\sqrt{ \sum_{i=1}^N (y_i - \bar{y})^2 }
#' \sqrt{ \sum_{i=1}^N (\mu_i - \bar{\mu})^2 }
#'  } \right)^2
#'  }
#'
#' # Left-censored data
#'
#' If the observed value is censored, and the predicted value is less than the
#' reported LOQ, then the observed value is (temporarily) set equal to the
#' predicted value, for an effective error of zero.
#'
#' If the observed value is censored, and the predicted value is greater than
#' the reported LOQ, the the observed value is (temporarily) set equal to the
#' reported LOQ, for an effective error of (LOQ - predicted).
#'
#' # Log10 transformation
#'
#' If `log10_trans %in% TRUE`, then both the observed and predicted values will
#' be (natural) log-transformed before calculating the RMSE. In the case where
#' observed values are reported in summary format, each sample mean and sample
#' SD (reported on the natural scale, i.e. the mean and SD of natural-scale
#' individual observations) are used to produce an estimate of the log-scale
#' sample mean and sample SD (i.e., the mean and SD of log-transformed
#' individual observations), using [convert_summary_to_log10()].
#'
#' The formulas are as follows. Again, \eqn{\bar{y}_i} is the sample mean for
#' group \eqn{i}. \eqn{s_i} is the sample standard deviation for group \eqn{i}.
#'
#' \deqn{\textrm{log-scale sample mean}_i = \log
#' \left(\frac{\bar{y}_i^2}{\sqrt{\bar{y}_i^2 + s_i^2}} \right)}
#'
#' \deqn{\textrm{log-scale sample SD}_i = \sqrt{\log \left(1 +
#' \frac{s_i^2}{\bar{y}_i^2} \right)}}
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions and compute RMSsE. If NULL (the default), then R-squared will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Time.Units`, `Dose`, `Route`,
#'   `Media`, `Conc`, `Conc_SD`, `N_Subjects`, `Detect`. If variable
#'   `Time_trans` is not present, then `Time` will be transformed according to
#'   the transformation in `obj$scales$time` before making predictions;
#'   otherwise, `Time_trans` will be used to make predictions.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions and calculate R-squared. If NULL (the default), R-squared will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions and calculate R-squared. If NULL (the default),
#'   RMSEs will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @param exclude Logical: `TRUE` to compute the R-squared excluding any observations
#'   in the data marked for exclusion (if there is a variable `exclude` in the
#'   data, an observation is marked for exclusion when `exclude %in% TRUE`).
#'   `FALSE` to include all observations, regardless of exclusion status.
#'   Default `TRUE`.
#' @param use_scale_conc Possible values: `TRUE`, `FALSE`, or a named list with
#'   elements `dose_norm` and `log10_trans` which themselves should be either
#'   `TRUE` or `FALSE`. If `use_scale_conc = TRUE` (the default for this
#'   function), then the concentration scaling/transformations in `obj` will be
#'   applied to both predicted and observed concentrations before the
#'   R-squared is computed. If `use_scale_conc = FALSE`, then no
#'   concentration scaling or transformation will be applied before the
#'   R-squared is computed. If `use_scale_conc = list(dose_norm = ...,
#'   log10_trans = ...)`, then the specified dose normalization and/or
#'   log10-transformation will be applied.
#' @param ... Additional arguments. Not currently in use.
#' @return  A dataframe with one row for each `data_group`, `model` and `method`.
#'   The final column contains the R-squared of the model fitted by the corresponding
#'   method, using the data in `newdata`.
#' @export
#' @author Caroline Ring
#' @family fit evaluation metrics
#' @family methods for fitted pk objects
rsq.pk <- function(obj,
                   newdata = NULL,
                   model = NULL,
                   method = NULL,
                   exclude = TRUE,
                   use_scale_conc = TRUE,
                   ...){
  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  if (is.null(model)) model <- names(obj$stat_model)
  if (is.null(method)) method <- obj$optimx_settings$method
  if (is.null(newdata)) newdata <- obj$data

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
                                           "Conc_SD",
                                           "N_Subjects",
                                           "Detect"),
                              exclude = exclude)


  # Conc_trans columns will contain transformed values,
  conc_scale <- conc_scale_use(obj = obj,
                               use_scale_conc = use_scale_conc)
  message("Transformations used: \n",
          "Dose-normalization ", conc_scale$dose_norm, "\n",
          "log-transformation ", conc_scale$log10_trans)

  #Get predictions
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
      newdata <- newdata %>%
        dplyr::filter(exclude %in% FALSE)
    }
  }

  req_vars <- c(names(preds),
                "Conc",
                "Conc_SD",
                "N_Subjects",
                "Detect",
                "exclude")


  new_preds <- suppressMessages(dplyr::left_join(preds, newdata) %>%
    dplyr::select(dplyr::all_of(req_vars)) %>%
    dplyr::ungroup())


  #apply dose-normalization if specified
  # conditional mutate ifelse
  rsq_df <- new_preds %>%
    # needs to be rowwise first then by
    dplyr::rowwise() %>%
    dplyr::mutate(
      Conc_set = ifelse(conc_scale$dose_norm,
                        Conc / Dose,
                        Conc),
      Conc_set_SD = ifelse(conc_scale$dose_norm,
                           Conc_SD / Dose,
                           Conc_SD)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!!obj$data_group,
                    model, method) %>%
    dplyr::summarize(
      Rsq = calc_rsq(obs = Conc_set,
                       obs_sd = Conc_set_SD,
                       pred = Conc_est,
                       n_subj = N_Subjects,
                       detect = Detect,
                       log10_trans = conc_scale$log10_trans)) %>%
    dplyr::distinct() %>%
    dplyr::ungroup()

  message("Groups: \n",
          paste(sapply(unlist(obj$data_group), rlang::as_label),
                collapse = ", "),
          ", ", " method, model")

  return(rsq_df)
}
