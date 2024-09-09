#' Root mean squared error
#'
#' Extract root mean squared error of a fitted `pk` object
#'
#' # Formula for RMSE
#'
#' RMSE is calculated using the following formula, to properly handle summary
#' data:
#'
#' \deqn{
#' \sqrt{
#' \frac{1}{N}
#'  \sum_{i=1}^G \left( (n_i - 1) s_i^2 +
#'   n_i \bar{y}_i ^2 - 2 n_i \bar{y}_i \mu_i + \mu_i^2 \right)
#'    }
#' }
#'
#' In this formula, there are \eqn{G} observations, each of which may be for one
#' subject or for multiple subjects.
#'
#' - \eqn{n_i} is the number of subjects for observation \eqn{i}.
#' - \eqn{\bar{y}_i} is the sample mean concentration for observation \eqn{i}, with no transformations applied.
#' - \eqn{s_i} is the sample standard deviation of concentrations for observation \eqn{i}, with no transformations applied.
#' - \eqn{\mu_i} is the model-predicted concentration for observation \eqn{i}, with no transformations applied.
#'
#' \eqn{N} is the grand total of subjects across observations:
#'
#' \deqn{N = \sum_{i=1}^G n_i}
#'
#' For the non-summary case (\eqn{N} single-subject observations, with all
#' \eqn{n_i = 1}, \eqn{s_i = 0}, and \eqn{\bar{y}_i = y_i}), this formula
#' reduces to the familiar RMSE formula
#'
#' \deqn{\sqrt{\frac{1}{N} \sum_{i=1}^N (y_i - \mu_i)^2}}
#'
#' # Left-censored data
#'
#' If the observed value is censored, and the predicted value is less than the
#' reported LOQ, then the predicted value is (temporarily) set equal to the LOQ,
#' for an effective error of zero.
#'
#' If the observed value is censored, and the predicted value is greater than
#' the reported LOQ, the the observed value is treated as the reported LOQ (so
#' that the effective error is the difference between the LOQ and the predicted
#' value).
#'
#' # Log10 transformation
#'
#' If `log10_trans %in% TRUE`, then both the observed and predicted values will be
#' log10-transformed before calculating the RMSE. In the case where
#' observed values are reported in summary format, each sample mean and sample
#' SD (reported on the natural scale, i.e. the mean and SD of natural-scale
#' individual observations) are used to produce an estimate of the log10-scale
#' sample mean and sample SD (i.e., the mean and SD of log10-transformed
#' individual observations), using [convert_summary_to_log10()].
#'
#' The formulas are as follows. Again, \eqn{\bar{y}_i} is the sample mean for
#' group \eqn{i}. \eqn{s_i} is the sample standard deviation for group \eqn{i}.
#'
#' \deqn{\textrm{log10-scale sample mean}_i = \log_{10}
#' \left(\frac{\bar{y}_i^2}{\sqrt{\bar{y}_i^2 + s_i^2}} \right)}
#'
#' \deqn{\textrm{log10-scale sample SD}_i = \sqrt{\log_{10} \left(1 +
#' \frac{s_i^2}{\bar{y}_i^2} \right)}}
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions and compute RMSEs. If NULL (the default), then RMSEs will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Time.Units`, `Dose`, `Route`,
#'   `Media`, `Conc`, `Conc_SD`, `N_Subjects`, `Detect`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions and calculate RMSEs. If NULL (the default), RMSEs will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions and calculate RMSEs. If NULL (the default),
#'   RMSEs will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @param exclude Logical: `TRUE` to compute the RMSE excluding any observations
#'   in the data marked for exclusion (if there is a variable `exclude` in the
#'   data, an observation is marked for exclusion when `exclude %in% TRUE`).
#'   `FALSE` to include all observations, regardless of exclusion status.
#'   Default `TRUE`.
#' @param use_scale_conc Possible values: `FALSE` (default, `TRUE`, or a named list with
#'   elements `dose_norm` and `log10_trans` which themselves should be either
#'   `TRUE` or `FALSE`.  If `use_scale_conc = FALSE` (the default for this
#'   function), then no concentration scaling or transformation will be applied
#'   when the RMSE is computed. If `use_scale_conc = TRUE, then the
#'   concentration scaling/transformations in `obj` will be applied to both
#'   predicted and observed concentrations when the RMSE is computed (see
#'   [calc_rmse()] for details).If `use_scale_conc = list(dose_norm = ...,
#'   log10_trans = ...)`, then the specified dose normalization and/or
#'   log10-transformation will be applied when the RMSE is computed.
#' @param rmse_group A list of quosures provided in the format `vars(...)` that
#'   determines the data groupings for which RMSE is calculated. Default NULL,
#'   in which case RMSE is calculated for each data group defined in the
#'   object's `data_group` element (use [get_data_group.pk()] to access the
#'   object's `data_group`).
#' @param ... Additional arguments. Not currently used.
#' @return A `data.frame` with calculated RMSE as the final column. There is one row per
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data, each [optimx::optimx()] methods (specified in
#'   [settings_optimx()]), `rmse_group` specified.
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @family fit evaluation metrics
#' @family methods for fitted pk objects
#' @seealso [calc_rmse()]
rmse.pk <- function(obj,
                    newdata = NULL,
                    model = NULL,
                    method = NULL,
                    exclude = TRUE,
                    use_scale_conc = FALSE,
                    rmse_group = NULL,
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
  if(is.null(rmse_group)) rmse_group <- obj$data_group

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
  message("rmse.pk(): Computing RMSE on transformed concentration scale. Transformations used: \n",
          "Dose-normalization ", conc_scale$dose_norm, "\n",
          "log10-transformation ", conc_scale$log10_trans)

  #Get predictions
  #do NOT apply transformations at this stage
  preds <- predict(obj,
                   newdata = newdata,
                   model = model,
                   method = method,
                   type = "conc",
                   exclude = exclude,
                   use_scale_conc = FALSE)


  #remove any excluded observations & corresponding predictions, if so specified
  if (exclude %in% TRUE) {
    if ("exclude" %in% names(newdata)) {
      newdata <- newdata %>%
        dplyr::filter(exclude %in% FALSE)
    }
  }

  #Requested variables
  #Note that we take the NON-transformed concentrations.
  #Any dose-normalization will be done in the next step.
  #Any log10 transformations will be handled within the calc_rmse() function.
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
  rmse_df <- new_preds %>%
    dplyr::mutate(
      Conc_set = ifelse(rep(conc_scale$dose_norm,
                            NROW(Dose)),
                        Conc / Dose,
                        Conc),
      Conc_set_SD = ifelse(rep(conc_scale$dose_norm,
                               NROW(Dose)),
                           Conc_SD / Dose,
                           Conc_SD),
      Conc_est = ifelse(rep(conc_scale$dose_norm,
                            NROW(Dose)),
                        Conc_est / Dose,
                        Conc_est)
                        ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!!rmse_group,
                    model, method) %>%
    dplyr::summarize(
      RMSE = calc_rmse(obs = Conc_set,
                       obs_sd = Conc_set_SD,
                       pred = Conc_est,
                       n_subj = N_Subjects,
                       detect = Detect,
                       log10_trans = conc_scale$log10_trans)) %>%
    dplyr::ungroup()

  if (conc_scale$log10_trans == FALSE) {
    message("rmse.pk(): RMSE calculated by groups: \n",
            paste(sapply(unlist(rmse_group), rlang::as_label),
                  collapse = ", "),
            ", method, model")
  } else {
    message("rmse.pk(): RMSLE calculated by groups: \n",
            paste(sapply(unlist(rmse_group), rlang::as_label),
                  collapse = ", "),
            ", method, model")
    rmse_df <- rmse_df %>% rename(RMSLE = RMSE)
  }

  return(rmse_df)
}
