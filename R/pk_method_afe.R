#' Calculate average fold error
#'
#' Calculate average fold error
#'
#' Average fold error is calculated as
#'
#' \deqn{
#' 10^{\frac{1}{N}\sum{log_{10} \left( \frac{\textrm{predicted}}{\textrm{observed}} \right)}}
#' }
#'
#' # Left-censored data
#'
#' If the observed value is censored, and the predicted value is less than the
#' reported LOQ, then the observed value is (temporarily) set equal to the
#' predicted value, for an effective error of zero.
#'
#' If the predicted value is less than the reported LOQ, then the user may
#' choose whether to (temporarily) set the predicted value equal to LOQ, using
#' argument `sub_pLOQ`).
#'
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions and compute AFE. If NULL (the default), then AFE will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Time.Units`, `Dose`, `Route`,
#'   `Media`, `Conc`, `Conc_SD`, `N_Subjects`, `Detect`, `pLOQ`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions and calculate AFE. If NULL (the default), AFE will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions and calculate AFE. If NULL (the default),
#'   RMSEs will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @param exclude Logical: `TRUE` to compute the AFE excluding any observations
#'   in the data marked for exclusion (if there is a variable `exclude` in the
#'   data, an observation is marked for exclusion when `exclude %in% TRUE`).
#'   `FALSE` to include all observations, regardless of exclusion status.
#'   Default `TRUE`.
#' @param AFE_group Default: Chemical, Species. Determines what the data
#' grouping that is used to calculate average fold error (AFE). Should be set to lowest number
#' of variables that still would return unique experimental conditions.
#' Input in the form of `ggplot2::vars(Chemical, Species, Route, Media, Dose)`.
#' @param sub_pLOQ TRUE (default): Substitute all predictions below the LOQ with
#'   the LOQ before computing AFE. FALSE: do not.
#' @param ... Additional arguments. Not currently in use.
#' @return  A dataframe with one row for each `data_group`, `model` and `method`.
#'   The final column contains the AFE of the model fitted by the corresponding
#'   method, using the data in `newdata`.
#' @export
#' @author Caroline Ring
#' @family fit evaluation metrics
#' @family methods for fitted pk objects
AFE.pk <- function(obj,
                   newdata = NULL,
                   model = NULL,
                   method = NULL,
                   exclude = TRUE,
                   use_scale_conc = FALSE,
                   AFE_group = NULL,
                   sub_pLOQ = TRUE,
                   ...) {
  # ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  if (is.null(model)) model <- names(obj$stat_model)
  if (is.null(method)) method <- obj$optimx_settings$method
  if (is.null(newdata)) newdata <- obj$data
  if (is.null(AFE_group)) AFE_group <- obj$data_group

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
                                           "LOQ",
                                           "Detect",
                                           "pLOQ"),
                              exclude = exclude)



  AFE_group_char <- sapply(AFE_group, rlang::as_label)


  # Get predictions
  # Do NOT apply any conc transformations at this stage
  # (conc transformations will be handled later)
  preds <- predict(obj,
                   newdata = newdata,
                   model = model,
                   method = method,
                   type = "conc",
                   exclude = exclude,
                   use_scale_conc = FALSE)


  # remove any excluded observations & corresponding predictions, if so specified
  if (exclude %in% TRUE && "exclude" %in% names(newdata)) {
    newdata <- subset(newdata, exclude %in% FALSE)
  }

  req_vars <- unique(c(names(preds),
                       AFE_group_char,
                       "Conc",
                       "Conc_SD",
                       "N_Subjects",
                       "Detect",
                       "exclude",
                       "pLOQ"))


  new_preds <- suppressMessages(dplyr::left_join(preds, newdata) %>%
                                  dplyr::select(dplyr::all_of(req_vars)) %>%
                                  dplyr::ungroup())

  #replace below-LOQ preds with pLOQ if specified
  if(sub_pLOQ %in% TRUE){
    message("AFE.pk(): Predicted conc below pLOQ substituted with pLOQ")
    new_preds <- new_preds %>%
      dplyr::mutate(Conc_est_tmp = dplyr::if_else(Conc_est < pLOQ,
                                              pLOQ,
                                              Conc_est))
  }

  #if Conc censored and Conc_est_tmp < LOQ, make fold error 1
  new_preds <- new_preds %>%
    dplyr::mutate(Conc_tmp = dplyr::if_else(Detect %in% FALSE &
                                              Conc_est_tmp < LOQ,
                                            Conc_est_tmp,
                                            Conc))

  # apply dose-normalization if specified
  # conditional mutate ifelse
  AFE_df <- new_preds %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!!AFE_group,
                    model, method) %>%
    dplyr::summarize(
      AFE = 10^mean(log10(Conc_est_tmp/Conc_tmp))) %>%
    # dplyr::distinct() %>%
    dplyr::ungroup()

  message("AFE.pk)(): Groups: \n",
          toString(AFE_group_char),
          ", method, model")

  return(AFE_df)
}
