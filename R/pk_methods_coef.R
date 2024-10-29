#' Get coefficients
#'
#' Extract coefficients from a fitted [pk()] object
#'
#' This function extracts fitted model parameter values from a fitted [pk()]
#' object.
#'
#' @param obj A [pk()] object
#' @param model Optional: Specify (as a `character` vector) one or more of the
#'   fitted models whose coefficients to return. If `NULL` (the default),
#'   coefficients will be returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify (as a `character` vector)one or more of the
#'   [optimx::optimx()] methods whose coefficients to return. If `NULL` (the
#'   default), coefficients will be returned for all of the models in
#'   `obj$settings_optimx$method`.
#' @param drop_sigma Logical: `FALSE` by default. Determines whether to include
#'  sigma in the output.
#' @param include_NAs Logical: `FALSE` by default. Determines whether to include
#'  aborted fits which have NAs as coefficients.
#' @param ... Additional arguments currently not in use.
#' @return A data.frame with a row for each `data_group` x `method` x `model` combination
#'  in a fitted [pk()] object. When `drop_sigma = TRUE` there is also a row for each
#'  unique standard deviation hyper-parameter defined by `error_group` in the fitted [pk()] object.
#'  There is a column for all parameter estimates given each model in `model`.
#'  A list-column `coefs_vector` summarizes all estimated parameters into a named vector.
#'  This named vector is used in functions that call upon the model functions, such as [predict()].
#' @import glue
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @family methods for fitted pk objects
coef.pk <- function(obj,
                    model = NULL,
                    method = NULL,
                    drop_sigma = FALSE,
                    include_NAs = FALSE,
                    ...) {
  # Check fit status
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)

  if (check %in% FALSE) stop(attr(check, "msg")) # Stop if not fitted

  if (is.null(model))
    model <- names(obj$stat_model)
  if (is.null(method))
    method <- obj$settings_optimx$method

  method_ok <- check_method(obj = obj, method = method)
  model_ok <- check_model(obj = obj, model = model)

  data_group_vars <- sapply(obj$data_group, rlang::as_label)

  # Get a unique list of possible parameters for each model used
  possible_model_params <- sapply(obj$stat_model, `[[`, "params") %>%
    unlist() %>%
    unique()

  # Get the parameters that were held constant from par_DF
  # And only keep the name and starting value for the parameter
  # along with the unique identifying columns model and data_group
  parDF <- subset(
    obj$prefit$par_DF,
    use_param == TRUE & optimize_param == FALSE &
      param_name %in% possible_model_params) %>%
    dplyr::select(model, !!!obj$data_group, param_name, start)

  coefs <- subset(
    obj$fit,
    subset = (use_param == TRUE)) %>%
    dplyr::select(model, method,
                  !!!obj$data_group,
                  param_name,
                  estimate,
                  convcode)

  parDF <- coefs %>%
    dplyr::distinct(model, method, !!!obj$data_group, convcode) %>%
    dplyr::inner_join(parDF, by = c(data_group_vars, "model")) %>%
    dplyr::rename(estimate = "start")

  coefs <- dplyr::bind_rows(coefs, parDF)

  # drop the sigma parameters (not used in some functions)
  if (drop_sigma == TRUE) {
    coefs <- coefs %>%
      dplyr::filter(
        stringr::str_detect(param_name, pattern = "^sigma_",
                            negate = TRUE))
  }

  # include NA values from aborted fits
  if (!include_NAs) {
    coefs <- coefs %>%
      dplyr::filter(!(convcode %in% 9999),
                    !(convcode %in% -9999))
  }

  # Get the columns describing time units and their (possibly) transformed units
  time_group <- get_data(obj = obj) %>%
    dplyr::select(!!!obj$data_group, Time.Units, Time_trans.Units) %>%
    dplyr::distinct()

  # Create the coefs vector
  coefs_tidy <- coefs %>%
    dplyr::group_by(model, method, !!!obj$data_group) %>%
    dplyr::reframe(coefs_vector = purrr::map2(
      estimate, param_name,
      .f = \(x,y){
        setNames(estimate, param_name)
      }
    )) %>%
    dplyr::left_join(time_group,
                     by = data_group_vars)

  # Various filtering steps and checks
  # By optimization method
  if (is.character(method)) {
    method_vector <- method
    message("coef.pk(): Filtering by method(s): ", paste(method, collapse = " "))
    coefs_tidy <- coefs_tidy %>% dplyr::filter(method %in% method_vector)
  }
  # By models used
  if (is.character(model)) {
    model_vector <- model
    message("coef.pk(): Filtering by model(s): ", paste(model, collapse = " "))
    coefs_tidy <- coefs_tidy %>% dplyr::filter(model %in% model_vector)
  }

  return(coefs_tidy)
}
