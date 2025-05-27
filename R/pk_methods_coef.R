#' Get coefficients
#'
#' Extract coefficients from a fitted [pk()] object
#'
#' This function extracts fitted model parameter values from a fitted [pk()]
#' object.
#'
#' @param object A [pk] object.
#' @param model Optional: Specify one or more of the fitted models whose
#'   coefficients to return. If NULL (the default), coefficients will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   whose coefficients to return. If NULL (the default), coefficients will be
#'   returned for all of the models in `obj$settings_optimx$method`.
#' @param drop_sigma Logical: `FALSE` by default. Determines whether to include
#'   sigma in the output.
#' @param include_NAs Logical: `FALSE` by default. Determines whether to include
#'   aborted fits which have NAs as coefficients.
#' @param include_type Character: `"use"` (default) will return all parameters
#'   used in evaluating the model, including those that were held constant.
#'   `"optimize"` will return only parameters that were optimized, dropping all
#'   that were held constant. `"constant"` will return *only* parameters that
#'   were held constant (used, but not optimized). (`"optimize"` and
#'   `"constant"` are useful, for example, when evaluating the Hessian of the
#'   log-likelihood function, which requires differentiating between parameters
#'   that were optimized and those that were held constant.) Any value other
#'   than `"use"`, `"optim"`, or `"const"` will return an error.
#' @param suppress.messages Logical: `NULL` by default to use the setting in
#'   `object$settings_preprocess$suppress.messages`. Determines whether to
#'   display messages.
#' @param ... Additional arguments currently not in use.
#' @return A data.frame with a row for each `data_group` x `method` x `model`
#'   combination in a fitted [pk()] object. When `drop_sigma = TRUE` there is
#'   also a row for each unique standard deviation hyper-parameter defined by
#'   `error_group` in the fitted [pk()] object. There is a column for all
#'   parameter estimates given each model in `model`. A list-column
#'   `coefs_vector` summarizes all estimated parameters into a named vector.
#'   This named vector is used in functions that call upon the model functions,
#'   such as [predict()].
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @family methods for fitted pk objects
coef.pk <- function(object,
                    model = NULL,
                    method = NULL,
                    drop_sigma = FALSE,
                    include_NAs = FALSE,
                    include_type = "use",
                    suppress.messages = NULL,
                    ...) {

  if (is.null(suppress.messages)) {
    suppress.messages <- object$settings_preprocess$suppress.messages
  }
  # Check fit status
  check <- check_required_status(obj = object,
                                 required_status = status_fit)

  if (check %in% FALSE) stop(attr(check, "msg")) # Stop if not fitted

  if (is.null(model))
    model <- names(object$stat_model)
  if (is.null(method))
    method <- object$settings_optimx$method

  method_ok <- check_method(obj = object, method = method)
  model_ok <- check_model(obj = object, model = model)

  if(!(include_type %in% c("use",
                           "const",
                           "optim"))){
    stop(paste0("Error in coef.pk(): `include_type` is\n",
                deparse(substitute(include_type)),
                "\n",
                "but it must be one of ",
                " 'use', 'optim', or 'const'."))
  }

  data_group_vars <- sapply(object$data_group, rlang::as_label)

  # Get a unique list of possible parameters for each model used
  possible_model_params <- sapply(object$stat_model, `[[`, "params") %>%
    unlist() %>%
    unique()

  coefs <- object$fit %>%
    dplyr::select(model, method,
                  !!!object$data_group,
                  param_name,
                  estimate,
                  convcode,
                  optimize_param,
                  use_param)

  # drop the sigma parameters (not used in some functions)
  if (drop_sigma %in% TRUE) {
    coefs <- coefs %>%
      dplyr::filter(
        !startsWith(param_name, "sigma_")
      )
  }

  # include NA values from aborted fits?
  if (include_NAs %in% FALSE) {
    coefs <- coefs %>%
      dplyr::filter(!(convcode %in% 9999),
                    !(convcode %in% -9999))
  }

  # Get the columns describing time units and their (possibly) transformed units
  time_group <- get_data(obj = object) %>%
    dplyr::select(!!!object$data_group, Time.Units, Time_trans.Units) %>%
    dplyr::distinct()

  # Create the coefs vector
  coefs_tidy <- coefs %>%
    dplyr::group_by(model, method, !!!object$data_group) %>%
    dplyr::summarise(
      coefs_vector = {
          outval <- setNames(estimate, param_name)
          #return only the selected "include_type"
          #this will be an empty vector if there are no params of the selected type
          if(include_type %in% "use"){
            list(outval[use_param %in% TRUE])
          }else if(include_type %in% "optim"){
            list(outval[optimize_param %in% TRUE &
                               use_param %in% TRUE])
          }else if(include_type %in% "const"){
            list(outval[optimize_param %in% FALSE &
                               use_param %in% TRUE])
          }
        }
      ) %>%
    dplyr::distinct() %>%
    dplyr::left_join(time_group,
                     by = data_group_vars) %>%
    dplyr::ungroup()

  # Various filtering steps and checks
  # By optimization method
  if (is.character(method)) {
    method_vector <- method
    if (suppress.messages %in% FALSE) {
    message("coef.pk(): Filtering by method(s): ", paste(method, collapse = " "))
    }
    coefs_tidy <- coefs_tidy %>% dplyr::filter(method %in% method_vector)
  }
  # By models used
  if (is.character(model)) {
    model_vector <- model
    if (suppress.messages %in% FALSE) {
      message("coef.pk(): Filtering by model(s): ", paste(model, collapse = " "))
    }
    coefs_tidy <- coefs_tidy %>% dplyr::filter(model %in% model_vector)
  }

  return(coefs_tidy)
}
