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
#' @return A named list of numeric matrixes. There is one list element named for
#'   each model in `model`. Each list element is a matrix with as many rows as
#'   items in `method`. The row names are the method names. The matrix column
#'   names are the names of the fitted parameters, including any error standard
#'   deviation hyperparameters (whose names begin with "sigma"). The matrix
#'   elements are the values of the corresponding model parameters as fitted by
#'   the corresponding method.
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @family methods for fitted pk objects
coef.pk <- function(obj,
                    model = NULL,
                    method = NULL,
                    drop_sigma = FALSE,
                    ...) {
  # Check fit status
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)

  if (check %in% FALSE) stop(attr(check, "msg")) # Stop if not fitted

  const_pars <- subset(obj$prefit$par_DF,
                       optimize_param %in% FALSE &
                         use_param %in% TRUE) %>% dplyr::select(model,
                                                                !!!obj$data_group,
                                                                param_name,
                                                                start)
  const_pars <- const_pars %>%
    tidyr::pivot_wider(names_from = param_name,
                values_from = start)

  possible_model_params <- sapply(obj$stat_model, `[[`, "params") %>%
    unlist() %>%
    unique()

  coefs <- obj$fit %>% tidyr::unnest(cols = fit) %>%
    dplyr::group_by(model, method, !!!obj$data_group) %>%
    dplyr::select(dplyr::starts_with("sigma_"),
                  dplyr::any_of(possible_model_params)) %>%
    # If user wants sigmas
    # this pivot_longer provides the sigmas per error group
    tidyr::pivot_longer(cols = starts_with("sigma_"),
                        names_to = "error_group",
                        values_to = "sigma_value") %>%
    dplyr::filter(!is.na(sigma_value))

  coefs <- dplyr::left_join(coefs, const_pars)

  coefs_tidy <- coefs %>%
    tidyr::nest(coefs_tibble = dplyr::any_of(possible_model_params)) %>%
    dplyr::mutate(coefs_vector = map(coefs_tibble,
                              .f = \(x){
                                as.data.frame(x %>%
                                                dplyr::select(!where(is.na))) %>%
                                  unlist()
                              })) %>%
    tidyr::unnest(coefs_tibble)



  if (is.character(method)) {
    method_vector <- method
    message("Filtering by method(s): ", paste(method, collapse = " "))
    coefs_tidy <- coefs_tidy %>% dplyr::filter(method %in% method_vector)
  }
  if (is.character(model)) {
    model_vector <- model
    message("Filtering by model(s): ", paste(model, collapse = " "))
    coefs_tidy <- coefs_tidy %>% dplyr::filter(model %in% model_vector)
  }

  if (drop_sigma) {
    coefs_tidy <- coefs_tidy %>%
      dplyr::select(!c(sigma_value, error_group)) %>%
      dplyr::distinct()
  }

  return(coefs_tidy)
}

