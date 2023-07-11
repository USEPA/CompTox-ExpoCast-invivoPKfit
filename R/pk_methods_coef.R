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
#' @return A data.frame with a row for each `data_group` x `method` x `model` combination
#'  in a fitted [pk()] object. When `drop_sigma = TRUE` there is also a row for each
#'  unique standard deviation hyper-parameter defined by `error_group` in the fitted [pk()] object.
#'  There is a column for all parameter estimates given each model in `model`.
#'  A list-column `coefs_vector` summarizes all estimated parameters into a named vector.
#'  This named vector is used in functions that call upon the model functions, such as [predict()].
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

  for (this_param in intersect(names(const_pars), possible_model_params)) {
    coefs[this_param] <- tidyr::replace_na(coefs[[this_param]],
                                           replace = unique(const_pars[[this_param]]))
  }

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

