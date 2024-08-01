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

#### Adding sigmas and constant parameters

  coefs <- suppressMessages(obj$fit %>% tidyr::unnest(cols = fit) %>%
                              dplyr::group_by(model,method, !!!obj$data_group) %>%
                              dplyr::select(dplyr::starts_with("sigma_"),
                                            tidyselect::any_of(possible_model_params)) %>%
                              tidyr::pivot_longer(cols = starts_with("sigma_"),
                                                  names_to = "error_group",
                                                  values_to = "sigma_value") %>%
                              dplyr::filter(stringr::str_detect(error_group,
                                                                paste(Chemical,
                                                                      Species,
                                                                      sep = ".")))
                            )




  coefs <- suppressMessages(dplyr::left_join(coefs, const_pars,
                                             by = c("model",
                                                    sapply(obj$data_group,
                                                           rlang::as_label))))

  if (any(stringr::str_detect(names(coefs), pattern = "\\.(x|y)$"))) {
    message("Coalescing parameter values...")
    # This next thing should replace NAs with the constant value...
    for (this_param in intersect(names(const_pars), possible_model_params)) {
      if (sum(stringr::str_detect(names(coefs), pattern = this_param)) > 1) {
        coefs <- dplyr::mutate(coefs,
                               {{this_param}} := coalesce(
                                 .data[[paste(this_param, "x",
                                              sep = ".")]],
                                 .data[[paste(this_param, "y",
                                              sep = ".")]]),
                               .keep = "unused")
      }
    }

  }

####

  time_group <- get_data(obj = obj) %>%
    dplyr::select(!!!obj$data_group, Time.Units, Time_trans.Units) %>%
    dplyr::distinct()

  coefs_tidy <- suppressMessages(coefs %>%
    tidyr::nest(coefs_tibble = dplyr::any_of(possible_model_params)) %>%
    dplyr::mutate(coefs_vector = purrr::map(coefs_tibble,
                              .f = \(x){
                                as.data.frame(x %>%
                                                dplyr::select(!where(is.na))) %>%
                                  unlist()
                              })) %>%
    tidyr::unnest(coefs_tibble) %>%
    dplyr::left_join(time_group))



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

  if (!include_NAs) {
    no_fits <- obj$prefit$fit_check %>%
      dplyr::filter(fit_decision %in% "abort") %>%
      dplyr::select(model, Chemical, Species)
    coefs_tidy <- suppressMessages(dplyr::anti_join(coefs_tidy, no_fits))
  }

  return(coefs_tidy)
}

