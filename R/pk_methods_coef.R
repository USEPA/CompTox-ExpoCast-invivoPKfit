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
#' @author Caroline Ring
#' @family methods for fitted pk objects
coef.pk <- function(obj,
                    data_group = NULL,
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
    pivot_wider(names_from = param_name,
                values_from = start)

  possible_model_params <- sapply(obj$stat_model, `[[`, "params") %>%
    unlist() %>%
    unique()

  coefs <- obj$fit %>% unnest(cols = fit) %>%
    group_by(model, method, !!!obj$data_group) %>%
    dplyr::select(starts_with("sigma_"),
                  any_of(possible_model_params)) %>%
    unite(starts_with("sigma_"),
          col = "sigma",
          sep = "",
          na.rm = TRUE)

  coefs <- left_join(coefs, const_pars)

  coefs_tidy <- coefs %>%
    nest(coefs_tibble = any_of(possible_model_params)) %>%
    mutate(coefs_vector = map(coefs_tibble,
                              .f = \(x){
                                as.data.frame(x %>%
                                                dplyr::select(!where(is.na))) %>%
                                  unlist()
                              })) %>%
    unnest(coefs_tibble)
  return(coefs_tidy)
}
