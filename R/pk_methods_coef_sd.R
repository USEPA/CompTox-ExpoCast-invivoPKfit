#' Get coefficient standard deviations
#'
#' Extract coefficient/parameter standard deviations from a fitted `pk` object
#'
#' The coefficient standard deviations are estimated by computing a numerical
#' approximation to the model Hessian (the matrix of second derivatives of the
#' model objective function with respect to each model parameter) and then
#' attempting to invert it. This procedure yields a variance/covariance matrix
#' for the model parameters. The square root of the diagonal elements of this
#' matrix represent the parameter standard deviations.
#'
#' If the Hessian is not invertible, an attempt is made to calculate a
#' pseudovariance matrix, following the procedure outlined in Gill & King
#' (2004). First, the generalized inverse of the Hessian is calculated using
#' [MASS::ginv()]. Then, a generalized Cholesky decomposition (to ensure
#' positive-definiteness) is calculated using [base::chol()] with argument
#' `pivot = TRUE`. The square root of the diagonal elements of this matrix
#' represent the parameter standard deviations.
#'
#' If neither of these procedures is successful, then `NA_real_` is returned for
#' all coefficient standard deviations.
#'
#' @param obj A [pk] object
#' @param model Optional: Specify one or more of the fitted models whose
#'   coefficients to return. If NULL (the default), coefficients will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   whose coefficients to return. If NULL (the default), coefficients will be
#'   returned for all of the models in `obj$settings_optimx$method`.
#' @param suppress.messages Logical. `TRUE` (the default) to suppress
#'   informative messages. `FALSE` to see them.
#' @param ... Additional arguments. Not in use right now.
#' @return A dataframe with one row for each `data_group`, `model` and `method`.
#'   The remaining columns include the parameters & hyperparameters as returned by
#'   [coef.pk()], as well as their calculated standard deviations.
#' @export
#' @importFrom MASS ginv
#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import numDeriv
#' @author Caroline Ring and Gilberto Padilla Mercado
#' @family methods for fitted pk objects
#' @references Gill J, King G. (2004) What to Do When Your Hessian is Not
#'   Invertible: Alternatives to Model Respecification in Nonlinear Estimation.
#'   Sociological Methods & Research 33(1):54-87. DOI: 10.1177/0049124103262681
coef_sd.pk <- function(obj,
                       model = NULL,
                       method = NULL,
                       suppress.messages = TRUE, ...){

  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  # get coefs data.frame for each model and method
  # but exclude non-optimized parameters and sigma values for error_group
  # coef does method & model checks
  coefs <- coef(
    obj = obj,
    model = model,
    method = method,
    drop_sigma = FALSE)

  other_vars <- ggplot2::vars(
    Value,
    Value.Units,
    Time_trans.Units,
    Conc_trans,
    Conc_trans.Units,
    data_sigma_group,
    exclude
  )

  data_grp_vars <- sapply(obj$data_group, rlang::as_label)

  # Get required variables for log_likelihood()
  req_vars <- ggplot2::vars(Time,
                            Time.Units,
                            Time_trans,
                            Dose,
                            Route,
                            Media,
                            Conc,
                            Conc_SD,
                            N_Subjects,
                            Detect,
                            LOQ)

  # Convert Time_trans to hours
  newdata <- obj$data %>%
    dplyr::select(!!!union(obj$data_group, req_vars),
                  !!!other_vars) %>%
    # log_likelihood() takes Time_trans so this must be converted to hours
    # so it is in concordance with coef()
    dplyr::mutate(data_sigma_group = factor(data_sigma_group),
                  Time_trans = convert_time(x = Time_trans,
                                            from = Time_trans.Units,
                                            to = "hours"),
                  Time_trans.Units = "hours") %>%
    dplyr::group_by(!!!obj$data_group) %>%
    tidyr::nest(.key = "observations") %>%
    dplyr::ungroup()

  # This setup allows for a more stable call to the model functions later on
  fun_models <- data.frame(
    model_name = unname(sapply(my_pk$stat_model, \(x) {x$name})),
    model_fun = unname(sapply(my_pk$stat_model, \(x) {x$conc_fun}))
  )

  newdata <- dplyr::left_join(coefs, newdata,
                              by = data_grp_vars) %>%
    dplyr::left_join(fun_models, join_by(model == model_name))


  newdata <- newdata %>%
    dplyr::rowwise() %>%
    dplyr::filter(!is.null(observations)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

  newdata <- suppressWarnings(newdata %>%
    dplyr::rowwise() %>%
    dplyr::mutate(hessian_mat = list(numDeriv::hessian(func = function(x){
      log_likelihood(par = x,
                     data = observations,
                     data_sigma_group = observations$data_sigma_group,
                     modelfun = model_fun,
                     dose_norm = obj$scales$conc$dose_norm,
                     log10_trans = obj$scales$conc$log10_trans,
                     negative = TRUE,
                     force_finite = TRUE)
    },
    x = coefs_vector,
    method = 'Richardson'))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sds = purrr::map2(hessian_mat, coefs_vector, \(x, y) {
      tryCatch(diag(solve(x))^(1/2) %>% as.numeric(),
               error = function(err){
                 if (!suppress.messages) {
                   message(paste0("Hessian can't be inverted, ",
                                  "using pseudovariance matrix ",
                                  "to estimate parameter uncertainty."))
                 }
                 # pseudovariance matrix
                 # see http://gking.harvard.edu/files/help.pdf
                 tryCatch(
                   suppressWarnings(diag(chol(MASS::ginv(x),
                             pivot = TRUE))^(1/2)),
                   error = function(err){
                     if (!suppress.messages) {
                       message(paste0("Pseudovariance matrix failed,",
                                      " returning NAs"))
                     }
                     rep(NA_real_, nrow(x))
                   })
               }) %>%
        set_names(nm = paste0(names(y), "_sd"))
    }),
    alerts = purrr::map2(hessian_mat, coefs_vector, \(x, y) {
      tryCatch({
        diag(solve(x))^(1/2) %>% as.numeric()
        return(paste0("Hessian successfully inverted"))
      },
      error = function(err){
        tryCatch(
          suppressWarnings(diag(chol(MASS::ginv(x),
                                     pivot = TRUE))^(1/2)),
          error = function(err){
            paste0("Pseudovariance matrix failed,",
                   " returning NAs")
          })
        paste0("Hessian can't be inverted, ",
               "using pseudovariance matrix ",
               "to estimate parameter uncertainty.")
      })
    })))

  # Limit data.frame columns to those that unique identify needed values
  newdata <- newdata %>%
    dplyr::select(model, method, !!!obj$data_group, coefs_vector, sds, alerts)

  newdata <- newdata %>%
    dplyr::mutate(sds_tibble = purrr::map(sds,
                                          \(x) as.list(x) %>%
                                            as.data.frame() %>%
                                            tidyr::pivot_longer(
                                              cols = tidyselect::everything(),
                                              names_to = "param_name_sd",
                                              values_to = "param_sd"))) %>%
    tidyr::unnest(sds_tibble) %>%
    dplyr::mutate(coefs_tibble = purrr::map(coefs_vector, \(x) as.list(x) %>%
                                              as.data.frame() %>%
                                              tidyr::pivot_longer(
                                                cols = tidyselect::everything(),
                                                names_to = "param_name",
                                                values_to = "param_value"
                                              ))) %>%
    tidyr::unnest(coefs_tibble) %>%
    dplyr::filter(stringr::str_detect(param_name_sd, param_name)) %>%
    dplyr::select(-c(coefs_vector, sds, param_name_sd)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(param_value = dplyr::na_if(param_value, NaN),
                  param_sd = dplyr::na_if(param_sd, NaN)) %>%
    dplyr::relocate(param_name, param_value, param_sd, .after = alerts)

  return(newdata)

}
