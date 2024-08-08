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
#' @param table_format Logical. `FALSE` by default, this determines whether the
#'   the output is in an unnested format.
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
                       table_format = FALSE,
                       suppress.messages = TRUE, ...){

  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$settings_optimx$method

  method_ok <- check_method(obj = obj, method = method)
  model_ok <- check_model(obj = obj, model = model)

  noptim_params <- obj$prefit$par_DF %>%
    dplyr::filter(optimize_param == FALSE,
                  use_param == TRUE) %>%
    dplyr::pull(param_name) %>%
    unique()

  other_vars <- ggplot2::vars(
    Value,
    Value.Units,
    Time_trans.Units,
    Conc_trans,
    Conc_trans.Units,
    data_sigma_group,
    exclude
  )

  # get coefs data.frame for each model and method
  # but exclude non-optimized parameters and sigma values for error_group
  coefs <- suppressMessages(
    coef(
      obj = obj,
      model = model,
      method = method,
      drop_sigma = FALSE
    ) %>%
      dplyr::select(coefs_vector, sigma_value, error_group) %>%
      dplyr::mutate(coefs_vector = purrr::map(coefs_vector,
                                              \(x) {
                                                sigma_transfer <- sigma_value
                                                names(sigma_transfer) <- error_group
                                                c(x[!(names(x) %in% noptim_params)], sigma_transfer)
                                              })) %>%
      dplyr::select(-sigma_value, -error_group)
  )


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
                            Detect)

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

  newdata <- tidyr::expand_grid(expand_grid(model, method),
                                newdata)

  newdata <- suppressMessages(dplyr::left_join(coefs, newdata))


  newdata <- suppressMessages(
    newdata %>%
      dplyr::rowwise() %>%
      dplyr::filter(!is.null(observations)) %>%
      dplyr::mutate(model_fun = obj$stat_model[[model]]$conc_fun) %>%
      dplyr::left_join(obj$prefit$par_DF %>%
                  dplyr::filter(param_name %in% noptim_params) %>%
                  dplyr::select(!!!obj$data_group, param_name, start)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(const_pars = setNames(start, param_name)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct())

  newdata <- newdata %>%
    dplyr::rowwise() %>%
    dplyr::mutate(hessian_mat = list(numDeriv::hessian(func = function(x){
      log_likelihood(par = x,
                     const_params = const_pars,
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
    }))

  # Limit data.frame columns to those that unique identify needed values
  newdata <- newdata %>%
    dplyr::select(model, method, !!!obj$data_group, coefs_vector, sds, alerts)

  # Pivot data.frame so that each parameter has a row
  if (table_format) {
    newdata <- newdata %>%
    dplyr::mutate(sds_tibble = purrr::map(sds,
                                          \(x) as.list(x) %>% as.data.frame)) %>%
      tidyr::unnest(sds_tibble) %>%
      tidyr::pivot_longer(cols = tidyselect::starts_with("sigma_"),
                          names_to = "error_group",
                          values_to = "sigma.value_sd") %>%
      dplyr::filter(!is.na(sigma.value_sd)) %>%
      dplyr::mutate(coefs_tibble = purrr::map(coefs_vector, \(x) as.list(x) %>% as.data.frame)) %>%
      tidyr::unnest(coefs_tibble) %>%
      tidyr::pivot_longer(cols = tidyselect::starts_with("sigma_"),
                          names_to = "error_group_value",
                          values_to = "sigma.value") %>%
      dplyr::filter(!is.na(sigma.value)) %>%
      dplyr::filter(error_group_value == gsub(x = error_group,
                                              pattern = "_sd", "")) %>%
      dplyr::select(-coefs_vector, -sds, -error_group_value) %>%
      dplyr::relocate(error_group, sigma.value, sigma.value_sd, .after = method)
  }

  return(newdata)

}
