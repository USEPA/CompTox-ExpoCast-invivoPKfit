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
#' @return A named list of numeric matrixes. There is one list element named for
#'   each model in `model`. Each list element is a matrix with as many rows as
#'   items in `method`. The row names are the method names. The matrix column
#'   names are the names of the fitted parameters, including any error standard
#'   deviation hyperparameters (whose names begin with "sigma").
#' @export
#' @author Caroline Ring
#' @family methods for fitted pk objects
#' @references Gill J, King G. (2004) What to Do When Your Hessian is Not
#'   Invertible: Alternatives to Model Respecification in Nonlinear Estimation.
#'   Sociological Methods & Research 33(1):54-87. DOI: 10.1177/0049124103262681
coef_sd.pk <- function(obj,
                       model = NULL,
                       method = NULL,
                       suppress.messages = TRUE){

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
  coefs <- coef(
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

  newdata <- obj$data %>%
    dplyr::select(!!!union(obj$data_group, req_vars),
                  !!!other_vars) %>%
    dplyr::mutate(data_sigma_group = factor(data_sigma_group)) %>%
    dplyr::group_by(!!!obj$data_group) %>%
    tidyr::nest(.key = "observations") %>%
    dplyr::ungroup()

  newdata <- tidyr::expand_grid(expand_grid(model, method),
                                newdata)

  newdata <- dplyr::left_join(coefs, newdata)


  newdata <- newdata %>%
    rowwise() %>%
    filter(!is.null(observations)) %>%
    mutate(model_fun = obj$stat_model[[model]]$conc_fun) %>%
    left_join(obj$prefit$par_DF %>%
                filter(optimize_param == FALSE,
                       use_param == TRUE) %>%
                dplyr::select(!!!obj$data_group, param_name, start)) %>%
    distinct() %>%
    mutate(const_pars = setNames(start, param_name)) %>%
    ungroup() %>%
    distinct()

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
    ungroup() %>%
    mutate(sds = map2(hessian_mat, coefs_vector, \(x, y) {
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
                   diag(chol(MASS::ginv(x),
                             pivot = TRUE))^(1/2),
                   error = function(err){
                     if (!suppress.messages) {
                       message(paste0("Pseudovariance matrix failed,",
                                      " returning NAs"))
                     }
                     rep(NA_real_, nrow(x))
                   })
               }) %>%
        set_names(nm = names(y))
    }))

  newdata <- newdata %>%
    dplyr::select(model, method, !!!obj$data_group, coefs_vector, sds)

  return(newdata)

}
