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
#' A first attempt is made to invert the Hessian using [solve()] (see
#' [hess_sd1()]). If the Hessian is singular, an attempt is made to calculate a
#' pseudovariance matrix, following the procedure outlined in Gill & King (2004)
#' (see [hess_sd2()]). First, the generalized inverse of the Hessian is calculated using
#' [MASS::ginv()]. Then, a generalized Cholesky decomposition (to ensure
#' positive-definiteness) is calculated using [Matrix::Cholesky] with argument
#' `perm = TRUE`. The generalized inverse is reconstructed from the generalized
#' Cholesky factorization. The square root of the diagonal elements of this
#' matrix represent the parameter standard deviations.
#'
#' If neither of these procedures is successful, then `NA_real_` is returned for
#' all coefficient standard deviations.
#'
#' @param object A [pk] object.
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
#'   The remaining columns include the parameters & hyperparameters as returned
#'   by [coef.pk()], as well as their calculated standard deviations.
#' @export
#' @import dplyr
#' @import purrr
#' @import tidyr
#' @author Caroline Ring and Gilberto Padilla Mercado
#' @family methods for fitted pk objects
#' @references Gill J, King G. (2004) What to Do When Your Hessian is Not
#'   Invertible: Alternatives to Model Respecification in Nonlinear Estimation.
#'   Sociological Methods & Research 33(1):54-87. DOI: 10.1177/0049124103262681
coef_sd.pk <- function(obj,
                       model = NULL,
                       method = NULL,
                       suppress.messages = TRUE, ...) {

  # ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  data_grp_vars <- sapply(obj$data_group, rlang::as_label)

  # get coefs data.frame for each model and method
  # but exclude non-optimized parameters and sigma values for error_group
  # coef does method & model checks

  #get optimized parameter vectors
  coefs_opt <- coef(
    obj = obj,
    model = model,
    method = method,
    drop_sigma = FALSE,
    include_type = "optim",
    suppress.messages = suppress.messages) %>%
    dplyr::rename(coefs_opt_vector = coefs_vector)

  #get constant parameter vectors
  coefs_const <- coef(
    obj = obj,
    model = model,
    method = method,
    drop_sigma = FALSE,
    include_type = "const",
    suppress.messages = suppress.messages) %>%
    dplyr::rename(coefs_const_vector = coefs_vector)

  #get all params used
  coefs_use <- coef(
    obj = obj,
    model = model,
    method = method,
    drop_sigma = FALSE,
    include_type = "use",
    suppress.messages = suppress.messages)

  coefs <- coefs_opt %>%
    dplyr::left_join(coefs_const,
                     by = c("model", "method",
                            data_grp_vars,
                            "Time.Units",
                            "Time_trans.Units")) %>%
    dplyr::left_join(coefs_use,
                     by = c("model", "method",
                            data_grp_vars,
                            "Time.Units",
                            "Time_trans.Units"))

  other_vars <- ggplot2::vars(
    Value,
    Value.Units,
    Time_trans.Units,
    Conc_trans,
    Conc_trans.Units,
    data_sigma_group,
    exclude
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
                            Detect,
                            LOQ,
                            pLOQ)

  # Convert Time_trans to hours
  newdata <- obj$data %>%
    dplyr::select(!!!union(obj$data_group, req_vars), !!!other_vars) %>%
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
    model_name = unname(sapply(obj$stat_model, \(x) {x$name})),
    model_fun = unname(sapply(obj$stat_model, \(x) {x$conc_fun}))
  )

  newdata <- dplyr::left_join(coefs, newdata, by = data_grp_vars) %>%
    dplyr::left_join(fun_models, join_by(model == model_name))


  newdata <- newdata %>%
    dplyr::rowwise() %>%
    dplyr::filter(!is.null(observations)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()


  sds_alerts <- newdata %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      sd_tbl = list(
        calc_sds_alerts(pars_opt = coefs_opt_vector,
                    pars_const = coefs_const_vector,
                    observations = observations,
                    modelfun = model_fun,
                    dose_norm = obj$scales$conc$dose_norm,
                    log10_trans = obj$scales$conc$log10_trans)
      )
    )  %>%
    tidyr::unnest(sd_tbl) %>%
    dplyr::ungroup() %>%
    dplyr::select(model, method, !!!obj$data_group,
                  param_name, param_sd, sd_alert)


  coefs_long <- newdata %>%
    dplyr::select(model, method, !!!obj$data_group, coefs_vector) %>%
    dplyr::group_by(model, method, !!!obj$data_group) %>%
    dplyr::mutate(coefs_tibble = purrr::map(coefs_vector, \(x) as.list(x) %>%
                                              as.data.frame() %>%
                                              tidyr::pivot_longer(
                                                cols = tidyselect::everything(),
                                                names_to = "param_name",
                                                values_to = "param_value"
                                              ))) %>%
    tidyr::unnest(coefs_tibble) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(coefs_vector))

  output <- coefs_long %>%
    dplyr::left_join(sds_alerts,
                     by = c("model", "method",
                                 data_grp_vars,
                            "param_name")) %>%
  dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(param_value = dplyr::na_if(param_value, NaN),
                  param_sd = dplyr::na_if(param_sd, NaN)) %>%
    dplyr::select(model, method, !!!obj$data_group,
                  param_name, param_value, param_sd, sd_alert)

  #add information about whether each parameter was optimized or not
  #to explain why some params don't have an SD
  param_type <- obj$fit %>%
    dplyr::select(model, method,
                  !!!obj$data_group,
                  param_name,
                  optimize_param,
                  use_param)

  output <- output %>%
    dplyr::left_join(param_type,
                     by = c("model", "method",
                            data_grp_vars,
                            "param_name"))


  return(output)

}
