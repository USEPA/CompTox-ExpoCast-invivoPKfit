#'Get Hessian matrixes
#'
#'Extract Hessian matrixes from a fitted `pk` object
#'
#'This function computes a numerical approximation to the model Hessian for each
#'data group and each model in a fitted `pk` object. The Hessian is the matrix
#'of second derivatives of the model objective function with respect to each
#'model parameter. Here, the objective function is the negative log-likelihood
#'implemented in [log_likelihood()], evaluated jointly across the data that was
#'used to fit the model.
#'
#'
#'@param obj A [pk] object
#'@param model Optional: Specify one or more of the fitted models whose
#'  coefficients to return. If NULL (the default), coefficients will be returned
#'  for all of the models in `obj$stat_model`.
#'@param method Optional: Specify one or more of the [optimx::optimx()] methods
#'  whose coefficients to return. If NULL (the default), coefficients will be
#'  returned for all of the models in `obj$settings_optimx$method`.
#'@param suppress.messages Logical. `TRUE` (the default) to suppress informative
#'  messages. `FALSE` to see them.
#'@param ... Additional arguments. Not in use right now.
#'@return A dataframe with one row for each `data_group`, `model` and `method`.
#'  The remaining column is a `list` column containing the Hessian for each row.
#'@export
#'@import dplyr
#'@import purrr
#'@import tidyr
#'@import numDeriv
#'@author Caroline Ring, Gilberto Padilla Mercado
#'@family methods for fitted pk objects
#'@references Gill J, King G. (2004) What to Do When Your Hessian is Not
#'  Invertible: Alternatives to Model Respecification in Nonlinear Estimation.
#'  Sociological Methods & Research 33(1):54-87. DOI: 10.1177/0049124103262681
get_hessian.pk <- function(obj,
                           model = NULL,
                           method = NULL,
                           suppress.messages = TRUE, ...) {

  # ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  # get coefs data.frame for each model and method
  # but exclude non-optimized parameters and sigma values for error_group
  # coef does method & model checks

  #get optimized parameter vectors
  coefs_opt <- coef(
    object = obj,
    model = model,
    method = method,
    drop_sigma = FALSE,
    include_type = "optim",
    suppress.messages = suppress.messages) %>%
    dplyr::rename(coefs_opt_vector = coefs_vector)

  #get constant parameter vectors
  coefs_const <- coef(
    object = obj,
    model = model,
    method = method,
    drop_sigma = FALSE,
    include_type = "const",
    suppress.messages = suppress.messages) %>%
    dplyr::rename(coefs_const_vector = coefs_vector)

  #get all params used
  coefs_use <- coef(
    object = obj,
    model = model,
    method = method,
    drop_sigma = FALSE,
    include_type = "use",
    suppress.messages = suppress.messages)

  coefs <- coefs_opt %>%
    dplyr::left_join(coefs_const) %>%
    dplyr::left_join(coefs_use) %>%
    suppressMessages()

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
  fun_models <- get_stat_model(obj)

  newdata <- dplyr::left_join(coefs, newdata, by = data_grp_vars) %>%
    dplyr::left_join(fun_models, join_by(model))


  newdata <- newdata %>%
    dplyr::rowwise() %>%
    dplyr::filter(!is.null(observations)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()


  hess <- newdata %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      hessian = list(
        calc_hessian(pars_opt = coefs_opt_vector,
                        pars_const = coefs_const_vector,
                        observations = observations,
                        modelfun = modelfun,
                        dose_norm = obj$scales$conc$dose_norm,
                        log10_trans = obj$scales$conc$log10_trans)
      )
    )  %>%
    dplyr::ungroup() %>%
    dplyr::select(model, method, !!!obj$data_group,
                  hessian)


  return(hess)

}
