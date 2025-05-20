#' Bayesian information criterion
#'
#' Get the Bayesian information criterion (BIC) for a fitted `pk` object
#'
#' The BIC is calculated from the log-likelihood (LL) as follows:
#' \deqn{\textrm{BIC} = -2\textrm{LL} + \log(n_{obs}) n_{par}}
#'
#' where \eqn{n_{par}} is the number of parameters in the fitted model.
#'
#' Note that the BIC is just the AIC with \eqn{k = \log(n_{obs})}.
#'
#' @param object A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to compute
#'   log-likelihood. If NULL (the default), then BICs will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Time.Units`, `Dose`, `Route`,`Media`, `Conc`,
#'   `Detect`, `N_Subjects`. Before log-likelihood is calculated, `Time` will be
#'   transformed according to the transformation in `obj$scales$time` and `Conc`
#'   will be transformed according to the transformation in `obj$scales$conc`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   calculate BIC. If NULL (the default), log-likelihoods will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to calculate BICs. If NULL (the default),
#'   log-likelihoods will be returned for all of the methods in
#'   `obj$settings_optimx$method`.
#' @param exclude Logical: `TRUE` to compute the AIC after removing any
#'   observations in the data marked for exclusion (if there is a variable
#'   `exclude` in the data, an observation is marked for exclusion when `exclude
#'   %in% TRUE`). `FALSE` to include all observations, regardless of exclusion
#'   status. Default `TRUE`.
#' @param ... Additional arguments. Not in use.
#' @return A data.frame with log-likelihood values and calculated BIC using `newdata`.
#'   There is one row for each model in `obj`'s [stat_model()] element and
#'   each [optimx::optimx()] method (specified in [settings_optimx()]).
#' @family fit evaluation metrics
#' @family log likelihood functions
#' @family methods for fitted pk objects
#' @importFrom stats BIC
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
BIC.pk <- function(object,
                   newdata = NULL,
                   model = NULL,
                   method = NULL,
                   exclude = TRUE,
                   ...) {
  # ensure that the model has been fitted
  check <- check_required_status(obj = object,
                                 required_status = 5)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }
  if (is.null(model)) model <- names(object$stat_model)
  if (is.null(method)) method <- object$settings_optimx$method

  # get log-likelihoods
  ll <- logLik(object = object,
               newdata = newdata,
               model = model,
               method = method,
               negative = FALSE,
               force_finite = FALSE,
               exclude = exclude,
               drop_obs = FALSE)


  # obj$fit is a "long" data.frame with both parameters and sigma values
  params_df <- object$fit %>%
    dplyr::filter(optimize_param == TRUE) %>%
    group_by(!!!object$data_group, model, method) %>%
    dplyr::summarize(npar = dplyr::n())
  data_grp_vars <- sapply(object$data_group, rlang::as_label)

  # BIC requires knowing how many observations each group has
  ll <- ll %>% dplyr::rowwise() %>%
    dplyr::mutate(N_ROW = nrow(observations))

  # Combining log-likelihood table with parameters table
  ll <- ll %>%
    dplyr::select(!!!object$data_group,
                  model, method,
                  log_likelihood,
                  N_ROW) %>%
    dplyr::left_join(params_df, by = c(data_grp_vars, "model", "method"))

  # get number of parameters (excluding any constant, non-optimized parameters)

  BIC <- ll %>%
    dplyr::group_by(!!!object$data_group, model, method) %>%
    dplyr::mutate(BIC = log(.data$N_ROW) * .data$npar) - (2 * .data$log_likelihood)


  return(BIC)
}
