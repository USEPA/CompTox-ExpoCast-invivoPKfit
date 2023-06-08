#' Pre-fitting
#'
#' Do pre-fit calculations and checks
#'
#' This function does the following:
#'
#' - Based on the error model in `stat_error_model` and the pre-processed data, determines the number of residual standard deviations ("sigmas") hyperparameters to be estimated.
#' - Determines which "sigma" hyperparameter corresponds to each observation in the data.
#' - Calculates lower/upper bounds and starting guesses for each "sigma" hyperparameter
#' - For each model in `stat_model`, calls its `params_fun`, the function that, based on the data, determines whether to optimize each model parameter, and calculates lower/upper bounds and starting guesses for each model parameter to be optimized. Only non-excluded observations are passed to each model's `params_fun`.
#'
#'
#' Lower bounds for each "sigma" hyperparameter are set to
#' `sqrt(.Machine$double_eps)`.
#'
#' Upper bounds for each "sigma" hyperparameter are calculated as the standard
#' deviation of observations in the corresponding error SD group (see
#' [combined_sd()]). If the combined SD is non-finite or less than the sigma
#' lower bound, then the combined SD of all non-excluded data is substituted. If
#' that is still non-finite or less then the sigma lower bound, then a constant
#' value of 100 is substituted.
#'
#' The starting guess for each "sigma" hyperparameter is one-tenth of the upper
#' bound.
#'
#' @param obj A `pk_faceted` object
#' @return The same `pk_faceted` object, with list column `pk_object` modified
#'   by applying [prefit()] to each item
#' @export
#' @author Caroline Ring
prefit.pk_faceted <- function(obj){

  obj <- obj %>%
    dplyr::transmute(pk_object = purrr::map(pk_object,
                                            prefit)
    )

  return(obj)

}
