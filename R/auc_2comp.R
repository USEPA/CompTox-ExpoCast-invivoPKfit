#' Analytical AUC for the 2-compartment model
#'
#' Calculate area under the plasma concentration vs. time curve for the
#' 2-compartment model, using an analytical equation (the integral of the
#' 2-compartment model equation with respect to time).
#'
#' @section Required params:
#'
#' `params` must include the following named items:
#'   \describe{
#'   \item{k12}{Rate at which the compound moves from the central to peripheral compartment, 1/h.}
#'   \item{k21}{Rate at which the compound moves from peripheral to central compartment, 1/h.}
#'   \item{kelim}{Elimination rate, 1/h.}
#'   \item{V1}{Apparent volume of central compartment, L/kg BW.}
#'   }
#'
#' For oral administration (\code{route} FALSE), \code{params} must also include:
#' \describe{
#' \item{Fgutabs}{Oral bioavailability, unitless fraction.}
#' \item{kgutabs}{rate of absorption from gut, 1/h.}
#' }
#'
#' For oral administration, in lieu of "V1" and "Fgutabs", you may instead
#' provide "Fgutabs_V1", the ratio of Fgutabs to V1 (1/L). This is an
#' alternate parameterization for situations where "Fgutabs" and "V1" are not
#' identifiable separately (i.e. when oral data are available, but IV data are
#' not). If "Fgutabs" and "V1" are provided, then "Fgutabs_V1" will not be
#' used.
#'
#' @param params A named list of parameter values.
#'
#'
#' @param time A numeric vector of time values, in hours
#' @param dose A numeric vector of doses in mg/kg
#' @param route A logical vector: TRUE for single IV bolus dose, FALSE for single oral dose
#' @param medium A character string that determines the measured media. Default: "plasma".
#' @return A vector of plasma AUC values, evaluated at each time point in `time`.
#' @export auc_2comp
#' @author Caroline Ring, John Wambaugh
#' @family built-in model functions
#' @family 2-compartment model functions
#' @family model AUC functions
auc_2comp <- function(params, time, dose, route, medium = "plasma") {
  # fill any missing parameters with NAs, and impute Fgutabs_v1 from Fgutabs and
  # V1 if necessary
  params <- fill_params_2comp(params)

  # check that required params are present
  check_msg <- check_params_2comp(params,
                                  route,
                                  medium)

  if (check_msg != "Parameters OK") {
    cli::cli_abort(check_msg)
  }

  # get transformed parameters for 2-comp model
  trans_params <- transformed_params_2comp(params)

  # params used (assigned NULL to prevent global variable hiccup)
  A_iv_unit = B_iv_unit = NULL
  A_oral_unit = B_oral_unit = NULL
  alpha = beta = NULL
  Rblood2plasma = kgutabs =  NULL

  list2env(as.list(params), envir = as.environment(-1))

  # for readability, assign transformed params to variables inside this function
  list2env(as.list(trans_params), envir = as.environment(-1))

  auc <- rep(NA_real_, length(time))
  iv_vec <- (route == "iv")
  or_vec <- (route == "oral")
  blood_vec <- (medium == "blood")

  auc[iv_vec] <- dose[iv_vec] * (
    (A_iv_unit / alpha) - (A_iv_unit * exp(-time * alpha) / alpha) +
      (B_iv_unit / beta) - (B_iv_unit * exp(-time * beta) / beta)
  )

  auc[or_vec] <- dose[or_vec] * (
    (A_oral_unit / alpha) - (A_oral_unit * exp(-time[or_vec] * alpha) / alpha) +
      (B_oral_unit / beta) - (B_oral_unit * exp(-time[or_vec] * beta) / beta) +
      ((-A_oral_unit - B_oral_unit) / kgutabs) -
      ((-A_oral_unit - B_oral_unit) * exp(-time[or_vec] * kgutabs) / kgutabs)
  )

  auc[blood_vec] <- auc[blood_vec] * Rblood2plasma

  return(auc)
}
