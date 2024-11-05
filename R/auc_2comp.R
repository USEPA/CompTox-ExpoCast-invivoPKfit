#' Analytical AUC for the 2-compartment model
#'
#' Calculate area under the plasma concentration vs. time curve for the
#' 1-compartment model, using an analytical equation (the integral of the
#' 1-compartment model equation with respect to time).
#'
#' @md
#' @param params A named list of parameter values including the following:
#'   * k12: Rate at which compound moves from central to peripheral
#'   compartment, 1/h.
#'   * k21: Rate at which compound moves from peripheral to central compartment, 1/h.
#'   * kelim: Elimination rate, 1/h.
#'   * V1: Apparent volume of central compartment, L/kg BW. Or see below for "Fgutabs_V1"
#'
#'   For oral administration (\code{route} FALSE), \code{params} must also include:
#'   * Fgutabs: Oral bioavailability, unitless fraction. Or see below for "Fgutabs_V1"
#'   * kgutabs: Rate of absorption from gut, 1/h.
#'
#'   For oral administration, in lieu of "V1" and "Fgutabs", you may instead
#'   provide "Fgutabs_V1", the ratio of Fgutabs to V1 (1/L). This is an
#'   alternate parameterization for situations where "Fgutabs" and "V1" are not
#'   identifiable separately (i.e. when oral data are available, but IV data are
#'   not). If "Fgutabs" and "V1" are provided, then "Fgutabs_V1" will not be
#'   used.
#'
#' @param time A numeric vector of time values, in hours
#' @param dose A numeric vector of doses in mg/kg
#' @param route A logical vector: TRUE for single IV bolus dose, FALSE for single oral dose
#' @param medium A character string that determines the measured media. Default: "plasma".
#'@return A vector of plasma AUC values, evaluated at each time point in `time`.
#' @export auc_2comp
#' @author Caroline Ring, John Wambaugh
#' @family built-in model functions
#' @family 2-compartment model functions
#' @family model AUC functions
auc_2comp <- function(params, time, dose, route, medium = "plasma")
{
  #fill any missing parameters with NAs, and impute Fgutabs_v1 from Fgutabs and
  #V1 if necessary
  params <- fill_params_2comp(params)

  #check that required params are present
  check_msg <- check_params_2comp(params,
                                  route,
                                  medium)

  if (check_msg != "Parameters OK") {
    stop("auc_2comp():", check_msg)
  }

  #get transformed parameters for 2-comp model
  trans_params <- transformed_params_2comp(params)

  #for readability, assign params to variables inside this function
  for(x in names(params)){
    assign(x, unname(params[x]))
  }

  #for readability, assign transformed params to variables inside this function
  for(x in names(trans_params)){
    assign(x, unname(trans_params[x]))
  }

  auc <- dose * ifelse(route %in% "iv",
                A_iv_unit/alpha -
                  A_iv_unit*exp(-time*alpha)/
                  alpha +
                  B_iv_unit/beta -
                  B_iv_unit*exp(-time*beta)/
                  beta,
                A_oral_unit/alpha -
                  A_oral_unit*exp(-time*alpha)/
                  alpha +
                  B_oral_unit/beta -
                  B_oral_unit*exp(-time*beta)/
                  beta +
                  (-A_oral_unit -B_oral_unit)/
                  kgutabs -
                  (-A_oral_unit -B_oral_unit)*
                  exp(-time*kgutabs)/kgutabs
                )
auc <- ifelse(medium %in% "blood",
              Rblood2plasma * auc,
              auc
              )

  return(auc)
}
