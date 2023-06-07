#' Analytical AUC for the 2-compartment model
#'
#' Calculate area under the plasma concentration vs. time curve for the
#' 1-compartment model, using an analytical equation (the integral of the
#' 1-compartment model equation with respect to time).
#'
#' @param params A named list of parameter values including the following:
#'   \describe{
#'   \item{k12}{Rate at which compound moves from central to peripheral
#'   compartment, 1/h.}
#'   \item{k21}{Rate at which compound moves from peripheral to central
#'   compartment, 1/h.}
#'   \item{kelim}{Elimination rate, 1/h.}
#'   \item{V1}{Apparent volume of central compartment, L/kg BW. Or see below for
#'   "Fgutabs_V1"}
#'   }
#'
#'   For oral administration (\code{route} FALSE), \code{params} must also
#'   include:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   "Fgutabs_V1"}
#'   \item{kgutabs}{Rate of absorption from gut, 1/h.}
#'   }
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
#'@return A vector of plasma AUC values, evaluated at each time point in `time`.
#' @export auc_2comp
#' @author Caroline Ring, John Wambaugh
#' @family built-in model functions
#' @family 2-compartment model functions
#' @family model AUC functions
auc_2comp <- function(params, time, dose, route, medium = "plasma")
{

  params <- fill_params_2comp(params)

  #check that required params are present
  check_msg <- check_params_2comp(params,
                                  route,
                                  medium)

  if(!(check_msg %in% "Parameters OK")){
    stop(paste("auc_2comp():",
               check_msg))
  }

  #get transformed parameters for 2-comp model
  trans_params <- transformed_params_2comp(params,
                                           time,
                                           dose,
                                           route,
                                           medium)

  auc <- dose * ifelse(route %in% "iv",
                trans_params$A_iv_unit/trans_params$alpha -
                  trans_params$A_iv_unit*exp(-time*trans_params$alpha)/
                  trans_params$alpha +
                  trans_params$B_iv_unit/trans_params$beta -
                  trans_params$B_iv_unit*exp(-time*trans_params$beta)/
                  trans_params$beta,
                trans_params$A_oral_unit/trans_params$alpha -
                  trans_params$A_oral_unit*exp(-time*trans_params$alpha)/
                  trans_params$alpha +
                  trans_params$B_oral_unit/trans_params$beta -
                  trans_params$B_oral_unit*exp(-time*trans_params$beta)/
                  trans_params$beta +
                  (-trans_params$A_oral_unit -trans_params$B_oral_unit)/
                  params$kgutabs -
                  (-trans_params$A_oral_unit -trans_params$B_oral_unit)*
                  exp(-time*params$kgutabs)/params$kgutabs
                )
auc <- ifelse(medium %in% "blood",
              params$Rblood2plasma * auc,
              auc
              )

  return(auc)
}
