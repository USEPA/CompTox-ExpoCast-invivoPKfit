#' Time derivative of analytical 2-compartment model
#'
#' Calculates the time derivative (instantaneous rate of change) of plasma
#' concentration according to the analytical solution for the 2-compartment
#' model.
#'
#' This function is used by [postprocess_data()] to determine the time of peak
#' concentration for the 2-compartment model, by locating the point where the
#' time derivative of concentration crosses zero.
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
#'  For oral administration (`route %in% "oral"`), \code{params} must also
#'  include:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   "Fgutabs_V1"}
#'   \item{kgutabs}{Rate of absorption from gut, 1/h.}
#'   }
#'
#'  For oral administration, in lieu of "V1" and "Fgutabs", you may instead
#'  provide "Fgutabs_V1", the ratio of Fgutabs to V1 (1/L). This is an alternate
#'  parameterization for situations where "Fgutabs" and "V1" are not
#'  identifiable separately (i.e. when oral data are available, but IV data are
#'  not). If "Fgutabs" and "V1" are provided, then "Fgutabs_V1" will not be
#'  used.
#'
#' @param time A numeric vector of times in hours, reflecting the time points
#'  when concentration is measured after the corresponding single bolus dose.
#'  Must be same length as `dose` and `route`, or length 1.
#' @param dose A numeric vector of doses in mg/kg, reflecting single bolus doses
#'  administered at time 0. Must be same length as `time` and `route`, or length
#'  1.
#' @param route A character vector, reflecting the route of administration of
#'  each single bolus dose. Currently, only "iv" and "oral" are supported. Must
#'  be same length as `time` and `dose`, or length 1.
#' @param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as `time` and `dose`, or length 1.
#' @return A vector of instantaneous rates of change of plasma concentration
#'  values (mg/L/time) corresponding to each value in \code{time}
#' @author Caroline Ring, John Wambaugh
#' @export cp_2comp_dt
#' @family built-in model functions
#' @family 2-compartment model functions
cp_2comp_dt <- function(params, time, dose, route, medium) {
  # fill any missing parameters with NAs, and impute Fgutabs_v1 from Fgutabs and
  # V1 if necessary
  params <- fill_params_2comp(params)

  # check that required params are present
  check_msg <- check_params_2comp(params,
                                  route,
                                  medium)

  if (check_msg != "Parameters OK") {
    stop("cp_2comp_dt(): ", check_msg)
  }

  # get transformed parameters for 2-comp model
  trans_params <- transformed_params_2comp(params)

  A_iv_unit = B_iv_unit = B_oral_unit = A_oral_unit = NULL
  kgutabs = Rblood2plasma = NULL

  # for readability, assign params to variables inside this function
  list2env(as.list(params), envir = as.environment(-1))

  # for readability, assign transformed params to variables inside this function
  list2env(as.list(trans_params), envir = as.environment(-1))

  # get dcp/dt
  dcpdt <- dose * ifelse(route %in% "iv",
                 A_iv_unit * -alpha *
                    exp(-alpha * time) +
                   B_iv_unit * -beta *
                    exp(-beta * time),
                A_oral_unit * -alpha *
                   exp(-alpha * time) +
                   B_oral_unit * -beta *
                   exp(-beta * time) +
                   -(A_oral_unit + B_oral_unit) *
                   -kgutabs *
                   exp(-kgutabs * time)
                 )

  dcpdt <- ifelse(medium %in% "blood",
                  Rblood2plasma * dcpdt,
                  dcpdt)

  return(dcpdt)
}
