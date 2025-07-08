#' Analytical 2-compartment model
#'
#' Calculates plasma concentration according to the analytical solution for the
#' 2-compartment model.
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
#' For oral administration, in lieu of `V1` and `Fgutabs`, you may instead
#' provide `Fgutabs_V1`, the ratio of Fgutabs to V1 (1/volume). This is an alternate
#' parameterization for situations where `Fgutabs` and `V1` are not identifiable
#' separately (i.e. when oral TK data are available, but IV data are not). If
#' `Fgutabs` and `V1` are provided, they will override any value also provided
#' for `Fgutabs_V1`.
#'
#' If both oral and IV administration are specified (i.e., some `route %in% "iv"`
#' and some `route %in% "oral"`), then `V1` is required along with either
#' `Fgutabs` or `Fgutabs_V1`. (If `V1` and `Fgutabs_V1` are provided, but
#' `Fgutabs` is not provided, then `Fgutabs` will be calculated from `V1` and
#' `Fgutabs_V1`.)
#'
#' #'If `any(medium %in% 'blood')`, then `params` must also include
#' `Rblood2plasma`, the ratio of chemical concentration in whole blood to the
#' chemical concentration in blood plasma.
#'
#' @param params A named numeric vector of parameter values. See Details for requirements.
#' @param time A numeric vector of times, reflecting the time points
#'  when concentration is measured after the corresponding single bolus dose.
#'  Must be same length as other arguments, or length 1.
#' @param dose A numeric vector of doses, reflecting single bolus doses
#'  administered at time 0. Must be same length as other arguments, or
#'  length 1.
#' @param route A character vector, reflecting the route of administration of
#'  each single bolus dose: `'oral'` or `'iv'`.  Must be same length as `time`
#'  and `dose`, or length 1.
#' @param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as other arguments, or length 1.
#' @return A vector of blood or plasma concentration values (mass chemical/volume media) corresponding to each
#'  value in \code{time}
#' @export cp_2comp
#'
#' @author Caroline Ring, John Wambaugh
#' @family built-in model functions
#' @family 2-compartment model functions
#' @family model concentration functions
cp_2comp <- function(params,
                     time,
                     dose,
                     route,
                     medium = "plasma") {
  # fill any missing parameters with NAs, and impute Fgutabs_v1 from Fgutabs and
  # V1 if necessary
  params <- fill_params_2comp(params)

  # check that required params are present
  check_msg <- check_params_2comp(
    params,
    route,
    medium
  )
  if (check_msg != "Parameters OK") {
    stop("cp_2comp(): ", check_msg)
  }

  # get transformed parameters for 2-comp model
  trans_params <- transformed_params_2comp(params)

  A_iv_unit = B_iv_unit = alpha = beta = A_oral_unit = B_oral_unit = NULL
  kgutabs = Rblood2plasma = NULL

  list2env(as.list(params), envir = as.environment(-1))

  # for readability, assign transformed params to variables inside this function
  list2env(as.list(trans_params), envir = as.environment(-1))

  Cp <- rep(NA_real_, length(time))
  iv_vec <- (route == "iv")
  or_vec <- (route == "oral")

  # get predicted concentration
  if (any(iv_vec)) {
    Cp[iv_vec] <- dose[iv_vec] * (
      A_iv_unit * exp(-alpha * time[iv_vec]) +
        B_iv_unit * exp(-beta * time[iv_vec])
    )
  }

  if (any(or_vec)) {
    Cp[or_vec] <- dose[or_vec] * (
      (A_oral_unit * exp(-alpha * time[or_vec])) +
        (B_oral_unit * exp(-beta * time[or_vec])) -
        ((A_oral_unit + B_oral_unit) * exp(-kgutabs * time[or_vec]))
    )
  }

  Cp <- ifelse(medium %in% "blood", Cp * Rblood2plasma, Cp)

  return(Cp)
}
