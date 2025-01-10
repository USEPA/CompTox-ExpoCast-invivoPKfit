#' AUC for flat model
#'
#' Evaluates the area under the concentration-time curve for a "flat" model
#'
#' # Required parameters
#'
#' `params` must include the following named items:
#'   \describe{
#'   \item{Vdist}{Apparent volume of central compartment, L/kg BW. Or see below for
#'   `Fgutabs_Vdist`}
#'   }
#'
#' For oral administration (if any `route %in% "oral"`), `params` must also
#' include:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   `Fgutabs_Vdist`}
#'   }
#'
#' For oral administration, in lieu of `Vdist` and `Fgutabs`, you may instead
#' provide `Fgutabs_Vdist`, the ratio of Fgutabs to Vdist (1/L). This is an
#' alternate parameterization for situations where `Fgutabs` and `Vdist` are not
#' identifiable separately (i.e., when oral TK data are available, but IV data
#' are not). If `Fgutabs` and `Vdist` are provided, they will override any value
#' provided for `Fgutabs_Vdist`.
#'
#' If both oral and IV administration are specified (i.e., some `route %in% "iv"`
#' and some `route %in% "oral"`), then `Vdist` is required along with either
#' `Fgutabs` or `Fgutabs_Vdist`. (If `Vdist` and `Fgutabs_Vdist` are provided,
#' but `Fgutabs` is not provided, then `Fgutabs` will be calculated from `Vdist`
#' and `Fgutabs_Vdist`.)
#'
#' If `any(medium %in% 'blood')`, then `params` must also include
#' `Rblood2plasma`, the ratio of chemical concentration in whole blood to the
#' chemical concentration in blood plasma.
#'
#' @param params A named list of model parameter values. See Details for requirements.
#' @param time A numeric vector of times in hours.
#' @param dose A numeric vector of doses in mg/kg
#' @param route A logical vector: TRUE for single IV bolus dose; FALSE for
#'  single oral dose. Not used, but must be present for compatibility with other
#'  model functions.
#' @param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as other arguments, or length 1.
#'
#' @return A vector of plasma concentration values (mg/L) corresponding to
#'  \code{time}.
#'
#' @author Caroline Ring, John Wambaugh, Chris Cook
#' @export auc_flat
#' @family built-in model functions
#' @family flat model functions
#' @family model AUC functions
auc_flat <- function(time, params, dose, route, medium) {

  params <- fill_params_flat(params)

  check_msg <- check_params_flat(params = params,
                           route = route,
                           medium = medium)

  if (check_msg != "Parameters OK") {
    stop("cp_flat(): ", check_msg)
  }

  # for readability, assign params to variables inside this function
  list2env(as.list(params), envir = as.environment(-1))


  auc <- dose * ifelse(route %in% "iv",
                       time / Vdist,
                       time * Fgutabs_Vdist)

  auc <- ifelse(medium %in% "blood",
                Rblood2plasma * auc,
                auc)

  return(auc)
}
