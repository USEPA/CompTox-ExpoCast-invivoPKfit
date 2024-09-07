#'Flat model
#'
#'Evaluates a "flat" model for concentration vs. time
#'
#'This function is used for model comparison: does a 1- or 2-compartment TK
#'model fit the data any better than this naive "flat" model?
#'
#'# Required parameters
#'
#'`params` must include the following named items:
#'   \describe{
#'   \item{Vdist}{Apparent volume of central compartment, volume/unit BW. Or see below for
#'   `Fgutabs_Vdist`}
#'   }
#'
#'For oral administration (if any `route %in% "oral"`), `params` must also
#'include:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   `Fgutabs_Vdist`}
#'   }
#'
#'For oral administration, in lieu of `Vdist` and `Fgutabs`, you may instead
#'provide `Fgutabs_Vdist`, the ratio of Fgutabs to Vdist (1/volume). This is an
#'alternate parameterization for situations where `Fgutabs` and `Vdist` are not
#'identifiable separately (i.e., when oral TK data are available, but IV data
#'are not). If `Fgutabs` and `Vdist` are provided, they will override any value
#'provided for `Fgutabs_Vdist`.
#'
#'If both oral and IV administration are specified (i.e., some `route %in% "iv"`
#'and some `route %in% "oral"`), then `Vdist` is required along with either
#'`Fgutabs` or `Fgutabs_Vdist`. (If `Vdist` and `Fgutabs_Vdist` are provided,
#'but `Fgutabs` is not provided, then `Fgutabs` will be calculated from `Vdist`
#'and `Fgutabs_Vdist`.)
#'
#'If `any(medium %in% 'blood')`, then `params` must also include
#'`Rblood2plasma`, the ratio of chemical concentration in whole blood to the
#'chemical concentration in blood plasma.
#'
#' # Flat model equations
#'
#' ## IV administration
#'
#' \deqn{\textrm{Conc} = \frac{\textrm{Dose}}{V_{\textrm{dist}}}}
#'
#' ## Oral administration
#'
#' \deqn{\textrm{Conc} = \frac{F_{\textrm{gutabs}} \textrm{Dose}}{V_{\textrm{dist}}}}
#'
#'@param params A named list of parameter values. See Details for requirements.
#'@param time A numeric vector of times, reflecting the time points
#'  when concentration is measured after the corresponding single bolus dose.
#'  Must be same length as other arguments, or length 1.
#'@param dose A numeric vector of doses, reflecting single bolus doses
#'  administered at time 0. Must be same length as other arguments, or
#'  length 1.
#'@param route A character vector, reflecting the route of administration of
#'  each single bolus dose: `'oral'` or `'iv'`.  Must be same length as `time`
#'  and `dose`, or length 1.
#'@param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as other arguments, or length 1.
#'@param loq A numeric vector of LOQ values. For any predicted value greater
#'  than zero but less than half the LOQ, that value is set to half the LOQ.
#'@return A vector of plasma concentration values (mass chemical/volume) corresponding to
#'  \code{time}.
#'
#'@author Caroline Ring, John Wambaugh, Chris Cook
#'
#'@export cp_flat
#' @family built-in model functions
#' @family flat model functions
#' @family model concentration functions

cp_flat <- function(params, time, dose, route, medium, loq) {

  params <- fill_params_flat(params)

  check_msg <- check_params_flat(params = params,
                            route = route,
                            medium = medium)

  if(!(check_msg %in% "Parameters OK")){
    stop(paste("cp_flat():",
               check_msg))
  }

  #for readability, assign params to variables inside this function
  list2env(as.list(params), envir = as.environment(-1))

  cp <-  dose * ifelse(route %in% "iv",
               1/Vdist,
               Fgutabs_Vdist)

  cp <- ifelse(medium %in% "blood",
               cp * Rblood2plasma,
               cp)
  # Any value greater than zero but less than LOQ will be set to LOQ/2
  cp <- ifelse(cp > 0 & cp < loq/2,
               loq/2,
               cp)

  return(cp)
}
