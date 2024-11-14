#' Analytical 1-compartment model
#'
#' Calculates plasma concentrations vs. time according to the analytical solution
#' for the 1-compartment model, for single bolus doses (IV and/or oral).
#'
#'
#' # Required parameters
#'
#' `params` must include the following named items:
#'   \describe{
#'   \item{kelim}{Elimination rate, 1/time.}
#'   \item{Vdist}{Apparent volume of central compartment, volume/unit BW. Or see below for
#'   `Fgutabs_Vdist`}
#'   }
#'
#' For oral administration (if any `route %in% "oral"`), `params` must also
#' include:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   `Fgutabs_Vdist`}
#'   \item{kgutabs}{Rate of absorption from gut, 1/time.}
#'   }
#'
#' For oral administration, in lieu of `Vdist` and `Fgutabs`, you may instead
#' provide `Fgutabs_Vdist`, the ratio of Fgutabs to Vdist (1/volume). This is an
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
#' @param params A named numeric vector of model parameter values. See Details for
#'  requirements.
#' @param time A numeric vector of times, reflecting the time point when
#'  concentration is measured after the corresponding single bolus dose. Must be
#'  same length as `dose` and `route`, or length 1.
#' @param dose A numeric vector of doses, reflecting single bolus doses
#'  administered at time 0. Must be same length as `time` and `route`, or
#'  length 1.
#' @param route A character vector, reflecting the route of administration of
#'  each single bolus dose: `'oral'` or `'iv'`.  Must be same length as `time`
#'  and `dose`, or length 1.
#' @param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as `time` and `dose`, or length 1.
#' @return A vector of blood or plasma concentration values  corresponding
#'  to `time`.
#'
#' @author Caroline Ring, John Wambaugh
#'
#' @export cp_1comp
#' @family built-in model functions
#' @family 1-compartment model functions
#' @family model concentration functions
cp_1comp <- function(params,
                     time,
                     dose,
                     route,
                     medium = "plasma") {
  params <- fill_params_1comp(params)

  check_msg <- check_params_1comp(
    params = params,
    route = route,
    medium = medium
  )

  if (check_msg != "Parameters OK") {
    stop("cp_1comp(): ", check_msg)
  }

  list2env(as.list(params), envir = as.environment(-1))

  # compute plasma concentration
  cp <- dose * ifelse(route %in% "iv",
                      exp(-kelim * time) / Vdist,
                      # iv route
                      ifelse(
                        rep(
                          kelim != kgutabs, # oral route
                          length(time)
                        ),
                        # equation when kelim != kgutabs
                        (
                          (Fgutabs_Vdist * kgutabs)
                          * (exp(-kelim * time) - exp(-kgutabs * time))
                          / (kgutabs - kelim)
                        ),
                        # alternate equation when kelim == kgutabs
                        Fgutabs_Vdist * kelim * time * exp(-kelim * time)
                      )
  )

  cp <- ifelse(medium %in% "blood", Rblood2plasma * cp, cp)

  return(cp)
}
