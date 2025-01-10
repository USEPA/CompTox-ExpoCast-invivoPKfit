#' Analytical 1-compartment model with specific clearance
#'
#' Calculates plasma concentrations vs. time according to the analytical solution
#' for the 1-compartment model, for single bolus doses (IV and/or oral).
#'
#'
#' # Required parameters
#'
#' `params` must include the following named items:
#'   \describe{
#'   \item{Fup}{Fraction of compound unbound in plasma. Unitless.}
#'   \item{Clint}{Intrinsic clearance by hepatocytes. Units: 1/hr}
#'   \item{Q_totli}{Total blood flow through the liver. Units: L/h/kg body weight ^ (3/4)}
#'   \item{Q_gfr}{Glomerular filtration rate, how quickly do kidneys filter. Units: L/h/kg body weight ^ (3/4)}
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
#'  same length as `dose` and `iv.dose`, or length 1.
#' @param dose A numeric vector of doses, reflecting single bolus doses
#'  administered at time 0. Must be same length as `time` and `iv.dose`, or
#'  length 1.
#' @param route A character vector, reflecting the route of administration of
#'  each single bolus dose: `'oral'` or `'iv'`.  Must be same length as `time`
#'  and `dose`, or length 1.
#' @param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as `time` and `dose`, or length 1.
#' @param restrictive A logical value (TRUE or FALSE. Default: FALSE) that says whether the
#' assumption is that the clearance is restrictive or non-restrictive
#'
#' @return A vector of blood or plasma concentration values  corresponding
#'  to `time`.
#'
#' @author Caroline Ring, John Wambaugh
#'
#' @export cp_1comp
#' @family built-in model functions
#' @family 1-compartment model functions
#' @family model concentration functions
cp_1comp_cl <- function(params, time, dose, route, medium = 'plasma',
                        restrictive = FALSE) {
  params <- fill_params_1comp_cl(params)

  check_msg <- check_params_1comp_cl(params = params,
                                     route = route,
                                     medium = medium)

  if (check_msg != "Parameters OK") {
    stop("cp_1comp(): ", check_msg)
  }


  list2env(as.list(params), envir = as.environment(-1))

  # compute total clearance
  Clhep <- (Q_totli * Fup * Clint) / (Q_totli + (Fup * Clint / Rblood2plasma))
  Cltot <- (Fup * Q_gfr) + Clhep
  Cltot_Vdist <- Cltot / Vdist

  # compute plasma concentration
  cp <- dose * ifelse(
    route %in% "iv",
    exp(-Cltot_Vdist * time) / Vdist,
    # iv route
    ifelse(
      rep(Cltot_Vdist != kgutabs, # oral route
          length(route)),
      # equation when kelim != kgutabs
      (Fgutabs_Vdist *
         kgutabs) /
        (kgutabs - Cltot_Vdist) *
        (exp(-Cltot_Vdist * time) -
           exp(-kgutabs * time)),
      # alternate equation when kelim == kgutabs
      Fgutabs_Vdist * Cltot_Vdist *
        time *
        exp(-Cltot_Vdist * time)
    )
  )

  cp <- ifelse(medium %in% "blood",
               Rblood2plasma * cp,
               cp)

  return(cp)
}
