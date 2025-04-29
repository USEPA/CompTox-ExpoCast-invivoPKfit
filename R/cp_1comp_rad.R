#' Analytical 1-compartment model for total radiolabeled material
#'
#' Calculates plasma concentrations vs. time according to the analytical solution
#' for the 1-compartment model, for single bolus doses (IV and/or oral).
#'
#' @section Required parameters:
#'
#' `params` must include the following named items:
#'   \describe{
#'   \item{Fup}{Fraction of compound unbound in plasma. Unitless.}
#'   \item{Clint}{Intrinsic clearance by hepatocytes. Units: 1/hr}
#'   \item{Q_totli}{Total blood flow through the liver. Units: L/h/kg body weight ^ (3/4)}
#'   \item{Q_gfr}{Glomerular filtration rate, how quickly do kidneys filter. Units: L/h/kg body weight ^ (3/4)}
#'   \item{Vdist}{Apparent volume of central compartment, volume/unit BW.}
#'   }
#'
#' For oral administration (if any `route %in% "oral"`), `params` must also
#' include:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction.}
#'   \item{kgutabs}{Rate of absorption from gut, 1/time.}
#'   }
#'
#' For excreta (if any `media %in% "excreta"`), `params` must also include:
#'   \describe{
#'   \item{Frec}{Fraction recovered in excreta. This is a precalculated value.}
#'   }
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
#'  length 1. In this model, it is expected that this value represents a measurement
#'  of radioactive particles from a radiolabeling experiment.
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
#' @author Caroline Ring, John Wambaugh, Gilberto Padilla Mercado
#'
#' @export cp_1comp_rad
#' @family built-in model functions
#' @family 1-compartment radiation model functions
#' @family model concentration functions
cp_1comp_rad <- function(params, time, dose, route, medium = 'plasma',
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
  # note: Qgfr is based on plasma wheras Q_totli is based on blood
  Clhep <- (Q_totli * Fup * Clint) / (Q_totli + (Fup * Clint / Rblood2plasma))
  Clren <- Fup * Q_gfr
  Cltot <- CLren + Clhep
  Cltot_Vdist <- Cltot / Vdist

  A_t <- rep(NA_real_, length(time))
  iv_vec <- (route == "iv")
  oral_vec <- (route == "oral")
  circ_vec <- (media %in% c("blood", "plasma"))
  excr_vec <- (media %in% c("excreta"))

  # Setup short names for indices
  ivc <- (iv_vec & circ_vec)
  ive <- (iv_vec & excr_vec)
  orc <- (oral_vec & circ_vec)
  ore <- (oral_vec & circ_vec)

  # For IV, CLtot_vdist == kgutabs has no bearing
  A_t[ivc] <- dose[ivc] * exp(-CLtot_Vdist * time[ivc])
  A_t[ive] <- dose[ive] * Frec * exp(-CLtot_Vdist * time[ive])

  # Check "static conditions" for oral absorption
  if (CLtot_Vdist != kgutabs & any(orc & ore)) {
    A_t[orc] <- dose[orc] * (Fgutabs * kgutabs) / (kgutabs - Cltot_Vdist) *
      (exp(-Cltot_Vdist * time[orc]) - exp(-kgutabs * time[orc]))

    A_t[ore] <- dose[ore] * (Fgutabs * Frec) / (Cltot_Vdist - kgutabs) *
      (Cltot_Vdist * (1 - exp(-kgutabs * time[ore])) -
         kgutabs * (1 - exp(-Cltot_Vdist * time[ore])))
  } else {
    A_t[orc] <- dose[orc] * (Fgutabs * CLtot_Vdist) *
      (exp(-Cltot_Vdist * time[orc]))

    A_t[ore] <- dose[ore] * (Fgutabs * Frec) *
      (1 - exp(-Cltot_Vdist * time[ore]) - time[ore] *
         Cltot_Vdist * exp(-Cltot_Vdist * time[ore]))
  }

  A_t <- ifelse(medium %in% "blood",
               Rblood2plasma * A_t,
               A_t)

  return(A_t)
}
