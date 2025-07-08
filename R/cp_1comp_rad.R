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
cp_1comp_rad <- function(params, time, dose, route, medium = "plasma",
                        restrictive = TRUE) {

  params <- fill_params_1comp_cl(params)

  check_msg <- check_params_1comp_cl(params = params,
                                     route = route,
                                     medium = medium)

  if (check_msg != "Parameters OK") {
    stop("cp_1comp(): ", check_msg)
  }

  Q_totli = Q_gfr = Q_alv = Fup = Clint = Kblood2air = NULL
  Rblood2plasma = kelim = kgutabs = Fgutabs = Frec = Vdist = NULL


  list2env(as.list(params), envir = as.environment(-1))

  # compute total clearance
  # note: Qgfr is based on plasma wheras Q_totli is based on blood
  Clint_hep <- Clint * (6.6) / 1E6 # Convert to L/hr
  Clhep <- (Q_totli * Fup * Clint_hep) / (Q_totli + (Fup * Clint_hep / Rblood2plasma))
  Clren <- Fup * Q_gfr
  Clair <- (Rblood2plasma * Q_alv / Kblood2air)
  Cltot <- Clren + Clhep + Clair
  kelim <- Cltot / Vdist # Cltot is L/hr and Vdist is in L

  A_t <- rep(NA_real_, length(time))
  iv_vec <- (route == "iv")
  oral_vec <- (route == "oral")
  circ_vec <- (medium %in% c("blood", "plasma"))
  excr_vec <- (medium %in% c("excreta"))

  # Setup short names for indices
  ivc <- (iv_vec & circ_vec)
  ive <- (iv_vec & excr_vec)
  orc <- (oral_vec & circ_vec)
  ore <- (oral_vec & excr_vec)

  # For IV, kelim == kgutabs has no bearing
  A_t[ivc] <- dose * exp(-kelim * time[ivc])
  A_t[ive] <- dose * Frec * (1 - exp(-kelim * time[ive]))

  # Check "static conditions" for oral absorption
  if (kelim != kgutabs && any(orc & ore)) {
    A_t[orc] <- dose * (Fgutabs * kgutabs) / ((kgutabs - kelim)) *
      (exp(-kelim * time[orc]) - exp(-kgutabs * time[orc]))

    A_t[ore] <- dose * (Fgutabs * Frec) / ((kelim - kgutabs)) *
      (kelim * (1 - exp(-kgutabs * time[ore])) -
         kgutabs * (1 - exp(-kelim * time[ore])))
  } else {
    A_t[orc] <- dose * (Fgutabs * kelim * time[orc]) *
      (exp(-kelim * time[orc]))

    A_t[ore] <- dose * (Fgutabs * Frec) *
      (1 - exp(-kelim * time[ore]) - time[ore] *
         kelim * exp(-kelim * time[ore]))
  }

  A_t <- ifelse(medium %in% "blood",
               Rblood2plasma * A_t,
               A_t)

  return(A_t)
}
