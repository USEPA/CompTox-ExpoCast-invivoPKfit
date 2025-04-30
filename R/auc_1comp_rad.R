#' Analytic AUC for 1-compartment model with total radiolabeled compound
#'
#' Calculate area under the plasma concentration vs. time curve for the
#' 1-compartment model, using an analytical equation (the integral of the
#' 1-compartment model equation with respect to time).
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
#'
#'
#' @return A vector of plasma AUC values (concentration*time) corresponding to `time`.
#'
#' @author Caroline Ring, John Wambaugh, Gilberto Padilla Mercado
#' @export auc_1comp_cl
#' @family built-in model functions
#' @family 1-compartment model functions
#' @family model AUC functions

auc_1comp_rad <- function(params,
                         time,
                         dose,
                         route,
                         medium) {

  params <- fill_params_1comp_cl(params)

  check_msg <- check_params_1comp_cl(params = params,
                                     route = route,
                                     medium = medium)

  if (check_msg != "Parameters OK") {
    stop("cp_1comp_cl(): ", check_msg)
  }

  # for readability, assign params to variables inside this function
  list2env(as.list(params), envir = as.environment(-1))

  # note: Qgfr is based on plasma wheras Q_totli is based on blood
  Clhep <- (Q_totli * Fup * Clint) / (Q_totli + (Fup * Clint / Rblood2plasma))
  Clren <- Fup * Q_gfr
  Cltot <- Clren + Clhep
  kelim <- Cltot / Vdist

  auc <- rep(NA_real_, length(time))
  iv_vec <- (route == "iv")
  oral_vec <- (route == "oral")
  circ_vec <- (medium %in% c("blood", "plasma"))
  excr_vec <- (medium %in% c("excreta"))

  # Setup short names for indices
  ivc <- (iv_vec & circ_vec)
  ive <- (iv_vec & excr_vec)
  orc <- (oral_vec & circ_vec)
  ore <- (oral_vec & excr_vec)

  # IV
  auc[ivc] <- dose * (1 / (kelim) - exp(-time * kelim) / Cltot)

  # ORAL
  # kelim != kgutabs
  if (kelim != kgutabs) {
    auc[orc] <-  dose *
      (-1 * Fgutabs_Vdist * kgutabs * (1 / kgutabs - 1 / kelim) /
      ((-kelim + kgutabs)) +
      Fgutabs_Vdist * kgutabs * (exp(-time * kgutabs) / kgutabs -
         exp(-time * kelim) / kelim) /
      ((-kelim + kgutabs))
      )
  } else { # kelim == kgutabs
    auc[orc] <- dose *
      (Fgutabs_Vdist / (kelim) +
         (-Fgutabs_Vdist * time * kelim - Fgutabs) * exp(-time * kelim) / (kelim)
      )
  }

  # excreta are transformed to cumulative values therefore
  # AUC is simply proportional to `Frec` and `Fgutabs` and `Dose`.
  auc[ive] <- dose * Frec
  auc[ore] <- dose * Frec * Fgutabs

  auc <- ifelse(medium %in% "blood",
                Rblood2plasma * auc,
                auc)

  return(auc)
}
