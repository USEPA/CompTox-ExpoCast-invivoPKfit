#' Toxicokinetic statistics for 1-compartment model for radiolabelling experiments
#'
#' Calculate predicted toxicokinetic statistics for a 1-compartment model.
#' This does use the parameters for the `model_1comp_rad` that are taken from
#' [httk] estimates.
#'
#' @section Statistics computed:
#'
#' \subsection{Total clearance}{
#' \deqn{\textrm{CL}_{tot} = k_{elim} + V_{dist}}
#' }
#' \subsection{Steady-state plasma concentration for long-term daily dose of 1 mg/kg/day}{
#' The dosing interval \eqn{\tau = \frac{1}{\textrm{day}}} will be converted to
#' the same units as \eqn{k_{elim}}.
#'
#' To convert to steady-state *blood* concentration, multiply by the
#' blood-to-plasma ratio.
#' \subsection{Oral route}{
#' \deqn{C_{ss} = \frac{F_{gutabs} V_{dist}}{k_{elim} \tau}}
#' }
#' \subsection{Intravenous route}{
#' \deqn{C_{ss} = \frac{1}{24 * \textrm{CL}_{tot}}}
#' }
#' }
#' \subsection{Half-life of elimination}{
#' \deqn{\textrm{Halflife} = \frac{\log(2)}{k_{elim}}}
#' }
#' \subsection{Time of peak concentration}{
#' For oral route:
#'
#' \deqn{\frac{\log \left( \frac{k_{gutabs}}{k_{elim}} \right)}{k_{gutabs} -
#' k_{elim}}}
#'
#' For intravenous route, time of peak concentration is always 0.
#' }
#' \subsection{Peak concentration}{
#' Evaluate [cp_1comp_rad()] at the time of peak concentration.
#' }
#' \subsection{AUC evaluated at infinite time}{
#' Evaluate [auc_1comp_rad()] at time = `Inf`.
#' }
#' \subsection{AUC evaluated at the time of the last observation}{
#' Evaluate [auc_1comp_rad()] at time = `tlast`.
#' }
#'
#' @inheritParams tkstats_1comp_cl
#' @return A `data.frame` with two variables:
#' - `param_name` = `c("Cltot", "Cltot/Fgutabs", "Css", "halflife", "tmax", "Cmax", "AUC_infinity")`
#' - `param_value` = The corresponding values for each statistic (which may be NA if that statistic could not be computed).
#' @export
#' @author John Wambaugh, Caroline Ring, Gilberto Padilla Mercado
tkstats_1comp_rad <- function(pars,
                              route,
                              medium,
                              dose,
                              time_unit,
                              conc_unit,
                              vol_unit,
                              restrictive = FALSE,
                              ...) {

  params <- fill_params_1comp_cl(pars)

  Q_totli = Q_gfr = Q_alv = Fup = Clint = NULL
  Kblood2air = Rblood2plasma = NULL
  kgutabs = Vdist = Fgutabs_Vdist = NULL

  # for readability, assign params to variables inside this function
  list2env(as.list(params), envir = as.environment(-1))


  # Set a Fup specific to the liver for clearance
  if (!restrictive) {
    Fup_hep <- 1
  } else {
    Fup_hep <- Fup
  }

  # compute total clearance
  Clhep <- (Q_totli * Fup_hep * Clint) / (Q_totli + (Fup_hep * Clint / Rblood2plasma))
  Clren <- Fup * Q_gfr
  Clair <- (Rblood2plasma * Q_alv / Kblood2air)
  Cltot <- Clren + Clhep + Clair
  Cltot_Fgutabs <- (Cltot / Vdist) / Fgutabs_Vdist

  # convert dose interval of (1/day) into time units
  # this is now standardized because time_units will always be hours
  dose_int <- 1 / 24

  Css <- dose * ifelse(route %in% "oral",
                       Fgutabs_Vdist / (Cltot / Vdist) / dose_int,
                       1 / (Cltot * dose_int)) *
    ifelse(medium %in% "blood",
           Rblood2plasma,
           1)

  halflife <- log(2) / (Cltot / Vdist)

  tmax <- ifelse(route %in% "oral",
                 log(kgutabs / (Cltot / Vdist)) / (kgutabs - (Cltot / Vdist)),
                 0)

  Cmax <- cp_1comp_rad(params = pars,
                       time = tmax,
                       dose = dose,
                       route = route,
                       medium = medium)

  AUC_inf <- auc_1comp_rad(params = pars,
                           time = Inf,
                           dose = dose,
                           route = route,
                           medium = medium)

  return(data.frame(param_name = c("Cltot",
                                   "Cltot/Fgutabs",
                                   "Css",
                                   "halflife",
                                   "tmax",
                                   "Cmax",
                                   "AUC_infinity",
                                   "Vss",
                                   "Vss/Fgutabs"
  ),
  param_value = c(Cltot,
                  Cltot_Fgutabs,
                  Css,
                  halflife,
                  tmax,
                  Cmax,
                  AUC_inf,
                  Vdist,
                  1 / (Fgutabs_Vdist)
  ),
  param_units = c(paste0(vol_unit, "/", time_unit), # Cltot
                  paste0(vol_unit, "/", time_unit), # Cltot/Fgutabs
                  conc_unit, # Css
                  time_unit, # halflife
                  time_unit, # tmax
                  conc_unit, # Cmax
                  paste0(conc_unit, " * ", time_unit), # AUC_inf
                  vol_unit,
                  vol_unit)
  ))

}
