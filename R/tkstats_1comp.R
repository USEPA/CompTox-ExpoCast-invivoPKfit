#'Toxicokinetic statistics for 1-compartment model
#'
#'Calculate predicted toxicokinetic statistics for a 1-compartment model.
#'
#'# Statistics computed
#'
#'## Total clearance
#'
#'\deqn{\textrm{CL}_{tot} = k_{elim} + V_{dist}}
#'
#'## Steady-state plasma concentration for long-term daily dose of 1 mg/kg/day
#'
#'The dosing interval \eqn{\tau = \frac{1}{\textrm{day}}} will be converted to the same units as \eqn{k_{elim}}.
#'
#' To convert to steady-state *blood* concentration, multiply by the blood-to-plasma ratio.
#'
#'### Oral route
#'
#'\deqn{C_{ss} = \frac{F_{gutabs} V_{dist}}{k_{elim} \tau}}
#'
#'### Intravenous route
#'
#'\deqn{C_{ss} = \frac{1}{24 * \textrm{CL}_{tot}}}
#'
#' ## Half-life of elimination
#'
#' \deqn{\textrm{Halflife} = \frac{\log(2)}{k_{elim}}}
#'
#' ## Time of peak concentration
#'
#' For oral route:
#'
#' \deqn{\frac{\log \left( \frac{k_{gutabs}}{k_{elim}} \right)}{k_{gutabs} - k_{elim}}}
#'
#' For intravenous route, time of peak concentration is always 0.
#'
#'  ## Peak concentration
#'
#' Evaluate [cp_1comp()] at the time of peak concentration.
#'
#' ## AUC evaluated at infinite time
#'
#' Evaluate [auc_1comp()] at time = `Inf`.
#'
#' ## AUC evaluated at the time of the last observation
#'
#' Evaluate [auc_1comp()] at time = `tlast`.
#'
#'
#'
#'
#'@param pars A named vector of model parameters (e.g. from [coef.pk()]).
#'@param route Character: The route for which to compute TK stats. Currently only "oral" and "iv" are
#'  supported.
#'@param medium Character: the media (tissue) for which to compute TK stats. Currently only "blood" and
#'  "plasma" are supported.
#'@param dose Numeric: A dose for which to calculate TK stats.
#'@param time.units Character: the units of time used for the
#'  parameters `par`. For example, if `par["kelim"]` is in units of 1/weeks,
#'  then `time.units = "weeks"`. If `par["kelim"]` is in units of 1/hours, then
#'  `time.units = "hours"`. This is used to calculate the steady-state plasma/blood concentration for long-term daily dosing of 1 mg/kg/day.
#'@return A `data.frame` with two variables:
#' - `param_name` = `c("CLtot", "CLtot/Fgutabs", "Css", "halflife", "tmax", "Cmax", "AUC_infinity")`
#' - `param_value` = The corresponding values for each statistic (which may be NA if that statistic could not be computed).
#' @export
#' @author John Wambaugh, Caroline Ring
tkstats_1comp <- function(pars,
                          route,
                          medium,
                          dose,
                          # tlast,
                          time.units){

  missing_pars <- setdiff(`1comp`$params,
                          names(pars))
  pars[missing_pars] <- NA_real_

  kelim <- pars["kelim"]
  Fgutabs <- pars["Fgutabs"]
  Vdist <- pars["Vdist"]
  Fgutabs_Vdist <- pars["Fgutabs_Vdist"]
  kgutabs <- pars["kgutabs"]
  Rblood2plasma <- pars["Rblood2plasma"]

  if(is.na(Fgutabs_Vdist) &
     !is.na(Fgutabs) &
     !is.na(Vdist)){
    Fgutabs_Vdist <- Fgutabs/Vdist
  }

  CLtot <- kelim * Vdist

  CLtot_Fgutabs <- kelim / Fgutabs_Vdist

  #convert dose interval of (1/day) into time units
  dose_int <- convert_time(x = 1,
               from = "days",
               to = time.units,
               inverse = TRUE)

  Css_1mgkgday <- ifelse(route %in% "oral",
                      Fgutabs_Vdist / kelim / dose_int,
                      1/(kelim * Vdist * dose_int)) *
    ifelse(medium %in% "blood",
           Rblood2plasma,
           1)

  halflife <- log(2) / kelim

  tmax <- ifelse(route %in% "oral",
                 log(kgutabs / kelim) / (kgutabs - kelim),
                 0)

  Cmax <- cp_1comp(params = list(
    "kelim" = kelim,
    "Vdist" = Vdist,
    "Fgutabs" = Fgutabs,
    "Fgutabs_Vdist" = Fgutabs_Vdist,
    "kgutabs" = kgutabs,
    "Rblood2plasma" = Rblood2plasma
  ),
  time = tmax,
  dose = dose,
  route= route,
  medium = medium)

  AUC_inf <- auc_1comp(params = list(
    "kelim" = kelim,
    "Vdist" = Vdist,
    "Fgutabs" = Fgutabs,
    "Fgutabs_Vdist" = Fgutabs_Vdist,
    "kgutabs" = kgutabs,
    "Rblood2plasma" = Rblood2plasma
  ),
  time = Inf,
  dose = dose,
  route = route,
  medium = medium)

  # AUC_tlast <- auc_1comp(params = list(
  #   "kelim" = kelim,
  # "Vdist" = Vdist,
  # "Fgutabs" = Fgutabs,
  #   "Fgutabs_Vdist" = Fgutabs_Vdist,
  #   "kgutabs" = kgutabs,
  #   "Rblood2plasma" = Rblood2plasma
  # ),
  # time = tlast,
  # dose = dose,
  # route = route,
  # medium = medium)

  return(data.frame(param_name = c("CLtot",
                                   "CLtot/Fgutabs",
                                   "Css_1mgkgday",
                                   "halflife",
                                   "tmax",
                                   "Cmax",
                                   "AUC_infinity"
                                   # "AUC_tlast"
                                   ),
                    param_value = c(CLtot,
                                    CLtot_Fgutabs,
                                    Css_1mgkgday,
                                    halflife,
                                    tmax,
                                    Cmax,
                                    AUC_inf
                                    # AUC_tlast
                                    )))

}
