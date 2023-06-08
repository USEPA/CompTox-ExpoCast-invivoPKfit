#'Toxicokinetic statistics for 1-compartment model
#'
#'Calculate predicted toxicokinetic statistics for a 1-compartment model.
#'
#'# Statistics computed
#'
#'## Total clearance
#'
#'\deqn{\textrm{CL}_{tot} = k_{elim} + V_{1}}
#'
#'## Steady-state plasma concentration for long-term daily dose of 1 mg/kg/day
#'
#'To convert to steady-state *blood* concentration, multiply by the
#'blood-to-plasma ratio.
#'
#'The dosing interval \eqn{\tau = \frac{1}{\textrm{day}}} will be converted to
#'the same units as \eqn{k_{elim}}.
#'
#'### Oral route
#'
#'\deqn{C_{ss} = \frac{F_{gutabs} V_{1}}{k_{elim} \tau}}
#'
#'### Intravenous route
#'
#'\deqn{C_{ss} = \frac{1}{\textrm{CL}_{tot} \tau}}
#'
#'## Half-life of elimination
#'
#'\deqn{\textrm{Halflife} = \frac{\log(2)}{k_{elim}}}
#'
#'## Time of peak concentration
#'
#'For oral route:
#'
#'\deqn{\frac{\log \left( \frac{k_{gutabs}}{k_{elim}} \right)}{k_{gutabs} -
#'k_{elim}}}
#'
#'For intravenous route, time of peak concentration is always 0.
#'
#'## Peak concentration
#'
#'Evaluate [cp_1comp()] at the time of peak concentration.
#'
#'## AUC evaluated at infinite time
#'
#'Evaluate [auc_1comp()] at time = `Inf`.
#'
#'## AUC evaluated at the time of the last observation
#'
#'Evaluate [auc_1comp()] at time = `tlast`.
#'
#'
#'
#'
#'@param pars A named numeric vector of model parameters (e.g. from [coef.pk()]).
#'@param route Character: The route for which to compute TK stats. Currently
#'  only "oral" and "iv" are supported.
#'@param medium Character: the media (tissue) for which to compute TK stats..
#'  Currently only "blood" and "plasma" are supported.
#'@param dose Numeric: A dose for which to calculate TK stats.
#'@param time_unit Character: the units of time used for the parameters `par`.
#'  For example, if `par["kelim"]` is in units of 1/weeks, then `time_unit =
#'  "weeks"`. If `par["kelim"]` is in units of 1/hours, then `time_unit =
#'  "hours"`. This is used to calculate the steady-state plasma/blood
#'  concentration for long-term daily dosing of 1 mg/kg/day.
#'@param conc_unit Character: The units of concentration.
#'@param dose_unit Character: The units of dose.
#'@return A `data.frame` with two variables:
#' - `param_name` = `c("CLtot", "CLtot/Fgutabs", "Css", "halflife", "tmax", "Cmax", "AUC_infinity", "A", "B", "alpha", "beta", "Vbeta", "Vbeta_Fgutabs", "Vss", "Vss_Fgutabs")`
#' - `param_value` = The corresponding values for each statistic (which may be NA if that statistic could not be computed; e.g. all of the `"x_Fgutabs"` parameters can only be computed if `route = "oral"` ).
#'@export
#'@author John Wambaugh, Caroline Ring
#' @family built-in model functions
#' @family 2-compartment model functions
tkstats_2comp <- function(pars,
                          route,
                          medium,
                          dose,
                          time_unit,
                          conc_unit,
                          vol_unit,
                          ...){


  params <- fill_params_2comp(pars)

  #get transformed parameters for 2-comp model
  trans_params <- transformed_params_2comp(params = pars)

  #for readability, assign params to variables inside this function
  for(x in names(params)){
    assign(x, unname(params[x]))
  }

  #for readability, assign transformed params to variables inside this function
  for(x in names(trans_params)){
    assign(x, unname(trans_params[x]))
  }

  CLtot <- kelim * V1

  CLtot_Fgutabs <- kelim / Fgutabs_V1



Vbeta <-  V1 * kelim / beta
Vbeta_Fgutabs <- (1/Fgutabs_V1) * kelim / beta

Vss <- V1 * (k21 + k12) / k21
Vss_Fgutabs <- (1/Fgutabs_V1) * (k21 + k12) / k21


  #convert dose interval of (1/day) into time units
  dose_int <- convert_time(x = 1,
               from = "days",
               to = time_unit,
               inverse = TRUE)

  Css <- dose * ifelse(route %in% "oral",
                      Fgutabs_V1 / kelim / dose_int,
                      1/(kelim * V1 * dose_int)) *
    ifelse(medium %in% "blood",
           Rblood2plasma,
           1)

  halflife_terminal <- log(2) / beta

  tmax <- ifelse(route %in% "oral",
                 tryCatch(
                   uniroot( f = function(x){
                     cp_2comp_dt(params = pars,
                                 time = x,
                                 dose = dose,
                                 route = "oral",
                                 medium = medium)
                   },
                   lower = 0,
                   upper = 2 * log(kgutabs / kelim) / (kgutabs - kelim), #tmax for 1-compartment
                   extendInt = "downX", #function should be decreasing
                   maxiter = 1000,
                   tol = .Machine$double.eps)$root,
                   error = function(err) return(NA_real_)),
                 0)

  Cmax <- cp_2comp(params = pars,
  time = tmax,
  dose = dose,
  route= route,
  medium = medium)

  AUC_inf <- auc_2comp(params = pars,
  time = Inf,
  dose = dose,
  route = route,
  medium = medium)


  return(data.frame(param_name = c("CLtot",
                                   "CLtot/Fgutabs",
                                   "Css",
                                   "halflife",
                                   "tmax",
                                   "Cmax",
                                   "AUC_infinity",
                                   "Vdist_ss",
                                   "Vdist_ss/Fgutabs"),
                    param_value = c(CLtot,
                                    CLtot_Fgutabs,
                                    Css,
                                    halflife_terminal,
                                    tmax,
                                    Cmax,
                                    AUC_inf,
                                    Vss,
                                    Vss_Fgutabs),
                    param_units = c(paste0(vol_unit, "/", time_unit),
                                    paste0(vol_unit, "/", time_unit),
                                    conc_unit,
                                    time_unit,
                                    time_unit,
                                    conc_unit,
                                    paste0(conc_unit, " * ", time_unit),
                                    vol_unit,
                                    vol_unit)
                    )
         )
}
