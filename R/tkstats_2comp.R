#' Toxicokinetic statistics for 1-compartment model
#'
#' Calculate predicted toxicokinetic statistics for a 1-compartment model.
#'
#' @section Statistics computed:
#'
#' \subsection{Total clearance}{
#' \deqn{\textrm{CL}_{tot} = k_{elim} + V_{1}}
#' }
#' \subsection{Steady-state plasma concentration for long-term daily dose of 1 mg/kg/day}{
#' To convert to steady-state *blood* concentration, multiply by the
#' blood-to-plasma ratio.
#'
#' The dosing interval \eqn{\tau = \frac{1}{\textrm{day}}} will be converted to
#' the same units as \eqn{k_{elim}}.
#' }
#' \subsection{Oral route}{
#' \deqn{C_{ss} = \frac{F_{gutabs} V_{1}}{k_{elim} \tau}}
#' }
#' \subsection{Intravenous route}{
#' \deqn{C_{ss} = \frac{1}{\textrm{CL}_{tot} \tau}}
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
#' Evaluate [cp_1comp()] at the time of peak concentration.
#' }
#' \subsection{AUC evaluated at infinite time}{
#' Evaluate [auc_1comp()] at time = `Inf`.
#' }
#' \subsection{AUC evaluated at the time of the last observation}{
#' Evaluate [auc_1comp()] at time = `tlast`.
#' }
#'
#' @inheritParams tkstats_1comp
#' @return A `data.frame` with two variables:
#' - `param_name` = `c("CLtot", "CLtot/Fgutabs", "Css", "halflife", "tmax", "Cmax", "AUC_infinity", "A", "B", "alpha", "beta", "Vbeta", "Vbeta_Fgutabs", "Vss", "Vss_Fgutabs")`
#' - `param_value` = The corresponding values for each statistic (which may be NA if that statistic could not be computed; e.g. all of the `"x_Fgutabs"` parameters can only be computed if `route = "oral"` ).
#' @export
#' @author John Wambaugh, Caroline Ring
#' @family built-in model functions
#' @family 2-compartment model functions
tkstats_2comp <- function(pars,
                          route,
                          medium,
                          dose,
                          time_unit,
                          conc_unit,
                          vol_unit,
                          ...) {

  params <- fill_params_2comp(pars)

  # get transformed parameters for 2-comp model
  # these are parameters such alpha, beta
  trans_params <- transformed_params_2comp(params = params)

  Fgutabs_V1 = Rblood2plasma = Fgutabs = V1 = NULL
  kelim = kgutabs = k12 = k21 = NULL
  beta = alpha = NULL

  list2env(as.list(params), envir = as.environment(-1))
  list2env(as.list(trans_params), envir = as.environment(-1))

  CLtot <- kelim * V1

  CLtot_Fgutabs <- kelim / Fgutabs_V1



Vbeta <- V1 * kelim / beta
Vbeta_Fgutabs <- (1 / Fgutabs_V1) * kelim / beta

Vss <- V1 * (k21 + k12) / k21
Vss_Fgutabs <- (1 / Fgutabs_V1) * (k21 + k12) / k21


  # convert dose interval of (1/day) into time units
  # this is now standardized because time_units will always be hours
  dose_int <- 1 / 24

  Css <- dose * ifelse(route %in% "oral",
                      Fgutabs_V1 / kelim / dose_int,
                      1 / (kelim * V1 * dose_int)) *
    ifelse(medium %in% "blood",
           Rblood2plasma,
           1)

  halflife_terminal <- log(2) / beta

  tmax <- ifelse(route %in% "oral",
                 tryCatch(
                   uniroot(f = function(x) {
                     cp_2comp_dt(params = pars,
                                 time = x,
                                 dose = dose,
                                 route = "oral",
                                 medium = medium)
                   },
                   lower = 0,
                   upper = 2 * log(kgutabs / kelim) / (kgutabs - kelim), # tmax for 1-compartment
                   extendInt = "downX", # function should be decreasing
                   maxiter = 1000,
                   tol = .Machine$double.eps)$root,
                   error = function(err) return(NA_real_)),
                 0)

  Cmax <- cp_2comp(params = pars,
  time = tmax,
  dose = dose,
  route = route,
  medium = medium,
  ...)

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
                                   "Vss",
                                   "Vss/Fgutabs"),
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
