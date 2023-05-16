#'TK stats for flat model
#'
#'@param pars A named vector of model parameters (e.g. from [coef.pk()]).
#'@param route Character: The route for which to compute TK stats. Currently
#'  only "oral" and "iv" are supported.
#'@param medium Character: the media (tissue) for which to compute TK stats.
#'  Currently only "blood" and "plasma" are supported.
#'@param dose Numeric: A dose for which to calculate TK stats.
#'@param time.units Character: the units of time used for the parameters `par`.
#'  For example, if `par["kelim"]` is in units of 1/weeks, then `time.units =
#'  "weeks"`. If `par["kelim"]` is in units of 1/hours, then `time.units =
#'  "hours"`. This is used to calculate the steady-state plasma/blood
#'  concentration for long-term daily dosing of 1 mg/kg/day. Not used for flat model.
#'@return A `data.frame` with two variables:
#' - `param_name` = `c("CLtot", "CLtot/Fgutabs", "Css_1mgkgday", "halflife", "Cmax", "AUC_infinity")`
#' - `param_value` = The corresponding values for each statistic (which may be NA if that statistic could not be computed).
#'@export
#'@author John Wambaugh, Caroline Ring
tkstats_flat <- function(pars,
                         route,
                         medium,
                         dose,
                         time.units){
  #the only TK stat for a flat model is Css, since effectively it is "always" at Css
  missing_pars <- setdiff(model_flat$params,
                          names(pars))
  pars[missing_pars] <- NA_real_

  Fgutabs <- pars["Fgutabs"]
  Vdist <- pars["Vdist"]
  Fgutabs_Vdist <- pars["Fgutabs_Vdist"]
  Rblood2plasma <- pars["Rblood2plasma"]
  if(is.na(Fgutabs_Vdist) &
     !is.na(Fgutabs) &
     !is.na(Vdist)){
    Fgutabs_Vdist <- Fgutabs/Vdist
  }

  Css_1mgkgday <- cp_flat(params = as.list(pars[!is.na(pars)]),
                          time = Inf,
                          dose = dose,
                          route = route,
                          medium = medium)

 AUC_inf <- auc_flat(params = as.list(pars[!is.na(pars)]),
                          time = Inf,
                          dose = dose,
                          route = route,
                          medium = medium)

  return(data.frame("param_name" = c(
    "CLtot",
    "CLtot/Fgutabs",
    "Css_1mgkgday",
    "halflife",
    "Cmax",
    "AUC_infinity",
    "Vdist_ss",
    "Vdist_ss/Fgutabs"),
                    "param_value" = c(
                      0,
                      0,
                      Css_1mgkgday,
                      Inf,
                      Css_1mgkgday,
                      AUC_inf,
                      Vdist,
                      1/(Fgutabs_Vdist))
    )
    )
}
