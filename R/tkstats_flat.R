#'TK stats for flat model
#'
#'@param pars A named numeric vector of model parameters (e.g. from [coef.pk()]).
#'@param route Character: The route for which to compute TK stats. Currently
#'  only "oral" and "iv" are supported.
#'@param medium Character: the media (tissue) for which to compute TK stats.
#'  Currently only "blood" and "plasma" are supported.
#'@param dose Numeric: A dose for which to calculate TK stats.
#'@param time_unit Character: the units of time.
#'@param conc_unit Character: The units of concentration.
#'@param vol_unit Character: The units of dose.
#'@return A `data.frame` with two variables:
#' - `param_name` = `c("CLtot", "CLtot/Fgutabs", "Css_1mgkgday", "halflife", "Cmax", "AUC_infinity")`
#' - `param_value` = The corresponding values for each statistic (which may be NA if that statistic could not be computed).
#'@export
#'@author John Wambaugh, Caroline Ring
tkstats_flat <- function(pars,
                         route,
                         medium,
                         dose,
                         time_unit,
                         conc_unit,
                         vol_unit,
                         ...){
  params <- fill_params_flat(pars)

  check_msg <- check_params_flat(params = params,
                                 route = route,
                                 medium = medium)

  if(!(check_msg %in% "Parameters OK")){
    stop(paste("cp_flat():",
               check_msg))
  }

  #for readability, assign params to variables inside this function
  list2env(as.list(params), envir = as.environment(-1))


  Css <- cp_flat(params = pars,
                          time = Inf,
                          dose = dose,
                          route = route,
                          medium = medium)

 AUC_inf <- auc_flat(params = pars,
                          time = Inf,
                          dose = dose,
                          route = route,
                          medium = medium)

  return(data.frame("param_name" = c(
    "CLtot",
    "CLtot/Fgutabs",
    "Css",
    "halflife",
    "Cmax",
    "AUC_infinity",
    "Vss",
    "Vss/Fgutabs"),
                    "param_value" = c(
                      0,
                      0,
                      Css,
                      Inf,
                      Css,
                      AUC_inf,
                      Vdist,
                      1/(Fgutabs_Vdist)),
                      param_units = c(paste0(vol_unit, "/", time_unit),
                                      paste0(vol_unit, "/", time_unit),
                                      conc_unit,
                                      time_unit,
                                      time_unit,
                                      paste0(conc_unit, " * ", time_unit),
                                      vol_unit,
                                      vol_unit)
    )
    )
}
