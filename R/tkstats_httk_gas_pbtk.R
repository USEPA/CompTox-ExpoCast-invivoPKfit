#' TK stats for gas_pbtk model fit with inivivoPKfit
#'
#' @param pars A named numeric vector of model parameters (e.g. from [coef.pk()]).
#' @param route Character: The route for which to compute TK stats. Currently
#'  only "oral" and "iv" are supported.
#' @param medium Character: the media (tissue) for which to compute TK stats.
#'  Currently only "blood" and "plasma" are supported.
#' @param dose Numeric: A dose for which to calculate TK stats.
#' @param time_unit Character: the units of time.
#' @param conc_unit Character: The units of concentration.
#' @param vol_unit Character: The units of dose.
#' @param this_chem Character: The DTXSID of a chemical.
#' @param this_species Character: The species of a subject.
#' @param ... Additional arguments not currently in use.
#' @return A `data.frame` with two variables:
#' - `param_name` = `c("CLtot", "CLtot/Fgutabs", "Css_1mgkgday", "halflife", "Cmax", "AUC_infinity")`
#' - `param_value` = The corresponding values for each statistic (which may be NA if that statistic could not be computed).
#' @export
#' @author Gilberto Padilla Mercado
tkstats_httk_gas_pbtk <- function(pars,
                                  route,
                                  medium,
                                  dose,
                                  time_unit,
                                  conc_unit,
                                  vol_unit,
                                  this_chem,
                                  this_species,
                                  ...) {
  this_chem <- unique(this_chem)
  this_species <- unique(this_species)

  Css <- httk::calc_analytic_css(parameters = recalculate_httk_pbtk_params(pars),
                          model = "pbtk",
                          species = this_species,
                          dose = dose,
                          output.units = "mg/L",
                          suppress.messages = TRUE)

  CLtot <- httk::calc_total_clearance(parameters = recalculate_httk_pbtk_params(pars),
                                      model = "pbtk",
                                      species = this_species,
                                      suppress.messages = TRUE)


  Vss <- httk::calc_vdist(parameters = recalculate_httk_pbtk_params(pars),
                          species = this_species,
                          suppress.messages = TRUE)
  halflife <- Vss * log(2) / CLtot

  # At about 20 half-lives, which accounts for 99.99% of original dose
  time_vct <- unique(signif(seq(0, halflife * 20, length.out = 100), 6))

  predictions <- cp_httk_gas_pbtk(
    params = pars,
    time = time_vct[1:40], # Only need two first 2 half-lives
    dose,
    route,
    medium,
    this_chem,
    this_species
  )
  Cmax <- max(predictions)

  tmax <- time_vct[which.max(predictions)]

  AUC_infinity <- max(auc_httk_gas_pbtk(
    params = pars,
    time = time_vct,
    dose,
    route,
    medium,
    this_chem,
    this_species
  ))

  pars <- recalculate_httk_pbtk_params(pars)
  Fabsgut <- NULL
  list2env(as.list(pars), envir = as.environment(-1))

  CLtot_Fgutabs <- CLtot / Fabsgut
  Vss_Fgutabs <- Vss / Fabsgut



  return(data.frame(
    "param_name" = c("CLtot",
                     "CLtot/Fgutabs",
                     "Css",
                     "halflife",
                     "tmax",
                     "Cmax",
                     "AUC_infinity",
                     "Vss",
                     "Vss/Fgutabs"),
    "param_value" = c(CLtot,
                      CLtot_Fgutabs,
                      Css,
                      halflife,
                      tmax,
                      Cmax,
                      AUC_infinity,
                      Vss,
                      Vss_Fgutabs),
    param_units = c("L/hour/kg BW",
                    "L/hour/kg BW",
                    "mg/L",
                    "hours",
                    "hour",
                    "mg/L",
                    "mg/L/hour",
                    "L",
                    "L")
  )
  )
}
