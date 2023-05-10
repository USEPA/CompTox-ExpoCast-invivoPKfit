#' Evaluate model predictions of toxicokinetic (TK) summary statistics
#'
#' Evaluate how well a fitted toxicokinetic (TK) model predicts quantities that
#' summarize kinetics.
#'
#' Evaluate how well a fitted TK model predicts summary statistics such as:
#'
#' \eqn{t_{\textrm{max}}}: the time of peak concentration for a single oral
#' bolus dose
#'
#' \eqn{C_{\textrm{max}}(1 \textrm{mg/kg})}: the peak concentration for an oral
#' bolus dose of 1 mg/kg
#'
#' \eqn{AUC_{tlast}(1 \textrm{mg/kg})}: the area under the concentration-time
#' curve for a bolus dose of 1 mg/kg (IV or oral), evaluated at the last
#' observed time point
#'
#' These statistics are selected because they can be estimated from observed TK
#' data (observations of concentration vs. time in plasma and/or blood after a
#' single bolus dose, IV and/or oral).
#'
#' @param obs_data A `data.table` of observed concentration vs. time data, e.g.
#'   as produced by [preprocess_data()]. Must include variables `Time`, `Dose`,
#'   `Conc`, `Conc_SD`, `LOQ`, `N_Subjects`, `Route`, and `Media`.
#' @param fit_flat A `data.table` of fitting output for the flat model: the output of `fit_all()`
#'   with `model = "flat"`
#' @param fit_1comp A `data.table` of fitting output for the 1-compartment model: the output of
#'   `fit_all()` with `model = "1compartment"`
#' @param fit_2comp A `data.table` of fitting output for the 2-compartment model: the output of
#'   `fit_all()` with `model = "2compartment"`
#' @param group_by A character vector of variables to in `cvt_pre` to group data
#'   by. NCA stats will be computed for each group defined by a unique
#'   combination of the values of these variables. Default is `c("DTXSID",
#'   "Species", "Route", "Media")`. You could also use `c("DTXSID", "Species",
#'   "Route", "Media", "Dose")` to check dose-specific parameters.
#' @param dose_norm Logical: TRUE (the default) means to do NCA on dose-normalized
#'   concentrations (`Conc_Dose` in `cvt_pre`), and FALSE means to do NCA on
#'   non-dose-normalized concentrations (`Conc` in `cvt_pre`).
#' @return A `data.table` of goodness-of-fit measures for each dataset, each
#'   model, and each analysis (joint, separate, and pooled)
#' @author Caroline Ring
#' @export
#' @import data.table
evaluate_tkstats <- function(cvt_pre,
                             fit_flat,
                             fit_1comp,
                             fit_2comp,
                             group_by = c("DTXSID",
                                          "Species",
                                          "Route",
                                          "Media"),
                             dose_norm = TRUE){

  #get postprocessed parameter tables
  #these contain model-predicted tmax, Cmax for 1 mg/kg
  pk_fit <- merge_fits(fit_flat = fit_flat,
                       fit_1comp = fit_1comp,
                       fit_2comp = fit_2comp)


#get TK stats: tmax, Cmax, AUC
  nca_DT <- get_tkstats(cvt_pre = cvt_pre,
                        group_by = group_by,
                        dose_norm = dose_norm)

  #Merge
  nca_fit_DT <- nca_DT[pk_fit,
                       on = intersect(names(nca_DT),
                                              names(pk_fit)),
                       allow.cartesian = TRUE]

  return(nca_fit_DT)

}

#' Get non-compartmental TK statistics from observed data
#'
#' @param cvt_pre A `data.table` of concentration vs. time data, preprocessed as
#'   with [preprocess_data()]. Must contain variables Time, Conc_Dose, DTXSID,
#'   Species, Route, Media, Subject, Study_ID, and Series_ID.
#' @param group_by A character vector of variables to in `cvt_pre` to group data
#'   by. NCA stats will be computed for each group defined by a unique
#'   combination of the values of these variables. Default is `c("DTXSID",
#'   "Species", "Route", "Media")`. You could also use `c("DTXSID", "Species",
#'   "Route", "Media", "Dose")` to check dose-specific parameters.
#' @param dose_norm Logical: TRUE (the default) means to do NCA on dose-normalized
#'   concentrations (`Conc_Dose` in `cvt_pre`), and FALSE means to do NCA on
#'   non-dose-normalized concentrations (`Conc` in `cvt_pre`).
#' @return A `data.table` of non-compartmental statistics, with variables
#'   DTXSID, Species, Route, Media, tmax (time of peak concentration),
#'   Cmax_1mgkg (peak concentration at 1 mg/kg single bolus dose),
#'   AUC_tlast_1mgkg (AUC at last observed time point for 1 mg/kg single bolus
#'   dose), AUC_inf_1mgkg (AUC extrapolated to infinite time for 1 mg/kg single
#'   bolus dose), AUMC_inf_1mgkg (area under first moment curve extrapolated to
#'   infinite time for 1 mg/kg single bolus dose), MRT (maximum residence time,
#'   only for IV administration), halflife (terminal half-life, only for IV
#'   administration), Clearance (clearance rate), and Vss (volume of
#'   distribution at steady state, only for IV administration),
#'   Clearance_Fgutabs (clearance rate divided by oral bioavailability, only for
#'   oral administration), and MTT (maximum transit time, only for oral
#'   administration).
#' @author Caroline Ring
#'
get_tkstats <- function(cvt_pre,
                        group_by = c("DTXSID",
                                     "Species",
                                     "Route",
                                     "Media"),
                        dose_norm = TRUE){
  # observed tmax & Cmax values by chemical/species datasets that have oral data
  if(dose_norm %in% TRUE){
  max_DT <- cvt_pre[Route %in% "po",
                    get_peak(x = Time,
                             y = Conc_Dose),
                    by = group_by]
  setnames(max_DT,
           c("x", "y"),
           c("tmax", "Cmax_1mgkg"))
  }else{
    max_DT <- cvt_pre[Route %in% "po",
                      get_peak(x = Time,
                               y = Conc),
                      by = group_by]
    setnames(max_DT,
             c("x", "y"),
             c("tmax", "Cmax"))
  }

  #observed dose-normalized AUC at last time point
  #this needs to go by Route and Media as well
  nca_DT <- cvt_pre[, get_nca(obs_data = .SD,
                              dose_norm = dose_norm),
                    by = group_by,
                    .SDcols = names(cvt_pre)]

#halflife and Vss estimates are not valid for oral data
  nca_DT[Route %in% "po", c("halflife", "Vss") := NA_real_]
  #For oral data, clearance is actually Clearance/Fgutabs
  nca_DT[Route %in% "po", Clearance_Fgutabs := Clearance]
  nca_DT[Route %in% "po", Clearance := NA_real_]
  #For oral data, MRT is actually MTT
  nca_DT[Route %in% "po", MTT := MRT]
  nca_DT[Route %in% "po", MRT := NA_real_]

  out_DT <- merge(max_DT, nca_DT, by = group_by,
                  all = TRUE)

  return(out_DT)
}

