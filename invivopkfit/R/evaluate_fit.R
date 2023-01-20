#' Evaluate goodness of fit
#'
#' Compute multiple metrics to evaluate how well the various TK models fit the
#' data.
#'
#' @param cvt_pre Original concentration vs. time data, pre-processed using
#'   `preprocess_data()`.
#' @param fit_flat Fitting output for the flat model: the output of `fit_all()`
#'   with `model = "flat"`
#' @param fit_1comp Fitting output for the 1-compartment model: the output of
#'   `fit_all()` with `model = "1compartment"`
#' @param fit_2comp Fitting output for the 2-compartment model: the output of
#'   `fit_all()` with `model = "2compartment"`
#' @param fit_dose_norm Logical: Whether the model fits were done on
#'   dose-normalized data (TRUE), or not (FALSE). This needs to be the same as
#'   whatever you used when calling `fit_all()`, or the goodness-of-fit metrics
#'   will be incorrect!
#' @return A `data.table` of goodness-of-fit measures for each dataset, each
#'   model, and each analysis (joint, separate, and pooled)
#' @author Caroline Ring
#' @export
#' @import data.table patchwork
evaluate_fit <- function(cvt_pre,
                         fit_flat,
                         fit_1comp,
                         fit_2comp,
                         fit_conc_dose = TRUE,
                         optimx_args = list(
                           "method" = "bobyqa",
                           "itnmax" = 1e6,
                           "control" = list("kkt" = FALSE)
                         )){

 pred_DT <- get_all_predictions(cvt_pre =cvt_pre,
                               fit_flat = fit_flat,
                               fit_1comp = fit_1comp,
                               fit_2comp = fit_2comp,
                               fit_conc_dose = fit_conc_dose,
                               optimx_args = optimx_args)

 pred_DT[, resid := Conc - pred_conc]
 pred_DT[, resid_Dose := Conc_Dose - pred_conc/Dose]

 #calculate RMSE and R-squared by dataset, analysis, and model
rmse_rsq <- pred_DT[Dose > 0, .(RMSE = calc_rmse(group_mean = Conc,
                              group_sd = Conc_SD,
                              group_N = N_Subjects,
                              group_LOQ = LOQ,
                             group_pred = pred_conc),
            rsq = calc_rsq(group_mean = Conc,
                           group_sd = Conc_SD,
                           group_N = N_Subjects,
                           group_LOQ = LOQ,
                           group_pred = pred_conc),
            RMSE_Dose = calc_rmse(group_mean = Conc_Dose,
                                  group_sd = Conc_SD_Dose,
                                  group_N = N_Subjects,
                                  group_LOQ = LOQ_Dose,
                                  group_pred = pred_conc/Dose),
            rsq_Dose = calc_rsq(group_mean = Conc_Dose,
                                  group_sd = Conc_SD_Dose,
                                  group_N = N_Subjects,
                                  group_LOQ = LOQ_Dose,
                                  group_pred = pred_conc/Dose)),
        by = .(DTXSID, Species, Analysis_Type, Studies.Analyzed, model)]

 pk_fit <- merge_fits(fit_flat = fit_flat,
                      fit_1comp = fit_1comp,
                      fit_2comp = fit_2comp)

 #Enumerate combinations of chemical, species, analysis type, studies analyzed,
 #and model
 fit_all <- rbindlist(list("flat" = unique(fit_flat[, .(DTXSID,
                                                        Species,
                                                        Analysis_Type,
                                                        Studies.Analyzed,
                                                        model)]),
                           "1compartment" = unique(fit_1comp[, .(DTXSID,
                                                                 Species,
                                                                 Analysis_Type,
                                                                 Studies.Analyzed,
                                                                 model)]),
                           "2compartment" = unique(fit_2comp[, .(DTXSID,
                                                                 Species,
                                                                 Analysis_Type,
                                                                 Studies.Analyzed,
                                                                 model)])
 ),
 use.names = TRUE,
 fill = TRUE,
 idcol = "fit_model")

 fit_combs <- unique(fit_all[, .(DTXSID,
                                 Species,
                                 Analysis_Type,
                                 Studies.Analyzed,
                                 model)])
 setnames(fit_combs,
          names(fit_combs),
          paste(names(fit_combs),
                "in",
                sep = "_"))

  #How big are "overall" error SDs?

  #Get overall error SDs (pseudo-RMSE?) for various analysis types & models.
  #This is the best-fit standard deviation for the log-likelihood, evaluated
  #holding the PK model parameters constant.
  error_sd <- fit_combs[, {
    studies_list <- strsplit(Studies.Analyzed_in,
                             split = ", ")[[1]]
    #get data.table of newdata to be predicted
    #leave out DTXSID, Species columns
    #since these will be added back in from fit_combs
    newdata <- cvt_pre[DTXSID %in% DTXSID_in &
                         Species %in% Species_in &
                         Study %in% studies_list,
                       .SD,
                       .SDcols = setdiff(names(cvt_pre),
                                         c("DTXSID", "Species"))]
    get_error_sd(pk_fit = pk_fit,
                 newdata = newdata,
                 DTXSID_in = DTXSID_in,
                 Species_in = Species_in,
                 Analysis_Type_in = Analysis_Type_in,
                 Studies.Analyzed_in = Studies.Analyzed_in,
                 model_in = model_in,
                 fit_conc_dose = fit_conc_dose,
                 optimx_args = optimx_args)
  },
  by = .(DTXSID_in,
         Species_in,
         Analysis_Type_in,
         Studies.Analyzed_in,
         model_in)]



  #How well do we predict:
  #tpeak?
  #Cpeak? Dose-normalized or not?
  #AUC at the last timepoint? Dose-normalized or not?
  #AUC infinity? Dose-normalized or not?
  #Average concentration? Dose-normalized or not?

  #Should we evaluate these per chemical/species, per
  #chemical/species/route/media, per chemical/species/route/media/dose? All of the above?



}
