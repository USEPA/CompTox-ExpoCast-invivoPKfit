#' Make predictions for all model fits
#'
#' Predict concentrations for all chemicals, species, analysis types, and models
#'
#' @param cvt_pre Original concentration vs. time data, pre-processed using
#'   [preprocess_data()].
#' @param fit_flat Fitting output for the flat model: the output of [fit_all()]
#'   with `model = "flat"`
#' @param fit_1comp Fitting output for the 1-compartment model: the output of
#'   [fit_all()] with `model = "1compartment"`
#' @param fit_2comp Fitting output for the 2-compartment model: the output of
#'   [fit_all()] with `model = "2compartment"`
#'
#' @return A `data.table` with `nrow(cvt_pre) * 9` rows and `ncol(cvt_pre) + 5`
#'   variables. It contains model predictions for each observation in `cvt_pre`,
#'   for each analysis type ("Joint", "Separate", and "Pooled") and each model
#'   ("flat", "1compartment", and "2compartment"). It has the same variables as
#'   `cvt_pre`, plus additional variables `Analysis_Type`, `Studies.Analyzed`,
#'   `model`, `pred_conc`, and `pred_auc`. `Analysis_Type`, `Studies.Analyzed`,
#'   and `model` specify which analysis type, which set of studies, and which
#'   model were used to produce the model fit for which predictions are made.
#'   `pred_conc` contains the predicted concentration (in mg/L), and `pred_auc`
#'   contains the predicted cumulative AUC (area under the concentration-time
#'   curve evaluated from time 0 to the time in the corresponding row, units of
#'   mg/L*hour), for the dose, time, route, and medium corresponding to the
#'   observation in each row.
#' @author Caroline Ring
#' @export
#' @import data.table
get_all_predictions <-  function(cvt_pre,
           fit_flat,
           fit_1comp,
           fit_2comp,
           optimx_args = list(
             "method" = "bobyqa",
             "itnmax" = 1e6,
             "control" = list("kkt" = FALSE)
           )){

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

    #get model predictions
    cvt_pred <- fit_combs[, {
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
      get_predictions(pk_fit = pk_fit,
                      newdata = newdata,
                      DTXSID_in = DTXSID_in,
                      Species_in = Species_in,
                      Analysis_Type_in = Analysis_Type_in,
                      Studies.Analyzed_in = Studies.Analyzed_in,
                      model_in = model_in)
    },
    by = .(DTXSID_in,
           Species_in,
           Analysis_Type_in,
           Studies.Analyzed_in,
           model_in)]

setnames(cvt_pred,
         names(cvt_pred),
         gsub(x = names(cvt_pred),
              pattern = "_in$",
              replacement = ""))

return(cvt_pred)
}
