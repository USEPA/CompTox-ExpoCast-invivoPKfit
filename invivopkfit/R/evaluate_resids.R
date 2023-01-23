#' Evaluate goodness of fit based on residuals
#'
#' Evaluate multiple metrics of goodness of fit based on residuals
#'
#' For each dataset (chemical and species) and each fitted model, this function
#' first calculates the predictions and residuals, then evaluates the following:
#'
#' - RMSE (root mean squared error)
#' - R-squared of observations vs. predictions
#' - The p-value from the Breusch-Pagan test of heteroscedasticity
#' - Estimated "overall" error standard deviation (details below)
#'
#' Unique fitted models are defined by unique combination of analysis type
#' ("Joint", "Separate", or "Pooled"), studies analyzed, and model ("flat",
#' "1compartment", or "2compartment").
#'
#' # Root mean squared error (RMSE)
#'
#' RMSE is calculated in [calc_rmse()].
#'
#' # Estimated overall error standard deviation
#'
#' The "overall" error standard deviation is the value for residual error
#' standard deviation ("sigma") that minimizes the log-likelihood function when
#' the TK model parameters are fixed at their fitted values. This is essentially
#' a more rigorous RMSE, because it correctly models the likelihood of
#' non-detect (below-LOQ) observations as the cumulative probability that the
#' model prediction is below the reported LOQ, rather than making any
#' substitutions for non-detects.
#'
#' @param cvt_pre Original concentration vs. time data, pre-processed using
#'   `preprocess_data()`.
#' @param fit_flat Fitting output for the flat model: the output of `fit_all()`
#'   with `model = "flat"`
#' @param fit_1comp Fitting output for the 1-compartment model: the output of
#'   `fit_all()` with `model = "1compartment"`
#' @param fit_2comp Fitting output for the 2-compartment model: the output of
#'   `fit_all()` with `model = "2compartment"`
#' @param fit_conc_dose Logical: Whether the model fits were done on
#'   dose-normalized data (TRUE), or not (FALSE). This needs to be the same as
#'   whatever you used when calling `fit_all()`, or the goodness-of-fit metrics
#'   will be incorrect!
#' @param fit_log_conc Logical: Whether the model fits were done on
#'   log-transformed data (TRUE), or not (FALSE). This needs to be the same as
#'   whatever you used when calling `fit_all()`, or the goodness-of-fit metrics
#'   will be incorrect!
#' @return A `data.table` of goodness-of-fit measures for each dataset, each
#'   model, and each analysis (joint, separate, and pooled)
#' @author Caroline Ring
#' @export
#' @import data.table
evaluate_resids <- function(cvt_pre,
                         fit_flat,
                         fit_1comp,
                         fit_2comp,
                         fit_conc_dose = TRUE,
                         fit_log_conc = FALSE,
                         optimx_args = list(
                           "method" = "bobyqa",
                           "itnmax" = 1e6,
                           "control" = list("kkt" = FALSE)
                         )){


#Calculate RMSE, R-squared of obs vs pred,
  #and heteroscedasticity of residuals

#First: calculate residuals

#Get predicted concentrations & AUCs
 pred_DT <- get_all_predictions(cvt_pre =cvt_pre,
                               fit_flat = fit_flat,
                               fit_1comp = fit_1comp,
                               fit_2comp = fit_2comp,
                               fit_conc_dose = fit_conc_dose,
                               optimx_args = optimx_args)
#Calculate residuals, RMSE, Rsquared, and Breusch-Pagan p-value for
#heteroscedasticity
 if(fit_conc_dose %in% TRUE){
   if(fit_log_conc %in% TRUE){
     pred_DT[, resid := log(Conc_Dose) - log(pred_conc/Dose)]
     pred_DT[, log_Conc_Dose := log(Conc_Dose)]
     resid_DT <- pred_DT[, .(RMSE = calc_rmse(group_mean = Conc_Dose,
                                              group_sd = Conc_SD_Dose,
                                              group_N = N_Subjects,
                                              group_LOQ = LOQ_Dose,
                                              group_pred = pred_conc/Dose,
                                              log = TRUE),
                             rsq = calc_rsq(group_mean = Conc_Dose,
                                            group_sd = Conc_SD_Dose,
                                            group_N = N_Subjects,
                                            group_LOQ = LOQ_Dose,
                                            group_pred = pred_conc/Dose,
                                            log = TRUE),
                             BP_pval = bp_test(varformula = resid^2 ~ log_Conc_Dose,
                                               data = .SD)$p.value
     ),
     by = .(DTXSID, Species,
            Analysis_Type, Studies.Analyzed,
            model)]
   }else{ #if fit_conc_dose %in% TRUE & fit_log_conc %in% FALSE
     pred_DT[, resid := Conc_Dose - pred_conc/Dose]
     resid_DT <- pred_DT[, .(RMSE = calc_rmse(group_mean = Conc_Dose,
                                              group_sd = Conc_SD_Dose,
                                              group_N = N_Subjects,
                                              group_LOQ = LOQ_Dose,
                                              group_pred = pred_conc/Dose,
                                              log = FALSE),
                             rsq = calc_rsq(group_mean = Conc_Dose,
                                            group_sd = Conc_SD_Dose,
                                            group_N = N_Subjects,
                                            group_LOQ = LOQ_Dose,
                                            group_pred = pred_conc/Dose,
                                            log = FALSE),
                             BP_pval = bp_test(varformula = resid^2 ~ Conc_Dose,
                                               data = .SD)$p.value
                             ),
                         by = .(DTXSID, Species,
                                Analysis_Type, Studies.Analyzed,
                                model)]
   }
 }else{ #if fit_conc_dose %in% FALSE
   if(fit_log_conc %in% TRUE){
     pred_DT[, resid := log(Conc) - log(pred_conc)]
     pred_DT[, log_Conc := log(Conc)]
     resid_DT <- pred_DT[, .(RMSE = calc_rmse(group_mean = Conc,
                                              group_sd = Conc_SD,
                                              group_N = N_Subjects,
                                              group_LOQ = LOQ,
                                              group_pred = pred_conc,
                                              log = TRUE),
                             rsq = calc_rsq(group_mean = Conc,
                                            group_sd = Conc_SD,
                                            group_N = N_Subjects,
                                            group_LOQ = LOQ,
                                            group_pred = pred_conc,
                                            log = TRUE),
                             BP_pval = bp_test(varformula = resid^2 ~ log_Conc,
                                               data = .SD)$p.value
     ),
     by = .(DTXSID, Species,
            Analysis_Type, Studies.Analyzed,
            model)]
   }else{ #if fit_conc_dose %in% FALSE & fit_log_conc %in% FALSE
     pred_DT[, resid := Conc - pred_conc]
     resid_DT <- pred_DT[, .(RMSE = calc_rmse(group_mean = Conc,
                                              group_sd = Conc_SD,
                                              group_N = N_Subjects,
                                              group_LOQ = LOQ,
                                              group_pred = pred_conc,
                                              log = FALSE),
                             rsq = calc_rsq(group_mean = Conc,
                                            group_sd = Conc_SD,
                                            group_N = N_Subjects,
                                            group_LOQ = LOQ,
                                            group_pred = pred_conc,
                                            log = FALSE),
                             BP_pval = bp_test(varformula = resid^2 ~ Conc,
                                               data = .SD)$p.value
                             ),
                         by = .(DTXSID, Species,
                                Analysis_Type, Studies.Analyzed,
                                model)]
   }
 }


#merge together all tables of fitted parameters
 pk_fit <- merge_fits(fit_flat = fit_flat,
                      fit_1comp = fit_1comp,
                      fit_2comp = fit_2comp)

 #How big are "overall" error SDs?

 #Get overall error SDs (pseudo-RMSE?) for various analysis types & models.
 #This is the best-fit standard deviation for the log-likelihood, evaluated
 #holding the PK model parameters constant.

 #Enumerate combinations of chemical, species, analysis type, studies analyzed
 fit_combs <- unique(fit_flat[, .(DTXSID,
                                 Species,
                                 Analysis_Type,
                                 Studies.Analyzed)])
#add models to the table
 fit_combs <- fit_combs[, .(model = c("flat",
                                      "1compartment",
                                      "2compartment")),
                        by = names(fit_combs)]

 setnames(fit_combs,
          names(fit_combs),
          paste(names(fit_combs),
                "in",
                sep = "_"))

  error_sd_DT <- fit_combs[, {
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

  #remove the "_in" from variable names
  setnames(error_sd_DT,
           names(error_sd_DT),
           gsub(x = names(error_sd_DT),
                pattern = "_in",
                replacement = ""))


#merge together the diagnostics
  eval_DT <- merge(resid_DT, error_sd_DT,
                   by = c("DTXSID",
                          "Species",
                          "Analysis_Type",
                          "Studies.Analyzed",
                          "model"))

return(eval_DT)

}

#' Breusch-Pagan test of heteroscedasticity
#' Evaluate heteroscedasticity using Breusch-Pagan test
#' The Breusch-Pagan test performs a linear regression of squared residuals on one or more explanatory variables. The test statistic is formed by the R-squared value of this regression multiplied by the number of residuals. This test statistic obeys a chi-squared distribution with S-1 degrees of freedom, where S is the number of explanatory terms (plus an intercept). The p-value of the test is given accordingly.
#'
#' @param varformula A formula for the variance linear model. The left-hand side of this formula must be the squared residuals.
#' @param data A `data.frame` containing all variables referenced in `varformula`.
#' @return A list with the following elements:
#' - `statistic` : The test statistic
#' - `p.value`: The p-value of the test
#' - `r.squared`: The R-squared of the linear regression described in `varformula`
#' - `n`: The number of observations
#' - `df`: The chi-squared degrees of freedom
#' @author Caroline Ring
#' @references https://rpubs.com/cyobero/187387

bp_test <- function(varformula, data){
  bp_lm <- lm(formula = varformula,
              data = data)
  bp_rsq <- summary(bp_lm)$r.squared
  n_obs <- nrow(data)
  qval <- bp_rsq * n_obs
  #critical value
  n_xvar <- length(labels(terms(varformula))) + 1
  pval <- pchisq(q = qval, df = n_xvar - 1)
  return(list(statistic = qval,
              p.value = pval,
              r.squared = bp_rsq,
              n = n_obs,
              df = n_xvar - 1))
}
